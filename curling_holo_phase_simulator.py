import numpy as np
from scipy.integrate import solve_bvp, quad, cumulative_trapezoid
import matplotlib.pyplot as plt
from PIL import Image
from PIL.PngImagePlugin import PngInfo  # <--- Ajout pour gérer les métadonnées PNG

from scipy.constants import e, hbar, mu_0
import os

# =============================================================================
# 1. PARAMÈTRES PHYSIQUES & CONFIGURATION
# =============================================================================
# Constantes du matériau (Unités SI)
A  = 1.3e-11        # Constante d'échange (J/m)
Ku = -1e4           # Constante d'anisotropie uniaxiale (J/m^3)
Ms = 1.00531 / mu_0 # Aimantation (A/m)

# Géométrie du cylindre
R  = 100e-9       # Rayon extérieur du cylindre (m)
r0 = R * 1e-6     # Rayon de cœur (singularité numérique évitée), très petit devant R

# Paramètres adimensionnels et discrétisation pour solve_bvp
kappa = (Ku * R**2) / (2 * A) # Parametre adimensionnel
epsilon = r0 / R  # Limite inférieure du domaine adimensionnel
N_BVP = 600       # Nombre de points de discrétisation pour l'équation d'Euler-Lagrange

# Résolution de la carte de phase (Taille finale de l'image PNG : N_IMAGE x N_IMAGE)
N_IMAGE = 1000    

# =============================================================================
# 2. SYSTÈME DIFFÉRENTIEL (EULER-LAGRANGE) & CONDITIONS AUX LIMITES
# =============================================================================
def ode_curling(rho, y):
    """Définition de l'EDO d'Euler-Lagrange adimensionnée pour le mode curling."""
    omega, domega_drho = y
    d2omega_drho2 = -domega_drho / rho + (1.0 / (2 * rho**2) + kappa) * np.sin(2 * omega)
    return np.vstack([domega_drho, d2omega_drho2])

def bc_curling(ya, yb):
    """Conditions aux limites : omega(epsilon) = 0 et omega'(1) = 0."""
    return np.array([ya[0] - 0.0, yb[1] - 0.0])

# =============================================================================
# 3. FONCTIONS D'INTÉGRATION POUR LE PROFIL I(y)
# =============================================================================
def integrand_z(z, y, R, sol_bvp):
    """Fonction à intégrer selon l'axe z pour obtenir le profil I(y)."""
    r = np.sqrt(y**2 + z**2)
    rho = r / R
    
    # Sécurité pour respecter le domaine de définition de la spline [epsilon, 1.0]
    rho = np.clip(rho, epsilon, 1.0)
        
    # Évaluation de omega(rho) via l'objet spline 'sol_bvp' généré par solve_bvp
    omega = sol_bvp(rho)[0]
    
    if r == 0:
        return 1.0
    return np.cos(omega)

def compute_I(y, R, sol_bvp):
    """Calcule la projection intégrale I pour une coordonnée y donnée."""
    assert y>=0, 'y must be positive'
    if y>R: y=R
    bound = np.sqrt(R**2 - y**2)
    # Intégration numérique par rapport à z entre -limite et +limite
    res, _ = quad(integrand_z, -bound, bound, args=(y, R, sol_bvp), limit=5000, epsabs=1e-10,
    epsrel=1e-10)
    return -mu_0 * Ms * res

# =============================================================================
# 4. EXÉCUTION DU CALCUL
# =============================================================================
if __name__ == "__main__":
    
    # --- Détermination dynamique du répertoire du script ---
    # Cette variable pointe de manière absolue vers le dossier contenant ce fichier .py
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    # --- Étape A : Estimation de la solution et diagnostic de stabilité ---
    if Ku != 0:
        l_K = np.sqrt(A / abs(Ku))
        Rc = 1.8412 * l_K
        print(f"Rayon critique estimé R_c = {Rc*1e9:.2f} nm (R actuel = {R*1e9:.1f} nm)")
        if Ku > 0 or R < Rc:
            print(">>> Attention : Configuration stable attendue vers la solution triviale omega(rho) = 0.")
        else:
            print(">>> Configuration favorable : Profil de curling non trivial attendu.")

    # --- Étape B : Résolution de l'EDO (BVP) ---
    rho_grid = np.linspace(epsilon, 1.0, N_BVP)
    # Profil 'guess' linéaire pour forcer le solveur vers la solution non triviale
    omega_guess = (np.pi / 2) * (rho_grid - epsilon) / (1.0 - epsilon)
    y_guess = np.vstack([omega_guess, np.gradient(omega_guess, rho_grid)])

    sol = solve_bvp(ode_curling, bc_curling, rho_grid, y_guess, tol=1e-7, max_nodes=1000000)
    print("[Euler-Lagrange / BVP] Statut de convergence :", sol.message)
    omega_bvp = sol.sol(rho_grid)[0]

    # --- Étape C : Calcul du profil d'intégration I(y) ---
    # Vectorisation de la fonction pour une exécution optimisée sur tableau NumPy
    compute_I_vectorized = np.vectorize(lambda y: compute_I(y, R, sol.sol))
    
    # Génération du maillage fin direct (évite de calculer deux fois avec des tailles de grille différentes)
    bound = 2*R
    y_vals = np.linspace(0, bound, N_IMAGE//2+1)
    print(f"[Calcul] Évaluation de I(y) sur {N_IMAGE} points...")
    I_vals = compute_I_vectorized(y_vals)

    # --- Étape D : Calcul de la phase Phi(y) ---
    # Intégration par la méthode des trapèzes
    phase_vals = cumulative_trapezoid(I_vals, y_vals, initial=0)
    phase_vals *= -e / hbar

    # anti-symetrisation
    y_vals = np.concatenate((-y_vals[:1:-1], y_vals))
    phase_vals = np.concatenate((-phase_vals[:1:-1], phase_vals))

    # =============================================================================
    # 5. GÉNÉRATION DES GRAPHIQUES ET SAUVEGARDES
    # =============================================================================
    print("[Graphiques] Génération des tracés de contrôle...")
    
    # Fig 1 : Profil Curling adimensionné
    plt.figure(figsize=(6, 4))
    plt.plot(rho_grid, omega_bvp, '-', lw=2.5, color='tab:orange', label="BVP (Spline)")
    plt.xlabel(r'$\rho = r/R$')
    plt.ylabel(r'$\omega(\rho)$ [rad]')
    plt.title("Profil adimensionné du mode curling")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(SCRIPT_DIR, 'curling_mode_adimensionne.png'), dpi=150)
    plt.show()
    #plt.close()

    # Fig 2 : Profil de la projection intégrale I(y)
    # plt.figure(figsize=(6, 4))
    # plt.plot(y_vals * 1e9, I_vals, '-', lw=2.5, color='tab:blue', label=r'$I(y)$')
    # plt.xlabel(r'$y$ [nm]')
    # plt.ylabel(r'$I(y)$ [UA]')
    # plt.title(r"Profil de l'intégrale $I(y)$")
    # plt.grid(alpha=0.3)
    # plt.tight_layout()
    # plt.savefig(os.path.join(SCRIPT_DIR, 'integrale_I_y.png'), dpi=150)
    # plt.close()

    # Fig 3 : Profil de la Phase Phi(y)
    plt.figure(figsize=(6, 4))
    plt.plot(y_vals * 1e9, phase_vals, '-', lw=2.5, color='tab:purple', label=r'$\Phi(y)$')
    plt.axvline(0, color='black', linestyle='--', alpha=0.4)
    plt.axhline(0, color='black', linestyle='--', alpha=0.4)

    # --- AJOUT DE L'ENCADRÉ DES VALEURS EXTREMA ---
    phase_min, phase_max = np.min(phase_vals), np.max(phase_vals)
    
    # Texte à afficher (utilisant les expressions LaTeX pour les symboles grecs)
    info_text = f"$\\Phi_{{min}}$: {phase_min:.4f} rad\n$\\Phi_{{max}}$: {phase_max:.4f} rad"
    
    # Placement du texte en coordonnées relatives du graphique (0.05 = 5% du bord gauche, 0.80 = 80% du bas)
    plt.text(0.05, 0.80, info_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    # ----------------------------------------------

    plt.xlabel(r'$y$ [nm]')
    plt.ylabel(r'$\Phi(y)$ [rad]')
    plt.title(r'Profil de la phase cumulative $\Phi(y)$')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(SCRIPT_DIR, 'phase_profile_y.png'), dpi=150)
    plt.show()
    #plt.close()

    # =============================================================================
    # 6. EXPORT CARTE DE PHASE INTERFÉROMÉTRIQUE (PNG 16-BITS) + MÉTADONNÉES
    # =============================================================================
    print(f"[Image] Phase min : {phase_min:.4e} | Phase max : {phase_max:.4e}")

    # Normalisation linéaire stricte dans la dynamique de codage 16-bits [0, 65535]
    if phase_max != phase_min:
        phase_normalized = (phase_vals - phase_min) / (phase_max - phase_min) * 65535
    else:
        phase_normalized = np.zeros_like(phase_vals)

    # Mutation explicite du tableau vers le type d'entier non-signé 16 bits
    phase_16bit = phase_normalized.astype(np.uint16)

    # Création de la matrice bidimensionnelle 2D (invariance stricte de la phase selon l'axe z)
    phase_map_2d = np.tile(phase_16bit, (N_IMAGE, 1))

    # Génération de l'image grayscale 16-bits
    img = Image.fromarray(phase_map_2d, mode='I;16')
    
    # --- Création de l'objet Métadonnées ---
    metadata = PngInfo()
    
    # Ajout des paramètres de simulation (les valeurs doivent impérativement être des chaînes 'str')
    metadata.add_text("Simulation_Type", "Curling Mode Phase Map")
    metadata.add_text("Param_A", str(A))
    metadata.add_text("Param_Ku", str(Ku))
    metadata.add_text("Param_Ms", str(Ms))
    metadata.add_text("Param_R", str(R))
    metadata.add_text("Param_r0", str(r0))
    metadata.add_text("Param_kappa", str(kappa))
    metadata.add_text("Param_N_IMAGE", str(N_IMAGE))
    
    # Très important pour pouvoir refaire la simulation ou reconstruire les vraies valeurs de Phase :
    # On stocke les extrema qui ont servi à la normalisation [0, 65535]
    metadata.add_text("Phase_Min_SI", str(phase_min))
    metadata.add_text("Phase_Max_SI", str(phase_max))

    # Sauvegarde forcée dans SCRIPT_DIR en transmettant l'objet pnginfo
    output_image_path = os.path.join(SCRIPT_DIR, 'analytics_HOLO_PHASE.png')
    img.save(output_image_path, pnginfo=metadata)  # <--- Passage des métadonnées ici

    print(f"[Succès] Image 16-bits (avec métadonnées) sauvegardée dans :\n -> {output_image_path}")
