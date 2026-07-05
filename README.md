# Curling Holo Fresnel Suite : Simulation Magnétique & Imagerie Électronique

Ce projet implémente la simulation numérique de la distribution d'aimantation dans un cylindre ferromagnétique (mode curling). À partir de cette configuration micromagnétique, la suite logicielle calcule l'impact du champ magnétique sur la phase d'un faisceau d'électrons par **holographie électronique**, puis simule la propagation en espace libre pour générer le **contraste de Fresnel** (imagerie par défocalisation). 

L'ensemble des résultats est exporté sous forme d'images 2D au format PNG 16-bits incluant des métadonnées physiques détaillées.

## 📋 Table des Matières
- [Physique et Paramètres](#-physique-et-paramètres)
- [Structure Algorithmique](#-structure-algorithmique)
- [Modélisation Mathématique](#-modélisation-mathématique)
- [Dépendances](#-dépendances)

---

## 🔬 Physique et Paramètres

Le système modélise un cylindre ferromagnétique de rayon $R$. La configuration dépend de trois constantes physiques majeures :
* **$A$** : La constante d'échange.
* **$K_u$** : La constante d'anisotropie uniaxiale positive ou négative.
* **$M_s$** : L'aimantation à saturation.

Un paramètre d'anisotropie adimensionnel $\kappa$ est défini par :
$$\kappa = \frac{K_u R^2}{2A}$$

Le domaine spatial radial est discrétisé selon la variable adimensionnelle $\rho = r/R \in [\epsilon, 1]$, où $\epsilon = r_0/R$ introduit un rayon de cœur minimal pour éviter toute singularité numérique à l'origine.

La suite logicielle intègre également les paramètres de la microscopie électronique à transmission (TEM) pour la propagation :
* **$E_0$** : La tension d'accélération des électrons (en eV).
* **$\Delta z$** (defocus) : La distance de défocalisation (en mètres).
* **$C_s$** : La constante d'aberration sphérique de la lentille objectif (en mètres).

---

## 💻 Structure Algorithmique

Le script de simulation est articulé autour de 5 étapes clés :
1. **Résolution du Problème aux Limites (BVP) :** Détermination de l'angle d'aimantation $\omega(\rho)$ à l'aide de l'algorithme `scipy.integrate.solve_bvp`.
2. **Projection Magnétique :** Calcul de l'intégrale de projection $I(y)$ le long de la ligne de visée du faisceau d'électrons ($z$). L'évaluation est vectorisée via `np.vectorize` et calculée par quadrature numérique (`scipy.integrate.quad`).
3. **Intégration de la Phase :** Calcul de la phase cumulative $\Phi(y)$ par la méthode des trapèzes sur l'intervalle $[-2R, 2R]$. On utilise l'anti-symétrie de la phase, imposant $\Phi(0) = 0$ au centre.
4. **Simulation de Fresnel 1D (Optimisée) :** Propagation de l'onde électronique via une transformée de Fourier rapide (FFT) appliquée en 1D sur le profil de phase. Cette étape intègre un **padding miroir horizontal** à gauche afin d'assurer la continuité $C^0$ aux bords.
5. **Génération d'Images 2D & Métadonnées :** Duplication des profils 1D (phase et intensité de Fresnel) via `np.tile` pour générer des matrices 2D exploitant l'invariance stricte selon la direction perpendiculaire. Les deux cartes sont exportées au format PNG 16-bits (`uint16`) en y injectant les métadonnées physiques (`PngInfo`) nécessaires à leur décodage futur.

---

## 📐 Modélisation Mathématique

### 1. Équation d'Euler-Lagrange
La configuration micromagnétique est régie par l'équation différentielle du second ordre suivante :
$$\omega''(\rho) = -\frac{\omega'(\rho)}{\rho} + \left( \frac{1}{2\rho^2} + \kappa \right) \sin\big(2\omega(\rho)\big)$$

Soumise aux conditions aux limites :
$$\omega(\epsilon) = 0 \quad \text{et} \quad \omega'(1) = 0$$

### 2. Intégrale de Projection $I(y)$
Pour chaque abscisse $y \in [0, R]$, l'intégration le long de la ligne de visée géométrique $z_{\text{max}}(y) = \sqrt{R^2 - y^2}$ donne :
$$I(y) = -\mu_0 M_s \int_{-z_{\text{max}}(y)}^{z_{\text{max}}(y)} \cos\big(\omega(\rho)\big) \, \mathrm{d}z$$
où la position radiale locale est $\rho = \frac{\sqrt{y^2 + z^2}}{R}$.

### 3. Calcul de la Phase Cumulative $\Phi(y)$
En intégrant les constantes fondamentales issues de `scipy.constants` (charge élémentaire $e$ et constante de Planck réduite $\hbar$) et en imposant $\Phi(0) = 0$, la phase finale vaut :
$$\Phi(y) = -\frac{e}{\hbar} \int_{0}^{y} I(u) \, \mathrm{d}u$$

### 4. Propagation de Fresnel 1D
L'invariance directionnelle selon l'axe vertical ($ky = 0$) réduit l'équation de propagation. La fonction d'onde initiale paddée en miroir horizontal est $\psi_0(x) = \exp\big(1j \cdot \Phi_{\text{pad}}(x)\big)$.

La propagation dans l'espace des fréquences spatiales $k_x$ est régie par la fonction de transfert des aberrations $H(k_x)$ :
$$H(k_x) = \exp\big(-1j \cdot \chi(k_x)\big)$$
$$\chi(k_x) = \pi \lambda \Delta z k_x^2 + \frac{1}{2}\pi C_s \lambda^3 k_x^4$$
où $\lambda$ est la longueur d'onde relativiste de l'électron calculée à partir de $E_0$. L'intensité finale du contraste de Fresnel 1D est obtenue par :
$$I_{\text{Fresnel}}(x) = \left| \mathcal{TF}^{-1} \big( \mathcal{TF}(\psi_0) \cdot H \big) \right|^2$$

### 5. Images Finales 2D
Les images discrètes $M(j, i)$ de taille $N \times N$ (pour la phase et le contraste de Fresnel) subissent une normalisation linéaire stricte suivie d'une quantification sur 16 bits :
$$M(j, i) = \text{round}\left( 65535 \times \frac{V(i) - V_{\text{min}}}{V_{\text{max}} - V_{\text{min}}} \right)$$
où $V$ représente le vecteur 1D étendu en 2D par invariance verticale via `np.tile`.

---

## 🛠️ Dépendances

Le projet nécessite l'environnement Python standard pour le calcul scientifique :
* **NumPy** (gestion des tableaux, FFT 1D et duplication matricielle)
* **SciPy** (constantes physiques fondamentales, modules `integrate.solve_bvp` et `integrate.quad`)
* **Pillow (PIL)** (génération des fichiers PNG 16-bits et gestion des métadonnées structurelles `PngInfo`)
* **Matplotlib** (génération des graphiques de contrôle et affichage des résultats)
