# Curling Holo Phase Simulator : Simulation de la Phase Holographique Electronique

Ce projet implémente la simulation numérique de la distribution d'aimantation dans un cylindre ferromagnétique (mode curling). À partir de cette configuration micromagnétique, la suite logicielle calcule l'impact du champ magnétique sur la phase d'un faisceau d'électrons par **holographie électronique**. 

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

---

## 💻 Structure Algorithmique

Le script de simulation est articulé autour de 5 étapes clés :
1. **Résolution du Problème aux Limites (BVP) :** Détermination de l'angle d'aimantation $\omega(\rho)$ à l'aide de l'algorithme `scipy.integrate.solve_bvp`.
2. **Projection Magnétique :** Calcul de l'intégrale de projection $I(y)$ le long de la ligne de visée du faisceau d'électrons ($z$). L'évaluation est vectorisée via `np.vectorize` et calculée par quadrature numérique (`scipy.integrate.quad`).
3. **Intégration de la Phase :** Calcul de la phase cumulative $\Phi(y)$ par la méthode des trapèzes sur l'intervalle $[-2R, 2R]$. On utilise l'anti-symétrie de la phase, imposant $\Phi(0) = 0$ au centre.
5. **Génération d'Images 2D & Métadonnées :** Duplication des profils 1D (phase et épaisseur) via `np.tile` pour générer des matrices 2D exploitant l'invariance stricte selon la direction perpendiculaire. Les deux cartes sont exportées au format PNG 16-bits (`uint16`) en y injectant les métadonnées physiques (`PngInfo`) nécessaires à leur décodage futur.

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

### 3. Calcul de la Phase Holographique $\Phi(y)$
En intégrant les constantes fondamentales issues de `scipy.constants` (charge élémentaire $e$ et constante de Planck réduite $\hbar$) et en imposant $\Phi(0) = 0$, la phase finale vaut :
$$\Phi(y) = -\frac{e}{\hbar} \int_{0}^{y} I(u) \, \mathrm{d}u$$

### 4. Images Finales 2D
Les images discrètes $M(j, i)$ de taille $N \times N$ (pour la phase et l'épaisseur) subissent une normalisation linéaire stricte suivie d'une quantification sur 16 bits :
$$M(j, i) = \text{round}\left( 65535 \times \frac{V(i) - V_{\text{min}}}{V_{\text{max}} - V_{\text{min}}} \right)$$
où $V$ représente le vecteur 1D étendu en 2D par invariance verticale via `np.tile`.

## 🚀 Usage

Le script s'exécute en ligne de commande. Il prend deux arguments positionnels obligatoires (taille de l'image et taille du pixel) ainsi que des options pour configurer la physique du cylindre.

### Exemple de base
Pour générer les cartes de phase et d'épaisseur d'un cylindre de rayon $R = 100\text{ nm}$ avec une anisotropie négative ($K_u = -10^5\text{ J/m}^3$), sur une image de $1000 \times 1000\text{ px}$ (résolution de $1\text{ nm/px}$), exécutez :

```bash
./curling_holo_simulator 1000 1e-9 --R 100e-9 --Ku=-1e5
```
---

## 🛠️ Dépendances

Le projet nécessite l'environnement Python standard pour le calcul scientifique :
* **NumPy** (gestion des tableaux, FFT 1D et duplication matricielle)
* **SciPy** (constantes physiques fondamentales, modules `integrate.solve_bvp` et `integrate.quad`)
* **Pillow (PIL)** (génération des fichiers PNG 16-bits et gestion des métadonnées structurelles `PngInfo`)
* **Matplotlib** (génération des graphiques de contrôle et affichage des résultats)
