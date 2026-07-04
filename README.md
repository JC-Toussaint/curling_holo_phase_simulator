# Simulation Holographique : Mode de Curling Magnétique

Ce projet implémente la simulation numérique de la distribution d'aimantation dans un cylindre ferromagnétique (mode curling) et calcule son impact sur la phase d'un faisceau d'électrons par holographie électronique. Il génère à terme une carte de phase bidimensionnelle normalisée.

## 📋 Table des Matières
- [Physique et Paramètres](#-physique-et-paramètres)
- [Structure Algorithmique](#-structure-algorithmique)
- [Modélisation Mathématique](#-modélisation-mathématique)
- [Dépendances](#-dépendances)

---

## 🔬 Physique et Paramètres

Le système modélise un cylindre ferromagnétique de rayon $R$. La configuration dépend de trois constantes physiques majeures:
* **$A$** : La constante d'échange.
* **$K_u$** : La constante d'anisotropie uniaxiale positive ou négative.
* **$M_s$** : L'aimantation à saturation.

Un paramètre d'anisotropie adimensionnel $\kappa$ est défini par:
$$\kappa = \frac{K_u R^2}{2A}$$

Le domaine spatial radial est discrétisé selon la variable adimensionnelle $\rho = r/R \in [\epsilon, 1]$, où $\epsilon = r_0/R$ introduit un rayon de cœur minimal pour éviter toute singularité numérique à l'origine.

---

## 💻 Structure Algorithmique

Le script de simulation est articulé autour de 4 étapes clés :
1. **Résolution du Problème aux Limites (BVP) :** Détermination de l'angle d'aimantation $\omega(\rho)$ à l'aide de l'algorithme `scipy.integrate.solve_bvp`.
2. **Projection Magnétique :** Calcul de l'intégrale de projection $I(y)$ le long de la ligne de visée du faisceau d'électrons ($z$). L'évaluation est vectorisée via `np.vectorize` et calculée par quadrature numérique (`scipy.integrate.quad`).
3. **Intégration de la Phase :** Calcul de la phase cumulative $\Phi(y)$ par la méthode des trapèzes. On utilise l'anti-symétie de la phase, imposant $\Phi(0) = 0$ au centre.
4. **Génération 2D :** Duplication du profil 1D (via `np.tile`) pour générer une image 2D invariante selon l'axe longitudinal, convertie et enregistrée au format `uint16` (16 bits non signés).

---

## 📐 Modélisation Mathématique

### 1. Équation d'Euler-Lagrange
La configuration micromagnétique est régie par l'équation différentielle du second ordre suivante:
$$\omega''(\rho) = -\frac{\omega'(\rho)}{\rho} + \left( \frac{1}{2\rho^2} + \kappa \right) \sin\big(2\omega(\rho)\big)$$

Soumise aux conditions aux limites:
$$\omega(\epsilon) = 0 \quad \text{et} \quad \omega'(1) = 0$$

### 2. Intégrale de Projection $I(y)$
Pour chaque abscisse $y \in [0, R]$, l'intégration le long de la ligne de visée géométrique $z_{\text{max}}(y) = \sqrt{R^2 - y^2}$ donne:
$$I(y) = -\mu_0 M_s \int_{-z_{\text{max}}(y)}^{z_{\text{max}}(y)} \cos\big(\omega(\rho)\big) \, \mathrm{d}z$$
où la position radiale locale est $\rho = \frac{\sqrt{y^2 + z^2}}{R}$.

### 3. Calcul de la Phase Cumulative $\Phi(y)$
En intégrant les constantes fondamentales (charge élémentaire $e$ et constante de Planck réduite $\hbar$) et en imposant $\Phi(0) = 0$, la phase finale vaut:
$$\Phi(y) = -\frac{e}{\hbar} \int_{0}^{y} I(u) \, \mathrm{d}u$$

### 4. Image Finale 2D
L'image discrète $M(j, i)$ de taille $N \times N$ est calculée par la normalisation et la quantification suivante :
$$M(j, i) = \text{round}\left( 65535 \times \frac{\Phi(y_i) - \Phi_{\text{min}}}{\Phi_{\text{max}} - \Phi_{\text{min}}} \right)$$

---

## 🛠️ Dépendances

Le projet nécessite l'environnement Python standard pour le calcul scientifique :
* **NumPy** (pour la vectorisation et la duplication de matrices) 
* **SciPy** (modules `integrate.solve_bvp` et `integrate.quad`) 
