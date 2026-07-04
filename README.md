# Simulation Holographique : Mode de Curling Magnétique

[cite_start]Ce projet implémente la simulation numérique de la distribution d'aimantation dans un cylindre ferromagnétique (mode de curling) et calcule son impact sur la phase d'un faisceau d'électrons par holographie électronique[cite: 1]. [cite_start]Il génère à terme une carte de phase bidimensionnelle normalisée[cite: 13].

## 📋 Table des Matières
- [Physique et Paramètres](#-physique-et-paramètres)
- [Structure Algorithmique](#-structure-algorithmique)
- [Modélisation Mathématique](#-modélisation-mathématique)
- [Dépendances](#-dépendances)

---

## 🔬 Physique et Paramètres

[cite_start]Le système modélise un cylindre ferromagnétique de rayon $R$[cite: 1]. [cite_start]La configuration dépend de trois constantes physiques majeures[cite: 2]:
* [cite_start]**$A$** : La constante d'échange[cite: 2].
* [cite_start]**$K_u$** : La constante d'anisotropie effective[cite: 2].
* **$M_s$** : L'aimantation à saturation[cite: 2].

[cite_start]Un paramètre d'anisotropie adimensionnel $\kappa$ est défini par[cite: 3]:
$$\kappa = \frac{K_u R^2}{2A}$$

[cite_start]Le domaine spatial radial est discrétisé selon la variable adimensionnelle $\rho = r/R \in [\epsilon, 1]$, où $\epsilon = r_0/R$ introduit un rayon de cœur minimal pour éviter toute singularité numérique à l'origine[cite: 3].

---

## 💻 Structure Algorithmique

Le script de simulation est articulé autour de 4 étapes clés :
1. **Résolution du Problème aux Limites (BVP) :** Détermination de l'angle d'aimantation $\omega(\rho)$ à l'aide de l'algorithme `scipy.integrate.solve_bvp`[cite: 5].
2. **Projection Magnétique :** Calcul de l'intégrale de projection $I(y)$ le long de la ligne de visée du faisceau d'électrons ($z$)[cite: 6, 7]. L'évaluation est vectorisée via `np.vectorize` et calculée par quadrature numérique (`scipy.integrate.quad`)[cite: 9].
3. **Intégration de la Phase :** Calcul de la phase cumulative $\Phi(y)$ par la méthode des trapèzes[cite: 11]. Une correction est appliquée pour garantir une phase nulle au centre ($\Phi(0) = 0$)[cite: 11].
4. **Génération 2D :** Duplication du profil 1D (via `np.tile`) pour générer une image 2D invariante selon l'axe longitudinal, convertie et enregistrée au format `uint16` (16 bits non signés)[cite: 13].

---

## 📐 Modélisation Mathématique

### 1. Équation d'Euler-Lagrange
La configuration micromagnétique est régie par l'équation différentielle du second ordre suivante[cite: 5]:
$$\omega''(\rho) = -\frac{\omega'(\rho)}{\rho} + \left( \frac{1}{2\rho^2} + \kappa \right) \sin\big(2\omega(\rho)\big)$$

Soumise aux conditions aux limites[cite: 5]:
$$\omega(\epsilon) = 0 \quad \text{et} \quad \omega'(1) = 0$$

### 2. Intégrale de Projection $I(y)$
Pour chaque coordonnée d'impact transverse $y \in [-R, R]$, l'intégration le long de la ligne de visée géométrique $z_{\text{max}}(y) = \sqrt{R^2 - y^2}$ donne[cite: 7, 8, 9]:
$$I(y) = -\mu_0 M_s \int_{-z_{\text{max}}(y)}^{z_{\text{max}}(y)} \cos\big(\omega(\rho)\big) \, \mathrm{d}z$$
où la position radiale locale est $\rho = \frac{\sqrt{y^2 + z^2}}{R}$[cite: 8].

### 3. Calcul de la Phase Cumulative $\Phi(y)$
En intégrant les constantes fondamentales (charge élémentaire $e$ et constante de Planck réduite $\hbar$) et en imposant $\Phi(0) = 0$, la phase finale vaut[cite: 11, 12]:
$$\Phi(y) = -\frac{e}{\hbar} \int_{0}^{y} I(u) \, \mathrm{d}u$$

### 4. Image Finale 2D
L'image discrète $M(j, i)$ de taille $N \times N$ est calculée par la normalisation et la quantification suivante[cite: 13]:
$$M(j, i) = \text{round}\left( 65535 \times \frac{\Phi(y_i) - \Phi_{\text{min}}}{\Phi_{\text{max}} - \Phi_{\text{min}}} \right)$$

---

## 🛠️ Dépendances

Le projet nécessite l'environnement Python standard pour le calcul scientifique :
* **NumPy** (pour la vectorisation et la duplication de matrices) [cite: 9, 13]
* **SciPy** (modules `integrate.solve_bvp` et `integrate.quad`) [cite: 5, 9]
