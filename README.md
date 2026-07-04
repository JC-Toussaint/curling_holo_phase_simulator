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