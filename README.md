# Curling Holo Phase Simulator: Electron Holographic Phase Simulation

This project implements the numerical simulation of the magnetization distribution in a ferromagnetic cylinder (curling mode). Starting from this micromagnetic configuration, the software suite computes the impact of the magnetic field on the phase of an electron beam via **electron holography**.

All results are exported as 2D images in 16-bit PNG format, including detailed physical metadata, and as 1D profiles in CSV format (see [Output files](#output-files)).

## 📋 Table of Contents
- [Physics and Parameters](#-physics-and-parameters)
- [Algorithmic Structure](#-algorithmic-structure)
- [Mathematical Modeling](#-mathematical-modeling)
- [Dependencies](#-dependencies)

---

## 🔬 Physics and Parameters

The system models a ferromagnetic cylinder of radius $R$. The configuration depends on three major physical constants:
* **$A$**: The exchange constant.
* **$K_u$**: The uniaxial anisotropy constant, positive or negative.
* **$M_s$**: The saturation magnetization.

A dimensionless anisotropy parameter $\kappa$ is defined by:
$$\kappa = \frac{K_u R^2}{2A}$$

The radial spatial domain is discretized using the dimensionless variable $\rho = r/R \in [\epsilon, 1]$, where $\epsilon = r_0/R$ introduces a minimal core radius to avoid any numerical singularity at the origin.

**Axis convention.** The wire axis is $Oz$, and the magnetization reads $\mathbf{m} = \sin\omega\,\mathbf{e}_\varphi + \cos\omega\,\mathbf{e}_z$, so that $m_z = \cos\omega$ is the axial component. The electron (or X-ray) beam propagates along $s$, perpendicular to the wire; $y$ is the transverse coordinate in the image, and $(y, s, z)$ is a right-handed frame. A point of the cross-section is $(y, s)$, at radius $r = \sqrt{y^2 + s^2}$, and $\mathbf{e}_\varphi = \mathbf{e}_z \times \mathbf{e}_r$. In the 2D maps, the wire axis is vertical (invariance direction).

---

## 💻 Algorithmic Structure

The simulation script is organized around 4 key steps:
1. **Solving the Boundary Value Problem (BVP):** Determining the magnetization angle $\omega(\rho)$ using the `scipy.integrate.solve_bvp` algorithm. For $K_u < 0$, the initial guess $\omega_0(\rho) = \frac{\pi}{2}\tanh\big(\sqrt{2|\kappa|}\,\rho\big)$ varies on the physical length scale $\ell_K/R = 1/\sqrt{2|\kappa|}$, and the grid is refined quadratically near the core, so that the solver converges even for $|\kappa| \gg 1$. The solution is then checked: `solve_bvp` must have converged and $\omega(\rho)$ must be monotonically increasing, otherwise the simulation stops (see [Failure handling](#failure-handling)).
2. **Magnetic Projection:** Computing the projection integral $I(y)$ along the electron beam's line of sight ($s$). The evaluation is vectorized via `np.vectorize` and computed by numerical quadrature (`scipy.integrate.quad`).
3. **Phase Integration:** Computing the cumulative phase $\Phi(y)$ using the trapezoidal method over the interval $[-2R, 2R]$. The anti-symmetry of the phase is used, enforcing $\Phi(0) = 0$ at the center.
4. **2D Image Generation & Metadata:** Duplicating the 1D profiles (phase and thickness) via `np.tile` to generate 2D matrices exploiting strict invariance along the perpendicular direction. Both maps are exported in 16-bit PNG format (`uint16`), embedding the physical metadata (`PngInfo`) needed for their future decoding.

---

## 📐 Mathematical Modeling

### 1. Euler-Lagrange Equation
The micromagnetic configuration is governed by the following second-order differential equation:
$$\omega''(\rho) = -\frac{\omega'(\rho)}{\rho} + \left( \frac{1}{2\rho^2} + \kappa \right) \sin\big(2\omega(\rho)\big)$$

Subject to the boundary conditions:
$$\omega(\epsilon) = 0 \quad \text{and} \quad \omega'(1) = 0$$

### 2. Projection Integral $I(y)$
For each abscissa $y \in [0, R]$, integrating along the geometric line of sight $s \in [-s_{\text{max}}(y), s_{\text{max}}(y)]$, with $s_{\text{max}}(y) = \sqrt{R^2 - y^2}$, gives:
$$I(y) = -\mu_0 M_s \int_{-s_{\text{max}}(y)}^{s_{\text{max}}(y)} \cos\big(\omega(\rho)\big) \, \mathrm{d}s$$
where the local radial position is $\rho = \frac{\sqrt{y^2 + s^2}}{R}$. $I(y)$ is the integral of the axial component $m_z$ along the beam; it is only used to compute the phase.

### 3. Computing the Holographic Phase $\Phi(y)$
By incorporating the fundamental constants from `scipy.constants` (elementary charge $e$ and reduced Planck constant $\hbar$) and enforcing $\Phi(0) = 0$, the final phase is:
$$\Phi(y) = -\frac{e}{\hbar} \int_{0}^{y} I(u) \, \mathrm{d}u$$

### 4. Final 2D Images
The discrete images $M(j, i)$ of size $N \times N$ (for the phase and the thickness) undergo strict linear normalization followed by 16-bit quantization:
$$M(j, i) = \text{round}\left( 65535 \times \frac{V(i) - V_{\text{min}}}{V_{\text{max}} - V_{\text{min}}} \right)$$
where $V$ represents the 1D vector extended to 2D by vertical invariance via `np.tile`.

## 🚀 Usage

The script runs from the command line. It takes two required positional arguments (image size and pixel size) as well as options to configure the physics of the cylinder.

### Basic example
To generate the phase and thickness maps of a cylinder with radius $R = 100\text{ nm}$ and negative anisotropy ($K_u = -10^5\text{ J/m}^3$), on a $1000 \times 1000\text{ px}$ image (resolution of $1\text{ nm/px}$), run:

```bash
./curling_holo_simulator 1000 1e-9 --R 100e-9 --Ku=-1e5
```

### Output files
The files are written next to the script:
* `analytics_MBEAM_INTEGRAL.png`, `analytics_HOLO_PHASE.png`, `analytics_PATH_LENGTH.png`, `analytics_STXM_XMCD.png`: 16-bit maps of $M_{\text{beam}}(y)$, $\Phi(y)$, the thickness and the STXM contrast, with their metadata;
* `profile_MBEAM_INTEGRAL.csv`, `profile_HOLO_PHASE.csv`, `profile_PATH_LENGTH.csv`, `profile_STXM_XMCD.csv`: the corresponding profiles along Ox (comma-separated, one header line, $x$ in nm).

$M_{\text{beam}}(y)$ is the integral, along the beam, of the component of $\mu_0\mathbf{M}$ along the beam (in T·m), $M_{\text{beam}}(y) = \mu_0 M_s\,J(y)$, with $J(y)$ defined below. It is odd in $y$.

The STXM (XMCD) contrast uses the absorption coefficients $\lambda_+$ (`--lbp`, default $0.01\text{ nm}^{-1}$) and $\lambda_-$ (`--lbm`, default $0.018\text{ nm}^{-1}$). For the helicity $\sigma = \pm 1$, the absorption coefficient depends on the component of the magnetization along the beam, $\mathbf{m}\cdot\mathbf{e}_s = \sin\omega\,\frac{y}{r}$ (only the azimuthal part contributes, since the beam is perpendicular to the wire): $\frac{1}{2}\left(\lambda_+ + \lambda_- + \sigma(\lambda_+ - \lambda_-)\,\mathbf{m}\cdot\mathbf{e}_s\right)$. With the path length $t(y)$ and
$$J(y) = \int_{-s_{\text{max}}(y)}^{s_{\text{max}}(y)} \sin\big(\omega(\rho)\big)\,\frac{y}{\sqrt{y^2 + s^2}}\,\mathrm{d}s,$$
this gives:
$$a_\sigma(y) = \frac{1}{2}\left[(\lambda_+ + \lambda_-)\,t(y) + \sigma(\lambda_+ - \lambda_-)\,J(y)\right], \qquad C(y) = \frac{e^{-a_+(y)} - e^{-a_-(y)}}{e^{-a_+(y)} + e^{-a_-(y)}}$$
$J$, hence the contrast, is odd in $y$: it vanishes on the wire axis and has opposite signs on either side of it. Its global sign depends on the vortex chirality (sign of $\omega$) and on the beam direction ($+\mathbf{e}_s$ here).

The CSV profiles can be plotted directly with gnuplot:

```gnuplot
set datafile separator comma
set key autotitle columnhead      # the first line becomes the curve title
plot 'profile_HOLO_PHASE.csv' using 1:2 with lines
```

### Failure handling
There is no fallback solver. If the BVP solution is not physically admissible, the script prints an error on `stderr`, recalls the value of $\kappa$, and exits with status `1` without producing any image. This happens in two cases:
* `solve_bvp` did not converge, for example:
  ```
  [Euler-Lagrange / BVP] ERROR: solve_bvp did not converge (status 1: The maximum number of mesh nodes is exceeded.)
  kappa = -533.3: no physically admissible curling profile was found.
  ```
* $\omega(\rho)$ is not monotonically increasing; the message then gives the minimum of $\mathrm{d}\omega/\mathrm{d}\rho$.

On success, the script exits with status `0`. The exit status can therefore be checked in batch runs. On a sweep of 224 parameter sets ($|\kappa|$ from $2.5\times10^{-3}$ to $3.2\times10^{4}$), no failure occurred. See [`dev-notes/convergence_bvp_curling.tex`](dev-notes/convergence_bvp_curling.tex) for details.

---

## 🛠️ Dependencies

The project requires the standard Python environment for scientific computing:
* **NumPy** (array handling, 1D FFT, and matrix duplication)
* **SciPy** (fundamental physical constants, `integrate.solve_bvp` and `integrate.quad` modules)
* **Pillow (PIL)** (generation of 16-bit PNG files and management of structural metadata `PngInfo`)
* **Matplotlib** (generation of control plots and display of results)
