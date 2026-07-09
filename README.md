# Curling Holo Phase Simulator: Electron Holographic Phase Simulation

This project implements the numerical simulation of the magnetization distribution in a ferromagnetic cylinder (curling mode). Starting from this micromagnetic configuration, the software suite computes the impact of the magnetic field on the phase of an electron beam via **electron holography**.

All results are exported as 2D images in 16-bit PNG format, including detailed physical metadata.

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

---

## 💻 Algorithmic Structure

The simulation script is organized around 5 key steps:
1. **Solving the Boundary Value Problem (BVP):** Determining the magnetization angle $\omega(\rho)$ using the `scipy.integrate.solve_bvp` algorithm.
2. **Magnetic Projection:** Computing the projection integral $I(y)$ along the electron beam's line of sight ($z$). The evaluation is vectorized via `np.vectorize` and computed by numerical quadrature (`scipy.integrate.quad`).
3. **Phase Integration:** Computing the cumulative phase $\Phi(y)$ using the trapezoidal method over the interval $[-2R, 2R]$. The anti-symmetry of the phase is used, enforcing $\Phi(0) = 0$ at the center.
5. **2D Image Generation & Metadata:** Duplicating the 1D profiles (phase and thickness) via `np.tile` to generate 2D matrices exploiting strict invariance along the perpendicular direction. Both maps are exported in 16-bit PNG format (`uint16`), embedding the physical metadata (`PngInfo`) needed for their future decoding.

---

## 📐 Mathematical Modeling

### 1. Euler-Lagrange Equation
The micromagnetic configuration is governed by the following second-order differential equation:
$$\omega''(\rho) = -\frac{\omega'(\rho)}{\rho} + \left( \frac{1}{2\rho^2} + \kappa \right) \sin\big(2\omega(\rho)\big)$$

Subject to the boundary conditions:
$$\omega(\epsilon) = 0 \quad \text{and} \quad \omega'(1) = 0$$

### 2. Projection Integral $I(y)$
For each abscissa $y \in [0, R]$, integrating along the geometric line of sight $z_{\text{max}}(y) = \sqrt{R^2 - y^2}$ gives:
$$I(y) = -\mu_0 M_s \int_{-z_{\text{max}}(y)}^{z_{\text{max}}(y)} \cos\big(\omega(\rho)\big) \, \mathrm{d}z$$
where the local radial position is $\rho = \frac{\sqrt{y^2 + z^2}}{R}$.

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
---

## 🛠️ Dependencies

The project requires the standard Python environment for scientific computing:
* **NumPy** (array handling, 1D FFT, and matrix duplication)
* **SciPy** (fundamental physical constants, `integrate.solve_bvp` and `integrate.quad` modules)
* **Pillow (PIL)** (generation of 16-bit PNG files and management of structural metadata `PngInfo`)
* **Matplotlib** (generation of control plots and display of results)
