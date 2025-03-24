# CARPP: Core Analysis via Radiative Transfer and Profile Parameters

**CARPP** is a powerful tool for analyzing dense molecular cores using multi-wavelength dust continuum observations. It leverages 3D radiative transfer under spherical symmetry to derive accurate density, temperature, and spectral index profiles, enabling researchers to distinguish between hydrostatic equilibrium (Bonnor-Ebert spheres) and dynamically collapsing cores (power-law profiles).

---

## 🌟 Features

- **3D Radiative Transfer Modeling**: Accurately accounts for temperature gradients and optical depth effects.
- **Multi-Wavelength Data Fitting**: Optimizes density and temperature profiles using dust continuum images at different wavelengths.
- **Flexible Profile Models**:
  - Generalized Plummer-like density profiles (supports Bonnor-Ebert spheres and power-law regimes).
  - Continuous temperature profiles with customizable central and ambient parameters.
- **Robust Performance**: Achieves <10% average parameter error under optimal noise/resolution conditions.


## 📥 Installation and Usage

The CARPP package in this repo is the IDL version of CARPP. To use CARPP, you should: 

1. Download the repo
2. Make sure you have IDL installed on your device
3. Add the path of CARPP to your IDL at `Preferences - IDL - Paths`
4. Do `idl run_carpp.pro` in `Test_Data` folder, it should start to calculate the test core
5. You may now edit and run `run_carpp_template.pro`


## 📜 Citation

TBC
