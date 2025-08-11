<p align="center">
  <img src="logo.png" alt="LontMKM Logo" width="500">
</p>

# LontMKM

## English
LontMKM (Light-Oriented Nano-catalysis Toolkit for MicroKinetic Modeling) is a Python package for simulating microkinetic reactions involving light-driven catalysis.
It supports rate constant calculations via Marcus theory, Eyring theory, and Transient Absorption Spectroscopy (k = 1/τ).

Features:
- YAML-based input
- Energy units in eV
- Reversible and irreversible reactions
- Automatic stiffness detection for ODE solvers (RK45 or BDF)
- Output as PNG graph and CSV file
- Customizable output names via YAML

Installation:
pip install git+https://github.com/mhputra94/LontMKM.git

Usage example:
1. Create a YAML file, for example example.yaml
2. Run: lontmkm example.yaml
3. Outputs: buaya.png and buaya.csv


