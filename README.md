# Braking Simulation & Visualization Toolkit:

This repository contains a small simulation and visualization toolkit for exploring vehicle braking performance under varying physical and control parameters.
It allows users to sweep over two independent variables, visualize stopping distance as a 3D surface, and then reduce that surface to meaningful 2D curves using statistical aggregation.
The project is designed for exploratory analysis, intuition-building, and parameter sensitivity studies rather than production-grade vehicle modeling.

# Features:

- Physics-based braking simulation (via simulate_stopping)
- Automated 2-parameter sweeps with validation
- 3D visualization of stopping distance surfaces
- 2D reduction using:
- Mean stopping distance
- Minimum stopping distance
- Standard deviation
- Interactive example runner
- Perceptually uniform colormaps for accurate visual interpretation

# Repository Structure:


├── simulate_stopping.py   # Core braking physics simulation

├── plots.py               # 2D and 3D plotting utilities

├── get_data.py            # Parameter sweeping, data generation, and reduction

└── examples.py            # Interactive demo script with predefined scenarios


# File Overview - simulate_stopping.py:

-   Implements the core braking simulation.
-   Given vehicle parameters (mass, wheelbase, friction coefficients, torque distribution, etc.), this module computes the total stopping distance from an initial velocity.
-   This file contains no visualization code and is intended to be called programmatically.

# File Overview - plots.py:

Provides reusable plotting utilities:
- graph3d() — renders a 3D scatter plot of stopping distance data
- graph2d() — renders 2D line plots for reduced datasets

Includes:
- Nonlinear sigmoid-based color scaling
- Validation of perceptually uniform Matplotlib colormaps

# File Overview - get_data.py

Acts as the analysis engine of the project.
Responsibilities:
- Validates that exactly two inputs are parameter sweeps (range)
- Runs full simulation grids using simulate_stopping
- Returns or visualizes 3D stopping-distance surfaces
- Handles axis labeling and parameter ordering automatically
- This module is intended to be imported by both scripts and notebooks.
- ## Reduces 3D data to 2D using:
- Minimum
- Mean
- Standard deviation

# File Overview - examples.py

An interactive command-line demo showcasing common analyses.
Users can choose from predefined scenarios such as:
- Front brake force distribution vs. total braking torque
- Effects of aggregation method (mean vs. minimum)
- Influence of wheelbase on braking performance
## Each example:
- Runs a 3D parameter sweep
- Displays the 3D surface
- Collapses the surface into a 2D curve

This file is meant for exploration and intuition-building, not reuse.

# Installation & Requirements

## Required Python packages:
- numpy
- matplotlib

# Usage

Run the interactive examples:

>> python examples.py

Follow the terminal prompt to select an example (1–4).
Plots are displayed immediately and are not saved automatically.

# Design Notes

- Exactly two parameters must be specified as range objects for 3D sweeps.
- All other parameters are treated as constants.
- Standard deviation uses the sample (n − 1) formulation.
- Returning 0 for undefined statistics (e.g., SD with n < 2) is a deliberate design choice.

# Intended Audience

This project is suited for:
- Engineering students
- Controls or vehicle dynamics exploration
- Parameter sensitivity studies
- Visualization-driven intuition development

It is not intended to replace high-fidelity vehicle dynamics solvers.

# License

This project is provided as-is for educational and exploratory use.
No warranty or guarantees are provided.
