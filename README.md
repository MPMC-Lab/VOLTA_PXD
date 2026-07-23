# PXD: Multi-Dimensional Electrochemical-Thermal Model Simulation  
*(Supports P2D, P3D, and P4D configurations)*

This project implements a simulation framework for electrochemical-thermal models of battery systems under various dimensional settings (P2D, P3D, and P4D). The code simulates cell voltage, temperature evolution, and internal electrochemical variables using finite difference methods with adaptive time stepping.

## Features

- **Multi-Dimensional Simulation:**  
  Supports P2D, P3D, and P4D configurations for battery modeling.

- **Thermal Modeling:**  
  Incorporates heat generation and conduction effects, including ohmic, entropic, and reaction heat contributions.

- **Intercalation Kinetics:**  
  Implements Butler-Volmer kinetics (adapted for pore-wall flux) for intercalation processes.

- **Adaptive Time Stepping:**  
  Adjusts the time step automatically based on error tolerance to maintain simulation stability.

- **Modular Structure:**  
  Organized into separate modules for parameter initialization, simulation options, state updates, and numerical operators.

## Files

- **ButlerVolmer.m**  
  Calculates pore-wall flux using a Butler-Volmer-like formulation and computes its derivatives.

- **init_Parameters.m**  
  Initializes physical constants, model parameters (e.g., diffusion coefficients, thermal properties), and initial conditions for temperature, concentrations, and potentials.  
  *Modify parameters in this file to adjust the model setup.*

- **init_Simulation.m**  
  Sets up the simulation by creating the computational grid, numerical operators, valid index sets, and initial simulation variables.  
  **Note:** For P4D simulations, you may need to modify this file to account for the tab configuration.

- **main.m**  
  Main driver script that runs simulations for various C-rates using different configurations (P2D, P3D, P4D). It plots cell voltage and temperature profiles.

- **Numerical_operators.m**  
  Provides finite difference operators (harmonic mean, gradient, and diffusion) used for spatial discretization.

- **opt_Simulation.m**  
  Contains simulation options including time step sizes, grid discretization, applied current settings, solid phase diffusion model selection, and adaptive time stepping parameters.  
  *Modify options in this file as needed.*

- **parameter_functions.m**  
  Defines functions to compute model-specific parameters such as diffusivities, conductivities, reaction rates, transference number, open-circuit potential, and their derivatives.

- **sys1_update_parameter.m**  
  Updates system 1 parameters (e.g., diffusivities, reaction rates, potentials) based on the current state of the simulation.

- **sys1_update_variable.m**  
  Updates state variables for system 1 by computing time derivatives, diffusion terms, and charge conservation equations.  
  The updated cell voltage is stored in `mysim.Vcell`, and simulation time in `mysim.time`.  
  Internal electrochemical variables are saved in `mysim.result.sys1`.

- **sys2_update_parameter.m**  
  Updates system 2 parameters (e.g., temperature-related properties) based on the current state variables.

- **sys2_update_variable.m**  
  Updates the temperature field (stored in `mysim.result.sys1T`) and the mean ion flux in system 2 by solving the thermal model and incorporating heat generation terms.

## Requirements

- **MATLAB (R2020 or later recommended)**
- No additional toolboxes are used.

## Usage

1. **Parameter and Option Configuration:**  
   Edit `init_Parameters.m` and `opt_Simulation.m` to adjust model parameters and simulation options.
   - **Caution:** Ensure that `nx` is set proportionally to `Lx` to maintain a uniform grid and stable simulation.
   - For P4D configurations, you may need to modify `init_Simulation.m` based on the tab configuration.

2. **Run Simulation:**  
   Execute the `main.m` script in MATLAB. The simulation will update over time, storing:
   - Time and cell voltage in `mysim.time` and `mysim.Vcell`.
   - Temperature in `mysim.result.sys1T`.
   - Internal electrochemical variables in `mysim.result.sys1`.

3. **Examine Results:**  
   Results are plotted and stored for further analysis.

## Notes

- **Grid Consistency:**  
  It is critical to set the grid (i.e., set `nx` in proportion to `Lx`) correctly to ensure a uniform level of discretization and stable simulation performance.

- **P4D Specifics:**  
  When running P4D simulations, consider the tab configuration requirements and modify `init_Parameters.m`, `opt_Simulation.m`, and `init_Simulation.m` accordingly.

## Authors

- **Sanghyun Kim** – [shnkim@yonsei.ac.kr](mailto:shnkim@yonsei.ac.kr)
- **Jung-Il Choi** – [jic@yonsei.ac.kr](mailto:jic@yonsei.ac.kr)

## Version

- **v2 (July 2026).** This release adds a staggered initialization of the
  thermal/mean-flux subsystem (system 2 advances half a time step ahead of the
  electrochemical subsystem, so every inter-system coupling is evaluated at the
  midpoint of its interval), a delta-form linear solve, and a streamlined
  Butler--Volmer Jacobian assembly. The default P4D in-plane resolution is
  `ny = 11`, matching the baseline grid of the manuscript. The
  numerical-verification results reported in Appendix A of the revised
  manuscript were produced with this version.
- **v1 (March 2025).** Original release accompanying the initial submission.

## License

This project is released under the terms of the MIT License.
