\<immersive type="text/markdown" id="ricker\_model\_readme"\>

# Coupled Ricker Map Simulation with 4D Parameter Sweep

This repository contains a C++ application for running large-scale simulations of a coupled Ricker map model. The simulation performs a comprehensive 4-dimensional parameter sweep to study the system's stability and behavior under various conditions of noise and coupling. The accompanying Jupyter Notebook is used for in-depth analysis and visualization of the generated data.

The core of this project is a robust C++ script designed for high-performance computing, featuring parallel processing with OpenMP, real-time progress tracking, and a resume-from-crash capability to handle long-running simulations.

## Features

  - **High-Performance Simulation**: Written in C++ and parallelized using OpenMP to efficiently explore the vast parameter space.
  - **4D Parameter Sweep**: Systematically investigates the model's behavior across ranges of four key parameters: `epsilon`, `local_sigma`, `global_scale`, and `shock_prob`.
  - **Resume Capability**: The simulation can be stopped and resumed, automatically skipping previously completed parameter sets by reading the existing output CSV file. This is crucial for long-running sweeps.
  - **Live Progress Tracking**: A dedicated thread provides real-time updates on the simulation's progress, showing the percentage of completed parameter sets.
  - **Graceful Exit**: Captures `Ctrl+C` signals to shut down gracefully, ensuring that all progress up to that point is saved.
  - **Data Analysis**: Includes a Jupyter Notebook (`RickerMapsNew.ipynb`) with Python code for loading, analyzing, and creating insightful visualizations from the simulation output.

## The Model

The simulation is based on a coupled Ricker map, a model often used in theoretical ecology to describe the dynamics of two interacting populations. The core equations are:

`X1(t+1) = η((1 - ε) * R(X1(t)) + ε * R(X2(t))) + g`
`X2(t+1) = η((1 - ε) * R(X2(t)) + ε * R(X1(t))) + g`

Where:

  - `X1(t)` and `X2(t)` are the population sizes at time `t`.
  - `R(x) = x * exp(r * (1 - x))` is the Ricker map function.
  - `η` represents local demographic noise, controlled by `local_sigma`.
  - `g` represents global environmental noise, controlled by `global_scale` and `shock_prob`.
  - `ε` (epsilon) is the coupling strength between the two populations.
  - `r` is the intrinsic growth rate.

The simulation calculates the characteristic timescale `τ` (tau), a measure of the system's stability, for each combination of parameters.

## Getting Started

### Prerequisites

  - A C++ compiler that supports OpenMP (e.g., GCC, Clang).
  - For data analysis: Python 3 with Jupyter Notebook and the following libraries: `pandas`, `numpy`, `matplotlib`, `scipy`.

### C++ Simulation Script (`main.cpp`)

1.  **Save the Code**: Save the C++ code into a file named `main.cpp`.
2.  **Compile**: Open a terminal and compile the code with OpenMP support.
    ```bash
    g++ -o ricker_sim -fopenmp main.cpp
    ```
3.  **Run**: Execute the compiled program.
    ```bash
    ./ricker_sim
    ```
    The simulation will start, creating or appending to the output file `RickerSim_...csv` in the same directory.

### Jupyter Notebook for Analysis (`RickerMapsNew.ipynb`)

1.  **Ensure Dependencies**: Make sure you have all the required Python libraries installed.
    ```bash
    pip install pandas numpy matplotlib scipy
    ```
2.  **Run the Notebook**: After the C++ simulation has generated the `.csv` file, open the Jupyter Notebook and run the cells to analyze and visualize the results.

## Output

The C++ script generates a CSV file named `RickerSim_YYYY-MM-DD_HH-MM-SS.csv`. Each row in this file represents the result of the simulation for one set of parameters, averaged over `N_RUNS_PER_POINT` runs.

  - **Columns**: `epsilon`, `local_sigma`, `global_scale`, `shock_prob`, `mean_tau`, `std_tau`.
  - `mean_tau`: The average characteristic timescale `τ` for the given parameter set.
  - `std_tau`: The standard deviation of `τ`.

## How to Contribute

Contributions are welcome\! If you have ideas for improvements or find any issues, please feel free to:

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature/YourFeature`).
3.  Commit your changes (`git commit -m 'Add some feature'`).
4.  Push to the branch (`git push origin feature/YourFeature`).
5.  Open a Pull Request.

