# TMV Mechanical Analysis Repository

This repository contains the source code and data necessary to reproduce the computational results presented in the paper:

**"The tubular cavity of tobacco mosaic virus shields mechanical stress and regulates disassembly"**
*A. Díez-Martínez et al.*

The work investigates the biomechanical properties of Tobacco Mosaic Virus (TMV) using a combination of atomic force microscopy (AFM) experiments, coarse-grained simulations, and finite elements modeling. The repository provides the Python scripts used for data analysis, simulation setup, and visualization of the results.

---

## Repository Structure

- **analysis/**
  Contains scripts for processing and analyzing indentation and stress data (e.g., `indentationAnalysis.py`, `stressComputation.py`, `processIndentation.py`, `tmvAnalysis.py`, and `visualize.py`).

- **data/**
  Contains simulation input files, JSON configuration files for different TMV sizes, and Protein Data Bank (PDB) files representing TMV structures.

- **experiments/**
  Provides scripts to set up and run various AFM indentation simulations, including force–velocity experiments and stress tests (e.g., `afm_tmv_eps.py`, `afm_tmv_size.py`, `afm_tmv_stress.py`, and `afm_tmv_vel.py`).

- **README.md**
  This file.

---

## Requirements

- Python 3.x
- Required Python libraries: NumPy, Matplotlib, MDAnalysis, PyMOL (with its Python API), and VLMP (for high-throughput AFM simulation).
- UAMMD library for GPU-based coarse-grained simulations.

Please refer to the respective documentation of each dependency for installation instructions.

---

## Usage

1. **Simulation Setup:**
   Use the scripts in the `experiments/` directory to configure and run the AFM indentation simulations. Adjust the parameters (e.g., TMV size, interaction energy, tip velocity) in the scripts as needed.

2. **Data Analysis:**
   The `analysis/` scripts process raw simulation outputs, generate force-indentation curves, compute stress distributions, and produce visualizations. For example, run:
   ```bash
   python analysis/indentationAnalysis.py <path_to_afm_data_file>
   ```
   to analyze indentation data.

3. **Visualization:**
   The `visualize.py` script leverages PyMOL to generate images from simulation output files. This script can be customized to create snapshots or movies of TMV deformation.

---

## Reproducing the Results

To replicate the computational results:

1. Generate simulation input files using the provided JSON configurations in the `data/jsons/` folder.
2. Run the AFM simulation scripts located in the `experiments/` folder.
3. Analyze the generated output using the tools in the `analysis/` folder.
4. Visualize the simulation outcomes with the `visualize.py` script.

Detailed instructions are provided as comments within each script.

---

## Contact

For questions or further information, please contact [p.ibanez.fre@gmail.com](mailto:p.ibanez.fre@gmail.com).

---
