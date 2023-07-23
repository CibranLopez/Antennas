# Antenna Array Simulation and Radiation Pattern Synthesis

Antenna arrays are an essential component of modern wireless communication systems, and their radiation patterns play a crucial role in determining their performance. This repository provides a set of Python scripts and programs that allow researchers and engineers to simulate antenna arrays, analyze their radiation patterns, and optimize their design for specific applications.

This repository contains a collection of programs for simulating antenna arrays and studying their radiation patterns. These programs have been used in the research paper titled "An improved pattern synthesis iterative method in planar arrays for obtaining efficient footprints with arbitrary boundaries," co-authored by Cibrán López Álvarez (2022).

## Contents

The repository is organized into the following directories:

1. `src`: This directory contains the source code for the antenna array simulation and radiation pattern analysis programs. All programs can be executed directly from the Jupyter Notebook `Antennas.ipynb`.

2. `Masks`: This directory contains the definition of arbitrary contours that can be used as input for the antenna simulations. Users can define their custom shapes in this folder.

3. `Results`: The results of the simulations and radiation pattern studies are stored in this directory. The program automatically saves the resulting radiation pattern and the necessary shape of the antenna with the intensity of each radiating element in this folder.

## Installation

To download the repository and install the dependencies:

```bash
git clone https://github.com/CibranLopez/Antennas.git
cd Antennas
pip3 install -r requirements.txt
```

## Execution

1. Open the Jupyter Notebook `Antennas.ipynb` to access all the programs for antenna array simulation and radiation pattern analysis.

2. In the notebook, you can specify the contour to be used for the antenna simulation by providing the path to the contour file in the `Masks` folder.

3. The program will automatically optimize the antenna design by removing low-excited elements.

4. After the simulation is complete, the resulting radiation pattern and the shape of the antenna with the intensity of each radiating element will be saved in the `Results` folder.

## Paper Reference

If you use the programs from this repository in your research or work, we kindly request that you cite the following paper:

- Title: An improved pattern synthesis iterative method in planar arrays for obtaining efficient footprints with arbitrary boundaries
- Authors: A.A. Salas-Sánchez, C. López-Álvarez, J.A. Rodríguez-González, M.E. López-Martín, F.J. Ares-Pena
- Journal: Sensors
- Publication Date: March 2022
- DOI: [https://www.mdpi.com/1424-8220/21/7/2358](https://www.mdpi.com/1424-8220/21/7/2358)

## Authors

This project is being developed by:

 - Cibrán López Álvarez

## Contact, questions and contributing

If you have questions, please don't hesitate to reach out at: cibran.lopez@upc.edu
