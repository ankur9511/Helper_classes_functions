# Helper_classes_functions
Contains some helper classes and functions that can be directly imported or patched for use

### USMD_plt.py

#### Function: axplt

Simplifies plotting using matplotlib

#### Class: US_1d

Requirements:

- Python 2
- PyMBAR

Creates a class for 1D umbrella histogram analysis using pymbar (based on pymbar's example from Prof. Shirt's & Chodera's work):
- Set histogram analysis settings like temperature, order parameter, production region
- Extract information from plumed input and output file
- Use Pymbar for computing PMF
- Class abstraction: User access is limited to final results and histogram



### pipeline

bash based helper script to do molecular simulation with gromacs
