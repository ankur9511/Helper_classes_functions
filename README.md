# Helper_classes_functions
Contains some helper classes and functions that can be directly imported or patched for use

--------------------------------------------------------------

### USMD_plt.py

#### Function: axplt

Simplifies plotting using matplotlib

#### Class: US_1d

Requirements:

- Python 2
- PyMBAR https://github.com/choderalab/pymbar

Creates a class for 1D umbrella histogram analysis using pymbar (based on pymbar's example from Prof. Shirt's & Chodera's work):
- Set histogram analysis settings like temperature, order parameter, production region
- Extract information from plumed input and output file
- Use Pymbar for computing PMF
- Class abstraction: User access is limited to final results and histogram

References:

Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states. J. Chem. Phys. 129:124105 (2008). DOI

Chodera JD, Swope WC, Pitera JW, Seok C, and Dill KA. Use of the weighted histogram analysis method for the analysis of simulated and parallel tempering simulations. J. Chem. Theor. Comput. 3(1):26-41 (2007). DOI

----------------------------------------------------------------------------------

### pipeline.sh

bash based helper script to do molecular simulation with gromacs
