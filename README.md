# Helper_classes_functions
Contains some helper classes and functions that can be directly imported or patched for use

#### USMD_plt.py
Simplifies plotting using matplotlib \n
Creates a class for 1D umbrella histogram analysis using pymbar (based on pymbar's example from Prof. Shirt's & Chodera's group):
- Set histogram analysis settings like temperature, order parameter, production region
- Extract information from plumed input and output file
- Use Pymbar for computing PMF
- Class abstraction: User access is limited to final results and histogram
