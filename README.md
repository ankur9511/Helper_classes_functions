# Helper_classes_functions
Contains some helper classes and functions that can be directly imported or patched for use

---

### USMD_plt.py

#### Function: axplt

Simplifies plotting using matplotlib

#### Class: US_1d 
Usage shown in file example_USMD_plt.ipynb for reproducing bulk water solvation energy using indirect umbrella sampling ( https://github.com/seanmarks/indus )
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

---

### pipeline.sh

bash based helper script to do molecular simulation with gromacs


---

### zeolite_connectivity_addH.py

Zeolite are composed of Si (Silicon) and O (Oxygen) atoms

Input: 
- Coordinate file (eg test.gro) of a first guess cleaved zeolite slab exposed to vacuum
  - At present, this can be directly used for orthorhombic crystal systems
  - For rhombic or other systems, preliminary system transformation is necessary
  - In future, the script will be applicable for other systems by inclusion of periodic boundary conditions for other systems

- Direction of surface normal, where the protons will be added to partially-coordinated O atoms
  - 0 : x direction
  - 1 : y direction
  - 2 : z direction

Usage example:

`
python zeolite_connectivity_addH.py test.gro 2
`

Python2 based script to:
- create topology files compatible with TrappeFF (http://siepmann.chem.umn.edu/166) for zeolite
- fine tune cleaved 2D interface with:
  - Cleaving 2D slabs of zeolite can leave unphysical and un-coordinated O atoms
  - This script iteratively checks for coordination number completion of Si, and removes unphysical O atoms
  - Once a good slab is obtained, it populates partially-coordinated surface O-atoms with protons
  - The script can be tweaked to include protonation constant to account for an equilibrium : (Zeo-OH <--> Zeo-O + H)
  - Save the new coordinate file
  - Create topology files: interaction parameters, bonds, angles, dihedrals parameters for the zeolite as dictated by FF
