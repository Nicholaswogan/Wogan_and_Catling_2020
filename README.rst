Code used in Wogan and Catling (2019) (link here)
==========
This Github directory contains all the Python, Matlab, and Fortan code used in Wogan and Catling (2019), which is an article about the co-evolution of chemical disequilibrium and life on the early Earth. To reproduce our results, please follow these steps. Most of the code was developed using macOS 10.14.

Requirements:

- gcc (i.e. fortran compiler)
- MATLAB 2017 or newer
- Python 3 with packages numpy, matplotlib, os, sys, and subprocesses
- Jupyter notebook

1. Calculate volcanic outgassing fluxes on the Hadean Earth
--------
Use Jupyter Notebook to open the file titled "Volcanism_outgassing_fluxes.ipynb" and run every cell within the notebook. This will generate ranges of volcanic outgassing fluxes on the prebiotic Earth, and save these fluxes in the file "Atmos/Data/volc_photochem_inputs.txt". This portion of the code is described in Wogan and Catling (2019) in Appendix A.

2. Photochemical modeling of the prebiotic and chemotrophic Earth
--------
Open a terminal and navigate to the "Atmos" directory. The Atmos photochemical code, which we use here, was developed by the Virtual Planetary Laboratory (https://github.com/VirtualPlanetaryLaboratory/atmos). Compile the Atmos photochemical code with the commands

.. code-block:: bash

   make -f PhotoMake clean
   make -f PhotoMake
   
In my experience, in order for the code to compile without any issues, you should be using gcc-4.9. To install gcc-4.9 on MacOS use the Homebrew command:

.. code-block:: bash

   brew install gcc@4.9
   
Once the code is compiled, run the Python command.

.. code-block:: bash

   python volcan_iter.py

This will run photochemical models for every outgassing scenario produced in step 1. Two outputs are produced for each outgassing scenario: The first is for the prebiotic Earth, and the second is for the Earth influenced by a chemotrophic ecosystem (for details of the ecosystem see Wogan and Catling (2019)). Note that modeling the chemotrophic Earth requires running multiple photochemical models which eventually converge on a solution. See the function "solve_eco" in the file "volc_func.py" for details of this algorithm.

The code will output results in the directory "Atmos/Output/disequilibrium/Volc_iter". Approximate run time is 30 minutes. 

3. Calculate the chemical disequilibrium of each modeled atmosphere
--------
Open a terminal and navigate to the same directory as this README file. Copy the results from the photochemical modeling to a new directory with the command

.. code-block:: bash

   cp -r Atmos/Output/disequilibrium/Volc_iter Gibbs_minimization/
   
Now navigate to the directory "Gibbs_minimization" and run the command

.. code-block:: bash

   python volc_gibbs.py
   
This calculates the atomsphere-ocean chemical disequilibrium, in terms of avaliable Gibbs energy, of each atmosphere-ocean system produced with the Atmos photochemical model. The code works by reacting all the molecules and atoms in the atmosphere and ocean system to a state of chemical equilibrium. The chemical disequilibrium is then defined by the Gibbs energy difference between the initial/observed and equilibrium state.

For this calculation we use an open-source Gibbs minimzation code that was produced and described by `Krissansen-Totton et al. (2016)
<https://www.liebertpub.com/doi/full/10.1089/ast.2015.1327?casa_token=WCllthEfTOEAAAAA%3AhqcOavOnfVYJItKyfv6Xq-TGNveOG2S9RQNkfj65iGV1EVcwnlTl4wSjK4DTXBi26hVF0AJOAX3t>`_. The source code can be downloaded at http://www.krisstott.com/publications.html .

Note: The atmosphere-ocean gibbs minimization is a optimization problem that has multiple local minima. The code tries to find a global minimum by attempting the minimization from many different random starting points. The number of attempted minimizations can be changed in the file "Gibbs_minimization/Main_script_iterate.m". More iterations will ensure a global minimum is found, although it will also slow down the calculation. I have had good luck with ~30 iterations, although, if the plots produced in the final step (step 4) do not look smooth, then you should re-run volc_gibbs.py with more iterations.

Approximate run-time is ~24 hours.

"Gibbs_minimization/Volc_iter_completed" contains the results from completing steps 1 through 3. 

4. Plot the results
--------
To reproduce Figure 1 and 2 (with the exception of the atmosphere-only disequilibrium calculation) in Wogan and Catling 2019, run the following Python scripts in the root directory.

.. code-block:: bash

   python plot_figure1.py
   python plot_figure2.py
   
If you have questions, please contact me at wogan@uw.edu.
