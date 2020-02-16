This set of Matlab scripts and databases calculates the multiphase equilibrium state for the Earth through time, and returns the Gibbs free energy difference between the initial (observed) and equilibrium states. The science is described in J. Krissansen-Totton, S. Olson and D.C. Catling (2018) "Disequilibrium biosignatures over Earth history and implications for detecting exoplanet life", Science Advances, 4, eaao5747 and J. Krissansen-Totton, D. S. Bergsman and D.C. Catling (2016) "On Detecting Biospheres from Thermodynamic Disequilibrium in Planetary Atmospheres", Astrobiology 16, 39-67.

As a matter of courtesy, we request that people using this code please cite Krissansen-Totton et al. (2018). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author (joshkt@uw.edu).

Examples are provided for the Archean and Proterozoic Earth. These examples reproduce the maximum and minimum disequilibria results reported in Krissansen-Totton et al. (2018). The code can, in principle, be modified to calculate multiphase equilibria for other aqueous systems, but users should be aware that aqueous chemistry approximations make break down if the desired system is qualitatively different to the Earth's atmosphere-ocean system. We cannot guarantee the accuracy of the code for other aqueous systems.

This work was supported by the NASA Astrobiology Institute's Virtual Planetary Laboratory (grant NNA13AA93A) and by the NASA Exobiology Program (grant NNX15AL23G) awarded to D.C.C. J.K.-T. was supported by NASA Headquarters under the NASA Earth and Space Science Fellowship program (grant NNX15AR63H).

REQUIREMENTS: Matlab 2014a or later (earlier versions are untested).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HOW TO RUN CODE:
(1) Put all the Matlab scripts and databases in the same directory, and ensure Matlab is working in this directory.
(2) Open and run the Matlab scripts loadDatabase, loadDatabaseB, loadDatabaseC, loadDatabaseD, and load_Pitzer. These scripts load the Gibbs free energies of formation, fugacity coefficients, and various aqueous species databases.
(3) Open the Matlab script load_vectors and uncomment the input file you wish to use. Provided examples include the Archean maximum disequilibrium, Archean minimum disequilibrium, Proterozoic maximum disequilibrium, and Proterozoic minimum disequilibrium. Run the load_vectors script once an input file has been selected.
(4) Open and run the Matlab script Main_script_iterate. Initial abundance, final abundances, and the Gibbs free energy change are displayed, along with a histogram comparing initial and equilibrium states.
(5) Important note on interpreting outputs: When completed, Main_script_iterate will display an array containing the Gibbs energy difference for each iteration e.g. 

array_of_Gibbs =

   1.0e+06 *

   0.003684320334054
   0.011704380516955
  -0.000153006969957
   0.003195400484238
   1.018585217903734
  -0.000233933835597
  -0.000234039836910
  -0.000233847837638
   1.883360942323555
   0.009415257485166

It is necessary to run multiple iterations because the fmincon optimizer often gets stuck at local minima or does not find a solution at all, which results in fill values of 9999.99 J/mol. This is because the multiphase equilibrium problem is not convex and may have many local minima (see Krissansen-Totton et al. (2016) for further explanation). To get around this problem, the code is automatically run multiple times starting at different initial conditions. There is no rigorous criterion for deciding how many iterations are enough to reliably obtain the global minimum. However, in our experience, if the same (negative) change in Gibbs energy is obtained multiple times starting from different initial conditions, then this is almost certainly the global minimum. All the global minima reported in Krissansen-Totton et al. (2018) were independently validated using the commercial computational chemistry package Aspen Plus. ~10 iterations will usually be enough to obtain the approximate global minimum, whereas ~100 iterations will reliably return a precise global minimum. The number of iterations is specified by the variable ‘num’ in Main_script_iterate. The script Main_script_iterate will output the largest, negative Gibbs energy change as the global minimum (‘deltaG_value’), but this should always be checked against the array of Gibbs energy differences to ensure that it is a reliable global minimum.

In the example above (Archean maximum disequilibrium), the global minimum is 234 J/mol. This is shown in the following output:

Global minimum (J/mol):

deltaG_value =

    -2.340398369102261e+02

Observe that if only 5 iterations had been run, then the global minimum reported would have been 153 J/mol, which is clearly not the true global minimum.

(6) Finally, note that the code may take a while to run. For simple systems with only a few molecular species, the code should be virtually instantaneous. However, for the complex Archean and Proterozoic examples provided, each iteration may take several minutes or even longer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPLANATION OF CODE STRUCTURE:

%% Main_script_iterate
This script loads the input file (see below) and performs the global optimization routine by calling the Gibbs energy minimization script for multiple, randomized initial conditions. The minimum value (largest negative change in Gibbs energy) is selected and displayed, and the plotting script (Plot_outputs) is called. Warnings are displayed if the number of iterations is too small to reliably recover the global minimum. The number of iterations is specified in this script by the variable ‘num’.


%% Gibbs_energy_minimization
This script calculates the multiphase equilibrium state given some initial (observed) out-of-equilibrium state, and calculates the Gibbs free energy difference between the initial and equilibrium states. The following inputs are required (the input section is clearly demarcated in the code):
- Pressure of the system (default P=1.013 bar)
- Temperature of the system (default T=288.15 K)
- Vectors containing true initial state species and abundances. These are automatically loaded from the load_vectors script.
- Random initial state vectors to be used as the arbitrary starting point for the optimization routine. These are automatically provided by Main_script_iterate.

Given these inputs, the "Gibbs_energy_minimization" script will adjust Na(+) abundances of the true initial state to ensure the system has zero charge. Atom and charge conservation conditions will be calculated from the true initial state. Specifically, the "makeAs" function is used to create an array that specifies the abundances of all elements. The total Gibbs free energy of the true initial (observed) state is calculated by calling the function totalGibbsInternal (see below).

Next, the random initial condition along with atom and charge constraints from the true initial state are fed into the fmincon optimization routine. fmincon is an interior points optimization function available in Matlab. This function minimizes the total Gibbs energy of the system, as defined by the function totalGibbsInternal, subject to the constraint that atoms and charge are conserved and that mixing ratios for all species are non-negative. In fact, we bound the mixing ratios between 1e-30 and 1 to avoid singularities from taking the natural logarithm of zero. The fmincon function returns the molar abundances in equilibrium, the Gibbs free energy of the final state, and an exitflag to ensure no errors occurred during minimization.

Finally, given the Gibbs energy of the final state, the difference in Gibbs energies between the initial and final state is calculated and displayed. 

Within "Gibbs_energy_minimization" there are three local functions that are called:
- totalGibbsInternal: This takes the molar abundances vector and calculates the Gibbs energy of the system according to equation 10 in the main text. This function calls various external functions such as fugCoef (calculates gas fugacity coefficients), fugCoefAQ (calculates activity coefficients for aqueous species), and Pitzer_activity_diseq (calculates aqueous activity coefficients using the Pitzer equations). The function also calculates the analytic gradient of Gibbs energy with respect to abundances, which ensures the optimization routine converges more reliably. Both the total Gibbs energy of the system and its derivative are returned.
- onlyPhase: This function returns a form of some original matrix with only species in it of the specified phase.
- mycon: Used for calculating charge balance.


%% Plot_outputs:
This function takes the initial and final abundances as inputs and plots them as histograms. Two figures are created, one that plots all species, and a second that separates gaseous and aqueous species into two subplots. Note that the second figure (with subplots) is set up specifically for the Archean and Proterozoic examples, and if the species list is modified then this part of the plotting script will need to be modified accordingly.


%% format_ticks_v2:
A script for formatting x-ticks in histogram plots, created by Alex Hayes. Called by Plot_ouputs.


%% gibbs:
The gibbs function takes the following inputs:
- list of names of gaseous molecular species (or the name of a single species)
- Temperature
It computes and returns standard Gibbs free energies of formation for gaseous species by calling the NEWNASA database. The NEWNASA database actually contains coefficients for computing Berman-Brown Gibbs free energies, but the "gibbs" function converts these to standard Gibbs free energies of formation.
The "gibbs" function calls the external functions "makeAs" and "convert2standard". The former creates a coefficient matrix whereas the latter converts certain elements to their standard forms for the purposes of Gibbs free energy calculations (see below).

Within "gibbs" there is one local function that is called:
"gibbsHelper": This function checks to make sure all species are available and that the temperature specified is within the appropriate range for the semi-empirical Gibbs free energy calculations. Errors/warnings are returned if this is not the case. The "gibbsHelper" function actually computes Berman-Brown Gibbs free energies by loading the NEWNASA database and calculating enthalpies and entropies. The "gibbsHelper" function also calls the external function "searchData" which searches the NEWNASA database for molecules that match the named inputs (returns errors if there is an ambiguity or if no match is found).


%% gibbsBB:
Identical to the function gibbs except the Berman-Brown Gibbs free energies are retained.
 

%% gibbsAQ:
This function takes the following inputs:
- list of names of aqueous molecular species (or the name of a single species)
- Temperature
- Pressure
It computes and returns Gibbs free energies of formation for aqueous species by calling the SUPCRT database (see loadDatabaseC). See Appendix C in Krissansen-Totton et al. 2016 for a complete description of how Gibbs energies are calculated for aqueous species.


%% makeAs:
The makeAs takes a cell array of species names and creates the coefficient matrix used to define the atom conservation constraint. For example, given the molecular names for Archean, makeAs(names) returns:
ans = 
    'H'    'O'    'N'    'C'    'S'    'Na'    'Cl'
    [2]    [1]    [0]    [0]    [0]    [ 0]    [ 0]
    [0]    [2]    [0]    [0]    [0]    [ 0]    [ 0]
    [0]    [0]    [2]    [0]    [0]    [ 0]    [ 0]
    [2]    [1]    [0]    [0]    [0]    [ 0]    [ 0]
    [0]    [2]    [0]    [1]    [0]    [ 0]    [ 0]
    [3]    [0]    [1]    [0]    [0]    [ 0]    [ 0]
    [4]    [0]    [0]    [1]    [0]    [ 0]    [ 0]
    [0]    [1]    [0]    [1]    [0]    [ 0]    [ 0]
    [2]    [0]    [0]    [0]    [0]    [ 0]    [ 0]
    [2]    [0]    [0]    [0]    [1]    [ 0]    [ 0]
    [1]    [0]    [0]    [0]    [0]    [ 0]    [ 0]
    [0]    [0]    [0]    [0]    [0]    [ 1]    [ 0]
    [0]    [0]    [0]    [0]    [0]    [ 0]    [ 1]
    [1]    [3]    [0]    [1]    [0]    [ 0]    [ 0]
    [0]    [2]    [0]    [1]    [0]    [ 0]    [ 0]
    [0]    [3]    [0]    [1]    [0]    [ 0]    [ 0]
    [1]    [1]    [0]    [0]    [0]    [ 0]    [ 0]
    [3]    [0]    [1]    [0]    [0]    [ 0]    [ 0]
    [4]    [0]    [1]    [0]    [0]    [ 0]    [ 0]
    [0]    [0]    [2]    [0]    [0]    [ 0]    [ 0]
    [0]    [2]    [0]    [0]    [0]    [ 0]    [ 0]
    [4]    [0]    [0]    [1]    [0]    [ 0]    [ 0]
    [0]    [1]    [0]    [1]    [0]    [ 0]    [ 0]
    [2]    [0]    [0]    [0]    [1]    [ 0]    [ 0]
    [0]    [4]    [0]    [0]    [1]    [ 0]    [ 0]
    [2]    [0]    [0]    [0]    [0]    [ 0]    [ 0]
This is the transpose of the 'a' matrix defined in Gibbs_energy_minimization. "makeAs" calls the external functions "breakdown" and "findElement"


%% convert2standard:
The "convert2standard" function ensures that molecular species are converted to their reference state for the purposes of Gibbs free energy calculations (e.g. O present as O2, carbon is present as graphite). It takes an array of molecular species as an input, and returns the same species in their standard form.


%% searchData:
The "searchData" function takes a molecular species as a (string) input. It then searches the NEWNASA database for a match with this molecule and returns the index that corresponds to its position. Errors/warnings are returned if there is an ambiguity or if no matches are found.


%% searchDataB:
The "searchDataB" function takes a molecular species as a (string) input. It then searches the fugacity database (databaseB) for a match with this molecule and returns the index that corresponds to its position. Errors/warnings are returned if there is an ambiguity or if no matches are found.


%% searchDataC:
The "searchDataC” function takes a molecular species as a (string) input. It then searches the sprons96 SUPCRT database (databaseC) for a match with this molecule and returns the index that corresponds to its position. Errors/warnings are returned if there is an ambiguity or if no matches are found.


%% searchDataD:
The "searchDataD” function takes a molecular species as a (string) input. It then searches the Aqueous coefficients database (databaseD) for a match with this molecule and returns the indices that corresponds to its position. This is used by fugCoefAQ to calculate aqueous activity coefficients using the Truesdell-Jones equation. This method has been largely superseded by calculating activity coefficients using the more accurate Pitzer equations (see below). 


%% breakdown:
The "breakdown" function returns a cell array of all of the elements in a given species in the first row and the number of each element in the second row. For example breakdown('H2O') returns the following:
ans = 
    'H'    'O'
    [2]    [1]


%% findElement:
This is a simple function that returns the location of a specified string in a larger array.


%% fugCoef:
This function returns fugacity coefficient values using the Soave equation (see Appendix A in Krissansen-Totton et al. 2016 for a complete description). It requires the following inputs:
- Temperature and pressure of the system
- List of names of molecular species for which fugacity coefficient is to be calculated.
- Array of molar abundances of these species
The function calls databaseB which is created by the function "loadDatabaseB", described below. databaseB contains the critical temperatures, critical pressures and acentric factors required to calculate fugacity coefficients. The "fugCoef" function returns a cell array containing log(fugacity coefficients) for each of the inputted species. Errors will be returned if the designated species cannot be located in the database or if the computer has trouble solving the cubic equation for the compressibility factor.


%% fugCoefAQ:
This function calculates aqueous activity coefficients using the Truesdell-Jones equation (see Krissansen-Totton et al. 2016 for details). This method has been largely superseded by calculating activity coefficients using the more accurate Pitzer equations (see below). 


%% Pitzer_activity_diseq:
This function calculates aqueous activity coefficients and water activity using the Pitzer equations. The methodology for calculating water activity is described in Krissansen-Totton et al. 2016, whereas the methodology for calculating activity coefficients for other aqueous species is described in the 2018 paper accompanying this code (equation 11). In short, this function takes as inputs the abundances, charges, and names of aqueous species and outputs their activity coefficients. Activity coefficients for anions and cations are calculated separately using the database of Pitzer coefficients “pitzer_coef.txt” (see load_Pitzer). The temperature dependence of aqueous activity coefficients is neglected. This function is called by totalGibbsInternal when calculating the total Gibbs energy of a system.


%% dielectric.m:
Calculates the dielectric constant of water, which is used to calculate the Gibbs energies of aqueous species in gibbsAQ.


%% load_vectors:
This script opens the designated input file containing species names, abundances, charges, and phase. The input file is read and each line is converted to an array or cell, as needed.


%% loadDatabase:
The "loadDatabase" function opens the NEWNASA.txt file that contains coefficients for calculating Gibbs free energies of formation. The coefficients for all species are stored as a Matlab array called "database" which is called by the "gibbs" function described above. The "loadDatabase" function must be run before "Main_script_iterate" can be successfully executed.


%% loadDatabaseB:
The "loadDatabaseB" function opens the fugacityCoefficientVariables.txt file that contains thermodynamic variables for calculating gas fugacity coefficients. The variables for all species are stored as a Matlab array called "databaseB" which is called by the "fugCoef" function described above. The "loadDatabaseB" function must be run before "Main_script_iterate" can be successfully executed.


%% loadDatabaseC:
This function opens the sprons96_edited2.dat file that contains coefficients for calculating aqueous Gibbs free energies of formation. The variables for all species are stored as the Matlab array called “databaseC” which is called by the “gibbsAQ” function described above. The "loadDatabaseC” function must be run before "Main_script_iterate" can be successfully executed.

%% loadDatabaseD:
This function opens the AQ_coefficients.txt file that contains coefficients for calculating aqueous activity coefficients using the Truesdell-Jones equation. The variables for all species are stored as the Matlab array called “databaseD” which is called by the “fugCoefAQ” function described above. Note, however, that this approach has been superseded by the Pitzer equation method. The "loadDatabaseD” function must be run before "Main_script_iterate" can be successfully executed.

%% load_Pitzer:
This function opens the pitzer_coef.txt file that contains parameters for calculating aqueous activity coefficients using a simplified version of the Pitzer equations. The variables for all species are stored in the Matlab arrays B0, B1, B2, and C0 (binary interaction parameters), which are called by function Pitzer_activity_diseq. The “load_Pitzer” function must be run before "Main_script_iterate" can be successfully executed.


END EXPLANATION OF CODE STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPLANATION OF INPUT FILES:

Four examples of input files are provided: inputs_Archean_max.txt, inputs_Archean_min.txt, inputs_Proterozoic_max.txt, and inputs_Proterozoic_min.txt. They all have the following structure:

- First line: Vector containing charge of all species. All gaseous species must have zero charge.
- Second line: Vector containing phase identifiers for each species.
    % Identifiers:
    % 0 - solvent (liquid water must be the first species in all vectors)
    % 1 - gases 
    % 2 - liquids (databases for other liquids not included, but this is included for a possible future expansion of the code)
    % 3 - solids (databases for solid species not included, but this is included for a possible future expansion of the code)
    % 4 - aqueous electrolytes
- Third line: List of molecular species. For gaseous species and liquid water, ensure there is one space after each species to unambiguously locate the string in the Gibbs energy database. Aqueous species do not require spacing.
- Fourth line: Vector containing molar abundances for each species. The molar abundances of gaseous species should be mixing ratios to ensure the final Gibbs free energy change has units of joules per mole of atmosphere. That is, the molar abundances of gaseous species should sum to 1. Aqueous species have units of moles per kg of water. Liquid water (first entry) has units of moles per mole of atmosphere. The default value, 436.77, ensures the correct water-to-gas ratio for the modern Earth’s ocean and atmosphere.
- Fifth line: Vector containing the scalars ocean volume (om), salinity scalar (ams), and a scaling factor (scale_factor). The ocean volume can be changed to perform sensitivity tests using different ocean volumes (but note that initial carbon chemistry must be re-equilibrated for a fair comparison). The salinity scalar is no longer used. The scaling factor helps ensure more rapid convergence of optimization algorithms. It should not need to be changed for any of the examples provided.

END EXPLANATION OF INPUT FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
DATABASE EXPLANATION

%% NEWNASA
The NEWNASA database contains the coefficients required to calculated the Gibbs free energy of formation as a function of temperature for a large number of molecular species. Appendix A in Krissansen-Totton et al. (2016) describes how these calculations are performed. The format of the data entries in the NASA text file is described here: http://www.grc.nasa.gov/WWW/CEAWeb/def_formats.htm. The NASA database itself can be accessed here: http://www.grc.nasa.gov/WWW/CEAWeb/ whilst an updated version of the text file is available here: "http://garfield.chem.elte.hu/Burcat/NEWNASA.TXT"


%% fugacityCoefficientsVariables
This text file contains critical temperatures, critical pressures, and acentric factors for a large number of gaseous species. There are used to calculate fugacity coefficients using the Soave equation as described in Appendix A of Krissansen-Totton et al. 2016. These thermodynamic parameters were obtained from Knovelís online database (see Perry, R. H., D. W. Green, and Knovel (Firm) (2008), Perry's chemical engineers' handbook, 8th ed., 1 v. (various pagings) pp., McGraw-Hill, New York.)


%% AQ_coefficients.txt
This text file contains coefficients for calculating aqueous activity coefficients using the Truesdell-Jones equation, as described in Appendix C of Krissansen-Totton et al. (2016). The coefficients were sourced from Langmuir (1997). However, this methodology has been superseded by calculating aqueous activity coefficients using the Pitzer equations, which is a more accurate approach.


%% pitzer_coef.txt
This text file contains the coefficients used for calculating aqueous activity coefficients using the a simplified version of the Pitzer equations, namely the binary interaction parameters B0, B1, B2, and C0. These interaction parameters were obtained from Appelo and Postma (2005) and Marion (2002). See Appendix C in Krissansen-Totton et al. (2016) and the Materials and Methods section in Krissansen-Totton et al. (2018) for a detailed description of how these parameters are used to calculate water activity and aqueous activity coefficients, respectively.


%% sprons96_edited2.dat
This file contains the sprongs96 SUPCRT database for calculating aqueous Gibbs energies of formation. The aqueous species are clearly demarcated under the header ‘aqueous species’. See Appendix C in Krissansen-Totton et al. (2016) for equations describing how Gibbs energies of formation are calculated using the species-specific coefficients provided in this database. 


END EXPLANATION OF DATABASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
COMMON PITFALLS (aka mistakes we made along the way that you should try to avoid)

- Carbonate chemistry equilibrium. Currently, the initial carbon-bearing species (pCO2, CO2(aq), CO3(-2), and HCO3(-)) in the examples provided are in equilibrium. Note that if you modify the input files, you will potentially change the carbonate equilibria. Therefore, any available energy you calculate will include an unwanted contribution from the re-equilibration of carbon-bearing species. To properly calculate the available Gibbs free, you will need to calculate the carbonate equilibria in isolation and use these final abundances as the new initial abundances for the full calculation. Additionally, if you want to change the pH of the system it is not sufficient to merely modify the initial H+ molality - carbonate equilibria and alkalinity must be adjusted accordingly, as described in Krissansen-Totton et al. (2018).

- Similarly, the initial abundances for other dissolved aqueous species in the examples provided are calculated from Henry’s law, and are therefore already in equilibrium with their gaseous counterparts. If the system parameters are changed dramatically, then these initial molalities may need to be recalculated from Henry’s law to avoid an unwanted phase-change contribution to the overall Gibbs energy difference.

- Be careful with aqueous species units. The input vector requires units of mmol/kg, whereas the output vectors have units of moles per mole of atmosphere. Be sure to convert back to modalities if you want to take output abundances and use them as input molalities for a new calculation.

- Pressure changes. Be aware that if you wish to change the pressure of the system, you may need to change more than merely the pressure parameter (P) depending on the system you are interested in modeling. If you wish to model a system with the same atmospheric mass as the Earth but a different pressure (e.g. because of higher gravity), then changing the pressure parameter is sufficient. However, if you wish to model a system like the early Earth with lower (higher) pressure because of smaller (greater) atmospheric mass, then you will also need to change the ocean mass (om) parameter in the input file to scale the ocean-to-atmosphere mass ratio accordingly. For example, to simulate a P = 0.5 bar early Earth, the P parameter must be changed, and om=2 in the input file.

- If you are only interested in calculating gas phase equilibria, do not use this version of the code. Instead, use the code associated with Krissansen-Totton et al. (2016), which is much easier to use and also available on the website of the lead author.

- Spacings in input ‘names’ cell. When modifying species names in input files, make sure to include spaces after gaseous species because otherwise the search scripts may have trouble uniquely identifying them in their respective databases. This is not necessary for aqueous species.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

-------------------------
Contact e-mail: joshkt@uw.edu

