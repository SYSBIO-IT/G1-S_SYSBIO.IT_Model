# G1-S_SYSBIO.IT_Model
DESCRIPTION:
The model focuses on the core mechanism of cell size control in steady state growth, no consideration being given yet on how changes in nutrients or growth rate affect cell size. A novel and conceptually simple mathematical model of the G1/S network is proposed, centered on the phosphorylation of Whi5, an inhibitor of G1/S-specific transcription. Although our model includes all known molecular players of the G1/S transition, only few critical parameters (that affect Whi5 phosphorylation) determine its most important functions. Model predictions have been validated by analyzing cell cycle entry and transcriptional activation in mutants expressing a Whi5 protein with constitutive (pseudo)phosphorylation. The model accounts for different data sets that were not considered in model design and parameter optimization, thereby providing a new unifying, comprehensive, molecular mechanism for the core critical cell size control and for the mitosis/mating switch.

FILES:
The Matlab Code to evaluate a single cell time evolution is composed by 7 files.

SingleCell.m
In the script the user can set different Whi5 phosphorylation configurations.
Moreover, we call:
-	the function to generate the matrices for the probability distributions of the different phosphorylation states (generate_matrices.m);
-	the function to start the simulation (d1_one_cell.m);
-	the function to generate the figures (Figures.m).
Besides the time course of all the molecular players involved, as well as of the gene activation curve, the outputs that we provide are:
-	T1, T2 and TG1 duration values
-	Hill coefficient and Median point values for the Hill function that best fits the fraction of the activated genes
-	Protein content at the end of G1 (critical size Ps) and volume at the end of T1 (Vs)
-	Linear growth rate α

d1_one_cell.m
In this function, we define:
-	all the model parameters;
-	the initial conditions;
-	the while loop that updates Euler integration step until the condition to exit G1-phase is fullfilled
-	the control to define the cell phase
The while loop makes the call to der_plugin.m, a function implementing the ODEs

der_plugin.m
In this function, we define:
-	ordinary differential equations of the model that describe the evolution of a cell in the G1/S period.

Figures.m
This function allows to generate the most significant figures presented in the paper. In particular, we plot the evolution of:
-	Swi6, Swi4, free unbound Whi5, Swi6Swi4, Swi6Swi4Whi5
-	Swi6, Mbp1, Swi6Mbp1
-	free Whi5 in the nucleus and Whi5 in the cytoplasm
-	Clb5 plus Clb6, Sic1 in the nucleus, Clb5Sic1 +Clb6Sic1, Clb5Sic1p +Clb6Sic1p, Sic1p, Sic1 in the cytoplasm
-	Clb5, Clb6, Nrm1, Cln1, Cln2
-	Percentage of activated genes and percentage of nuclear Whi5
-	Cln2, free Sic1, Clb5, percentage of activated genes

generate_matrices.m
Function that provides the matrices for the probability distributions of the different phosphorylation states.

Best_Hill.m and Loss_Hill.m
Functions to evaluate the hill coefficient and the median point for the fraction of G1/S regulon activated genes

RUN THE CODE
To execute the Matlab code, the SingleCell.m script must be run.
>>SingleCell

