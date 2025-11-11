Hello!
For this automation to work, a series of criteria must be met:

First: poconv, psrmin, psh0 and psh1 must be compiled and working inside the 'arquivos' folder

Second: MSCD-ATA also needs to be compiled and working inside the 'arquivos' folder

Third: The input file must always be named "input_cluster.txt" as it is the reference for creating cluster.i

Fourth: The file with the experimental data must be inside the "arquivos" folder and be in the MSCD format

Fifth: For new.py, which plots the theoretical pattern and calculates the R-factor, just leave the output file name as saida.out.
Or go into new.py and change the name of the file it will read in the last few lines.

Sixth: Due to a limitation of the program that calculates the phase shift, input_cluster.txt can only have one vector basis, but nothing stops you from using this automation just to create the ps and rm files and then take the input file to use in your simulations.

