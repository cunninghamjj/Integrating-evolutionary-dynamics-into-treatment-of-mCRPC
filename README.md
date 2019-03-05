Software for simulations in manuscript and user's manual.

Created by Jessica Cunningham (jessica.cunningham at moffitt.org) in the department of Imaging at the Moffitt Cancer Research Institute, Tampa, FL.

Publication: Zhang, J. Z., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. Nature Communications , 8, 1816

Table of Contents

Matlab Setup to run simulation
	1.	System requirements
	2.	Downloading simulation
	3.	Matlab path requirements

Files Included
	1. 	MatrixPermutations
	2. 	ESSAnalysis
	3. 	Abiraterone_Main

Running simulation
	1. 	Pre-calculated parameters
	2.	Choosing a patient
	3. 	Choosing Abiraterone treatment type
	4. 	Setting simulation parameters
	5.	Simulation output

Copyright and disclaimer

Acknowledgements


Matlab Setup to run simulation

1. System requirements
To setup and run the simulation you will need a computer running version Matlab 2010a or higher. No special toolboxes are required.

2. Downloading simulation
Download the ZIP of all contents from https://github.com/cunninghamjj/Applying-evolutionary-principles-to-optimize-control-of-mcrpc into a folder of your choice. Once downloaded, unzip the contents and discard the .zip file.

3. Matlab path requirements
Open Matlab and change path the folder with the downloaded files. Then go to File -> Set Path -> Add with Subfolders and select the folder with the downloaded files. 


Files Included

1. MatrixPermutations.m
This file creates all possible orderings of the six competition coefficients, then identifies the 22 ordering that satisfy the inequalities presented in the manuscript. Values from the set of [0.4, 0.5, 0.6, 0.7, 0.8, 0.9] are assigned to these orderings. For "Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer", matrixCoefficients row 7 is used for 'Representative patient #1' and matrixCoefficients row 5 is used for 'Representative patient #2.' The matrix 'matrixCoefficients' is used in both ESSAnalysis.m and Abiraterone_Main.m. 

2. ESSAnalysis.m
This file runs competitive Lotka-Volterra equations for each of the 22 coefficient sets defined in 'matrixCoefficients' in order to calculate the ESS frequencies and population densities. This is done using very fast growth rates and long simulation times in order to have confidence that convergence upon ESS is reached. The matrix 'ESS' is used to set initial conditions in Abiraterone_Main.m and the matrix 'ESSFrequencies' is reported in the supplemental material. 

3. Abiraterone_Main.m
This file applies various Abiraterone treatments to any of the 22 coefficient sets. This will recreate the figures shown in the published manuscript and provide results. 


Running Simulation

1. Pre-calculated parameters
The ESS information and the coefficient sets that are calculated in MatrixPermutations.m and ESSAnalysis.m are hardcoded in lines 5 and 6 in order to avoid running both previous scripts every time a new simulation is performed. Uncommenting line 4 %ESSAnalysis; will run the previous two scripts and provide the same values for 'ESS' and 'matrixCoefficients.' Leaving it commented and using the hardcoded is much faster. 

2. Choosing a patient
Setting matrixIndex to any value between 1-22 will choose which of the 22 coefficient sets will be used. matrixIndex = 7 is used for 'Representative patient #1' and matrixIndex = 5 is used for 'Representative patient #2'. 

3. Choosing an Abiraterone treatment type. 
There are four simulated treatment schedules. By uncommenting one of the four lines, the simulation will apply that treatment schedule. 'None' will give no abiraterone for the entire simulation time. 'SOC' will give maximum tolerated dose. 'Metronomic' will apply a predetermined schedule of application and breaks of abiraterone. 'Adaptive' will use the adaptive therapy protocol of administering abiraterone once the maximum PSA value is reached and removing abiraterone once the PSA value is cut in half. These treatments are created in detail in the 'Design treatment' section of the file. 

4. Setting simulation parameters
Here we set the competition coefficients from the 'matrixCoefficients' matrix, growth rates, PSA decay parameter, initial population densities are set to 40% of the ESS population densities, and maximum simulation time is chosen.

5. Once the simulation is complete, the features discussed in the manuscript are available. The time for the T- population to reach 5% of the total population is available as timeTo5percent. The time for the T- population to reach 90% of the total population is available at timeTo90Percent. Also the total amount of abiraterone given as a percentage of standard of care is available as percentageOfSOC. Two figures are created. One showing the population dynamics and one showing the PSA dynamics. Each are highlighted in tan where the abiraterone was actually administered. 

*** Commenting out the code under the label 'Highlight when Abi is given. - This is slow. Be patient.' in both figures will drastically speed up plotting.

COPYRIGHT AND DISCLAIMER
This software and documentation contained with it are copyright © 2017 by Jessica Cunningham and the Moffitt Cancer Research Institute. All rights are reserved.
This software is supplied 'as is' without any warranty or guarantee of support. The Moffitt Cancer Research Institute is not responsible for its use, misuse, or functionality. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability arising from, out of, or in connection with this software.

ACKNOWLEDGEMENTS
Thanks to both Joel Brown, Robert Gatenby, and Yannick Viossat for their patience and support through the many iterations of this simulation.

		 
