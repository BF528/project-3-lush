# Project Description

This respistory contains the code for project 3 Concordance of microarray and RNA-Seq differential gene expression. The code is a combination of Linux and R languages. The premise of the project was reproduce the results from Figure 2a and 3b+c and compare the pathway enrichment results reported in Wang et al. 

# Contributors

Athena Mustakis, Evan Holmes, Navin Ramanan 

# Repository Contents
Files: 

analyst_run_limma.R

-contains the functions necessary for parts 5 and 6 of the project (data analyst role)

analyst.RMD

-- contains markdown chunks used to run the functions in analyst_run_limma.R and create tables/plots (data analyst role)

programmer.R

--contains code necessary for programmer role

multiqc.qsub 

-contains the code necessary to run multiqc, run on the command line in terminal using qsub 

runSTAR.qsub 

-contains the code necessary to run the STAR alignment, run on the command line in terminal using qsub and specify the two samples you will be using 


