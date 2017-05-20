# Cancer Modeling Project (CMSC703 @ UMD) 

A repository includes the three projects from the *Network Analysis and Modeling of Biological Systems* class at the *University of Maryland, College Park.* This class was organized and taught by [Prof. Eytan Ruppin](https://sites.umiacs.umd.edu/ruppinlab/).

------

### Datasets (cancer data) 
 - [Achiles](https://software.broadinstitute.org/software/cprg/?q=node/10)
 - ESSbrca: cancer cell model.
 - HER2: cancer cell model.
 - recon1, recon2: the cell models provided in the class

### Source codes (data cleaning)
 - *findESSGenesFromCancerData.m*: a script to extract a list of essential genes from a given cancer data (using threshold score value)
 - *getHER2fromCancerData.m*: a script to extract HER2 cancer type from a given cancer data (assuming the ESSbrca is given)
 - *findESSGenesFromMModel.m*: a script to extract a list of essential genes from a given metabolic model (using FBA)
 - *getLowExprGenes.m*: a script to extract a list of genes that have low expression data by using p-value
 - *findRxnsWithGenesFromMModel*: a script to extract a list of reactions that are related to a given genes but not blocked and essential
 - *defineHumanMediaNCI60.m*: the NCI60 media definition 

### Source codes
 - *finalScript.m*: the main script---i.e., used for the class project
 - *runProjectScript.m*: a script that contains our algorithm implementation
 - *evaluations.py*: a Pypthon script to compute the accuracy of a wile-type model and all modified models

License
-------
This code is licensed under the MIT license. Please see the LICENSE file for more information.
