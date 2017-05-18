# CancerModeling - CMSC703 @ UMD

## Functions related to cancer data
- findESSGenesFromCancerData.m: get a list of essential genes from a given cancer data (using threshold score value)
- getHER2fromCancerData.m: extract HER2 cancer type from a given cancer data (assuming ESSbrca is given)

## Functions related to metabolic model
- findESSGenesFromMModel.m: get a list of essential genes from a given metabolic model (using FBA)

## Functions to choose reactions for removal
- getLowExprGenes.m: get a list of genes that have low expression data by using p-value
- findRxnsWithGenesFromMModel: get a list of reactions that are related to a given genes but not blocked and essential

## Other functions
- defineHumanMediaNCI60.m: the media definition function used for our project work

## Scripts
- finalScript.m: a script for our project work (you only need to execute this file)
- runProjectScript.m: a function/script that contains the procedures of our algorithm (it is called by finalScript.m)
