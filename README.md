# PDD-shiny
A shiny app to analyse the report.tsv, extract of results of DIA-NN software (https://github.com/vdemichev/DiaNN)
## Requirements
R4.1 or higher
RStudio 2021.09.0 Build 351 or higher
## Required data
1. report.tsv from DIA-NN software
The format of the name of run must not contain "." or "_", but for the same biological replicat, you must add "_increment"
example: you have 3 biologicals replicats for sampleA et sampleB, the name of run must be : sampleA_1,sampleA_2,sampleA_3 and sampleB_1,sampleB_2, sampleB_3
