# PDD-shiny
A shiny app to analyse the report.tsv, extract of results of DIA-NN 1.8 version software (https://github.com/vdemichev/DiaNN)
## Requirements
R4.1 or higher
RStudio 2021.09.0 Build 351 or higher
## Required data
1. report.tsv from DIA-NN software
The format of the name of run must not contain "." or "_", but for the same biological replicat, you must add "_increment"
example: you have 3 biologicals replicats for sampleA et sampleB, the name of run must be : sampleA_1,sampleA_2,sampleA_3 and sampleB_1,sampleB_2, sampleB_3
## What types of analyses can I do in PDD-shiny?

* __Data__

Create a table of proteins for each run, with differents informations extract in report.tsv : number of peptides, MaxLFQ,Protein.Q.Value and the peptide sequence used for the MaxLFQ. 

* __Visualization__

The plots currently present are essentially centered on the global analysis of proteins by condition. For the Upset plot ((like a Venn diagramm with bar, an app that i recommend just for this plot : https://github.com/hms-dbmi/UpSetR), the intersection data are downloadable.


- `Boxplot ` 
- `Barplot`
- `Multiscatterplot`
- `Upset plot` 
- `Violin plot`


