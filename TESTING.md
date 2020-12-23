# How to Test 



## Script, Libraries, and directories

**Script**

Please download the GOFIG.r script from the GitHub repository. If using RStudio, please be sure to source the script after opening. If using R, please run the script prior to using. 

**Directory**

The functions will automatically search and save the files to the working directory. Please be sure to set your working directory to the directory containing the input files using the `setwd()` function. 

**Library**

GOFIG is dependent on `dplyr`, `VennDiagram`, and `ggplot2` from the CRAN repository. To install and load these packages, please type the following into the terminal: 

```
install.packages("dplyr")
install.packages("VennDiagram") 
install.packages("ggplot2") 

library(dplyr)
library(VennDiagram)
library(ggplot2) 
```

GOFIG is also dependent on `org.Hs.eg.db` and `limma `from the BioConductor repository. To install and load these packages, please type the following into the terminal: 

```
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")

library(limma)
library(org.Hs.eg.db) 
```

After installing and loading the following, the script should have all of the necessary libraries loaded.

## Test Files 

With in the repository, there are 2 input files: InputFile1.csv and InputFile2.csv. 

These are the files which will be crucial for the analysis of the function in the GOFIG package. 

InputFile1.csv and InputFile2.csv must both be run through the `geneOnt()` function as both outputs will be needed for the `ontComp()` function. However, only the output from one of the other will be required to run the `ontBar()`, `ontb2bBar()`, and `ontBubble()` functions. 


## geneOnt() 

The `geneOnt()` function has arguments: the input file name, the gene ID type, and the output file name. 

Using InputFile1.csv, the `geneOnt()` function can be executed through typing the following: 

geneOnt("InputFile1.csv", "SYMBOL", "OutputFile1.csv") 

 - As this input file ID type is symbol, we can use "SYMBOL" as an argument. For the full list of potential ID type conversion please refer to the README .md file in the repository. 
 
 - The output name is trivial, however in the remainder of this guide it will be referred to as "OutputFile1.csv". 
 
 - This function in particular is prone to taking a fair amount of time. The genes are all being sorted through over 22,000 ontologies each, and as such, run times can be long depending on the machine it is being run on. 

As mentioned previously, we should run this function on InputFile2.csv as well as the output from this file will be used for other downstream analysis. This can be run the same way as InputFile1.csv, only the argument specifying the input file should be changed to InputFile2.csv and the argument specifying the output file should be changed to OutputFile2.csv. 

The out put file will contain the following columns. 

1. Gene Ontology ID
2. Enriched Term
3. Total number of genes associated with the ontology term
4. Number of genes in the ontology term from the input list
5. P-value of enrichment
6. Number of up regulated genes
7. Number of down regulated genes
8. Difference between the up and down regulated genes (a z-score like metric)
9. Option to include the term for Visualization.


## ontBar()

The `ontBar()` function has four arguments: the input file, the ontology cutoff (which is defaulted to 10), the color of the bars representing the number of genes in the ontology, and the color of the bars representing number of genes in the ontology (which is defaulted to light green), and the color of the point representing the P-value (which is defaulted to dark blue).

Using OutputFile1.csv, the `ontBar()` function can be executed through typing the following: 

ontBar("OutputFile1.csv")

 - The ontology cutoff can be changed using any number as the second argument. 
 
 - The colors of the bars and the points can be changed using strings representing the recognized R color values (ex. "red", "lightyellow") in their respective argument locations. 
 
 A PDF containing a bar plot for cellular component, molecular function, and biological processes will be saved to the working directory. 
 
## ontb2bBar() 

The `ontb2bBar()` function has five arguments: the input file, the ontology cutoff (which is defaulted to 10), the color of the bars representing the number of up regulated genes in the ontology (which is defaulted to light green), the color of the bars representing the number of down regulated genes in the ontology (which is defaulted to light blue), and the color of the point representing the P-value (which is defaulted to dark blue).

 
Using OutputFile1.csv, the `ontb2bBar()` function can be executed through typing the following: 

ontb2bBar("OutputFile1.csv")

 - The ontology cutoff can be changed using any number as the second argument. 

 - The colors of the bars and the points can be changed using strings representing the recognized R color values (ex. "red", "lightyellow") in their respective argument locations. 

 A PDF containing a back-to-back bar plots for cellular component, molecular function, and biological processes will be saved to the working directory. 


## ontBubble()

The `ontBubble()` function has four arguments: the input file, the ontology cutoff (which is defaulted to 10), the bubble outline color (which is defaulted to dark blue) and the bubble fill which is defaulted to light blue. 

Using OutputFile1.csv, the `ontBubble()` function can be executed through typing the following: 

- The ontology cutoff can be changed using any number as the second argument. 

 - The colors of the bubble outline and the bubble fill can be changed using strings representing the recognized R color values (ex. "red", "lightyellow") in their respective argument locations. 

 A PDF containing bubble plots for cellular component, molecular function, and biological processes will be saved to the working directory. 

## ontComp() 

The `ontComp()` function has two arguments: the name of the first input file (which should be an output of the `geneOnt()` function) and a second input file which should be a different output of the `geneOnt()` function. 

Using OutputFile1.csv and OutputFile2.csv the `ontComp()` function can be executed through typing the following: 

ontComp("OutputFile1.csv", "OutputFile2.csv") 

A CSV file will be produced of the name "out.csv". This output will have four columns in the following order: 

1. Gene Ontology Identification 
2. Gene Ontology Term 
3. The ontology category (BP for biological processes, CC for cellular component, and MF for molecular function) 
4. Whether the ontology is a component of the first input, the second input, or both.

## ontVenn() 

The `ontVenn()` function has six inputs: the input file (which is an output of the `ontComp()` function, the name of the first input file (which will be used to label that side of the Venn Diagram and is defaulted to Input 1), the name of the second input file (which will be used to label that side of the Venn Diagram and is defaulted to Input 2), the color of the first input in the diagram (which is defaulted to light green), the color of the second input (which is defaulted to blue), and the title of the graph (which is defaulted to Gene Ontology Comparison). 

Using the out.csv file which is the output of the `ontComp()` function, the `ontVenn()` function can be executed through typing the following: 

ontComp("out.csv") 

 - The colors of the circles in the Venn diagram can be changed using strings representing the recognized R color values (ex. "red", "lightyellow") in their respective argument locations. 
 
 The output of this function will be saved as a PNG file to the working directory. 

## ontVennSub() 

The `ontVennSub()` function has five inputs: the input file (which is an output of the `ontComp()` function, the name of the first input file (which will be used to label that side of the Venn Diagram and is defaulted to Input 1), the name of the second input file (which will be used to label that side of the Venn Diagram and is defaulted to Input 2), the color of the first input in the diagram (which is defaulted to light green), and the color of the second input (which is defaulted to blue).

Using the out.csv file which is the output of the `ontComp()` function, the `ontVennSub()` function can be executed through typing the following: 

 - The colors of the circles in the Venn diagram can be changed using strings representing the recognized R color values (ex. "red", "lightyellow") in their respective argument locations. 
 
 The outputs of this function will be saved as three PNG files to the working directory, one for cellular component, one for molecular function, one for biological processes. 

## Contact Information 
For any further clarification please do not hesitate to contact me for support at adevarakonda@augusta.edu
