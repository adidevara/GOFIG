
# GOFIG: A tool for Gene Ontology Enrichment Analysis and Visualization

Aditya K. Devarakonda, Eric G. Rafalovsky, Tae Lee Jin Ph.D., Ashok Sharma Ph.D.


# Introduction

GOFIG is a tool for gene ontology enrichment analysis and visualization. This tool is very helpful in generating publication quality figures for the manuscripts. A list of differentially expressed genes and log fold-changes are used as input. Various ontology terms over-represented in the gene list are identified and p-values are assigned, which represent the significance of enrichment of each ontology term. Along with the ontology terms and the p-values, the analysis output includes the total number of genes, the number of upregulated genes, and downregulated genes in each ontology term. Ontology terms will only be included in the output if they contain at least 10% of the total number of differentially expressed genes provided in the input list. All p-values are calculated based on the total number of genes in the ontology rather than just the number of upregulated or downregulated genes. After analysis, the output is saved to the working directory as a CSV file. This file can be used for image production using “ontBar” function available in the package.

## Required Packages

The required packages which must be installed and loaded for the functions in GOFIG to work include: ggplot2, dplyr, VennDiagram, limma, and org.Hs.eg.db. 

ggplot2, dplyr, and VennDiagram can all be loaded in CRAN using the install.packages function. 

limma and org.Hs.eg.db must be loaded through Bioconductor. To do this, use the install.packages function to install BiocManager after which you can use the BiocManager::install package to load limma and org.Hs.eg.db.

## Directory Setting 

The package automatically searches for input files and saves output files to the working directory. In order to ensure that the package works properly, make sure to set your working directory to the directory which contains your input file.

## geneOnt : A function for gene ontology enrichment analysis



**geneOnt(InputFile = “”, IDType = “SYMBOL”, OutputName = “”)**

InputFile = “” : The input file name which should be formatted as a CSV file with two columns. The first column should contain a list of gene ID’s and the second column should contain the log fold changes for each of the genes.

IDType = “”: The following types of Gene IDs can be used. If there is no ID type specified in the second argument, the default ID will be set to gene symbol (“SYMBOL”).

1. GenBank Accession Numbers  : “ACCNUM”
2. Ensembl ID  : “ENSEMBL”
3. Ensembl Protein Accession Numbers : “ENSEMBLPROT”
4. Ensembl Transcript Accession Numbers:“ENSEMBLTRANS”
5. Entrez Identification  : “ENTREZID”
6. Enzyme Commission Numbers  : “ENZYME”
7. Gene Name  :  “GENENAME”
8. Prosite Identifiers  :  “PROSITE”
9. RefSeq Identifiers  :  “REFSEQ”
10. Gene Symbol  : “SYMBOL”
11. Uniprot Accession Numbers  : “UNIPROT”

The third argument is the name of the output file and should include the “.csv” extension. The output of the analysis will be saved in the working directory. The output file will contain the following 11 columns:

1. Gene Ontology ID
2. Enriched Term
3. Total number of genes associated with the ontology term
4. Number of genes in the ontology term from the input list
5. P-value of enrichment
6. Number of up regulated genes
7. Number of down regulated genes
8. Difference between the up and down regulated genes (a z-score like metric)*
9. Option to include the term for Visualization.**

** The z-score like metric is simply a means to display the difference between the up and down regulated genes in relation to the total number of genes. The formula we used to produce this metric is: 

 $\frac{(up \ regulated \ genes - down \ regulated \ genes}{\sqrt[2]{(total \ genes)}}$  

**Option to include the term for Visualization allows for the user to select which gene ontology terms will be included to generate the final Figure. As some of the terms are very similar it may be redundant to include them in the visuals. If the user chooses to change the default entry from Y to N, the gene ontology will not be included in the visualization.

## ontBar(), ontb2bBar(), ontBubble: Functions for visualization of the gene ontology enrichment analysis results

**ontBar(InputFile = "", OntologyCutoff = "10", BarColor = "lightgreen", PointColor = "darkblue")**

This function produces three bar plots (cellular components, molecular functions, and biological processes), containing the most significant ontologies (by p-value) in each category. A line depicting the p-value of each ontology is also overlayed on the barplots.

InputFile = "": The name of the csv file generated from the _geneOnt_ function of the package as previously described.

OntologyCutoff: The number of ontologies to be included in each graph. The default value is 10, and if this number exceeds the number of ontologies present in the analysis, all will be included.

BarColor: The color of the bar. The default is “lightgreen”.

PointColor: The color of the p-value line. The default is “darkblue”.

**ontb2bBar****(InputFile = "", OntologyCutoff = "10", UpBarColor = "lightgreen", DownBarColor = "lightblue", PointColor = "darkblue")**

The above function production produces a three back to back bar plots, each containing the most significant ontologies (by p-value) in the cellular component, molecular function, and biological processes category. The bar on the right side represents the number of upregulated genes in an ontology and the bar on the left side represents the number of downregulated genes in an ontology. Overlayed on the bars on the right is a series of points which depict the p-value of each ontology.

InputFile = "": The name of the csv file generated from the _geneOnt_ function of the package as previously described.

OntologyCutoff: The number of ontologies to be included in each graph. The default value is 10, and if this number exceeds the number of ontologies present in the analysis, all will be included.

UpBarColor: The color of the bar representing the upregulated genes in the ontology. The default is “lightgreen”.

DownBarColor: The color of the bar representing the downregulated genes in the ontology. The default is “lightblue”.

PointColor: The color of the p-value line. The default is “darkblue”.

**ontBubble(****InputFile = "", OntologyCutoff = "10", BubbleOutline = “darkblue”, BubbleFill = “lightblue”)**

The above function production produces a three bubble plots each containing the most significant ontologies (by p-value) in the cellular component, molecular function, and biological processes category. The they are distributed across the x-axis based on their pseudo-z-score ( the calculation of which was described above, as well as the log of their adjusted p-values across the y-axis.

InputFile = "": The name of the csv file generated from the _geneOnt_ function of the package as previously described.

OntologyCutoff: The number of ontologies to be included in each graph. The default value is 10, and if this number exceeds the number of ontologies present in the analysis, all will be included.

BubbleOutline: The outline color of the bubbles produced on the plot. The default is “darkblue”.

BubbleFill: The fill color of the bubbles produced on the plot. The default is “lightblue”.

## ontComp(): A function for gene ontology enrichment analysis comparison

**ontComp(InputOne = "", InputTwo = "")**

The above function compares the values in two separate outputs of the geneOnt() functions by gene ontology (GO) IDs and writes a CSV to the working directory containing 4 columns: the GO IDs, the GO terms, the major ontology within which they are contained (cellular component, molecular function, or biological processes), and whether the ontology is present in the first input file, the second input file, or both.

InputOne = “”: The name of a CSV file produced from the geneOnt() function

InputTwo = “”: The name of a second CSV file produced from the geneOnt() function to compare against the first


# ontVenn() and ontVennSub(): A function for gene ontology enrichment analysis comparison

**ontVenn(InputFile, FirstInputFileName = "Input 1", SecondInputFileName = "Input 2", FirstInputColor = 'lightgreen', SecondInputColor = 'blue', Title = "Gene Ontology Comparison")**

The function above creates a Venn diagram of the ontologies listed from the output of the ontComp() function. The output file is saved as a PNG file to the working directory.

InputFile = “”: The name of the CSV file generated from the aformentioned ontComp() function

FirstInputFile: The name of the first input file to label the Venn diagram. The default is “Input 1”.

SecondInputFile: The name of the second input file to label the Venn diagram. The default is “Input 2”.

FirstInputColor: The color of the Venn diagram section for the first input. The default is ‘lightgreen’.

SecondInputColor: The color of the Venn diagram section for the first input. The default is ‘blue’.

Title: The title for the Venn diagram. The default in “Gene Ontology Comparison”.

**ontVennSub(InputFile, FirstInputFileName = "Input 1", SecondInputFileName = "Input 2", FirstInputColor = 'lightgreen', SecondInputColor = 'blue')**

The above function creates three Venn diagrams of the ontologies listed from the output of the ontComp() function. The diagram represents the ontologies present in the cellular component, molecular function, and biological processes subgroups. The output files are saved as PNG files to the working directory.

InputFile = “”: The name of the CSV file generated from the aformentioned ontComp() function

FirstInputFile: The name of the first input file to label the Venn diagram. The default is “Input 1”.

SecondInputFile: The name of the second input file to label the Venn diagram. The default is “Input 2”.

FirstInputColor: The color of the Venn diagram section for the first input. The default is ‘lightgreen’.

SecondInputColor: The color of the Venn diagram section for the first input. The default is ‘blue’.

Title: The title for the Venn diagram. The default in “Gene Ontology Comparison”.
