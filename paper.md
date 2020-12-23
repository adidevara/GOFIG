---
title: 'GOFIG: A tool for Gene Ontology Enrichment Analysis and Visualization'
tags:
  - R
  - Gene Ontology 
  - Enrichment Analysis 
  - Data visualization 
authors:
  - name: Aditya K. Devarakonda
    orcid: 0000-0002-7014-3029
    affiliation: "1"
  - name: Eric G. Rafalovsky 
    affiliation: "2"
  - name: Tae-Jin Lee Ph.D. 
    affiliation: "1"
  - name: Ashok Sharma Ph.D. 
    orcid: 0000-0001-9597-4374
    affiliation: "1"
affiliations:
 - name: Center for Biotechnology and Genomic Medicine, Augusta University, Augusta, Georgia, United States
   index: 1
 - name: Georgia Institute of Technology, Atlanta, Georgia, United States
   index: 2
date: 22 December 2020
bibliography: paper.bib
---


  

# Summary

The advent of high throughput sequencing has fundamentally changed both basic science practices and clinical medicine. RNA sequencing can be used to identify change genetic changes and epigenetic regulation in a variety of cell types and conditions. A key tool in any cellular biologists arsenal to parse out these effects is differential expression analysis. While this is a means of evaluating the changes in cellular function, the output lists containing hundreds of gene identifiers can be esoteric. A contextualized understanding of these results can be attained through a process known as Gene Ontology Enrichment Analysis. This process hierarchically categorizes genes into categories with terms that are much easier to understand and display.

# Statement of need

`GOFIG` is an R package which allows for Gene Ontology enrichment analysis along with downstream ontology comparison and data visualization. While there are a few existing programs which allow for Gene Ontology analysis, few are as easy to use while providing as much functionality as `GOFIG`.

The pipeline begin with a simple .csv file containing a list of gene identifications and fold change values. Using the `org.Hs.eg.db` package [@Carlson2020] in combination with the `limma` package [@Ritchie2015], the package is able to produce an enrichment analysis from various different gene identification types along with data regarding the number of up regulated and down regulated genes in each category. This feature only exists in a few packages for enrichment analysis currently. As the output is saved as a .csv, users can easily view the output while also marking ontologies which may be very similar to others and may not warrant display in the visualization.

For downstream visualization, a combination of the `dplyr` package [@Wickham2020] and `ggplot2` [@Wickham2016] were used. There are two different types of bar plots which can be produced: A single bar plot showing the number of genes in the most significant ontologies along with points to represent the respective p-value of these ontologies as well as a back-to-back bar plot which provides the same information along with the subdivisions of up regulated and down regulated genes in each ontology. Another visualization which can be produced is a bubble plot which is plotted based on a normalization to depict whether the ontology contains more up regulated or down regulated genes against the log of the p-values. The size of the bubbles is a representation of the number of genes in the ontology.

Aside from ease of use, `GOFIG` offers the unique functionality of comparing the overlapping ontologies between two different enrichment analyses and visualizing them through a Venn Diagram. This is done through the use of upstream enrichment analysis outputs along with `VennDiagram` package [@Chen2018].

The functionality along with the ease of use makes `GOFIG` and optimal tool for investigators looking to create data visuals to better represent the meaning of differential expression analysis data from RNA sequencing.


# References
