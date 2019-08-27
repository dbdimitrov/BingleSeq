# BingleSeq
BingleSeq - A user-friendly R package for Bulk and Single-cell RNA-Seq data analyses


## About
BingleSeq provides a comprehesnive solution to both Bulk and scRNA-Seq analyses, as such it is best to look at BingleSeq as two separate parts.




### Bulk RNA-Seq
The Bulk RNA-Seq part follows a typical pipeline used for the DE analysis of Bulk RNA-Seq count data and it makes use of differential expression packages DESeq2 (Love, Huber, and Anders, 2014), edgeR (Robinson et al., 2010), and Limma (Ritchie et al., 2015).

![BingleSeq Bulk RNA-Seq pipeline](/figures/Bulk.jpg)

#### Input

BingleSeq's Bulk RNA-Seq pipeline accepts count tables in the following format:

![BingleSeq Bulk RNA-Seq format](/figures/Bulk_format.jpg)


### scRNA-Seq
The scRNA-Seq part is based on Seurat’s pipeline (Satija et al., 2015) and follows a typical scRNA-Seq pipeline. Nonetheless, clustering can also be performed with monocle and SC3 packages (Trapnell et al., 2014; Kiselev et al., 2017). 

![BingleSeq Single-Cell RNA-Seq pipeline](/figures/sc.jpg)



#### Input

BingleSeq's scRNA-Seq pipeline accepts 10x genomics data as well as count tables in the following format:

![BingleSeq Bulk RNA-Seq format](/figures/sc_format.jpg)


### Prerequisites

BingleSeq requires R>= 3.5.3


### Installing

To install BingleSeq on your machine simply type the following R code:

```
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
setwd("*Folder Containing BingleSeq's zip file*")
d <- getwd()

untar(file.path(d, "BingleSeq_0.2.0.tar.gz"), exdir=d)
devtools::install(file.path(d, "BingleSeq"), dependencies=TRUE,
                  repos=BiocManager::repositories())


library(BingleSeq)  # Load BingleSeq
BingleSeq::startBingleSeq()  # Starts the application
```

## Typical Workflow

### Bulk RNA-Seq Data analysis typical workflow with BingleSeq.

#### 0. Data used
For the purpose of representation, a simulated two-group dataset with 4 replicates was generated with compcodeR package (Soneson, 2014).

#### 1.	Load Count Table
To begin DE analysis of Bulk RNA-Seq data, first a count table with a specific format (shown above) must be loaded. BingleSeq allows some flexibility in terms of the ‘separator’ used. Once a count table is loaded, it is displayed to allow the user to check whether it was loaded appropriately. 

![BingleSeq Bulk RNA-Seq Load Data](/figures/bulk_loadData.PNG)



#### 2.	Quality Control
When a count table is uploaded the ‘Quality Control’ tab is generated. This tab enables the user to filter genes below certain counts per million (CPM), Max, or Median thresholds. The results of the filtering and the raw data are displayed as summary tables, alongside histograms.

![BingleSeq Bulk RNA-Seq QC Data](/figures/bulk_qcData.PNG)


#### 3.	Differential Expression 
Subsequent to filtering is the ‘Differential Expression’ tab, where the user is given the option to run DESeq2, edgeR, and limma pipelines. Upon DE pipeline completion, the results are displayed as a table that contains the log2 expression fold-change (logFC), package specific test statistics, p-value, multiple-testing adjusted p-value (FDR).

![BingleSeq Bulk RNA-Seq de Data](/figures/bulk_deData.PNG)

*Note that normalization between samples is done automatically with DE analysis using the package-specific methods for DESeq2 and edgeR, while limma also uses edgeR’s TMM normalization approach.*



#### 4.	Visualization
In a typical DE analysis workflow, the next stage following DE analysis would be to proceed to the various visualization techniques. BingleSeq offers this functionality within the ‘Visualize Data’ tab, which allows users to pick from several key plotting techniques including a PCA plot (A), Scree plot (B), Barchart (C), Volcano plot (D), MA plot (E), Venn Diagram (F), and a Heatmap (G). *Note that the Venn Diagram was generated using a three-group simulated dataset.*

![BingleSeq Bulk RNA-Seq de Data](/figures/bulk_visData.PNG)

BingleSeq's visualization techniques were implemented with customization in mind, as users can specify parameters such as p-value threshold, fold-change threshold, and contrast of interest. Due to their broad usability, heatmaps (G) were designed as BingleSeq’s most customizable plotting component. Furthemore, users can download the genes displayed in the heatmap.

![BingleSeq Bulk RNA-Seq heat Data](/figures/bulk_heatmap.PNG)


#### 5.	Functional Annotation
Following DE analysis, BingleSeq enables the Functional Annotation analysis of DE results within the ‘Functional Annotation’ tab via GOseq package (Young et al., 2010). The GOseq pipeline enables users to obtain results from KEGG pathway analysis and three types of GO categories, including ‘Cellular Component’, ‘Molecular Function’, and ‘Biological Function’. To run the pipeline users are first prompted to filter the DEGs according to logFC and adjusted p-value (FDR). Users can then select several parameters before running the pipeline with the previously obtained subset of DEGs. These parameters include the GO category, multiple-testing corrected or uncorrected p-value, gene symbol type, and genome of interest.

![BingleSeq Bulk RNA-Seq bulkGO Tab](/figures/bulk_GOtab.PNG)

Once the GOseq pipeline is run and completed, a table with results is returned.

![BingleSeq Bulk RNA-Seq bulkGO Results](/figures/bulk_GOresults.PNG)
*Note that these results were obtained with real data taken from McFarlane et al., 2019*


Users can also generate GO term histograms with the top 10 GO terms and to choose whether to display their GO identifiers (GO:IDs) or their corresponding terms.

![BingleSeq Bulk RNA-Seq bulkGO hist](/figures/bulk_GOhist.PNG)

Moreover, users can obtain further information about a given GO term by querying its GO:ID using the ‘GO.db’ package (Carlson et al., 2019). Note that in the current state of BingleSeq, only Mouse and Human genomes are supported (Carlson, 2019A; Carlson, 2019B).

![BingleSeq Bulk RNA-Seq bulkGO query](/figures/bulk_GOquery.PNG)



#### 6.	DE Package Comparison
BingleSeq supplies users with an option to assess the agreement between the different DE analysis packages. This is done using a Venn diagram which represents the overlap of DE analysis results obtained using DESeq2, edgeR, and limma on the same dataset. Moreover, users can download the genes from the different intersects of the Venn Diagram. 

![BingleSeq Bulk RNA-Seq compVenn Data](/figures/bulk_compVenn.PNG)

Furthermore, the DE results from the same packages are used to generate a Rank-based consesus. The Rank-based consesus is displayed as a table, alongisde adjusted p-values and the ranks for each gene, as calculated by the different packages.

![BingleSeq Bulk RNA-Seq compRank Data](/figures/bulk_compRank.PNG)


## Built With

* [shiny](https://shiny.rstudio.com/)
* [RStudio](https://www.rstudio.com/)


## References
Kiselev, V.Y., Kirschner, K., Schaub, M.T., Andrews, T., Yiu, A., Chandra, T., Natarajan, K.N., Reik, W., Barahona, M., Green, A.R. and Hemberg, M., 2017. SC3: consensus clustering of single-cell RNA-seq data. Nature methods, 14(5), pp. 483–486.

Love, M.I., Huber, W. and Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12). doi:10.1186/s13059-014-0550-8

McFarlane, S., Orr, A., Roberts, A.P., Conn, K.L., Iliev, V., Loney, C., da Silva Filipe, A., Smollett, K., Gu, Q., Robertson, N. and Adams, P.D., 2019. The histone chaperone HIRA promotes the induction of host innate immune defences in response to HSV-1 infection. PLoS pathogens, 15(3). doi:10.1371/journal.ppat.1007667

Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. and Smyth, G.K., 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), pp.e47-e47.

Robinson, M.D., McCarthy, D.J. and Smyth, G.K., 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), pp.139-140.

Satija, R., Farrell, J.A., Gennert, D., Schier, A.F. and Regev, A., 2015. Spatial reconstruction of single-cell gene expression data. Nature biotechnology, 33(5), pp.495-502.

Trapnell, C., Cacchiarelli, D., Grimsby, J., Pokharel, P., Li, S., Morse, M., Lennon, N.J., Livak, K.J., Mikkelsen, T.S. and Rinn, J.L., 2014. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nature biotechnology, 32(4), pp. 381–386.

Young, M.D., Wakefield, M.J., Smyth, G.K. and Oshlack, A., 2010. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome biology, 11(R14). doi: https://doi.org/10.1186/gb-2010-11-2-r14.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
