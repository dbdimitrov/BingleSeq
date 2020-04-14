# BingleSeq
BingleSeq - A user-friendly R package for Bulk and Single-cell RNA-Seq data analyses


## About
BingleSeq was developed as part of my MSc Bioinformatics thesis at the University of Glasgow under the supervision of Dr Quan Gu. The application provides a comprehesnive solution to both Bulk and scRNA-Seq analyses, as such it is best to look at BingleSeq as two separate parts.


### Bulk RNA-Seq
The Bulk RNA-Seq part follows the structure of a typical pipeline used for the DE analysis of Bulk RNA-Seq count data and it makes use of differential expression packages DESeq2 (Love, Huber, and Anders, 2014), edgeR (Robinson et al., 2010), and Limma (Ritchie et al., 2015).

![BingleSeq Bulk RNA-Seq pipeline](/figures/Bulk.png)


#### Input

BingleSeq's Bulk RNA-Seq pipeline accepts count tables in the following format:

![BingleSeq Bulk RNA-Seq format](/figures/Bulk_format.jpg)

**Note that a metadata table must also be provided for the Bulk RNA-Seq pipeline.**

#### Metadata table format
A metadata table is also required and its appropriate formating is key for the acquisition of correct results using the Bulk RNA-Seq part of the application. Metadata tables must be in this specific format:

![BingleSeq Bulk meta format](/figures/meta_format.jpg)

*The use of metadata tables was insipired by similar applications preceding BingleSeq - DEapp and DEBrowser (Li and Andrade, 2017; Kucukural et al., 2019).*



### scRNA-Seq
The scRNA-Seq part is based on Seurat’s pipeline (Satija et al., 2015) and follows a typical scRNA-Seq count analysis structure. Furthermore, clustering can be performed with monocle and SC3 packages (Trapnell et al., 2014; Kiselev et al., 2017). 

![BingleSeq Single-Cell RNA-Seq pipeline](/figures/sc.jpg)


#### Input

BingleSeq's scRNA-Seq pipeline accepts 10x genomics data as well as count tables in the following format:

![BingleSeq Bulk RNA-Seq format](/figures/sc_format.jpg)


## Typical Workflow

### Bulk RNA-Seq Data analysis typical workflow.

#### 0. Data used
For the purpose of representation, a simulated two-group dataset with 4 replicates was generated with compcodeR package (Soneson, 2014).

#### 1.	Load Count Table and Metadata Table
To begin DE analysis of Bulk RNA-Seq data, first a count table with a specific format (shown above) must be loaded. BingleSeq allows some flexibility in terms of the ‘separator’ used. Once a count table is loaded, it is displayed to allow the user to check whether it was loaded appropriately. A metadata table is also required in a specific format (shown above).

![BingleSeq Bulk RNA-Seq Load Data](/figures/bulk_loadData.PNG)




#### 2.	Quality Control
When a count table is uploaded the ‘Quality Control’ tab is generated. This tab enables the user to filter genes below certain counts per million (CPM), Max, or Median thresholds. The results of the filtering and the raw data are displayed as summary tables, alongside histograms.

![BingleSeq Bulk RNA-Seq QC Data](/figures/bulk_qcData.PNG)

Moreover, if required Batch effect correction is available with Harman and ComBat packages (Leek et al. 2016; Oytam et al., 2019).

![BingleSeq Bulk RNA-Seq Batch Correction](/figures/batchCorrected.jpg)
*Batch effect corrected using Harman package on example data acquired from HarmanData package (Bowden, Ross, Oytam, 2019).*


#### 3.	Differential Expression 
Subsequent to filtering is the ‘Differential Expression’ tab, where the user is given the option to run DESeq2, edgeR, and limma pipelines. Upon DE pipeline completion, the results are displayed as a table that contains the log2 expression fold-change (logFC), package specific test statistics, p-value, and multiple-testing adjusted p-value (FDR).

![BingleSeq Bulk RNA-Seq de Data](/figures/bulk_deData.PNG)

*Note that normalization between samples is done automatically with DE analysis using the package-specific methods for DESeq2 and edgeR, while limma also uses edgeR’s TMM normalization approach.*



#### 4.	Visualization
In a typical DE analysis workflow, the next stage following DE analysis would be to proceed to the various visualization techniques. BingleSeq offers this functionality within the ‘Visualize Data’ tab, which allows users to pick from several key plotting techniques including a PCA plot (A), Scree plot (B), Barchart (C), Volcano plot (D), MA plot (E), Venn Diagram (F), and a Heatmap (shown below). 

![BingleSeq Bulk RNA-Seq de Data](/figures/bulk_visData.PNG)


BingleSeq's visualization techniques were implemented with customization in mind, as users can specify parameters such as p-value threshold, fold-change threshold, and contrast of interest. Due to their versatility, heatmaps were designed as BingleSeq’s most customizable plotting component. Furthemore, users can download the genes displayed in the heatmap.

![BingleSeq Bulk RNA-Seq heat Data](/figures/bulk_heatmap.PNG)


#### 5.	Functional Annotation
Following DE analysis, BingleSeq enables the Functional Annotation analysis of DE results within the ‘Functional Annotation’ tab via GOseq package (Young et al., 2010). The GOseq pipeline enables users to obtain results from KEGG pathway analysis and three types of GO categories, including ‘Cellular Component’, ‘Molecular Function’, and ‘Biological Function’. To run the pipeline users are first prompted to filter the DEGs according to logFC and adjusted p-value (FDR). Users can then select several parameters before running the pipeline with the previously obtained subset of DEGs. These parameters include the GO category, multiple-testing corrected or uncorrected p-value, gene symbol type, and genome of interest.

![BingleSeq Bulk RNA-Seq bulkGO Tab](/figures/bulk_GOtab.PNG)

Once the GOseq pipeline is run and completed, a table with results is returned.

![BingleSeq Bulk RNA-Seq bulkGO Results](/figures/bulk_GOresults.PNG)
*Note that these results were generated using real data taken from McFarlane et al., 2019.*


Users can also generate GO term histograms with the top 10 GO terms and to choose whether to display their GO identifiers (GO:IDs) or their corresponding terms.

![BingleSeq Bulk RNA-Seq bulkGO hist](/figures/bulk_GOhist.jpg)

Moreover, users can obtain further information about a given GO term by querying its GO:ID using the ‘GO.db’ package (Carlson et al., 2019). 

![BingleSeq Bulk RNA-Seq bulkGO query](/figures/bulk_GOquery.PNG)
  
Note that in the current state of BingleSeq, only Mouse and Human genomes are supported.
  
  
#### 6.	DE Package Comparison
BingleSeq supplies users with an option to assess the agreement between the different DE analysis packages. This is done using a Venn diagram which represents the overlap of DE analysis results obtained using DESeq2, edgeR, and limma on the same dataset. Moreover, users can download the genes from the different intersects of the Venn Diagram. 

![BingleSeq Bulk RNA-Seq compVenn Data](/figures/bulk_compVenn.PNG)

Furthermore, the DE results from the same packages are used to generate a **Rank-based consesus**. The Rank-based consesus is displayed as a table, alongisde adjusted p-values and the ranks for each gene, as calculated by the different packages.

![BingleSeq Bulk RNA-Seq compRank Data](/figures/bulk_compRank.PNG)




### Single-Cell RNA-Seq Data analysis typical workflow.

#### 0. Data used
For the purpose of representation, a 10x genomics dataset was utilized (https://bit.ly/2Z3IUUk), which is also used in the Seurat package tutorial (https://bit.ly/2HlBfKx).
  
  
#### 1.	Load Count Data
To begin the analysis of scRNA-Seq data, users can supply gene count data in two input types. The first input type is ‘10x Genomics data’ in the form of a directory containing 10x Genomics protocol output files. These files include a matrix.mtx, barcodes.tsv, and genes.tsv files which represent the expression matrix, cell barcodes, and gene symbols, respectively. The second input type is a count table which must follow a specific format (shown above).
  
  
#### 2.	Quality Control
Once the data is loaded, the ‘Quality Control’ tab is generated which enables users to filter unwanted cells and features. Users can filter genes detected below a certain number of cells and cells with less than a certain number of expressed genes. These parameters are used when creating the initial Seurat object. Cell outliers can then be filtered according to the number of expressed features (i.e. genes) per cell. Visual aid is provided for filtering in terms of Violin plots which represent the number of genes (nFeature) and unique molecules (nCount_RNA) per cell.

![BingleSeq Bulk RNA-Seq sc qcData](/figures/sc_qcData.PNG)
  
  
#### 3.	Normalization
After excluding unwanted cells and features from the dataset, the next step is to normalize the data. BingleSeq provides two Seurat-based global-scaling normalization options. The first one is the “LogNormalize” method in which gene counts for each cell are divided by the total counts for that cell, multiplied by a ‘scale factor’, and then natural-log transformed. The second method is “Relative Counts” which follows the same procedure excluding the log transformation. Seurat's authors recommend using the former method when working with integer counts and the latter when working with relative counts. 10e4 is the recommended and default scale factor option, but when using CPM values the scale factor should be set to 10e6.

Simultaneously with normalization, the highly variable features within the dataset are identified and these features are later used when clustering with Seurat.

![BingleSeq Bulk RNA-Seq sc normData](/figures/sc_normData.PNG)

Once normalization and feature selection methods are complete, a Variable Features plot is returned and displayed. Seurat’s ‘Feature Selection’ methods include “VST”, “Mean Variance Plot”, and “Dispersion”.

![BingleSeq Bulk RNA-Seq sc varPlot](/figures/sc_variancePlot.PNG)
*This plot was genered using the recommended/default settings and the "VST" variance estimation method.*
  
*Note that ‘Feature selection’ does not apply to monocle and SC3 clustering approaches and hence their inbuilt pre-clustering filter procedures were implemented. These procedures have a similar purpose to Seurat’s ‘Feature Selection’, as they can be used to filter out unwanted noise.*
  
  
#### 4.	Clustering
Following normalization, the ‘Clustering’ tab is generated which implements functionality for scaling of the data, dimensionality reduction with PCA, PC selection, and unsupervised clustering. The former three are done simultaneously and Seurat’s ‘PCElbowPlot’ is used to generate and return an elbow plot. The returned elbow plot provides an ad hoc method for determining the true dimensionality of the dataset (i.e. PC Selection). Selecting which PCs to include in Seurat and monocle clustering methods is an essential step as it enables a large portion of technical noise to be excluded.

![BingleSeq Bulk RNA-Seq sc clustEblow](/figures/sc_clustElbow.PNG)
  
  
In addition to the Elbow plot, BingleSeq implements Seurat's PC heatmaps option - to be used as a complementary tool to the heuristic nature of the elbow plot.
  
![BingleSeq Bulk RNA-Seq sc clustHeat](/figures/sc_clustHeat.PNG)

*A) represents the 1st PC Heatmap with the top 10 most variable Genes and it is highly likely to represent part of the true dimensionality of the dataset. In contrast, B) represents the 15th PC Heatmap which likely represent mainly noise and not true signal.*



Once the count data is scaled and linear dimensionality reduction performed, users can proceed to unsupervised clustering with Seurat, SC3, and monocle. 

When using Seurat for unsupervised clustering, users must specify the number of PCs to be included in the analysis as well as the value of its ‘Resolution’ parameter. The latter parameter is used to set the ‘granularity’ of the clustering and as such it controls the number of clusters. The authors suggest that the optimal Resolution for datasets with ~3000 cells is 0.6-1.2 and it is typically higher for larger datasets. Users can also choose from Seurat’s inbuilt clustering algorithms including Louvain and SLM algorithms.

![BingleSeq Bulk RNA-Seq sc clustSeurat](/figures/sc_clustSeurat.PNG)

*tSNE plot produced using 0.5 as granularity parameter and the first 10 PCs.*
  
  
When clustering with monocle, users are requested to specify the number of PCs to be included in the analysis. Also, if required users can further minimize noise by filtering the gene counts according to the minimum expression level via the ‘Lower Detection Parameter’. Users can also pick from monocle’s inbuilt clustering algorithms, which include Density Peak and Louvain algorithms. Furthermore, Monocle enables the number of clusters to be explicitly specified as well as to be estimated.
  
  
![BingleSeq Bulk RNA-Seq sc clustMono](/figures/sc_clustMonocle.PNG)

*tSNE plot produced by explicitly setting the number of clusters to 9 and using the first 10 PCs (without any additional filtering).*


Unsupervised clustering with SC3 in BingleSeq utilizes SC3's k-means-based clusering approach. Users must specify the number of random initial centroid selections (sets). A larger number of initial centroid configurations (nStart) is likely to produce a better clustering result, but has a high toll on computational time. By default, this parameter is set to 1000 when working with less than 2000 cells and to 50 when working with more than 2000 cells, in accordance to the authors recommendantions.
Users can also use SC3’s inbuilt filtering options to further reduce noise by filtering out genes below and above certain dropout (zero value) percentage thresholds. Similarly to monocle, the number of clusters can be supplied by user or estimated with SC3.
It is worth noting that the k-means approach of SC3 is likely too computantionally demanding when working with large datasets (e.g. when N=2000, computational time is ~20 mins), hence future updates of BingleSeq are likely to also implement SC3's SVM-hybrid approach as an alternative solution. 

![BingleSeq Bulk RNA-Seq sc clustSC3](/figures/sc_clustSC3.PNG)

*tSNE plot produced by explicitly setting the number of clusters to 9 and 50 random initial centroid sets (without any additional filtering).*
  
  
*Each tSNE plot in BingleSeq is generated using the package with which clustering was performed.
Also, when performing clustering with SC3 or monocle, the data used to create the required objects to run these pipelines is the same data that was previously filtered, normalized, and scaled using Seurat’s pipeline.*


#### 5.	Differential Expression
Following clustering, DE analysis can be conducted using Seurat’s inbuilt functionality to identify marker genes. Users can perform marker gene identification using the following inbuilt DE testing methods: Student’s T test, Wilcoxon Rank Sum test, and Logistic regression. Additionally, DE analysis can also be performed with DESEq2 and MAST packages (Love, Huber, and Anders, 2014; Finak et al., 2015). It is worth noting that DESeq2 is magnitudes slower than other DE tests; thus, making it impractical when working with large datasets.
Prior to running the DE analysis, users are prompted to enter the following filtering parameters: genes expressed in a minimum fraction of cells, fold-change, and adjusted p-value.

![BingleSeq Bulk RNA-Seq sc deTab](/figures/sc_deTab.PNG)

Furthermore, by using Seurat’s inbuilt visualization options, BingleSeq provides tools for the exploration of DE results. These tools include cluster heatmap with user-specified gene number as well as exploration of specific genes via Violin, Feature, and Ridge plots.

![BingleSeq Bulk RNA-Seq sc deFigs](/figures/sc_deFigs.PNG)

*A) Heatmap showing the top 10 genes for each cluster in the 2700 PBMCs dataset, while Violin B), Feature C), and Ridge D) plots are shown for MS4A1 gene – a biomarker of B lymphocytes.*


#### 6.	Functional Annotation
The scRNA-Seq part of BingleSeq incorporates functional annotation in an analogous manner to its Bulk RNA-Seq counterpart. The only difference is that the subsets of DEGs to be used in the GOseq pipeline can be filtered according to the cluster they belong to; thus, allowing users to assess each cluster independently. 
  
  
#### 7.	DE Method Comparison
The scRNA-Seq part also implements a ‘DE Method Comparison’ tab analogous to the ‘DE Package Comparison’ tab in Bulk RNA-Seq. The only difference is that scRNA-Seq Overlap functionality enables filtering according to the same parameters used in marker gene identification. Furthermore, rather than comparing the different packages, it compares the DE Methods implemented within Seurat. These include: DE testing with MAST, Wilcoxon Rank Sum Test, and Student’s T test.

*Also, note that Rank-based consensus is yet to be implemented for the scRNA-Seq pipeline.*



### Prerequisites

BingleSeq requires R>= 3.5.3



### Installing

To install BingleSeq on your machine simply copy the following R code:

```
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
setwd("*Folder Containing BingleSeq's zip file*")
d <- getwd()

untar(file.path(d, "BingleSeq_0.3.0.tar.gz"), exdir=d)
devtools::install(file.path(d, "BingleSeq"), dependencies=TRUE,
                  repos=BiocManager::repositories())


library(BingleSeq)  # Load BingleSeq
BingleSeq::startBingleSeq()  # Starts the application
```


## Built With

* [shiny](https://shiny.rstudio.com/)
* [RStudio](https://www.rstudio.com/)


## References
Bowden J, Ross J, Oytam Y (2019). HarmanData: Data for the Harman package. R package version 1.12.0, http://www.bioinformatics.csiro.au/harman/.

Carlson, M. Falcon, S., Pages, H., Li, N., 2019. GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.8.2.

Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A.K., Slichter, C.K., Miller, H.W., McElrath, M.J., Prlic, M. and Linsley, P.S., 2015. MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome biology, 16(278) doi: 10.1186/s13059-015-0844-5.

Kiselev, V.Y., Kirschner, K., Schaub, M.T., Andrews, T., Yiu, A., Chandra, T., Natarajan, K.N., Reik, W., Barahona, M., Green, A.R. and Hemberg, M., 2017. SC3: consensus clustering of single-cell RNA-seq data. Nature methods, 14(5), pp. 483–486.

Kucukural, A., Yukselen, O., Ozata, D.M., Moore, M.J. and Garber, M., 2019. DEBrowser: interactive differential expression analysis and visualization tool for count data. BMC genomics, 20(6). DOI: https://doi.org/10.1186/s12864-018-5362-x

Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Storey JD, Zhang Y, Torres LC (2019). sva: Surrogate Variable Analysis. R package version 3.32.1.

Li, Y. and Andrade, J., 2017. DEApp: an interactive web interface for differential expression analysis of next generation sequence data. Source code for biology and medicine, 12(2), doi: 10.1186/s13029-017-0063-4

Love, M.I., Huber, W. and Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12). doi:10.1186/s13059-014-0550-8

McFarlane, S., Orr, A., Roberts, A.P., Conn, K.L., Iliev, V., Loney, C., da Silva Filipe, A., Smollett, K., Gu, Q., Robertson, N. and Adams, P.D., 2019. The histone chaperone HIRA promotes the induction of host innate immune defences in response to HSV-1 infection. PLoS pathogens, 15(3). doi:10.1371/journal.ppat.1007667

Oytam Y, Sobhanmanesh F, Duesing K, Bowden JC, Osmond-McLeod M, Ross J (2016). “Risk-conscious correction of batch effects: maximising information extraction from high-throughput genomic datasets.” BMC Bioinformatics, 17(1), 1–17. doi: 10.1186/s12859-016-1212-5, http://dx.doi.org/10.1186/s12859-016-1212-5.

Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. and Smyth, G.K., 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), pp.e47-e47.

Robinson, M.D., McCarthy, D.J. and Smyth, G.K., 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), pp.139-140.

Satija, R., Farrell, J.A., Gennert, D., Schier, A.F. and Regev, A., 2015. Spatial reconstruction of single-cell gene expression data. Nature biotechnology, 33(5), pp.495-502.

Soneson, C., 2014. compcodeR—an R package for benchmarking differential expression methods for RNA-seq data, Bioinformatics, 30(17), pp.2517–2518.

Trapnell, C., Cacchiarelli, D., Grimsby, J., Pokharel, P., Li, S., Morse, M., Lennon, N.J., Livak, K.J., Mikkelsen, T.S. and Rinn, J.L., 2014. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nature biotechnology, 32(4), pp. 381–386.

Young, M.D., Wakefield, M.J., Smyth, G.K. and Oshlack, A., 2010. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome biology, 11(R14). doi: https://doi.org/10.1186/gb-2010-11-2-r14.



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
