# BingleSeq
BingleSeq - A user-friendly R package for Bulk and Single-cell RNA-Seq data analyses.

Manuscript available as a [PeerJ publication](https://doi.org/10.7717/peerj.10469).


### Installation

BingleSeq can be installed directly from GitHub using the following code:

```
library("devtools")
install_github("dbdimitrov/BingleSeq")

# Start the application
library(BingleSeq)
startBingleSeq()
```


### Prerequisites

BingleSeq requires R>= 3.6.3


## About
BingleSeq provides a comprehesnive solution to both Bulk and scRNA-Seq analyses, as such it is best to look at BingleSeq as two separate parts.


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
In a typical DE analysis workflow, the next stage following DE analysis would be to proceed to the various visualization techniques. BingleSeq offers this functionality within the ‘Visualize Data’ tab, which allows users to pick from several key plotting techniques. Including A) Barchart plot presenting the number of up- and downregulated genes B) PCA plot, C) Volcano plot, and D) MA plot. Note that these results were generated with limma using real data taken from McFarlane et al., 2019.

![BingleSeq Bulk RNA-Seq de Data](/figures/bulk_visData.PNG)


BingleSeq's visualization techniques were implemented with customization in mind, as users can specify parameters such as p-value threshold, fold-change threshold, and contrast of interest. Due to their versatility, heatmaps were designed as BingleSeq’s most customizable plotting component. Furthemore, users can download the genes displayed in the heatmap.

![BingleSeq Bulk RNA-Seq heat Data](/figures/bulk_heatmap.PNG)

#### Functional Annotation
##### 5a.	Over-representation Analysis
Following DE analysis, BingleSeq enables the Functional Annotation analysis of DE results within the ‘Functional Annotation’ tab via GOseq package (Young et al., 2010). The GOseq pipeline enables users to obtain results from KEGG pathway analysis and three types of GO categories, including ‘Cellular Component’, ‘Molecular Function’, and ‘Biological Function’. To run the pipeline users are first prompted to filter the DEGs according to logFC and adjusted p-value (FDR). Users can then select several parameters before running the pipeline with the previously obtained subset of DEGs. These parameters include the GO category, multiple-testing corrected or uncorrected p-value, gene symbol type, and genome of interest.

Users can also generate GO term histograms with the top 10 GO terms and to choose whether to display their GO identifiers (GO:IDs) or their corresponding terms.

![BingleSeq Bulk RNA-Seq bulkGO hist](/figures/bulk_GOhist.jpg)

Moreover, users can obtain further information about a given GO term by querying its GO:ID using the ‘GO.db’ package (Carlson et al., 2019). 

![BingleSeq Bulk RNA-Seq bulkGO query](/figures/bulk_GOquery.PNG)

##### As of v3.6, BingleSeq implements the following GOSeq Model organisms:
Homo sapiens, Mus musculus, Danio rerio, Drosophila melanogaster, E. coli K12.
Please do not hesitate to contact us or open a GitHub issue if you want us to include additional model organism.
However, it should be noted that GOSeq requires gene lengths for its usual functions and not all model organisms' genomes have such available.
Refer to the GOSeq manual for further information as well as a way to obtain 

![BingleSeq Bulk RNA-Seq bulkGO Tab](/figures/bulk_GOtab.PNG)

Once the GOseq pipeline is run and completed, a table with results is returned.

![BingleSeq Bulk RNA-Seq bulkGO Results](/figures/bulk_GOresults.PNG)

##### 5b.	Footprint Analysis
BingleSeq also enables the use of footprint analysis tools DoRothEA and PROGENy
(Schubert et al., 2018; Garcia-Alonso, et al., 2019; Holland et al., 2019,
Holland et al., 2020). These tools are used infer the activity of TFs and
pathways, respectively. Footprint-based strategies such as the aforementioned
packages infer TF/pathway activity from the expression of molecules considered
to be downstream of a given TF or pathway (in the case of these tools).

DoRothEA is a gene set resource containing signed TF-target interactions that 
can be coupled with different statistical methods to estimate TF activity.
In BingleSeq, DoRothEA is coupled to the statistical method VIPER
(Alvarez et al., 2016).
PROGENy is based on downstream gene signatures observed to be consistently
deregulated in pertrubation experiments. PROGENy estimates the activity of 14
signalling pathways from gene expression using a linear model.
For more information, please refer to the cited publications.
DoRothEA and PROGENy are available for mouse and human data.

![BingleSeq Bulk RNA-Seq Footprint Analysis Results](/figures/footprint.PNG)
A) TF activities (using DE method statistic - e.g. t-value) and B) TF activities per sample (using normalized gene counts) estimated with DoRothEA and viper.
Pathway activity estimation with PROGENy showing C) the Normalized Enrichment Scores (NES) for each pathway and D) PROGENy pathway scores per sample.

*Note that these results were generated using real data taken from McFarlane et al., 2019.*  

#### 6.	DE Package Comparison
BingleSeq supplies users with an option to assess the agreement between the different DE analysis packages. This is done using a Venn diagram which represents the overlap of DE analysis results obtained using DESeq2, edgeR, and limma on the same dataset. Moreover, users can download the genes from the different intersects of the Venn Diagram. 

![BingleSeq Bulk RNA-Seq compVenn Data](/figures/bulk_compVenn.PNG)

Furthermore, the DE results from the same packages are used to generate a **Rank-based consensus**. The Rank-based consensus is displayed as a table, alongside adjusted p-values and the ranks for each gene as calculated by the different packages.

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
Following normalization, the ‘Clustering’ tab is generated which implements functionality for scaling of the data, dimensionality reduction with PCA, PC selection, and unsupervised clustering. The former three are done simultaneously and Seurat’s ‘PCElbowPlot’ is used to generate and return an elbow plot. The returned elbow plot serves as a heuristic method for determining the true dimensionality of the dataset (i.e. PC Selection). Selecting which PCs to include in Seurat and monocle clustering methods is an essential step as it enables a large portion of technical noise to be excluded.

![BingleSeq Bulk RNA-Seq sc clustEblow](/figures/sc_clustElbow.PNG)
  
  
In addition to the Elbow plot, BingleSeq implements Seurat's PC heatmaps option - to be used as a complementary tool to the elbow plot.
  
![BingleSeq Bulk RNA-Seq sc clustHeat](/figures/sc_clustHeat.PNG)

*A) represents the 1st PC Heatmap with the top 10 most variable Genes and it is highly likely to represent part of the true dimensionality of the dataset. In contrast, B) represents the 15th PC Heatmap which likely represent mainly noise and not true signal.*



Once the count data is scaled and linear dimensionality reduction performed, users can proceed to unsupervised clustering with Seurat, SC3, and monocle. 

When using Seurat for unsupervised clustering, users must specify the number of PCs to be included in the analysis as well as the value of its ‘Resolution’ parameter. The latter parameter is used to set the ‘granularity’ of the clustering and as such it controls the number of clusters. The authors suggest that the optimal Resolution for datasets with ~3000 cells is 0.6-1.2 and it is typically higher for larger datasets. Users can also choose from Seurat’s inbuilt algorithms including Louvain and SLM algorithms.

![BingleSeq Bulk RNA-Seq sc clustSeurat](/figures/sc_clustSeurat.PNG)

*tSNE plot produced using 0.5 as granularity parameter and the first 10 PCs.*
  
  
When clustering with monocle, users are requested to specify the number of PCs to be included in the analysis. Also, if required users can further minimize noise by filtering the gene counts according to the minimum expression level via the ‘Lower Detection Parameter’. Users can also pick from monocle’s inbuilt algorithms, which include Density Peak and Louvain algorithms. Furthermore, Monocle enables the number of clusters to be explicitly specified as well as to be estimated.
  
  
![BingleSeq Bulk RNA-Seq sc clustMono](/figures/sc_clustMonocle.PNG)

*tSNE plot produced by explicitly setting the number of clusters to 9 and using the first 10 PCs (without any additional filtering).*


Unsupervised clustering with SC3 in BingleSeq utilizes SC3's k-means-based clusering approach. Users must specify the number of random initial centroid selections (sets). A larger number of initial centroid configurations (nStart) is likely to produce a better clustering result, but has a high toll on computational time. By default, this parameter is set to 1000 when working with less than 2000 cells and to 50 when working with more than 2000 cells, in accordance to the authors recommendantions.
Users can also use SC3’s inbuilt filtering options to further reduce noise by filtering out genes below and above certain dropout (zero value) percentage thresholds. Similarly to monocle, the number of clusters can be supplied by user or estimated with SC3.
It is worth noting that the k-means approach of SC3 is likely too computationally demanding when working with large datasets (e.g. when N=2000, computational time is ~20 mins), hence future updates of BingleSeq are likely to also implement SC3's SVM-hybrid approach as an alternative solution. 

![BingleSeq Bulk RNA-Seq sc clustSC3](/figures/sc_clustSC3.PNG)

*tSNE plot produced by explicitly setting the number of clusters to 9 and 50 random initial centroid sets (without any additional filtering).*
  
  
*Each tSNE plot in BingleSeq is generated using the package with which clustering was performed.
Also, when performing clustering with SC3 or monocle, the data used to create the required objects to run these pipelines is the same data that was previously filtered, normalized, and scaled using Seurat’s pipeline.*


#### 5.	Differential Expression
Following clustering, DE analysis can be conducted using Seurat’s inbuilt functionality to identify marker genes. Users can perform marker gene identification using the following inbuilt DE testing methods: Student’s T test, Wilcoxon Rank Sum test, and Logistic regression. Additionally, DE analysis can also be performed with DESEq2 and MAST packages (Love, Huber, and Anders, 2014; Finak et al., 2015). Prior to running DE analysis, users are prompted to enter the following filtering parameters: genes expressed in a minimum fraction of cells, fold-change, and adjusted p-value.

![BingleSeq Bulk RNA-Seq sc deTab](/figures/sc_deTab.PNG)

Furthermore, by using Seurat’s inbuilt visualization options, BingleSeq provides tools for the exploration of DE results. These tools include cluster heatmap with user-specified gene number as well as exploration of specific genes via Violin, Feature, and Ridge plots.

![BingleSeq Bulk RNA-Seq sc deFigs](/figures/sc_deFigs.PNG)

*A) Heatmap showing the top 10 genes for each cluster in the 2700 PBMCs dataset, while Violin B), Feature C), and Ridge D) plots are shown for MS4A1 gene – a biomarker of B lymphocytes.*


#### 6.	Functional Annotation
The scRNA-Seq pipeline of BingleSeq incorporates functional annotation in an analogous manner to its Bulk RNA-Seq counterpart. The only difference is that the subsets of DEGs to be used in the GOseq pipeline can be filtered according to the cluster they belong to; thus, allowing users to assess each cluster independently. Accordingly, it also implements both PROGENy and Dorothea-Viper which were implemented to estimate Pathway and TF activity per cell cluster, respectively.
  
  
#### 7.	DE Method Comparison
The scRNA-Seq part also implements a ‘DE Method Comparison’ tab analogous to the ‘DE Package Comparison’ tab in Bulk RNA-Seq. The only difference is that scRNA-Seq Overlap functionality enables filtering according to the same parameters used in marker gene identification. Furthermore, rather than comparing the different packages, it compares the DE Methods implemented within Seurat. These include: DE testing with MAST, Wilcoxon Rank Sum Test, and Student’s T test.

*Also, note that Rank-based consensus is yet to be implemented for the scRNA-Seq pipeline.*

## Test Data
As of v0.3.2 BingleSeq features test data for both Bulk- and scRNA-Seq.
Bulk data - contrast between HSV-1 infected control and interferon B treatment (taken from McFarlane et al., 2019) 
Single cell data - Cell Ranger 10x Genomics public dataset looking at filtered data of 2700 PBMCs
[10x Genomics link](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_unsorted_3k)



## Built With

* [shiny](https://shiny.rstudio.com/)
* [RStudio](https://www.rstudio.com/)


## References
Alvarez M.J., Shen Y., Giorgi F.M., Lachmann A., Ding B.B., Ye B.H., Califano A. (2016). “Functional characterization of somatic mutations in cancer using network-based inference of protein activity.” Nature genetics, 48(8), 838–47.

Alvarez MJ, Shen Y, Giorgi FM, Lachmann A, Ding BB, Ye BH, Califano A (2016). “Functional characterization of somatic mutations in cancer using network-based inference of protein activity.” Nature genetics, 48(8), 838–47.

Bowden J, Ross J, Oytam Y (2019). HarmanData: Data for the Harman package. R package version 1.12.0, http://www.bioinformatics.csiro.au/harman/.

Carlson, M. Falcon, S., Pages, H., Li, N., 2019. GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.8.2.

Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A.K., Slichter, C.K., Miller, H.W., McElrath, M.J., Prlic, M. and Linsley, P.S., 2015. MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome biology, 16(278) doi: 10.1186/s13059-015-0844-5.

Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. "Benchmark and integration of resources for the estimation of human transcription factor activities." Genome Research. 2019. DOI: 

Holland CH, Szalai B, Saez-Rodriguez J. "Transfer of regulatory knowledge from human to mouse for functional genomics analysis." Biochimica et Biophysica Acta (BBA) - Gene Regulatory Mechanisms. 2019. DOI: 10.1016/j.bbagrm.2019.194431.

Holland CH, Tanevski J, Perales-Patón J, Gleixner J, Kumar MP, Mereu E, Joughin BA, Stegle O, Lauffenburger DA, Heyn H, Szalai B, Saez-Rodriguez, J. "Robustness and applicability of transcription factor and pathway analysis tools on single-cell RNA-seq data." Genome Biology. 2020. DOI: 10.1186/s13059-020-1949-z.

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

Schubert M, Klinger B, Klünemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Blüthgen N, Saez-Rodriguez J (2018). “Perturbation-response genes reveal signaling footprints in cancer gene expression.” Nature communications, 9(20).

Soneson, C., 2014. compcodeR—an R package for benchmarking differential expression methods for RNA-seq data, Bioinformatics, 30(17), pp.2517–2518.

Trapnell, C., Cacchiarelli, D., Grimsby, J., Pokharel, P., Li, S., Morse, M., Lennon, N.J., Livak, K.J., Mikkelsen, T.S. and Rinn, J.L., 2014. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nature biotechnology, 32(4), pp. 381–386.

Young, M.D., Wakefield, M.J., Smyth, G.K. and Oshlack, A., 2010. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome biology, 11(R14). doi: https://doi.org/10.1186/gb-2010-11-2-r14.



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



## R Session Info:
```r
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.3 (2020-10-10)
 os       Ubuntu 20.04.2 LTS          
 system   x86_64, linux-gnu           
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Berlin               
 date     2021-07-11                  

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 ! package              * version  date       lib source        
   abind                  1.4-5    2016-07-21 [1] CRAN (R 4.0.2)
   annotate               1.68.0   2020-10-27 [1] Bioconductor  
   AnnotationDbi        * 1.52.0   2020-10-27 [1] Bioconductor  
   askpass                1.1      2019-01-13 [1] CRAN (R 4.0.2)
   assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.0.2)
   backports              1.2.1    2020-12-09 [1] CRAN (R 4.0.2)
   bcellViper             1.26.0   2020-10-29 [1] Bioconductor  
   beachmat               2.6.4    2020-12-20 [1] Bioconductor  
   beeswarm               0.2.3    2016-04-25 [1] CRAN (R 4.0.2)
   BiasedUrn            * 1.07     2015-12-28 [1] CRAN (R 4.0.2)
 P BingleSeq            * 0.3.6    2021-02-06 [?] local         
   Biobase              * 2.50.0   2020-10-27 [1] Bioconductor  
   BiocFileCache          1.14.0   2020-10-27 [1] Bioconductor  
   BiocGenerics         * 0.36.0   2020-10-27 [1] Bioconductor  
   BiocNeighbors          1.8.2    2020-12-07 [1] Bioconductor  
   BiocParallel         * 1.24.1   2020-11-06 [1] Bioconductor  
   BiocSingular           1.6.0    2020-10-27 [1] Bioconductor  
   biomaRt                2.46.3   2021-02-09 [1] Bioconductor  
   Biostrings             2.58.0   2020-10-27 [1] Bioconductor  
   bit                    4.0.4    2020-08-04 [1] CRAN (R 4.0.2)
   bit64                  4.0.5    2020-08-30 [1] CRAN (R 4.0.2)
   bitops                 1.0-6    2013-08-17 [1] CRAN (R 4.0.2)
   blob                   1.2.1    2020-01-20 [1] CRAN (R 4.0.2)
   bluster                1.0.0    2020-10-27 [1] Bioconductor  
   broom                  0.7.8    2021-06-24 [1] CRAN (R 4.0.3)
   bslib                  0.2.4    2021-01-25 [1] CRAN (R 4.0.3)
   cachem                 1.0.4    2021-02-13 [1] CRAN (R 4.0.3)
   Cairo                  1.5-12.2 2020-07-07 [1] CRAN (R 4.0.3)
   callr                  3.5.1    2020-10-13 [1] CRAN (R 4.0.2)
   car                    3.0-10   2020-09-29 [1] CRAN (R 4.0.2)
   carData                3.0-4    2020-05-22 [1] CRAN (R 4.0.2)
   cellranger             1.1.0    2016-07-27 [1] CRAN (R 4.0.2)
   class                  7.3-17   2020-04-26 [4] CRAN (R 4.0.0)
   cli                    2.5.0    2021-04-26 [1] CRAN (R 4.0.3)
   cluster                2.1.0    2019-06-19 [4] CRAN (R 4.0.0)
   codetools              0.2-16   2018-12-24 [4] CRAN (R 4.0.0)
   colorspace             2.0-0    2020-11-11 [1] CRAN (R 4.0.2)
   combinat               0.0-8    2012-10-29 [1] CRAN (R 4.0.2)
   cowplot                1.1.1    2020-12-30 [1] CRAN (R 4.0.3)
   crayon                 1.4.1    2021-02-08 [1] CRAN (R 4.0.3)
   crosstalk              1.1.1    2021-01-12 [1] CRAN (R 4.0.3)
   curl                   4.3      2019-12-02 [1] CRAN (R 4.0.2)
   data.table             1.14.0   2021-02-21 [1] CRAN (R 4.0.3)
   DBI                    1.1.1    2021-01-15 [1] CRAN (R 4.0.3)
   dbplyr                 2.1.0    2021-02-03 [1] CRAN (R 4.0.3)
   DDRTree              * 0.1.5    2017-04-30 [1] CRAN (R 4.0.2)
   DelayedArray           0.16.1   2021-01-22 [1] Bioconductor  
   DelayedMatrixStats     1.12.3   2021-02-03 [1] Bioconductor  
   deldir                 0.2-10   2021-02-16 [1] CRAN (R 4.0.3)
   densityClust           0.3      2017-10-24 [1] CRAN (R 4.0.2)
   DEoptimR               1.0-8    2016-11-19 [1] CRAN (R 4.0.2)
   desc                   1.2.0    2018-05-01 [1] CRAN (R 4.0.2)
   DESeq2               * 1.30.1   2021-02-19 [1] Bioconductor  
   devtools               2.3.2    2020-09-18 [1] CRAN (R 4.0.3)
   digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.2)
   docopt                 0.7.1    2020-06-24 [1] CRAN (R 4.0.2)
   doParallel             1.0.16   2020-10-16 [1] CRAN (R 4.0.2)
   doRNG                  1.8.2    2020-01-27 [1] CRAN (R 4.0.2)
   dorothea             * 1.2.1    2021-02-11 [1] Bioconductor  
   dplyr                * 1.0.7    2021-06-18 [1] CRAN (R 4.0.3)
   dqrng                  0.2.1    2019-05-17 [1] CRAN (R 4.0.2)
   DT                   * 0.17     2021-01-06 [1] CRAN (R 4.0.3)
   e1071                  1.7-4    2020-10-14 [1] CRAN (R 4.0.2)
   edgeR                * 3.32.1   2021-01-14 [1] Bioconductor  
   ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.0.3)
   factoextra           * 1.0.7    2020-04-01 [1] CRAN (R 4.0.2)
   fansi                  0.4.2    2021-01-15 [1] CRAN (R 4.0.3)
   farver                 2.0.3    2020-01-16 [1] CRAN (R 4.0.2)
   fastcluster          * 1.1.25   2018-06-07 [1] CRAN (R 4.0.2)
   fastICA                1.2-2    2019-07-08 [1] CRAN (R 4.0.2)
   fastmap                1.1.0    2021-01-25 [1] CRAN (R 4.0.3)
   fitdistrplus           1.1-3    2020-12-05 [1] CRAN (R 4.0.2)
   FNN                    1.1.3    2019-02-15 [1] CRAN (R 4.0.2)
   forcats              * 0.5.1    2021-01-27 [1] CRAN (R 4.0.3)
   foreach                1.5.1    2020-10-15 [1] CRAN (R 4.0.2)
   foreign                0.8-80   2020-05-24 [1] CRAN (R 4.0.2)
   formatR                1.7      2019-06-11 [1] CRAN (R 4.0.2)
   fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.2)
   futile.logger        * 1.4.3    2016-07-10 [1] CRAN (R 4.0.2)
   futile.options         1.0.1    2018-04-20 [1] CRAN (R 4.0.2)
   future                 1.21.0   2020-12-10 [1] CRAN (R 4.0.3)
   future.apply           1.7.0    2021-01-04 [1] CRAN (R 4.0.3)
   genefilter           * 1.72.1   2021-01-21 [1] Bioconductor  
   geneLenDataBase      * 1.26.0   2020-10-29 [1] Bioconductor  
   geneplotter            1.68.0   2020-10-27 [1] Bioconductor  
   generics               0.1.0    2020-10-31 [1] CRAN (R 4.0.2)
   GenomeInfoDb         * 1.26.2   2020-12-08 [1] Bioconductor  
   GenomeInfoDbData       1.2.4    2020-12-14 [1] Bioconductor  
   GenomicAlignments      1.26.0   2020-10-27 [1] Bioconductor  
   GenomicFeatures        1.42.1   2020-11-12 [1] Bioconductor  
   GenomicRanges        * 1.42.0   2020-10-27 [1] Bioconductor  
   ggbeeswarm             0.6.0    2017-08-07 [1] CRAN (R 4.0.2)
   ggplot2              * 3.3.3    2020-12-30 [1] CRAN (R 4.0.3)
   ggpubr                 0.4.0    2020-06-27 [1] CRAN (R 4.0.2)
   ggrepel              * 0.9.1    2021-01-15 [1] CRAN (R 4.0.3)
   ggridges               0.5.3    2021-01-08 [1] CRAN (R 4.0.3)
   ggsignif               0.6.0    2019-08-08 [1] CRAN (R 4.0.2)
   globals                0.14.0   2020-11-22 [1] CRAN (R 4.0.2)
   glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.2)
   GO.db                * 3.12.1   2021-02-21 [1] Bioconductor  
   goftest                1.2-2    2019-12-02 [1] CRAN (R 4.0.2)
   goseq                * 1.42.0   2020-10-27 [1] Bioconductor  
   gridExtra            * 2.3      2017-09-09 [1] CRAN (R 4.0.2)
   gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.0.2)
   Harman               * 1.18.0   2020-10-27 [1] Bioconductor  
   haven                  2.3.1    2020-06-01 [1] CRAN (R 4.0.2)
   hms                    1.0.0    2021-01-13 [1] CRAN (R 4.0.3)
   HSMMSingleCell         1.10.0   2020-10-29 [1] Bioconductor  
   htmltools              0.5.1.1  2021-01-22 [1] CRAN (R 4.0.3)
   htmlwidgets            1.5.3    2020-12-10 [1] CRAN (R 4.0.2)
   httpuv                 1.5.5    2021-01-13 [1] CRAN (R 4.0.3)
   httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.2)
   ica                    1.0-2    2018-05-24 [1] CRAN (R 4.0.2)
   igraph                 1.2.6    2020-10-06 [1] CRAN (R 4.0.2)
   IRanges              * 2.24.1   2020-12-12 [1] Bioconductor  
   irlba                * 2.3.3    2019-02-05 [1] CRAN (R 4.0.2)
   iterators              1.0.13   2020-10-15 [1] CRAN (R 4.0.2)
   jquerylib              0.1.3    2020-12-17 [1] CRAN (R 4.0.3)
   jsonlite               1.7.2    2020-12-09 [1] CRAN (R 4.0.2)
   kernlab                0.9-29   2019-11-12 [1] CRAN (R 4.0.2)
   KernSmooth             2.23-17  2020-04-26 [4] CRAN (R 4.0.0)
   labeling               0.4.2    2020-10-20 [1] CRAN (R 4.0.2)
   lambda.r               1.2.4    2019-09-18 [1] CRAN (R 4.0.2)
   later                  1.1.0.1  2020-06-05 [1] CRAN (R 4.0.2)
   lattice                0.20-41  2020-04-02 [4] CRAN (R 4.0.0)
   lazyeval               0.2.2    2019-03-15 [1] CRAN (R 4.0.2)
   leiden                 0.3.7    2021-01-26 [1] CRAN (R 4.0.3)
   lifecycle              1.0.0    2021-02-15 [1] CRAN (R 4.0.3)
   limma                * 3.46.0   2020-10-27 [1] Bioconductor  
   listenv                0.8.0    2019-12-05 [1] CRAN (R 4.0.2)
   lmtest                 0.9-38   2020-09-09 [1] CRAN (R 4.0.2)
   locfit                 1.5-9.4  2020-03-25 [1] CRAN (R 4.0.2)
   lubridate              1.7.9.2  2020-11-13 [1] CRAN (R 4.0.2)
   magrittr               2.0.1    2020-11-17 [1] CRAN (R 4.0.2)
   MASS                   7.3-53   2020-09-09 [4] CRAN (R 4.0.2)
   MAST                 * 1.16.0   2020-10-27 [1] Bioconductor  
   Matrix               * 1.3-4    2021-06-01 [1] CRAN (R 4.0.3)
   MatrixGenerics       * 1.2.1    2021-01-30 [1] Bioconductor  
   matrixStats          * 0.58.0   2021-01-29 [1] CRAN (R 4.0.3)
   memoise                2.0.0    2021-01-26 [1] CRAN (R 4.0.3)
   mgcv                 * 1.8-33   2020-08-27 [4] CRAN (R 4.0.2)
   mime                   0.10     2021-02-13 [1] CRAN (R 4.0.3)
   miniUI                 0.1.1.1  2018-05-18 [1] CRAN (R 4.0.2)
   mixtools               1.2.0    2020-02-07 [1] CRAN (R 4.0.2)
   modelr                 0.1.8    2020-05-19 [1] CRAN (R 4.0.2)
   monocle              * 2.18.0   2020-10-27 [1] Bioconductor  
   munsell                0.5.0    2018-06-12 [1] CRAN (R 4.0.2)
   mvtnorm                1.1-1    2020-06-09 [1] CRAN (R 4.0.2)
   nlme                 * 3.1-149  2020-08-23 [4] CRAN (R 4.0.2)
   openssl                1.4.3    2020-09-18 [1] CRAN (R 4.0.2)
   openxlsx               4.2.3    2020-10-27 [1] CRAN (R 4.0.2)
   org.Hs.eg.db         * 3.12.0   2021-02-21 [1] Bioconductor  
   org.Mm.eg.db         * 3.12.0   2021-02-21 [1] Bioconductor  
   parallelly             1.24.0   2021-03-14 [1] CRAN (R 4.0.3)
   patchwork              1.1.1    2020-12-17 [1] CRAN (R 4.0.3)
   pbapply                1.4-3    2020-08-18 [1] CRAN (R 4.0.2)
   pcaPP                  1.9-73   2018-01-14 [1] CRAN (R 4.0.2)
   pheatmap             * 1.0.12   2019-01-04 [1] CRAN (R 4.0.2)
   pillar                 1.6.1    2021-05-16 [1] CRAN (R 4.0.3)
   pkgbuild               1.2.0    2020-12-15 [1] CRAN (R 4.0.3)
   pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.2)
   pkgload                1.1.0    2020-05-29 [1] CRAN (R 4.0.2)
   plotly               * 4.9.3    2021-01-10 [1] CRAN (R 4.0.3)
   plyr                   1.8.6    2020-03-03 [1] CRAN (R 4.0.2)
   png                    0.1-7    2013-12-03 [1] CRAN (R 4.0.2)
   polyclip               1.10-0   2019-03-14 [1] CRAN (R 4.0.2)
   prettyunits            1.1.1    2020-01-24 [1] CRAN (R 4.0.2)
   processx               3.4.5    2020-11-30 [1] CRAN (R 4.0.2)
   progeny              * 1.12.0   2020-10-27 [1] Bioconductor  
   progress               1.2.2    2019-05-16 [1] CRAN (R 4.0.2)
   promises               1.2.0.1  2021-02-11 [1] CRAN (R 4.0.3)
   ps                     1.5.0    2020-12-05 [1] CRAN (R 4.0.2)
   purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.2)
   qlcMatrix              0.9.7    2018-04-20 [1] CRAN (R 4.0.2)
   R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.2)
   RANN                   2.6.1    2019-01-08 [1] CRAN (R 4.0.2)
   rappdirs               0.3.3    2021-01-31 [1] CRAN (R 4.0.3)
   RColorBrewer           1.1-2    2014-12-07 [1] CRAN (R 4.0.2)
   Rcpp                   1.0.6    2021-01-15 [1] CRAN (R 4.0.3)
   RcppAnnoy              0.0.18   2020-12-15 [1] CRAN (R 4.0.3)
   RCurl                  1.98-1.3 2021-03-16 [1] CRAN (R 4.0.3)
   readr                * 1.4.0    2020-10-05 [1] CRAN (R 4.0.2)
   readxl                 1.3.1    2019-03-13 [1] CRAN (R 4.0.2)
   remotes                2.2.0    2020-07-21 [1] CRAN (R 4.0.2)
   reprex                 1.0.0    2021-01-27 [1] CRAN (R 4.0.3)
   reshape2             * 1.4.4    2020-04-09 [1] CRAN (R 4.0.2)
   reticulate             1.18     2020-10-25 [1] CRAN (R 4.0.2)
   rio                    0.5.16   2018-11-26 [1] CRAN (R 4.0.2)
   rlang                  0.4.10   2020-12-30 [1] CRAN (R 4.0.3)
   rngtools               1.5      2020-01-23 [1] CRAN (R 4.0.2)
   robustbase             0.93-7   2021-01-04 [1] CRAN (R 4.0.3)
   ROCR                   1.0-11   2020-05-02 [1] CRAN (R 4.0.2)
   rpart                  4.1-15   2019-04-12 [4] CRAN (R 4.0.0)
   rprojroot              2.0.2    2020-11-15 [1] CRAN (R 4.0.2)
   rrcov                  1.5-5    2020-08-03 [1] CRAN (R 4.0.2)
   Rsamtools              2.6.0    2020-10-27 [1] Bioconductor  
   RSQLite                2.2.3    2021-01-24 [1] CRAN (R 4.0.3)
   rstatix                0.7.0    2021-02-13 [1] CRAN (R 4.0.3)
   rstudioapi             0.13     2020-11-12 [1] CRAN (R 4.0.2)
   rsvd                   1.0.3    2020-02-17 [1] CRAN (R 4.0.2)
   rtracklayer            1.50.0   2020-10-27 [1] Bioconductor  
   Rtsne                  0.15     2018-11-10 [1] CRAN (R 4.0.2)
   rvest                  0.3.6    2020-07-25 [1] CRAN (R 4.0.2)
   S4Vectors            * 0.28.1   2020-12-09 [1] Bioconductor  
   sass                   0.3.1    2021-01-24 [1] CRAN (R 4.0.3)
   SC3                  * 1.18.0   2020-10-27 [1] Bioconductor  
   scales                 1.1.1    2020-05-11 [1] CRAN (R 4.0.2)
   scater                 1.18.5   2021-02-16 [1] Bioconductor  
   scattermore            0.7      2020-11-24 [1] CRAN (R 4.0.3)
   scran                * 1.18.5   2021-02-04 [1] Bioconductor  
   sctransform            0.3.2    2020-12-16 [1] CRAN (R 4.0.3)
   scuttle                1.0.4    2020-12-17 [1] Bioconductor  
   segmented              1.3-2    2021-02-09 [1] CRAN (R 4.0.3)
   sessioninfo            1.1.1    2018-11-05 [1] CRAN (R 4.0.2)
   Seurat               * 4.0.3    2021-06-10 [1] CRAN (R 4.0.3)
   SeuratObject         * 4.0.2    2021-06-09 [1] CRAN (R 4.0.3)
   shiny                * 1.6.0    2021-01-25 [1] CRAN (R 4.0.3)
   shinyFiles           * 0.9.0    2020-11-09 [1] CRAN (R 4.0.2)
   shinyjs              * 2.0.0    2020-09-09 [1] CRAN (R 4.0.2)
   shinyWidgets         * 0.5.7    2021-02-03 [1] CRAN (R 4.0.3)
   SingleCellExperiment * 1.12.0   2020-10-27 [1] Bioconductor  
   slam                   0.1-48   2020-12-03 [1] CRAN (R 4.0.2)
   sparseMatrixStats      1.2.1    2021-02-02 [1] Bioconductor  
   sparsesvd              0.2      2019-07-15 [1] CRAN (R 4.0.2)
   spatstat.core          2.0-0    2021-03-23 [1] CRAN (R 4.0.3)
   spatstat.data          2.1-0    2021-03-21 [1] CRAN (R 4.0.3)
   spatstat.geom          2.0-1    2021-03-22 [1] CRAN (R 4.0.3)
   spatstat.sparse        2.0-0    2021-03-16 [1] CRAN (R 4.0.3)
   spatstat.utils         2.1-0    2021-03-15 [1] CRAN (R 4.0.3)
   statmod                1.4.35   2020-10-19 [1] CRAN (R 4.0.2)
   stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.2)
   stringr              * 1.4.0    2019-02-10 [1] CRAN (R 4.0.2)
   SummarizedExperiment * 1.20.0   2020-10-27 [1] Bioconductor  
   survival               3.2-7    2020-09-28 [1] CRAN (R 4.0.2)
   sva                  * 3.38.0   2020-10-27 [1] Bioconductor  
   tensor                 1.5      2012-05-05 [1] CRAN (R 4.0.2)
   testthat               3.0.2    2021-02-14 [1] CRAN (R 4.0.3)
   tibble               * 3.1.2    2021-05-16 [1] CRAN (R 4.0.3)
   tidyr                * 1.1.3    2021-03-03 [1] CRAN (R 4.0.3)
   tidyselect             1.1.0    2020-05-11 [1] CRAN (R 4.0.2)
   tidyverse            * 1.3.0    2019-11-21 [1] CRAN (R 4.0.3)
   tinytex                0.29     2021-01-21 [1] CRAN (R 4.0.3)
   usethis                2.0.1    2021-02-10 [1] CRAN (R 4.0.3)
   utf8                   1.1.4    2018-05-24 [1] CRAN (R 4.0.2)
   uwot                   0.1.10   2020-12-15 [1] CRAN (R 4.0.3)
   vctrs                  0.3.8    2021-04-29 [1] CRAN (R 4.0.3)
   VennDiagram          * 1.6.20   2018-03-28 [1] CRAN (R 4.0.2)
   VGAM                 * 1.1-5    2021-01-14 [1] CRAN (R 4.0.3)
   viper                * 1.24.0   2020-10-27 [1] Bioconductor  
   vipor                  0.4.5    2017-03-22 [1] CRAN (R 4.0.2)
   viridis                0.5.1    2018-03-29 [1] CRAN (R 4.0.2)
   viridisLite            0.3.0    2018-02-01 [1] CRAN (R 4.0.2)
   waiter               * 0.2.0    2021-01-14 [1] CRAN (R 4.0.3)
   withr                  2.4.1    2021-01-26 [1] CRAN (R 4.0.3)
   WriteXLS               6.1.0    2020-11-23 [1] CRAN (R 4.0.2)
   xfun                   0.21     2021-02-10 [1] CRAN (R 4.0.3)
   XML                    3.99-0.5 2020-07-23 [1] CRAN (R 4.0.2)
   xml2                   1.3.2    2020-04-23 [1] CRAN (R 4.0.2)
   xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.0.2)
   XVector                0.30.0   2020-10-27 [1] Bioconductor  
   yaml                   2.2.1    2020-02-01 [1] CRAN (R 4.0.2)
   zip                    2.1.1    2020-08-27 [1] CRAN (R 4.0.2)
   zlibbioc               1.36.0   2020-10-27 [1] Bioconductor  
   zoo                    1.8-8    2020-05-02 [1] CRAN (R 4.0.2)

[1] /home/dbdimitrov/R/x86_64-pc-linux-gnu-library/4.0
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```