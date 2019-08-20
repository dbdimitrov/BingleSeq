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

To Run BingleSeq, you must have R>= 3.5.3


### Installing

To install BingleSeq on your machine you must type the following R code:

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


## Built With

* [shiny](https://shiny.rstudio.com/) - The R framework used
* [RStudio](https://www.rstudio.com/) - Dependency Management



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
