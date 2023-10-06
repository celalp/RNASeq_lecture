# CCM RNA-Seq lecture files

These are the necessary files and code for CCM 2023 RNA-Seq lecture demonstration. In addition to these files you will also need[R](https://www.r-project.org/) and [Rstudio](https://rstudio.com/) installed. You can find system specific instructions on the respective websites. 

We will be using a bunch of 3rd party packages, to speed things up you can install these packages using the code snippet below:

```R
install.packages("BiocManager")

BiocManager::install("ggplot2", "dplyr", "pheatmap", "DESeq2", "AnnotationHub", 
                     "clusterProfiler", "tximport")
```

The data for this excercise if from the paper:
Natural history of a mouse model of X-linked myotubular myopathy ([DOI](https://doi.org/10.1242/dmm.049342))

