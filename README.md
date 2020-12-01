# CCM RNA-Seq lecture files

These are the necessary files for RNA-Seq lecture demonstration. In addition to these files you will need [R](https://www.r-project.org/) and [Rstudio](https://rstudio.com/) installed. You can find instructions in their respective pages. 

After installation you will need to donwload and install some required packages. Please launch Rstudio and run 

```{r}
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "GenomicFeatures", "remotes", "patchwork"))

remotes::install_github("js229/Vennerable")
```
