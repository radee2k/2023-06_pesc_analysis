06-23_pesc_ev_analysis
================
2023-06-22

<body style="background-color:FloralWhite;">

This script identifies lineages of EVs among cells from pleural effusion
samples.

------------------------------------------------------------------------

# Prepare the workign enviornment:

-   install and load packages
-   set %notin% and %notlike%
    -   set ggplot’s theme
-   set the working directory
-   set a plan for multithreading

``` r
# install.packages("Seurat")
# install.packages("remotes")
# BiocManager::install(version = '3.16')
# BiocManager::install("glmGamPoi")
# remotes::install_github("stephenturner/annotables")
# install.packages("glmGamPoi")
# BiocManager::install("DESeq2")
# install.packages("data.filt")

library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(ggplot2)
library(scales) # to better customize Seurat's plots
library(ggpubr)
library(ggrepel)
library(patchwork)
library(ggbeeswarm) # for better positioning of arrows on log2 change plots - position_jitter doesn't allow to mark points
library(future)
library(annotables) # for turning Ensembl ID to symbol
library(sctransform) # for normalization  
library(glmGamPoi) # for SCTransform
# library(svglite) # for vectorized, lightweight plotting
library(systemfonts) # to set the font for svg outputs
library(DESeq2)
library(MAST)
# library(data.filt) # <---------------- GET IT

"%notin%" <- Negate("%in%")
"%notlike%" <- Negate("%like%")


# set the theme for plotting (Seurat uses ggplot's themes)
theme_set(new = theme_classic())
theme_update(
  axis.text.x = element_text(vjust = 0.5),
  strip.background = element_rect(fill = '#FFFFFF'),
  plot.title = element_text(hjust = 0.5, size = 25),
  plot.subtitle = element_text(size = 20, hjust = .5),
  axis.title = element_text(size = 23),
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.key.size = unit(2, 'line'),
  legend.title = element_text(size = 20, hjust = .5, vjust = .5)
  # text = element_text(family= "mono")
)

# That's not necessary (rmarkdown sets its directory as the one the .Rmd file is in.)
wd <- "/disk2/user/radgro/projects/2023-06_pesc_analysis"
knitr::opts_knit$set(root.dir = wd)



fonts <- list(
  mono = "Consolas",
  sans = "Lato"
)

# set svglite as a default for all the plots
# knitr::opts_chunk$set(knitr.chunk.dev = 'svglite')
# knitr::opts_chunk$set(dev = 'svglite', system_fonts = fonts)
knitr::opts_chunk$set(dev = 'svglite', dev.args = list(system_fonts = fonts),
                       cache.path = "pesc_analysis_1_cache/gfm/", cache = T) # cache of a github_document doesn't work if the path to the gfm folder is not provided!!!

# knitr::opts_chunk$set(cache.extra = 1) # RESETS CACHE

# plan("multicore", workers = 8) # Not allowed on the server
# plan()
```

## Load and prepare the data and metadata.

``` r
# for EVs
d_dc <- as.sparse(read.csv("data/DC_matrix.txt", sep = "\t"))

# for cells
d_c <- as.sparse(read.csv("data/PEsc_matrix.txt", sep = "\t"))

# Gene names
gt <- fread("data/genes_title.txt", sep = "\t")

# get gene symbols from annotables
gt_s <- setDT(grch38[, c("ensgene", "symbol")])
colnames(gt_s)[1] <- "Geneid"
```

Join gene names with gene symbols and change symbols for Geneids when
the grch38 doesn’t include any particular ID, mult set to ‘first’, to
exclude duplicates/synonyms.

``` r
gt_join <- gt_s[gt, on = .(Geneid), mult = 'first'][symbol %in% NA | symbol == "", symbol := Geneid]
gt_sym <- gt_join[, symbol]
```

``` r
gt_mito <- readLines("code_for_daniel/human_mitochondrial_genes_list.txt")[-1] # Could be used for Daniel's approach
gt_join[Geneid %in% gt_mito] # same genes are mitochondrial here and in Daniel's table - all is fine
```

Insert joined gene names and symbols into datasets.

``` r
rownames(d_c) <- gt_sym
rownames(d_dc) <- gt_sym

rm(list = c("gt", "gt_s", "gt_sym"))
```

# QC and filtering of cells.

``` r
ds_c <- CreateSeuratObject(count = d_c, min.cells = 0, min.features = 1, project = "cells")
```

    ## Warning: Non-unique features (rownames) present in the input matrix, making
    ## unique

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
rm(d_c)
suppressMessages(gc())
```

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  6679317 356.8   12136614  648.2   7478701  399.5
    ## Vcells 27639360 210.9  242108892 1847.2 299738434 2286.9

``` r
ds_c <- PercentageFeatureSet(ds_c, pattern = "^MT-", col.name = "percent_mt")
ds_c <- PercentageFeatureSet(ds_c, "^RP[SL]", col.name = "percent_ribo")
ds_c <- PercentageFeatureSet(ds_c, "^HB[^(P)]", col.name = "percent_hb")
ds_c <- PercentageFeatureSet(ds_c, "PECAM1|PF4", col.name = "percent_plat")
```

## Plotting

### Main QC plots

``` r
VlnPlot(ds_c, features = c('nCount_RNA','nFeature_RNA', 'percent_mt', 'percent_hb', "percent_ribo", "percent_plat"), pt.size = 1.3, ncol = 1) +  NoLegend()
```

![](pesc_analysis_1_files/figure-gfm/Main%20QC%20plots,%20-1.svg)<!-- -->

### UMI-gene count correlation

``` r
FeatureScatter(ds_c, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = T) + scale_x_continuous(labels = scales::scientific) + NoLegend()
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-4-1.svg)<!-- -->

### Hemoglobin content

``` r
c_hemo <- ds_c[[c("orig.ident", "nFeature_RNA", "percent_hb")]]

ggplot(c_hemo) +
  geom_jitter(aes(x = orig.ident, y = nFeature_RNA, size = percent_hb, color = percent_hb)) +
  scale_color_gradient2(midpoint = 30, mid = "orange", low="blue", high="red") +
  theme(axis.text.x = element_text(angle = 50)) +
  geom_hline(yintercept = 500)
```

![](pesc_analysis_1_files/figure-gfm/hemoglobin_content_cells-1.svg)<!-- -->

``` r
rm(c_hemo)
```

------------------------------------------------------------------------

## Filtering

How many cells have more than 800 RNAs or more than 1000 RNAs and genes
found.

``` r
length(WhichCells(ds_c, expression = nCount_RNA > 800))
```

    ## [1] 1765

``` r
length(WhichCells(ds_c, expression = nCount_RNA > 800 & nFeature_RNA > 500 & nCount_RNA > 1000 & nFeature_RNA > 600 & percent_mt < 25 & percent_ribo > 1.5
                & percent_hb < 1))
```

    ## [1] 1292

Filter out cells that have \<=800 features and 500 UMIs (subsetting on
the percentages of different transcript types can be done later).

The percent of ribosomal genes might be low due to low transcriptional
activity of the cells (5%, which was a threshold used before, removes
too many cells).

``` r
ds_cf <- subset(x = ds_c, subset = nCount_RNA > 1000 & nFeature_RNA > 600 & percent_mt < 20 & percent_ribo > 1.5
                & percent_hb < 1)
rm(ds_c)
suppressMessages(gc())
```

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  6907749 369.0   12136614  648.2  12136614  648.2
    ## Vcells 34065420 259.9  193687114 1477.8 299738434 2286.9

### Find which genes contribute to the nCount_RNA the most

``` r
counts_cells <- ds_cf@assays$RNA@counts
counts_cells <- Matrix::t(Matrix::t(counts_cells)/Matrix::colSums(counts_cells)) * 100

most_expr_cells <- order(apply(counts_cells, 1, median), decreasing = T)[20:1]
most_expr_counts_cells <- as.matrix(t(counts_cells[most_expr_cells,]))

rm(list = c("counts_cells", "most_expressed_cells"))
```

    ## Warning in rm(list = c("counts_cells", "most_expressed_cells")): object
    ## 'most_expressed_cells' not found

``` r
par(mar=c(5,10,1,1))
boxplot(most_expr_counts_cells, cex = 1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
```

![](pesc_analysis_1_files/figure-gfm/most_expressed_cells_filtered-1.svg)<!-- -->

### QC plotting of filtered reads

``` r
VlnPlot(ds_cf, features = c('nFeature_RNA','nCount_RNA'), log = T) ## UGLY
```

![](pesc_analysis_1_files/figure-gfm/qc_plots_cells_filtered-1.svg)<!-- -->

``` r
VlnPlot(ds_cf, features = c('percent_mt', "percent_ribo"))
```

![](pesc_analysis_1_files/figure-gfm/qc_plots_cells_filtered-2.svg)<!-- -->

### QC plotting of filtered reads with ggplot

``` r
p1 <- ggplot(as.data.table(ds_cf$nFeature_RNA)) +
  geom_violin(aes(x = "", y = V1)) +
  geom_jitter(aes(x = "", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('Genes') +
  theme(plot.title = element_text(size = 18)) +
  xlab("") +
  ylab("count")

p2 <- ggplot(as.data.table(ds_cf$nCount_RNA)) +
  geom_violin(aes(x = "", y = V1)) +
  geom_jitter(aes(x = "", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('RNAs') +
  theme(plot.title = element_text(size = 18)) +
  xlab("") +
  ylab("")

p3 <- ggplot(as.data.table(ds_cf$percent_mt)) +
  geom_violin(aes(x = "", y = V1)) +
  geom_jitter(aes(x = "", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of\nmitochondrial') +
  theme(plot.title = element_text(size = 18)) +
  xlab("") +
  ylab("percent")

p4 <- ggplot(as.data.table(ds_cf$percent_ribo)) +
  geom_violin(aes(x = "", y = V1)) +
  geom_jitter(aes(x = "", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of\nribosomal') +
  theme(plot.title = element_text(size = 18)) +
  xlab("") +
  ylab("")

p5 <- ggplot(as.data.table(ds_cf$percent_hb)) +
  geom_violin(aes(x = "", y = V1)) +
  geom_jitter(aes(x = "", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of\nhemoglobin') +
  theme(plot.title = element_text(size = 18)) +
  xlab("") +
  ylab("")

p6 <- ggplot(as.data.table(ds_cf$percent_plat)) +
  geom_violin(aes(x = "", y = V1)) +
  geom_jitter(aes(x = "", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of\nplatelet') +
  theme(plot.title = element_text(size = 18)) +
  xlab("") +
  ylab("")


ggarrange(p1, p2, p3, p4, p5, p6, nrow = 1)
```

![](pesc_analysis_1_files/figure-gfm/qc_plots_cells_filtered_ggplot-1.svg)<!-- -->

``` r
# ggarrange(ggarrange(p1, p2), 
          # ggarrange(p3, p4, p5, p6, ncol = 2, nrow = 2, widths = 2),
          # nrow = 2)

#ggsave("qc_violins_filt_p3.png", qcplots, device = 'png', dpi = "retina", width = 21, height = 12, bg = "white")
```

### UMI-gene count correlation

``` r
FeatureScatter(ds_cf, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = T) + scale_x_continuous(labels = scales::scientific) + NoLegend()
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-6-1.svg)<!-- -->

### Hemoglobin content

``` r
cf_hemo <- ds_cf[[c("orig.ident", "nFeature_RNA", "percent_hb")]]

ggplot(cf_hemo) +
  geom_jitter(aes(x = orig.ident, y = nFeature_RNA, size = percent_hb, color = percent_hb)) +
  scale_color_gradient2(midpoint = .5, mid = "orange", low="blue", high="red", ) +
  theme(axis.text.x = element_text(angle = 50))
```

![](pesc_analysis_1_files/figure-gfm/hemoglobin_content_cells_filtered-1.svg)<!-- -->

``` r
rm(cf_hemo)
```

------------------------------------------------------------------------

# QC and filtering of EVs.

``` r
ds_dc <- CreateSeuratObject(count = d_dc, min.cells = 0, min.features = 1, project = "ev")
```

    ## Warning: Non-unique features (rownames) present in the input matrix, making
    ## unique

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
rm(d_dc)
suppressMessages(gc())
```

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  6790216 362.7   12535812 669.5 10082339 538.5
    ## Vcells 59501028 454.0   98968996 755.1 97846734 746.6

``` r
ds_dc <- PercentageFeatureSet(ds_dc, pattern = "^MT-", col.name = "percent_mt")
ds_dc <- PercentageFeatureSet(ds_dc, "^RP[SL]", col.name = "percent_ribo")
ds_dc <- PercentageFeatureSet(ds_dc, "^HB[^(P)]", col.name = "percent_hb")
ds_dc <- PercentageFeatureSet(ds_dc, "PECAM1|PF4", col.name = "percent_plat")
```

``` r
VlnPlot(ds_dc, features = c('nCount_RNA','nFeature_RNA', 'percent_mt', 'percent_hb', "percent_ribo", "percent_plat"), pt.size = 1.3, ncol = 1) +  NoLegend()
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-8-1.svg)<!-- -->

### Other qc plots

``` r
FeatureScatter(ds_dc, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = T) + 
  scale_x_continuous(labels = scales::scientific) + 
  scale_y_continuous(labels = scales::scientific)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-9-1.svg)<!-- -->
\_\_\_ ## Filtering

``` r
length(WhichCells(ds_dc, expression = nFeature_RNA > 100))
```

    ## [1] 569

``` r
length(WhichCells(ds_dc, expression = nFeature_RNA > 100 & nCount_RNA > 500))
```

    ## [1] 541

``` r
# most common genes for cells:
sort(Matrix::rowSums(ds_cf), decreasing = T)[1:10]
```

    ## MT-RNR2  MALAT1    ACTB MT-RNR1    CD74  MT-ND4  MT-CYB     FTL     B2M    PSAP 
    ## 1020525  734223  713029  585686  554177  512719  463503  460913  459010  416967

``` r
# most common genes for EVs
sort(Matrix::rowSums(ds_dc), decreasing = T)[1:10]
```

    ##         HELLPAR         MT-RNR2 ENSG00000279010 ENSG00000280156 ENSG00000279184 
    ##           72147           57170           55177           47776           36823 
    ##         MT-RNR1 ENSG00000279080 ENSG00000279738          MT-ND5          MT-ND4 
    ##           34482           29941           28346           23232           23109

## filter out EV samples that have x features and UMIs

``` r
ds_dcf <- subset(x = ds_dc, subset = nCount_RNA > 500 & nFeature_RNA > 100 & percent_mt < 90 & percent_ribo < 90)
rm(ds_dc)
suppressMessages(gc())
```

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  6837112 365.2   12535812 669.5 12535812 669.5
    ## Vcells 65730114 501.5   98968996 755.1 98960718 755.1

## Find which genes contribute to the nCount_RNA the most

``` r
C <- ds_dcf@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

par(mar=c(5,10,1,1))
boxplot(as.matrix(t(C[most_expressed,])), cex = 1, las = 1, xlab = "% total count per fraction",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
```

![](pesc_analysis_1_files/figure-gfm/most_expressed_evs_filtered-1.svg)<!-- -->

------------------------------------------------------------------------

## QC plotting of filtered reads

**Logarithmic axes.**

``` r
VlnPlot(ds_dcf, features = 'nFeature_RNA', pt.size = 1.3, log = T) +  NoLegend()
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-12-1.svg)<!-- -->

``` r
VlnPlot(ds_dcf, features = 'nCount_RNA', pt.size = 1.3, log = T) +  NoLegend()
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-12-2.svg)<!-- -->

------------------------------------------------------------------------

### Other qc plots

``` r
FeatureScatter(ds_dcf, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = T) + 
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-13-1.svg)<!-- -->

# Ratio of well annotated genes (with symbols) to the rest

**On both graphs there’re two dots per cell or EV sample - one for genes
with and one for genes without a symbol.**

``` r
ratc <- as.data.table(ds_cf@assays$RNA@counts, keep.rownames = T)
ratc <- melt(ratc, id.vars = "rn")

ratc[, "code_symb" := ifelse(rn %like% "ENSG", "ens", "symb")]
ratc <- ratc[, sum(value), by = .(variable, code_symb)]
ratc[, "ratio" := V1/sum(V1), by = .(variable)]

ratc[code_symb == "ens", mean(ratio)]
```

    ## [1] 0.06248185

``` r
ggplot(ratc) +
  geom_point(aes(x = V1, y = ratio, color = code_symb)) +
  guides(color = guide_legend(title = "", size = 19, override.aes = list(size = 5))) +
  xlab("RNA count") +
  ylab("ratio to all genes per cell") +
  ggtitle("Ratios of RNA counts with and without a symbol - cells") +
  theme(plot.title = element_text(size = 18), axis.text = element_text(size = 15), axis.title = element_text(size = 17))
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-14-1.svg)<!-- -->

``` r
ratev <- as.data.table(ds_dcf@assays$RNA@counts, keep.rownames = T)
ratev <- melt(ratev, id.vars = "rn")

ratev[, "code_symb" := ifelse(rn %like% "ENSG", "ens", "symb")]
ratev <- ratev[, sum(value), by = .(variable, code_symb)]
ratev[, "ratio" := V1/sum(V1), by = .(variable)]
ratev[code_symb == "ens", mean(ratio)]
```

    ## [1] 0.1342733

``` r
ggplot(ratev) +
  geom_point(aes(x = V1, y = ratio, color = code_symb)) +
  guides(color = guide_legend(title = "", size = 19, override.aes = list(size = 5))) +
  xlab("RNA count") +
  ylab("ratio to all genes per EV sample") +
  ggtitle("Ratios of RNA counts with and without a symbol - EVs") +
  theme(plot.title = element_text(size = 18), axis.text = element_text(size = 15), axis.title = element_text(size = 17))
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-15-1.svg)<!-- -->

``` r
rm(list = c('ratev', 'ratc'))
suppressMessages(gc())
```

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  6871932 367.1   12535812  669.5  12535812  669.5
    ## Vcells 65872346 502.6  603832016 4606.9 754772198 5758.5

------------------------------------------------------------------------

# Metadata

## Cells

Assign cells to the patient, …

``` r
full_names <- colnames(ds_cf)
reg.pat <- regmatches(full_names, regexpr("PEsc[0-9]{1,6}", full_names))
reg.pat <- gsub("PEsc", "", reg.pat)
ds_cf@meta.data$pat <- reg.pat
```

…the replicate…

``` r
reg.rep <- regmatches(full_names, regexpr("PEsc[0-9]{1,6}(.|_)[0-9]", full_names))
reg.rep <- gsub("PEsc[0-9]{1,6}(.|_)", "", reg.rep)
ds_cf@meta.data$rep <- reg.rep
```

…and diagnosis and type of biopsy.

``` r
metd <- fread("PEpatients_2_edited.csv")
metd[, c("type", "diag") := lapply(.SD, factor), .SDcols = c("type", 'diag')]
metd[, pat_numb := regmatches(metd[, pat_numb], regexpr("[0-9]{3,5}", metd[, pat_numb]))]

# replace missing type values ("") with "unknown".
metd[type == "", type := "unknown"]

patd <- as.data.table(reg.pat)
colnames(patd) <- "pat_numb"
met.full <- metd[patd, on = .(pat_numb)]


ds_cf@meta.data$type <- met.full[, type]
ds_cf@meta.data$diag <- met.full[, diag]

rm(list = c("reg.rep", "reg.pat"))
```

## EVs

Assign EVs to the fraction they’re from,…

``` r
full_names <- colnames(ds_dcf)
reg.frac <- regmatches(full_names, regexpr("_[A-Z]{1,2}_", full_names))
reg.frac <- gsub("_", "", reg.frac)

reg.frac[reg.frac == "C"] <- "cells" # change C to cells

ds_dcf@meta.data$frac <- reg.frac
```

…the replicate…

``` r
reg.rep <- regmatches(full_names, regexpr("X[0-9]{1,6}_[A-Z]{1,2}_[0-9]", full_names))
reg.rep <- gsub("X[0-9]{1,6}_[A-Z]{1,2}_", "", reg.rep)
ds_dcf@meta.data$rep <- reg.rep
```

…the patient the type of biopsy.

``` r
reg.pat <- regmatches(full_names, regexpr("X[0-9]{1,6}", full_names))
reg.pat <- gsub("X", "", reg.pat)
reg.pat[reg.pat == 804] <- 604 # a patient number correction
ds_dcf@meta.data$pat <- reg.pat


patd <- as.data.table(reg.pat)
colnames(patd) <- "pat_numb"
met.full <- metd[patd, on = .(pat_numb)]

met.full[type %like% "Considered", type := "Considered benign but \nanother sample same \nday was adca (lung)"] # add a breaks in longest labels for EVs
met.full[type %like% "Wrong", type := "Wrong number, \nnot pleural effusion"]

ds_dcf@meta.data$type <- met.full[, type]
ds_dcf@meta.data$diag <- met.full[, diag]
```

------------------------------------------------------------------------

# Clustering analysis

## Cells

### Normalization, scaling and dimensionality reduction

**Normalization was done using the SC transform described here:
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1>
as it is supposed to be depth-independent, which justifies its use in
Smartseq3 EV sequencing.**

**During the analysis a change in mapping occurred, without
significantly affecting its overall profile, possibly due to package
update.**

``` r
ds_cf <- SCTransform(ds_cf, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
  # RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE, n.components = 2)
# table(ds_cf@meta.data[c("pat", "diag", "type")])
```

Data allows to stratify cells according to their origin and diagnosis,
however different populations seem to be too mixed.

``` r
DimPlot(ds_cf, group.by = c("pat", "diag", "type", "rep"))
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-22-1.svg)<!-- -->

### PCs

**Inspect PCs - Jackstraw doesn’t work with SCTransformed data.**

**Elbow Plot** Around 10 top PCs should be enough to obtain proper
clustering.

``` r
ElbowPlot(ds_cf)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-24-1.svg)<!-- -->

### Variable features

**Top 10 most variable genes:**

| Gene            | Protein                                        | Description (from genecards)                                                                                                                                                                         |
|-----------------|------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| COL1A1, COL1A2  | Collagen 1                                     |                                                                                                                                                                                                      |
| SFTPB           | pulmonary-associated surfactant protein B      |                                                                                                                                                                                                      |
| SPARC           | cysteine-rich acidic matrix-associated protein | involved in extracellular matrix synthesis                                                                                                                                                           |
| SPP1            | Secreted Phosphoprotein 1                      | bone-osteoclasts interaction and a cytokine that upregulates expression of interferon-gamma and interleukin-12                                                                                       |
| ENSG00000203396 | WDR45-like pseudogene                          |                                                                                                                                                                                                      |
| TAGLN           | Transgelin                                     | calponin, a shape change and transformation sensitive actin-binding protein, early marker of smooth muscle differentiation, a tumor suppressor                                                       |
| IGFBP4          | insulin-like growth factor binding protein 4   | binds both insulin-like growth factors I and II and circulates in the plasma prolonging their half-life and changing interactions                                                                    |
| IGF2            | Insulin Like Growth Factor 2                   | epigenetic changes at its locus are associated with Wilms tumour, Beckwith-Wiedemann syndrome, rhabdomyosarcoma, and Silver-Russell syndrome, involved in Apoptotic Pathways in Synovial Fibroblasts |
| COL6A1          | Collagen Type VI Alpha 1 Chain                 | its mutations result in Bethlem myopathy and Ullrich muscular dystrophy, involved in platelet-derived growth factor binding                                                                          |

``` r
top10_c <- head(VariableFeatures(ds_cf), 10) # No need to use FindVariableFeatures, as SCTransform does that

p_var_c <- VariableFeaturePlot(ds_cf)
LabelPoints(p_var_c, points = top10_c, repel = T)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-25-1.svg)<!-- -->

**Top 10 variable genes and gene and RNA counts.** Most of the variable
genes are upregulated in cells from patients 3133 and 3256.

``` r
FeaturePlot(ds_cf, features = c('nFeature_RNA','nCount_RNA', top10_c), pt.size = 1, reduction = 'umap', slot = "scale.data")
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-26-1.svg)<!-- -->

### Feature plot - Umap

**RNA and gene counts seem to not influence the clustering outcome**

``` r
FeaturePlot(ds_cf, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 2, reduction = 'umap')
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-27-1.svg)<!-- -->

### Clustering

``` r
ds_cf <- FindNeighbors(ds_cf, reduction = "umap", verbose = FALSE, dims = 1:2) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

p1 <- DimPlot(ds_cf, group.by = c("pat", "diag", "type"), pt.size = 1.5, reduction = "umap")

alpha_colors <- hue_pal()(27) #  number of colors equal to number of clusters - required to change the alpha

p2 <- DimPlot(ds_cf, pt.size = 2, label = T, cols = alpha(alpha_colors, .6), reduction = "umap") + 
  labs(title = "clusters") + 
  theme(plot.title = element_text(hjust = .5))


p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/clustering_cells_umap_3-1.svg)<!-- -->

**If clusters are found using PCA:**

``` r
ds_cf <- FindNeighbors(ds_cf, reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(ds_cf, group.by = c("pat", "diag", "type"), pt.size = 1.5, reduction = "umap")


alpha_colors <- hue_pal()(14) #  number of colors equal to number of clusters - required to change the alpha

p2 <- DimPlot(ds_cf, pt.size = 2, label = T, cols = alpha(alpha_colors, .6), reduction = "umap") + 
  labs(title = "clusters") + 
  theme(plot.title = element_text(hjust = .5))


p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/clustering_cells_pca-1.svg)<!-- -->

### Post-clustering QC

**Data was further explored to exclude possible confounding effects in
data quality:**

#### Inter-replicate variance

``` r
p1 <- DimPlot(ds_cf, split.by = "pat", group.by = "rep", reduction = "pca")
p2 <- DimPlot(ds_cf, split.by = "pat", group.by = "rep")

p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/check%20inter-replicate%20variance-1.svg)<!-- -->

#### Heatmaps of genes with highest PC scores.

##### All patients

``` r
dim_plot_pca <- DimHeatmap(ds_cf, dims = 1:6, balanced = TRUE, reduction = "pca", fast = F, combine = F)

for (i in 1:length(dim_plot_pca)) {
  
  pc <- paste0("PC_", i)
  
  dim_plot_pca[[i]] <- dim_plot_pca[[i]] + theme(legend.position = "none", axis.text.y = element_text(size = 11),
          plot.title = element_text(size = 14)) +
          ggtitle(label = pc)
  
}

ggarrange(plotlist = dim_plot_pca, align = "v", common.legend = T, legend = "right")
```

![](pesc_analysis_1_files/figure-gfm/dim_heatmap_cells-1.svg)<!-- -->

##### PC1 of each patient.

``` r
plot_pat_dim_heat <- function (data, patient) {

  DimHeatmap(subset(data, subset = pat == patient), nfeatures = 30,  dims = 1, balanced = TRUE, fast = F) +
    theme(legend.position = "none", axis.text.y = element_text(size = 11),
          plot.title = element_text(size = 14)) +
    ggtitle(label = patient)

  }

pat_dim_heat <- lapply(unique(ds_cf$pat), plot_pat_dim_heat, data = ds_cf)
ggarrange(plotlist  = pat_dim_heat, common.legend = T, legend = "right")
```

![](pesc_analysis_1_files/figure-gfm/dim_heatmap_cells_patients-1.svg)<!-- -->

#### Patients 3133 and 3256.

**Cells of patients 3133, 3256 have nearly identical expression profile
despite different diagnoses - Both patients will be removed from the
analysis.**

``` r
list_heatmaps <- DimHeatmap(subset(ds_cf, subset = pat %in% c(3133, 3256)), dims = 1:6, balanced = TRUE, fast = F, combine = F) 
  

plot_dim_heat <- function (x) {
  
  heatmap_x <- list_heatmaps[[x]] + 
    ggtitle(paste0("PC_", x)) +
    theme(axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 14),
          legend.position = "none")
  
  return(heatmap_x)
}

heatmaps <- lapply(1:length(list_heatmaps), plot_dim_heat)

ggarrange(plotlist = heatmaps, align = "v", common.legend = T, legend = "right")
```

![](pesc_analysis_1_files/figure-gfm/dim_heatmap_cells_3133_3256-1.svg)<!-- -->

#### All patients - most expressed and most variable genes.

``` r
list_heatmap <- names(sort(colSums(most_expr_counts_cells), decreasing = T))

p1 <- DoHeatmap(ds_cf, features = list_heatmap, group.by = "pat", slot = 'counts') + 
    guides(color = "none") +
    theme(axis.text.y = element_text(size = 11), title = element_text(size = 13)) +
    ggtitle(label = "top 10 most expressed")
p2 <- DoHeatmap(ds_cf, features = top10_c, group.by = "pat", slot = "scale.data") + 
    guides(color = "none") +
    theme(axis.text.y = element_text(size = 11), title = element_text(size = 13)) +
    ggtitle(label = "top 10 most variable (scaled data)")

p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/dim_heatmap_cells_most_expr_most_var-1.svg)<!-- -->

``` r
DimPlot(ds_cf, group.by = c("pat", "diag", "type", "rep"))
```

## Cells without patients 3133 and 3256.

### PCA

``` r
# Subset the dataset
pat2 <- unique(ds_cf$pat)
pat2 <- pat2[1:6]

ds_cf2 <- subset(ds_cf, subset = pat %in% pat2)

ds_cf2 <- SCTransform(ds_cf2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE or
    ## useNames = TRUE.

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

### Elbow Plot

Around 15 top PCAs should be enough to obtain proper clustering.

``` r
ElbowPlot(ds_cf2, ndims = 30)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-28-1.svg)<!-- -->

**10 PCAs was allowing to obtain lesser separation of cells with regard
to replicates, however, 15 provided with better resolution of
clusters.**

### Umap

``` r
title_list <- c("patient", "diagnosis", "type of biopsy", "replicate")
plot_list <- DimPlot(ds_cf2, group.by = c("pat", "diag", "type", "rep"), combine = F)

for (i in 1:length(plot_list)) {

  plot_list[[i]] <- plot_list[[i]] + ggtitle(label = title_list[i])
}

plot <- ggarrange(plotlist = plot_list, align = "v")

annotate_figure(plot, top = text_grob("UMAP of PCA results", face = "bold", size = 18))
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-29-1.svg)<!-- -->

### Correlation of metadata and PCs

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
covariates <- ds_cf2@meta.data

pcs <- ds_cf2@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

heatmap(r2, cexCol = 1)
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata-1.svg)<!-- -->

### PCs

``` r
DimPlot(ds_cf2, group.by = c("pat", "diag", "type", "rep"), reduction = "pca")
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-30-1.svg)<!-- -->

### Inter-replicate variance

``` r
p1 <- DimPlot(ds_cf2, split.by = "pat", group.by = "rep", reduction = "pca")
p2 <- DimPlot(ds_cf2, split.by = "pat", group.by = "rep")

p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/check%20inter-replicate%20variance2-1.svg)<!-- -->

### Heatmaps of genes with highest PC scores.

#### First 20 PCs of all patients

``` r
dim_plot_pca <- DimHeatmap(ds_cf2, dims = 1:20, nfeatures = 20, balanced = TRUE, fast = F, combine = F)

for (i in 1:length(dim_plot_pca)) {
  
  pc <- paste0("PC_", i)
  
  dim_plot_pca[[i]] <- dim_plot_pca[[i]] + theme(legend.position = "none", axis.text.y = element_text(size = 11),
          plot.title = element_text(size = 14)) +
          ggtitle(label = pc)
  
}

ggarrange(plotlist = dim_plot_pca, align = "v", ncol = 3, nrow = 7, common.legend = T, legend = "right")
```

![](pesc_analysis_1_files/figure-gfm/dim_heatmap_cells_2-1.svg)<!-- -->

#### PC1 of each patient.

``` r
plot_pat_dim_heat <- function (data, patient) {

  DimHeatmap(subset(data, subset = pat == patient), nfeatures = 30,  dims = 1, balanced = TRUE, fast = F) +
    theme(legend.position = "none", axis.text.y = element_text(size = 11),
          plot.title = element_text(size = 14)) +
    ggtitle(label = patient)

  }

pat_dim_heat <- lapply(unique(ds_cf2$pat), plot_pat_dim_heat, data = ds_cf)
ggarrange(plotlist  = pat_dim_heat, align = "v", common.legend = T, legend = "right")
```

![](pesc_analysis_1_files/figure-gfm/dim_heatmap_cells_patients_2-1.svg)<!-- -->

### Variable features

Top 22 variable genes are annotated.

``` r
top22_c <- head(VariableFeatures(ds_cf2), 22)

p_var_c <- VariableFeaturePlot(ds_cf2)
LabelPoints(p_var_c, points = top22_c, repel = T, xnudge = .2, ynudge = 1, max.overlaps = 20)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-31-1.svg)<!-- -->

Top 22 variable genes and

``` r
FeaturePlot(ds_cf2, features = c('nFeature_RNA','nCount_RNA', top22_c), pt.size = 1, reduction = 'umap', slot = "scale.data") & 
  theme(plot.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 10))
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-32-1.svg)<!-- -->

``` r
top22_c_an <- data.table(symbol = top22_c)
gt_s_an <- setDT(grch38[c("symbol", "description")])

top22_c_an <- gt_s_an[top22_c_an, on = .(symbol), mult = 'first']

knitr::kable(top22_c_an, caption = "**Top 22 most variable genes**")
```

| symbol          | description                                              |
|:----------------|:---------------------------------------------------------|
| SPP1            | secreted phosphoprotein 1                                |
| SFTPB           | surfactant protein B                                     |
| ENSG00000203396 | NA                                                       |
| KRT18           | keratin 18                                               |
| ZNF717          | zinc finger protein 717                                  |
| COL14A1         | collagen type XIV alpha 1 chain                          |
| COL6A1          | collagen type VI alpha 1 chain                           |
| KRT19           | keratin 19                                               |
| CRYAB           | crystallin alpha B                                       |
| COL3A1          | collagen type III alpha 1 chain                          |
| ISG15           | ISG15 ubiquitin like modifier                            |
| LTBP2           | latent transforming growth factor beta binding protein 2 |
| COL18A1         | collagen type XVIII alpha 1 chain                        |
| SAA1            | serum amyloid A1                                         |
| COL1A1          | collagen type I alpha 1 chain                            |
| LTB             | lymphotoxin beta                                         |
| IGF2            | insulin like growth factor 2                             |
| SAMD15          | sterile alpha motif domain containing 15                 |
| COL11A1         | collagen type XI alpha 1 chain                           |
| PTGDS           | prostaglandin D2 synthase                                |
| AHNAK2          | AHNAK nucleoprotein 2                                    |
| COL1A2          | collagen type I alpha 2 chain                            |

**Top 22 most variable genes**

### Clustering

``` r
ds_cf2 <- FindNeighbors(ds_cf2, verbose = FALSE, dims = 1:15) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
```

``` r
title_list <- c("patient", "diagnosis", "type of biopsy", "replicate")
plot_list <- DimPlot(ds_cf2, group.by = c("pat", "diag", "type", "rep"),  pt.size = 1, combine = F)

for (i in 1:length(plot_list)) {

  plot_list[[i]] <- plot_list[[i]] + ggtitle(label = title_list[i])
}

p1 <- ggarrange(plotlist = plot_list, align = "v")

alpha_colors <- hue_pal()(10) #  number of colors equal to number of clusters - required to change the alpha

p2 <- DimPlot(ds_cf2, pt.size = 2, label = T, cols = alpha(alpha_colors, .6), reduction = "umap") + 
  labs(title = "clusters") + 
  theme(plot.title = element_text(hjust = .5))


plot <- ggarrange(p1, p2, widths = c(3, 2))

annotate_figure(plot, top = text_grob("UMAP of PCA results and clustering analysis", face = "bold", size = 18))
```

![](pesc_analysis_1_files/figure-gfm/finding_clusters_cells2_plotting-1.svg)<!-- -->

## DE analysis

### DESeq2 - high minimum detection rate - .4 - not run

``` r
c_markers <- FindAllMarkers(ds_cf2, assay = "RNA", test.use = "DESeq2", min.pct = 0.4, verbose = F) # minimal feature frequency of 40%

for(i in levels(c_markers$cluster)) {print(head(c_markers[c_markers$cluster == i,], 10))}
```

Plot most differentially expressed genes

``` r
setDT(c_markers)
highest_de <- c_markers[order(avg_log2FC, decreasing = T), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]
lowest_de <- c_markers[order(avg_log2FC), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]

most_diff_c <- merge(highest_de, lowest_de, all = T)

ggplot(c_markers, aes(x = cluster, y = avg_log2FC, color = -log10(p_val))) +
  geom_point(position = position_quasirandom()) +
  geom_text_repel(most_diff_c, mapping = aes(label = gene),
                  size = 5, color = "red", fontface = "bold", 
                  force = 3, min.segment.length = 0) +
  geom_hline(yintercept = 0) +
  scale_color_continuous(trans = "log", breaks = c(1, 10, 50, 125)) + 
  ggtitle(label = "Fold change and p-values of differentially expressed genes among clusters", 
          subtitle = "most differentially expressed genes are labeled")
```

``` r
highest_de <- c_markers[order(avg_log2FC, decreasing = T), .SD[1:10] ,by = cluster][p_val < 5e-2][order(cluster)]


DoHeatmap(ds_cf2, group.by = "seurat_clusters", features = highest_de$gene) + 
  NoLegend() + 
  ggtitle(label = "Expression of top markers of clusters of cells found in the dataset")
```

Genes:
`DBR1, TMEM52B, SHPK, ENSG00000269044, MAP4K2, OTUD4, LINC02503, LINC00504, NCBP3, GIHCG, TRBV12-3, INTS6, CIITA, HLA-DQA1, ASB3, GPR82, TAP1, NAAA, NIBAN1, MT-TE, DSE, ASB11, LINC02242, ENSG00000225113, ENSG00000261177, DSTYK, MMS22L, RAB12, ANAPC5, RCN2, ENSG00000280138, HBA2, CD36, ENSG00000256618, F13A1, TMEM176B, TNFSF13, CMKLR1, NRP1, MRC1`
were present in the differentially expressed set but not among the most
variable genes which might indicate, that the .4 minimum detection rate
might have been too high and excluded too many genes from the DESeq
analysis.

### DESeq2 - low minimum detection rate - .25

``` r
c_markers2 <- FindAllMarkers(ds_cf2, assay = "RNA", test.use = "DESeq2", min.pct = .25, verbose = T) # minimal feature frequency of 25%

for(i in levels(c_markers$cluster)) {
  
  print(head(c_markers[c_markers$cluster == i,], 10))
  
  }
```

Plot most differentially expressed genes

``` r
setDT(c_markers2)
highest_de <- c_markers2[order(avg_log2FC, decreasing = T), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]
lowest_de <- c_markers2[order(avg_log2FC), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]

most_diff_c <- merge(highest_de, lowest_de, all = T)

ggplot(c_markers2, aes(x = cluster, y = avg_log2FC, color = -log10(p_val))) +
  geom_point(position = position_quasirandom()) +
  geom_text_repel(most_diff_c, mapping = aes(label = gene),
                  size = 5, color = "red", fontface = "bold", 
                  force = 3, min.segment.length = 0) +
  geom_hline(yintercept = 0) +
  scale_color_continuous(trans = "log", breaks = c(1, 10, 50, 125)) + 
  ggtitle(label = "Fold change and p-values of differentially expressed genes among clusters", 
          subtitle = "most differentially expressed genes are labeled")  
```

``` r
highest_de <- c_markers2[order(avg_log2FC, decreasing = T), .SD[1:10] ,by = cluster][p_val < 5e-2][order(cluster)]


DoHeatmap(ds_cf2, group.by = "seurat_clusters", features = highest_de$gene) + 
  NoLegend() + 
  ggtitle(label = "Expression of top markers of clusters of cells found in the dataset")
```

### DE with Wilcoxon test

``` r
c_markers3 <- FindAllMarkers(ds_cf2, verbose = T)

# for(i in levels(c_markers$cluster)) {print(head(c_markers[c_markers$cluster == i,], 10))}
```

Plot most differentially expressed genes

``` r
setDT(c_markers3)
highest_de <- c_markers3[order(avg_log2FC, decreasing = T), .SD[1:5] , by = cluster][p_val < 5e-2][order(cluster)]
lowest_de <- c_markers3[order(avg_log2FC), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]

most_diff_c <- merge(highest_de, lowest_de, all = T)

ggplot(c_markers3, aes(x = cluster, y = avg_log2FC, color = -log10(p_val))) +
  geom_point(position = position_quasirandom()) +
  geom_text_repel(most_diff_c, mapping = aes(label = gene),
                  size = 5, color = "red", fontface = "bold", 
                  force = 3, min.segment.length = 0) +
  geom_hline(yintercept = 0) +
  scale_color_continuous(trans = "log", breaks = c(1, 10, 50, 115)) + 
  ggtitle(label = "Fold change and p-values of differentially expressed genes among clusters", 
          subtitle = "most differentially expressed genes are labeled") +
  theme(legend.title.align = 1,
      legend.justification = c("center"))  
```

``` r
highest_de <- c_markers3[order(avg_log2FC, decreasing = T), .SD[1:10] ,by = cluster][p_val < 5e-2][order(cluster)]

heat_cols <- hcl.colors(3, palette = "gnbu")

DoHeatmap(ds_cf2, group.by = "seurat_clusters", features = highest_de$gene, angle = 0) + 
  guides(color = "none") +
  scale_fill_gradient2("Average log2 fold change", low = heat_cols[1], mid = heat_cols[2], high = heat_cols[3], guide = guide_colourbar(title.position = "top")) +
  ggtitle(label = "Expression of top markers of clusters of cells found in the dataset")
```

### DE with MAST

**Possibly the best algorithm for sc data.**

``` r
c_markers_mast <- FindAllMarkers(ds_cf2, verbose = T, test.use = "MAST")
```

Plot most differentially expressed genes

``` r
setDT(c_markers_mast)
highest_de <- c_markers_mast[order(avg_log2FC, decreasing = T), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]
lowest_de <- c_markers_mast[order(avg_log2FC), .SD[1:5] ,by = cluster][p_val < 5e-2][order(cluster)]

most_diff_c <- merge(highest_de, lowest_de, all = T)

ggplot(c_markers_mast, aes(x = cluster, y = avg_log2FC, color = -log10(p_val))) +
  geom_point(position = position_quasirandom()) +
  geom_text_repel(most_diff_c, mapping = aes(label = gene),
                  size = 5, color = "red", fontface = "bold", 
                  force = 3, min.segment.length = 0,
                  segment.size = 1) +
  geom_hline(yintercept = 0) +
  scale_color_continuous(trans = "log", breaks = c(1, 10, 50, 115)) + 
  ggtitle(label = "Fold change and p-values of differentially expressed genes among clusters", 
          subtitle = "most differentially expressed genes are labeled") +
  theme(legend.title.align = 1,
      legend.justification = c("center"))  
```

    ## Warning: ggrepel: 3 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](pesc_analysis_1_files/figure-gfm/most_diff_bar_mast-1.svg)<!-- -->

``` r
highest_de <- c_markers_mast[order(avg_log2FC, decreasing = T), .SD[1:10] ,by = cluster][p_val < 5e-2][order(cluster)]

heat_cols <- hcl.colors(3, palette = "gnbu")

DoHeatmap(ds_cf2, group.by = "seurat_clusters", features = highest_de$gene, angle = 0) + 
  guides(color = "none") +
  scale_fill_gradient2("Average log2 fold change", low = heat_cols[1], mid = heat_cols[2], high = heat_cols[3], guide = guide_colourbar(title.position = "top")) +
  ggtitle(label = "Expression of top markers of clusters of cells found in the dataset")
```

    ## Warning in DoHeatmap(ds_cf2, group.by = "seurat_clusters", features =
    ## highest_de$gene, : The following features were omitted as they were not found
    ## in the scale.data slot for the SCT assay: SLC5A3, ENSG00000268854, MASP1,
    ## ENSG00000271755, SDC2, GPAM, FARP1, STK24, COL8A1, HHIP-AS1, CENPW, PTTG1,
    ## STMN1, TBC1D20, ERICH2-DT, KDELR2, HBP1, KARS1, DSE, LMO2, EMC3, LDLR, HNRNPAB,
    ## MCRIP2, ZNF273, BICD2, TGFBR2, LRP1, TNFSF12-TNFSF13, MRC1, TNFSF13, ITGAM

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](pesc_analysis_1_files/figure-gfm/cluster_heatmap_mast-1.svg)<!-- -->

``` r
highest_de <- c_markers_mast[order(avg_log2FC, decreasing = T), .SD[1:10] ,by = cluster][p_val < 5e-2][order(cluster)]

heat_cols <- hcl.colors(3, palette = "gnbu")

DoHeatmap(ds_cf2, group.by = "diag", features = highest_de$gene, angle = 0) + 
  guides(color = "none") +
  scale_fill_gradient2("Average log2 fold change", low = heat_cols[1], mid = heat_cols[2], high = heat_cols[3], guide = guide_colourbar(title.position = "top")) +
  ggtitle(label = "Expression of top markers of clusters of cells across diagnoses")
```

    ## Warning in DoHeatmap(ds_cf2, group.by = "diag", features = highest_de$gene, :
    ## The following features were omitted as they were not found in the scale.data
    ## slot for the SCT assay: SLC5A3, ENSG00000268854, MASP1, ENSG00000271755, SDC2,
    ## GPAM, FARP1, STK24, COL8A1, HHIP-AS1, CENPW, PTTG1, STMN1, TBC1D20, ERICH2-DT,
    ## KDELR2, HBP1, KARS1, DSE, LMO2, EMC3, LDLR, HNRNPAB, MCRIP2, ZNF273, BICD2,
    ## TGFBR2, LRP1, TNFSF12-TNFSF13, MRC1, TNFSF13, ITGAM

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](pesc_analysis_1_files/figure-gfm/cluster_heatmap_mast_diag-1.svg)<!-- -->

### GENE ONTOLOGY - CELLS - 2DO!!!!

## Cells - diagnosis separated

### Cells

#### Normalization, scaling and dimensionality reduction

Cells were split depending on patient’s diagnosis to look at
inter-diagnosis variation.

``` r
ds_cf_diag <- SplitObject(ds_cf2, split.by = "diag")

map.list.seurat <- function(x) {
  SCTransform(x, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)  
}
  

ds_cf_diag <- lapply(ds_cf_diag, map.list.seurat)

# ds_cf_diag <- merge(ds_cf_diag_list[[1]], y = ds_cf_diag_list[-1])
```

#### PCs

**Elbow Plot** 13 top PCs should be enough to obtain proper clustering.

``` r
plots.names <- names(ds_cf_diag)

plot.list <- lapply(plots.names, FUN = function(x) { ElbowPlot(ds_cf_diag[[x]]) +
                                                               plot_annotation(title = x)
                                                             })

ggarrange(plotlist = plot.list, ncol = 1, nrow = 3)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-34-1.svg)<!-- -->

**Top variable genes**

``` r
plots.names <- names(ds_cf_diag)

top10_c_diag <- function(x) {
  
  top10_c <- head(VariableFeatures(ds_cf_diag[[x]]), 10) # No need to use FindVariableFeatures, as SCTransform does that
  p_var_c <- VariableFeaturePlot(ds_cf_diag[[x]]) + 
                  plot_annotation(title = x)
  
  LabelPoints(p_var_c, points = top10_c, repel = T)
  
  }

plot.list <- lapply(plots.names, top10_c_diag)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results
    ## When using repel, set xnudge and ynudge to 0 for optimal results
    ## When using repel, set xnudge and ynudge to 0 for optimal results

``` r
ggarrange(plotlist = plot.list, ncol = 1, nrow = 3)
```

    ## Warning: ggrepel: 3 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-35-1.svg)<!-- -->

#### UMAP

``` r
plots.names <- names(ds_cf_diag)

plot.list <- lapply(plots.names, FUN = function(x) { DimPlot(ds_cf_diag[[x]], group.by = c("pat", "rep")) +
                                                               plot_annotation(title = x)
                                                             })

ggarrange(plotlist = plot.list, ncol = 1, nrow = 3)
```

![](pesc_analysis_1_files/figure-gfm/umap_diag_sep_cells_meta-1.svg)<!-- -->

``` r
plots.names <- names(ds_cf_diag)

plot.list <- lapply(plots.names, FUN = function(x) { FeaturePlot(ds_cf_diag[[x]], features = c("nCount_RNA", "nFeature_RNA")) +
                                                               plot_annotation(title = x)
                                                             })

ggarrange(plotlist = plot.list, ncol = 1, nrow = 3)
```

![](pesc_analysis_1_files/figure-gfm/umap_diag_sep_cells_features-1.svg)<!-- -->

## EVs - all

### PCA

``` r
# remove the patients removed from the cells dataset and patient 5168 (Wrong number, not pleural effusion) and 5835 (Considered benign but another sample same day was adca (lung))
pat3 <- c("3133", "3256", "5168", "5835")

ds_dcf <- subset(ds_dcf, subset = pat %notin% pat3)

ds_dcf <- SCTransform(ds_dcf, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

Data doesn’t allow to stratify cells and EVs according to their
identities and diagnoses and the mapping is homogeneous with regard to
fractions.

``` r
# DimPlot(ds_dcf, group.by = c("pat", "diag", "type", "frac"), ncol = 2) +
  # scale_color_discrete(labels = ~ stringr::str_wrap(.x, width = 20))

title_list <- c("patient", "diagnosis", "type of biopsy", "fraction")
plot_list <- DimPlot(ds_dcf, group.by = c("pat", "diag", "type", "frac"), combine = F)

for (i in 1:length(plot_list)) {

  plot_list[[i]] <- plot_list[[i]] + ggtitle(label = title_list[i])
}

plot <- ggarrange(plotlist = plot_list, align = "v")

annotate_figure(plot, top = text_grob("UMAP of PCA results", face = "bold", size = 18))
```

![](pesc_analysis_1_files/figure-gfm/dimplot_evs_all-1.svg)<!-- -->

``` r
FeaturePlot(ds_dcf, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
```

![](pesc_analysis_1_files/figure-gfm/featureplot_evs_all-1.svg)<!-- -->

### Elbow Plot

Around 10 top PCs should be enough to obtain proper clustering.

``` r
ElbowPlot(ds_dcf, ndims = 30)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-36-1.svg)<!-- -->

## EVs - without cells

### PCA

``` r
ds_dcf2 <- subset(ds_dcf, subset = frac %notin% "cells")

ds_dcf2 <- SCTransform(ds_dcf2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

### UMAP

Data still doesn’t allow to stratify EVs according to their identities
and diagnoses and the mapping is homogenous with regard to fractions.

``` r
DimPlot(ds_dcf2, group.by = c("pat", "diag", "type", "frac"), ncol = 2) +
  scale_color_discrete(labels = ~ stringr::str_wrap(.x, width = 20)) 
```

![](pesc_analysis_1_files/figure-gfm/dimplot_evs_wo_cells-1.svg)<!-- -->

``` r
FeaturePlot(ds_dcf2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
```

![](pesc_analysis_1_files/figure-gfm/featureplot_evs_wo_cells-1.svg)<!-- -->

There’s no replicate-bias in the data.

``` r
DimPlot(ds_dcf2, group.by = "rep") +
  scale_color_discrete(labels = ~ stringr::str_wrap(.x, width = 20))
```

![](pesc_analysis_1_files/figure-gfm/check_inter_replicate_var1-1.svg)<!-- -->

``` r
p1 <- DimPlot(ds_dcf2, split.by = "pat", group.by = "rep", reduction = "pca", ncol = 11)
p2 <- DimPlot(ds_dcf2, split.by = "pat", group.by = "diag", ncol = 11)

p1 + p2
```

### Elbow plot

Around 10 top PCs should be enough to obtain proper clustering.

``` r
ElbowPlot(ds_dcf2, ndims = 30)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-37-1.svg)<!-- -->

### Correlation of metadata and PCs

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
covariates <- ds_dcf2@meta.data

pcs <- ds_dcf2@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

heatmap(r2, cexCol = 1)
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs-1.svg)<!-- -->

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}

r2.list <- list()

for (i in unique(ds_dcf2@meta.data$frac)) {
  
covariates <- setDT(subset(ds_dcf2@meta.data, subset = frac == i))
covariates <- covariates[, !"frac"] # remove frac as it will have only one type per loop

pcs <- subset(ds_dcf2, subset = frac == i)@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, get(cov)], pcs[, pc])
  }
}

r2.list[[i]] <- r2

}

lapply(names(r2.list), function(x) heatmap(r2.list[[x]], main = x, 
                                           Colv = NA, Rowv = NA, margins = c(5, 1)))
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction-1.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction-2.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction-3.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction-4.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction-5.svg)<!-- -->

    ## [[1]]
    ## [[1]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[1]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[1]]$Rowv
    ## NULL
    ## 
    ## [[1]]$Colv
    ## NULL
    ## 
    ## 
    ## [[2]]
    ## [[2]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[2]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[2]]$Rowv
    ## NULL
    ## 
    ## [[2]]$Colv
    ## NULL
    ## 
    ## 
    ## [[3]]
    ## [[3]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[3]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[3]]$Rowv
    ## NULL
    ## 
    ## [[3]]$Colv
    ## NULL
    ## 
    ## 
    ## [[4]]
    ## [[4]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[4]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[4]]$Rowv
    ## NULL
    ## 
    ## [[4]]$Colv
    ## NULL
    ## 
    ## 
    ## [[5]]
    ## [[5]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[5]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[5]]$Rowv
    ## NULL
    ## 
    ## [[5]]$Colv
    ## NULL

### Distribution of features among fractions

As percentages of different RNA subtypes seem to have a high impact on
PCs their distribution between different fractions was explored. To
decrease this and the bias caused by molecule and feature count FT will
be removed.

``` r
# fractions <- unique(ds_dcf2@meta.data$frac)
# 
# values.names <- c('nCount_RNA','nFeature_RNA', 'percent_mt', 'percent_hb', "percent_ribo", "percent_plat")
# 
# 
# plot.list <- list()
# 
# dist.frac <- function(x, y) {
#   
#   plot <- ggplot(data = y, aes(x = frac, y = get(x))) +
#     geom_boxplot() +
#     ggtitle(x)
#   
# }
# 
# plot.list <- lapply(values.names, FUN = dist.frac, y = ds_dcf2@meta.data)
# 
# ggarrange(plotlist = plot.list)

VlnPlot(ds_dcf2, features = c('nCount_RNA','nFeature_RNA', 'percent_mt', 'percent_hb', "percent_ribo", "percent_plat"), group.by = "frac")
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-38-1.svg)<!-- -->

## EVs - without cells and FT

### PCA

``` r
ds_dcf3 <- subset(ds_dcf2, subset = frac %notin% "FT")

ds_dcf3 <- SCTransform(ds_dcf3, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

### UMAP

Data still doesn’t allow to stratify cells and EVs according to their
identities and diagnoses and the mapping is homogenous with regard to
fractions.

``` r
DimPlot(ds_dcf3, group.by = c("pat", "diag", "type", "frac"), ncol = 2) +
  scale_color_discrete(labels = ~ stringr::str_wrap(.x, width = 20)) 
```

![](pesc_analysis_1_files/figure-gfm/dimplot_evs_wo_cells_FT-1.svg)<!-- -->

``` r
FeaturePlot(ds_dcf3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
```

![](pesc_analysis_1_files/figure-gfm/featureplot_evs_wo_cells_FT-1.svg)<!-- -->

### Elbow plot

Around 10 top PCs should be enough to obtain proper clustering.

``` r
ElbowPlot(ds_dcf3, ndims = 30)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-39-1.svg)<!-- -->

### Correlation of metadata and PCs

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
covariates <- ds_dcf3@meta.data

pcs <- ds_dcf3@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

heatmap(r2, cexCol = 1)
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs_wo_cells_FT-1.svg)<!-- -->

To further decrease the molecule and feature count bias SP will be
removed as it doesn’t correlate as much as other fractions with
diagnosis.

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}

r2.list <- list()

for (i in unique(ds_dcf3@meta.data$frac)) {
  
covariates <- setDT(subset(ds_dcf3@meta.data, subset = frac == i))
covariates <- covariates[, !"frac"] # remove frac as it will have only one type per loop

pcs <- subset(ds_dcf3, subset = frac == i)@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, get(cov)], pcs[, pc])
  }
}

r2.list[[i]] <- r2

}

lapply(names(r2.list), function(x) heatmap(r2.list[[x]], main = x, 
                                           Colv = NA, Rowv = NA, margins = c(5, 1)))
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT-1.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT-2.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT-3.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT-4.svg)<!-- -->

    ## [[1]]
    ## [[1]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[1]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[1]]$Rowv
    ## NULL
    ## 
    ## [[1]]$Colv
    ## NULL
    ## 
    ## 
    ## [[2]]
    ## [[2]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[2]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[2]]$Rowv
    ## NULL
    ## 
    ## [[2]]$Colv
    ## NULL
    ## 
    ## 
    ## [[3]]
    ## [[3]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[3]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[3]]$Rowv
    ## NULL
    ## 
    ## [[3]]$Colv
    ## NULL
    ## 
    ## 
    ## [[4]]
    ## [[4]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[4]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[4]]$Rowv
    ## NULL
    ## 
    ## [[4]]$Colv
    ## NULL

## EVs - without cells, FT and SP

### PCA

``` r
ds_dcf4 <- subset(ds_dcf3, subset = frac %notin% "SP")

ds_dcf4 <- SCTransform(ds_dcf4, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

### UMAP

Data still doesn’t allow to stratify cells and EVs according to their
identities and diagnoses and the mapping is homogenous with regard to
fractions.

``` r
DimPlot(ds_dcf4, group.by = c("pat", "diag", "type", "frac"), ncol = 2) +
  scale_color_discrete(labels = ~ stringr::str_wrap(.x, width = 20)) 
```

![](pesc_analysis_1_files/figure-gfm/dimplot_evs_wo_cells_FT_SP-1.svg)<!-- -->

``` r
FeaturePlot(ds_dcf4, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
```

![](pesc_analysis_1_files/figure-gfm/featureplot_evs_wo_cells_FT_SP-1.svg)<!-- -->

### Elbow plot

Around 10 top PCs should be enough to obtain proper clustering.

``` r
ElbowPlot(ds_dcf4, ndims = 30)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-40-1.svg)<!-- -->

### Correlation of metadata and PCs

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
covariates <- ds_dcf4@meta.data

pcs <- ds_dcf4@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

heatmap(r2, cexCol = 1)
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs_wo_cells_FT_SP-1.svg)<!-- -->

**PCs are still closely tied to the feature and molecule counts as well
as transcript of red cell, mitochondrial or platelet origin. That
indicates, that most other genes might be expressed homogeneously.**

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}

r2.list <- list()

for (i in unique(ds_dcf4@meta.data$frac)) {
  
covariates <- setDT(subset(ds_dcf4@meta.data, subset = frac == i))
covariates <- covariates[, !"frac"] # remove frac as it will have only one type per loop

pcs <- subset(ds_dcf4, subset = frac == i)@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, get(cov)], pcs[, pc])
  }
}

r2.list[[i]] <- r2

}

lapply(names(r2.list), function(x) heatmap(r2.list[[x]], main = x, 
                                           Colv = NA, Rowv = NA, margins = c(5, 1)))
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT_SP-1.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT_SP-2.svg)<!-- -->![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs%20per%20fraction_wo_cells_FT_SP-3.svg)<!-- -->

    ## [[1]]
    ## [[1]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[1]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[1]]$Rowv
    ## NULL
    ## 
    ## [[1]]$Colv
    ## NULL
    ## 
    ## 
    ## [[2]]
    ## [[2]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[2]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[2]]$Rowv
    ## NULL
    ## 
    ## [[2]]$Colv
    ## NULL
    ## 
    ## 
    ## [[3]]
    ## [[3]]$rowInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
    ## 
    ## [[3]]$colInd
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30
    ## 
    ## [[3]]$Rowv
    ## NULL
    ## 
    ## [[3]]$Colv
    ## NULL

## EVs - without cells, FT and SP and platelet, ribosomal, mitochondrial transcripts and hemoglobin

### PCA

``` r
genes_remove <- which(rownames(ds_dcf4@assays$RNA) %like% "^MT-" | 
                        rownames(ds_dcf4@assays$RNA) %like% "^RP[SL]" | 
                        rownames(ds_dcf4@assays$RNA) %like% "^HB[^(P)]" |
                        rownames(ds_dcf4@assays$RNA) %like% "PECAM1|PF4") 

genes_stay <- rownames(ds_dcf4@assays$RNA[-genes_remove])
ds_dcf5 <- subset(ds_dcf4, features = genes_stay)

# Check if all the genes were removed:
ds_dcf5 <- PercentageFeatureSet(ds_dcf5, "^MT-", col.name = "percent_mt")
ds_dcf5 <- PercentageFeatureSet(ds_dcf5, "^RP[SL]", col.name = "percent_ribo")
ds_dcf5 <- PercentageFeatureSet(ds_dcf5, "^HB[^(P)]", col.name = "percent_hb")
ds_dcf5 <- PercentageFeatureSet(ds_dcf5, "PECAM1|PF4", col.name = "percent_plat")
VlnPlot(ds_dcf5, features = c('nCount_RNA','nFeature_RNA', 'percent_mt', 'percent_hb', "percent_ribo", "percent_plat"), group.by = "frac")
```

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent_mt.

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent_hb.

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent_ribo.

    ## Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
    ## idents, : All cells have the same value of percent_plat.

![](pesc_analysis_1_files/figure-gfm/dim_reduc_evs_wo_cells_FT_SP_SP_genes-1.svg)<!-- -->

``` r
ds_dcf5 <- SCTransform(ds_dcf5, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
```

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
    ## or useNames = TRUE.

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

### UMAP

Data still doesn’t allow to stratify cells and EVs according to their
identities and diagnoses and the mapping is homogenous with regard to
fractions.

``` r
DimPlot(ds_dcf5, group.by = c("pat", "diag", "type", "frac"), ncol = 2) +
  scale_color_discrete(labels = ~ stringr::str_wrap(.x, width = 20)) 
```

![](pesc_analysis_1_files/figure-gfm/dimplot_evs_wo_cells_FT_SP_genes-1.svg)<!-- -->

``` r
FeaturePlot(ds_dcf5, features =c('nCount_RNA','nFeature_RNA'), ncol = 2)
```

![](pesc_analysis_1_files/figure-gfm/featureplot_evs_wo_cells_FT_SP_genes-1.svg)<!-- -->

### Elbow plot

Around 10 top PCs should be enough to obtain proper clustering.

``` r
ElbowPlot(ds_dcf5, ndims = 30)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-41-1.svg)<!-- -->

### Correlation of metadata and PCs

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
covariates <- setDT(ds_dcf5@meta.data)
covariates <- covariates[, !"percent_mt"][, !"percent_ribo"][, !"percent_hb"][, !"percent_plat"]
pcs <- ds_dcf5@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, ..cov][[1]], pcs[, pc])
  }
}

heatmap(r2, cexCol = 1)
```

![](pesc_analysis_1_files/figure-gfm/correlation%20of%20PCs%20with%20metadata%20EVs_wo_cells_FT_SP_genes-1.svg)<!-- -->

**PCs are still closely tied to the feature and molecule counts. That
indicates, again, that most other genes might be expressed
homogeneously.**

``` r
# adapted from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}

r2.list <- list()

for (i in unique(ds_dcf5@meta.data$frac)) {
  
covariates <- setDT(subset(ds_dcf5@meta.data, subset = frac == i))
covariates <- covariates[, !"frac"][,!"percent_mt"][, !"percent_ribo"][, !"percent_hb"][, !"percent_plat"] # remove frac as it will have only one type per loop and others removed from the dataset

pcs <- subset(ds_dcf5, subset = frac == i)@reductions$pca@cell.embeddings[,1:30]

r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, get(cov)], pcs[, pc])
  }
}

r2.list[[i]] <- r2

}

plot_list <- lapply(names(r2.list), function(x) heatmap(r2.list[[x]], main = x, 
                                           Colv = NA, Rowv = NA, margins = c(5, 1)))
```

### PCs

NOT FINAL - Copied for EVs

**Inspect PCs - Jackstraw doesn’t work with SCTransformed data.**

**Elbow Plot** Around 10 top PCs should be enough to obtain proper
clustering.

``` r
ElbowPlot(ds_dcf)
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-43-1.svg)<!-- -->

### Variable features

``` r
top10_c <- head(VariableFeatures(ds_cf), 10) # No need to use FindVariableFeatures, as SCTransform does that

p_var_c <- VariableFeaturePlot(ds_cf)
LabelPoints(p_var_c, points = top10_c, repel = T)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-44-1.svg)<!-- -->

**Top 10 variable genes and gene and RNA counts.** Most of the variable
genes are upregulated in cells from patients 3133 and 3256.

``` r
FeaturePlot(ds_cf, features = c('nFeature_RNA','nCount_RNA', top10_c), pt.size = 1, reduction = 'umap', slot = "scale.data")
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-45-1.svg)<!-- -->

### Feature plot - Umap

**RNA and gene counts seem to not influence the clustering outcome**

``` r
FeaturePlot(ds_cf, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 2, reduction = 'umap')
```

![](pesc_analysis_1_files/figure-gfm/unnamed-chunk-46-1.svg)<!-- -->

### Clustering

``` r
ds_cf <- FindNeighbors(ds_cf, reduction = "umap", verbose = FALSE, dims = 1:2) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

p1 <- DimPlot(ds_cf, group.by = c("pat", "diag", "type"), pt.size = 1.5, reduction = "umap")

alpha_colors <- hue_pal()(27) #  number of colors equal to number of clusters - required to change the alpha

p2 <- DimPlot(ds_cf, pt.size = 2, label = T, cols = alpha(alpha_colors, .6), reduction = "umap") + 
  labs(title = "clusters") + 
  theme(plot.title = element_text(hjust = .5))


p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/clustering_ev_umap_2-1.svg)<!-- -->

**If clusters are found using PCA:**

``` r
ds_cf <- FindNeighbors(ds_cf, reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(ds_cf, group.by = c("pat", "diag", "type"), pt.size = 1.5, reduction = "umap")


alpha_colors <- hue_pal()(14) #  number of colors equal to number of clusters - required to change the alpha

p2 <- DimPlot(ds_cf, pt.size = 2, label = T, cols = alpha(alpha_colors, .6), reduction = "umap") + 
  labs(title = "clusters") + 
  theme(plot.title = element_text(hjust = .5))


p1 + p2
```

![](pesc_analysis_1_files/figure-gfm/clustering_ev_pca-1.svg)<!-- -->

``` r
knitr::knit_exit()
```
