# 06.23
# This script looks for the lineage of EVs among cells, both from pleural effusion samples. 


# install.packages("Seurat")
# install.packages("remotes")
# BiocManager::install("glmGamPoi")
# remotes::install_github("stephenturner/annotables")
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(future)
library(annotables) # for turning Ensembl ID to symbol
library(sctransform) # for normalization
library(glmGamPoi) # for normalization

"%notin%" <- Negate("%in%")
"%notlike%" <- Negate("%like%")


theme_set(new = theme_classic())
theme_update(
  axis.text.x = element_text(vjust = 0.5),
  strip.background = element_rect(fill = '#FFFFFF'),
  plot.title = element_text(hjust = 0.5, size = 25),
  axis.title = element_text(size = 23),
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 18),
  legend.key.size = unit(2, 'line'),
  legend.title = element_text(size = 20)
)

wd <- "~/OneDrive - Karolinska Institutet/Dokument/nord_lab/seq/0423_pesc/"
setwd(wd)


plan("multicore", workers = 8)
plan()

### Load the data ----

# Commented-out is a demonstration of advantages of a sparse matrix ("d" for data).
# d_dc <- as.matrix(read.csv("data/DC_matrix.txt", sep = "\t"))
# s1 <- object.size(d_dc)
# d_dc <- as(d_dc, "sparseMatrix")
# s2 <- object.size(d_dc)
# s1-s2

# EVs
d_dc <- as.sparse(read.csv("data/DC_matrix.txt", sep = "\t"))
# Cells
d_c <- as.sparse(read.csv("data/PEsc_matrix.txt", sep = "\t"))
# gene names
gt <- fread("data/genes_title.txt", sep = "\t")

# get gene symbols from annotables
gt_s <- setDT(grch38[, c("ensgene", "symbol")])
colnames(gt_s)[1] <- "Geneid"

# join both and change symbols for Geneid when the grch38 doesn't include any particular ID, mult set to 'first', to exclude duplicates/synonyms
gt_join <- gt_s[gt, on = .(Geneid), mult = 'first'][symbol %in% NA | symbol == "", symbol := Geneid]
gt_sym <- gt_join[, symbol]

# gt_mito <- readLines("code_for_daniel/human_mitochondrial_genes_list.txt")[-1] # Could be used for Daniel's approach
# gt_join[Geneid %in% gt_mito] # same genes are mitochondrial here and in Daniel's table - all is fine

# insert gene names
rownames(d_c) <- gt_sym
rownames(d_dc) <- gt_sym


### Cells! ----

ds_c <- CreateSeuratObject(count = d_c, min.cells = 0, min.features = 1, project = "cells")

ds_c <- PercentageFeatureSet(ds_c, pattern = "^MT-", col.name = "percent_mt")
ds_c <- PercentageFeatureSet(ds_c, "^RP[SL]", col.name = "percent_ribo")
# ds_c <- PercentageFeatureSet(ds_c, "^HB[^(P)]", col.name = "percent_hb")
# ds_c <- PercentageFeatureSet(ds_c, "PECAM1|PF4", col.name = "percent_plat")


# plotting ----

unf <- VlnPlot(ds_c, features = c('percent_mt', "percent_ribo")) # unfiltered
unf_l <- VlnPlot(ds_c, features = c('nFeature_RNA','nCount_RNA'), log = T) ## UGLY
# ggsave('test.png', unf)
ggsave('qc_pesc_perc_unf.png', unf, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_counts_unf.png', unf_l, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")


# Other qc plots

CF <- FeatureScatter(ds_c, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = F) + scale_x_continuous(labels = scales::scientific)
ggsave("counts-features_unf.png", CF, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white") 


# Filtering ----

length(WhichCells(ds_c, expression = nCount_RNA > 800))
length(WhichCells(ds_c, expression = nFeature_RNA > 1000 & nCount_RNA > 1000))

rownames(ds_c)[Matrix::rowSums(ds_c) > 100000]

### filter out cells that have >1000 features and UMIs as they are the only ones You can work with (subsetting on the percentages of different transcript types can be done later)
# the percent of ribosomal genes might be low due to low transcriptional activity of the cells (5%, which was used before removes too many cells)
ds_cf <- subset(x = ds_c, subset = nCount_RNA > 1000 & nFeature_RNA > 600 & percent_mt < 25 & percent_ribo > .5)

# Find which genes contribute to the nCount_RNA the most

C <- ds_cf@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

png("most_expr_filt.png", width = 15, height = 12, bg = "white", units = "in", res = 320)
boxplot(as.matrix(t(C[most_expressed,])), cex = 1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

dev.off()

# QC plotting of filtered reads - LEAVE THE GGPLOT APPROACH FOR LATER ----


filt <- VlnPlot(ds_cf, features = c('percent_mt', "percent_ribo")) # unfiltered
filt_l <- VlnPlot(ds_cf, features = c('nFeature_RNA','nCount_RNA'), log = T) ## UGLY
# ggsave('test.png', unf)
ggsave('qc_pesc_perc_filt.png', filt, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_counts_filt.png', filt_l, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")



p1 <- ggplot(as.data.table(p3gf$nFeature_RNA)) +
  geom_violin(aes(x = "BM", y = V1)) +
  geom_jitter(aes(x = "BM", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('Genes') +
  xlab("") +
  ylab("count")

p2 <- ggplot(as.data.table(p3gf$nCount_RNA)) +
  geom_violin(aes(x = "BM", y = V1)) +
  geom_jitter(aes(x = "BM", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('RNAs') +
  xlab("") +
  ylab("")

p3 <- ggplot(as.data.table(p3gf$percent_mt)) +
  geom_violin(aes(x = "BM", y = V1)) +
  geom_jitter(aes(x = "BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of mitochondrial') +
  xlab("") +
  ylab("percent")

p4 <- ggplot(as.data.table(p3gf$percent_ribo)) +
  geom_violin(aes(x = "BM", y = V1)) +
  geom_jitter(aes(x = "BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of ribosomal') +
  xlab("") +
  ylab("")

p5 <- ggplot(as.data.table(p3gf$percent_hb)) +
  geom_violin(aes(x = "BM", y = V1)) +
  geom_jitter(aes(x = "BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of hemoglobin') +
  xlab("") +
  ylab("")

p6 <- ggplot(as.data.table(p3gf$percent_plat)) +
  geom_violin(aes(x = "BM", y = V1)) +
  geom_jitter(aes(x = "BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of platelet') +
  xlab("") +
  ylab("")

qcplots <- ggarrange(p1, p2, p3, p4, p5, p6, nrow = 1)
ggsave("qc_violins_filt_p3.png", qcplots, device = 'png', dpi = "retina", width = 21, height = 12, bg = "white")

### EVs! ----

ds_dc <- CreateSeuratObject(count = d_dc, min.cells = 0, min.features = 1, project = "ev")

ds_dc <- PercentageFeatureSet(ds_dc, pattern = "^MT-", col.name = "percent_mt")
ds_dc <- PercentageFeatureSet(ds_dc, "^RP[SL]", col.name = "percent_ribo")
# ds_c <- PercentageFeatureSet(ds_c, "^HB[^(P)]", col.name = "percent_hb")
# ds_c <- PercentageFeatureSet(ds_c, "PECAM1|PF4", col.name = "percent_plat")


## plotting ----

unf_ev <- VlnPlot(ds_dc, features = c('percent_mt', "percent_ribo")) # unfiltered
# unf_l_ev <- VlnPlot(ds_dc, features = c('nFeature_RNA','nCount_RNA'), log = T)
unf_ev_2 <- VlnPlot(ds_dc, features = c('nFeature_RNA','nCount_RNA'))
# ggsave('test.png', unf)
ggsave('qc_pesc_ev_perc_unf.png', unf_ev, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_ev_counts_unf.png', unf_ev_2, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")

## separate ----

# percentages
unf_ev_mit <- VlnPlot(ds_dc, features = 'percent_mt', pt.size = 1.3) +  NoLegend()
unf_ev_ribo <- VlnPlot(ds_dc, features = 'percent_ribo', pt.size = 1.3) +  NoLegend()
ggsave('qc_pesc_ev_perc1_unf.png', unf_ev_mit, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_ev_perc2_unf.png', unf_ev_ribo, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")

# counts
unf_ev_genes <- VlnPlot(ds_dc, features = 'nFeature_RNA', pt.size = 1.3) +  NoLegend()
unf_ev_rna <- VlnPlot(ds_dc, features = 'nCount_RNA', pt.size = 1.3) +  NoLegend()
ggsave('qc_pesc_ev_genes_unf.png', unf_ev_genes, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_ev_rna_unf.png', unf_ev_rna, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")

# log
unf_ev_genes <- VlnPlot(ds_dc, features = 'nFeature_RNA', pt.size = 1.3, log = T) +  NoLegend()
unf_ev_rna <- VlnPlot(ds_dc, features = 'nCount_RNA', pt.size = 1.3, log = T) +  NoLegend()
ggsave('qc_pesc_ev_genes_unf_log.png', unf_ev_genes, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_ev_rna_unf_log.png', unf_ev_rna, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")

# Other qc plots ----

CF <- FeatureScatter(ds_dc, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = F) + scale_x_continuous(labels = scales::scientific)
ggsave("counts-features_ev_unf.png", CF, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")


# Filtering  ----

length(WhichCells(ds_dc, expression = nFeature_RNA > 800))
length(WhichCells(ds_dc, expression = nFeature_RNA > 800 & nCount_RNA > 800))
rownames(ds_c)[Matrix::rowSums(ds_c) > 100000]

### filter out EVs that have x features and UMIs
ds_dcf <- subset(x = ds_dc, subset = nCount_RNA > 100 & nFeature_RNA > 500 & percent_mt < 90 & percent_ribo < 90)

# Find which genes contribute to the nCount_RNA the most
C <- ds_dcf@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

png("most_expr_ev_filt.png", width = 15, height = 12, bg = "white", units = "in", res = 320)
boxplot(as.matrix(t(C[most_expressed,])), cex = 1, las = 1, xlab = "% total count per fraction",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

dev.off()

# QC plotting of filtered reads ----
# percentages
filt_ev_mit <- VlnPlot(ds_dcf, features = 'percent_mt', pt.size = 1.3) +  NoLegend()
filt_ev_ribo <- VlnPlot(ds_dcf, features = 'percent_ribo', pt.size = 1.3) +  NoLegend()
ggsave('qc_pesc_ev_perc1_filt.png', filt_ev_mit, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_ev_perc2_filt.png', filt_ev_ribo, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")

# log
filt_ev_genes <- VlnPlot(ds_dcf, features = 'nFeature_RNA', pt.size = 1.3, log = T) +  NoLegend()
filt_ev_rna <- VlnPlot(ds_dcf, features = 'nCount_RNA', pt.size = 1.3, log = T) +  NoLegend()
ggsave('qc_pesc_ev_genes_filt_log.png', filt_ev_genes, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")
ggsave('qc_pesc_ev_rna_filt_log.png', filt_ev_rna, device = 'png', dpi = "retina", width = 16, height = 9, bg = "white")

# Other qc plots ----

CF <- FeatureScatter(ds_dcf, "nCount_RNA", "nFeature_RNA", pt.size = 1, plot.cor = F) + scale_x_continuous(labels = scales::scientific)
ggsave("counts-features_ev_filt.png", CF, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")


# Ratio of well annotated genes (with symbols) to the rest ----
# 
ratc <- ds_cf@assays$RNA@counts
ratev <- ds_dcf@assays$RNA@counts
# 
# ratc_s <- as.data.table(Matrix::rowSums(ratc), keep.rownames = T)
# ratev_s <- as.data.table(Matrix::rowSums(ratev), keep.rownames = T)
# 
# ratc_s[, "V3" := ifelse(V1 %like% "ENSG", "ens", "symb")]
# ratev_s[, "V3" := ifelse(V1 %like% "ENSG", "ens", "symb")]
# 
# ratc_s[, sum(V2), by = V3][, .SD[V3 == "ens", V1]/.SD[, sum(V1)]]
# ratev_s[, sum(V2), by = V3][, .SD[V3 == "ens", V1]/.SD[, sum(V1)]]

# How variable is that between samples?

ratc <- as.data.table(ratc, keep.rownames = T)
ratc <- melt(ratc, id.vars = "rn")

ratc[, "code_symb" := ifelse(rn %like% "ENSG", "ens", "symb")]
ratc <- ratc[, sum(value), by = .(variable, code_symb)]
ratc[, "ratio" := V1/sum(V1), by = .(variable)]
# the overall ratio
ratc[code_symb == "ens", mean(ratio)]

ggplot(ratc) +
  geom_point(aes(x = V1, y = ratio, color = code_symb)) +
  guides(color = guide_legend(title = "", size = 19, override.aes = list(size = 5))) +
  xlab("RNA count") +
  ylab("ratio to all genes per cell") +
  ggtitle("Ratios of RNA counts with and without a symbol - cells")



ratev <- as.data.table(ratev, keep.rownames = T)
ratev <- melt(ratev, id.vars = "rn")

ratev[, "code_symb" := ifelse(rn %like% "ENSG", "ens", "symb")]
ratev <- ratev[, sum(value), by = .(variable, code_symb)]
ratev[, "ratio" := V1/sum(V1), by = .(variable)]
# the overall ratio
ratev[code_symb == "ens", mean(ratio)]

ggplot(ratev) +
  geom_point(aes(x = V1, y = ratio, color = code_symb)) +
  guides(color = guide_legend(title = "", size = 19, override.aes = list(size = 5))) +
  xlab("RNA count") +
  ylab("ratio to all genes per EV sample") +
  ggtitle("Ratios of RNA counts with and without a symbol - EVs")

rm(list = c('ratev', 'ratc'))

save.image("pesc_analysis.RData")
### Normalization and scaling -----
## Normalization was done using the SC transform described here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1 as it is supposed to be depth-independent, which justifies its use in Smartseq3 EV sequencing.


### Metadata ----
save.image("pesc_analysis.RData")

# Assign cells to the patient...
full_names <- colnames(ds_cf)
reg.pat <- regmatches(full_names, regexpr("PEsc[0-9]{1,6}", full_names))
reg.pat <- gsub("PEsc", "", reg.pat)
ds_cf@meta.data$pat <- reg.pat

# ...and the replicate...
reg.rep <- regmatches(full_names, regexpr("PEsc[0-9]{1,6}(.|_)[0-9]", full_names))
reg.rep <- gsub("PEsc[0-9]{1,6}(.|_)", "", reg.rep)
ds_cf@meta.data$rep <- reg.rep

# ...and diagnosis and type.
metd <- fread("PEpatients_sorted.csv")
metd[, c("type", "diag") := lapply(.SD, factor), .SDcols = c("type", 'diag')]
metd[, pat_numb := regmatches(metd[, pat_numb], regexpr("[0-9]{3,5}", metd[, pat_numb]))]

patd <- as.data.table(reg.pat)
colnames(patd) <- "pat_numb"
met.full <- metd[patd, on = .(pat_numb)]

ds_cf@meta.data$type <- met.full[, type]
ds_cf@meta.data$diag <- met.full[, diag]

# Assign EVs to the fraction they're from,...
full_names <- colnames(ds_dcf)
reg.frac <- regmatches(full_names, regexpr("_[A-Z]{1,2}_", full_names))
reg.frac <- gsub("_", "", reg.frac)
ds_dcf@meta.data$frac <- reg.frac

# ...their patient...
reg.pat <- regmatches(full_names, regexpr("X[0-9]{1,6}", full_names))
reg.pat <- gsub("X", "", reg.pat)
ds_dcf@meta.data$pat <- reg.pat

# metd[pat_numb %in% reg.pat]

patd <- as.data.table(reg.pat)
colnames(patd) <- "pat_numb"
met.full <- metd[patd, on = .(pat_numb)]

ds_dcf@meta.data$type <- met.full[, type]
ds_dcf@meta.data$diag <- met.full[, diag]


### DimRed ----
## Cells

ds_cf <- SCTransform(ds_cf, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

table(ds_cf@meta.data[c("pat", "diag", "type")])

# VariableFeaturePlot(ds_cf)
FeaturePlot(ds_cf, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 2, reduction = 'umap')

ds_cf <- FindNeighbors(ds_cf, reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(ds_cf, group.by = c("pat", "diag", "type"), pt.size = 2)
p2 <- DimPlot(ds_cf, pt.size = 2)

p1 + p2

# Inspect PCs - Jackstraw doesn't work with SCTransformed data
  
# ds_cf <- JackStraw(ds_cf, num.replicate = 100)
# ds_cf <- ScoreJackStraw(ds_cf, dims = 1:20)
# JackStrawPlot(ds_cf, dims = 1:15)
ElbowPlot(ds_cf)
  
top10 <- head(VariableFeatures(ds_cf), 10)
plot1 <- VariableFeaturePlot(ds_cf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

FeaturePlot(ds_cf, features = c('nFeature_RNA','nCount_RNA', top10), pt.size = 2, reduction = 'umap')

clust_cell <- DimPlot(ds_cf, group.by = c("seurat_clusters", "pat", "diag", "type"))
ggsave("clust_cell_1.png", clust_cell, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")

# Differential expression analysis - cells ----

ds_cf <- PrepSCTFindMarkers(ds_cf)






  ## EVs

ds_dcf <- SCTransform(ds_dcf, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

DimPlot(ds_dcf, group.by = "pat", pt.size = 2)
DimPlot(ds_dcf, group.by = "frac", pt.size = 2)
DimPlot(ds_dcf, group.by = "diag", pt.size = 2)
DimPlot(ds_dcf, group.by = "type", pt.size = 2)
DimPlot(ds_dcf, group.by = c("frac", "pat", "diag", "type"), pt.size = 2)


FeaturePlot(ds_dcf, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 2, reduction = 'umap')
ElbowPlot(ds_dcf, ndims = 30)


top12 <- head(VariableFeatures(ds_dcf), 12)
plot1 <- VariableFeaturePlot(ds_dcf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
FeaturePlot(ds_dcf, features = top12, pt.size = 2, reduction = 'umap')


# histograms and violinplots of overall expression of genes
aggr_exp <- Matrix::rowSums(ds_cf@assays$SCT)
hist(aggr_exp[aggr_exp > 50], breaks = 500)
ggplot(as.data.table(aggr_exp)) +
  geom_histogram()

aggr_exp <- Matrix::rowSums(ds_dcf@assays$SCT)
hist(aggr_exp[aggr_exp > 50], breaks = 500)

frac_diff <-VlnPlot(ds_dcf, features = c('nFeature_RNA','nCount_RNA'), group.by = "frac", log = T)
ggsave("frac_diff_counts.png", frac_diff, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")
### Further EV clustering ----
wd2 <- "~/OneDrive - Karolinska Institutet/Dokument/nord_lab/seq/0423_pesc/plots/"
setwd(wd2)
save.image("../pesc_analysis.RData")

# without FT
ds_dcf1 <- subset(ds_dcf, subset = frac %notin% "FT") %>%
  SCTransform(vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

umap_ftless <- DimPlot(ds_dcf1, group.by = c("frac", "pat", "diag", "type"), pt.size = 2)
ggsave("ev_umap_wo_ft.png", umap_ftless, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")

top10 <- head(VariableFeatures(ds_dcf1), 10)
f_umap_ftless <- FeaturePlot(ds_dcf1, features = c('nFeature_RNA','nCount_RNA', top10), pt.size = 2, reduction = 'umap', ncol = 4)
ggsave("ev_umap_feat_wo_ft.png", f_umap_ftless, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")


# without FT & cells

ds_dcf2 <- subset(ds_dcf, subset = frac %notin% c("FT", "C")) %>%
  SCTransform(vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

umap_ftcless <- DimPlot(ds_dcf2, group.by = c("frac", "pat", "diag", "type"), 
                        pt.size = 2)
ggsave("ev_umap_wo_ftc.png", umap_ftcless, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")

top10 <- head(VariableFeatures(ds_dcf2), 10)
f_umap_ftcless <- FeaturePlot(ds_dcf2, features = c('nFeature_RNA','nCount_RNA', top10), 
                              pt.size = 2, reduction = 'umap', ncol = 3)
ggsave("ev_umap_feat_wo_ftc.png", f_umap_ftcless, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")


# without FT & cells & SP


ds_dcf3 <- subset(ds_dcf2, subset = frac %notin% c("FT", "C", "SP")) %>%
  SCTransform(vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

umap_ftcspless <- DimPlot(ds_dcf3, group.by = c("frac", "pat", "diag", "type"), 
                        pt.size = 2)
ggsave("ev_umap_wo_ftcsp.png", umap_ftcspless, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")

top10 <- head(VariableFeatures(ds_dcf3), 10)
f_umap_ftcspless <- FeaturePlot(ds_dcf3, features = c('nFeature_RNA','nCount_RNA', top10), 
                              pt.size = 2, reduction = 'umap', ncol = 3)
ggsave("ev_umap_feat_wo_ftcsp.png", f_umap_ftcspless, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")



### Dataset integration ----

wd3 <- "C:/Users/radgro/OneDrive - Karolinska Institutet/Dokument/nord_lab/seq/0423_pesc/plots/comb"
setwd(wd3)

ds_cf[["cev"]] <- "cell"
ds_cf[["frac"]] <- "cell"
ds_dcf3[["cev"]] <- "ev"


ds_comb_list <- list(cells = ds_cf, ev = ds_dcf3)
features <- SelectIntegrationFeatures(object.list = ds_comb_list, nfeatures = 6000)
ds_comb_list <- PrepSCTIntegration(object.list = ds_comb_list, anchor.features = features)

pe_anchors <- FindIntegrationAnchors(object.list = ds_comb_list, normalization.method = "SCT", anchor.features = features)
ds_comb <- IntegrateData(anchorset = pe_anchors, normalization.method = "SCT")


ds_comb <- SCTransform(ds_comb, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:3, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

ElbowPlot(ds_comb, ndims = 100)


DimPlot(ds_comb, pt.size = 1.5)
DimPlot(ds_comb, split.by = "cev")
DimPlot(ds_comb, group.by = "cev", pt.size = 2)
umap_comb <- DimPlot(ds_comb, group.by = c("cev", "frac", "diag", "seurat_clusters"), pt.size = 1.5)



# Inspect difference between EVs distant and close to the cells
umap_comb <- DimPlot(subset(ds_comb, seurat_clusters %in% c(6, 2, 4, 9)), 
                     group.by = c("cev", "frac", "diag", "seurat_clusters"), pt.size = 1.5)
ggsave("comb_cl2469_6k.png", umap_comb, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")

ds_comb <- PrepSCTFindMarkers(ds_comb)

# set idents to ev frac and unsupervised cluster
# ds_comb$frac.cl <- paste(ds_comb$) # not needed right now

head(Idents(ds_comb))
ds_comb.m <- FindMarkers(ds_comb, assay = "SCT", ident.1 = 9, ident.2 = 4,
                         verbose = FALSE, recorrect_umi = FALSE)
head(ds_comb.m, 15)
de_ev1 <- rownames(head(ds_comb.m, 6))

f_de_ev1 <- FeaturePlot(ds_comb, features = de_ev1, pt.size = 1.5)
ggsave("f_de_ev1_cl94.png",  f_de_ev1, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")

# Against all other cells/evs
ds_comb.m2 <- FindMarkers(ds_comb, assay = "SCT",   ident.1 = 9,
                         verbose = FALSE, recorrect_umi = FALSE)
head(ds_comb.m2, 15)
de_ev2 <- rownames(head(ds_comb.m2, 6))

f_de_ev2 <- FeaturePlot(ds_comb, features = de_ev2, pt.size = 1.5)
ggsave("f_de_ev1_cl9all.png",  f_de_ev2, device = 'png', dpi = "retina", width = 15, height = 12, bg = "white")



top10 <- head(VariableFeatures(ds_comb), 10)
f_umap_comb <- FeaturePlot(ds_comb, features = c('nFeature_RNA','nCount_RNA', top10),
                                pt.size = 1.5, reduction = 'umap', ncol = 3)

f_umap_comb <- FeaturePlot(ds_comb, features =  top10,
                           pt.size = 1.5, reduction = 'umap', ncol = 3)

f_umap_comb + DimPlot(ds_comb, pt.size = 1.5)

# optimize clustering

ds_comb.3 <- RunUMAP(ds_comb, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:3, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

DimPlot(ds_comb.3, pt.size = 1.5)




# Filter the cells with previously decided cutoffs ----

apf <- subset(x = p_all, subset = nCount_RNA > 1000 & nFeature_RNA > 1000)
apf <- PercentageFeatureSet(apf, pattern = "^MT-", col.name = "percent_mt")
apf <- PercentageFeatureSet(apf, "^HB[^(P)]", col.name = "percent_hb")
apf <- PercentageFeatureSet(apf, "PECAM1|apf4", col.name = "percent_plat")
apf <- PercentageFeatureSet(apf, "^RP[SL]", col.name = "percent_ribo")

# QC plotting of filtered reads ----

p1 <- ggplot(as.data.table(apf$nFeature_RNA)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('Genes') +
  xlab("") +
  ylab("count")

p2 <- ggplot(as.data.table(apf$nCount_RNA)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('RNAs') +
  xlab("") +
  ylab("")

p3 <- ggplot(as.data.table(apf$percent_mt)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of mitochondrial') +
  xlab("") +
  ylab("percent")

p4 <- ggplot(as.data.table(apf$percent_ribo)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of ribosomal') +
  xlab("") +
  ylab("")

p5 <- ggplot(as.data.table(apf$percent_hb)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of hemoglobin') +
  xlab("") +
  ylab("")

p6 <- ggplot(as.data.table(apf$percent_plat)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of platelet') +
  xlab("") +
  ylab("")

qcplots <- ggarrange(p1, p2, p3, p4, p5, p6, nrow = 1)
ggsave("qc_violins_filt_1.png", qcplots, device = 'png', dpi = "retina", width = 21, height = 12, bg = "white")

# Second filtering - ribo, mito, hem ----
apf <- subset(x = apf, subset = nCount_RNA > 1000 & nFeature_RNA > 1000
             & percent_hb < 5 & percent_mt < 25 & percent_ribo > 5 & percent_plat < 1)



# QC plotting of filtered reads ----
p1 <- ggplot(as.data.table(apf$nFeature_RNA)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('Genes') +
  xlab("") +
  ylab("count")

p2 <- ggplot(as.data.table(apf$nCount_RNA)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  scale_y_log10() +
  ggtitle('RNAs') +
  xlab("") +
  ylab("")

p3 <- ggplot(as.data.table(apf$percent_mt)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of mitochondrial') +
  xlab("") +
  ylab("percent")

p4 <- ggplot(as.data.table(apf$percent_ribo)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of ribosomal') +
  xlab("") +
  ylab("")

p5 <- ggplot(as.data.table(apf$percent_hb)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of hemoglobin') +
  xlab("") +
  ylab("")

p6 <- ggplot(as.data.table(apf$percent_plat)) +
  geom_violin(aes(x = "PBMC, BM", y = V1)) +
  geom_jitter(aes(x = "PBMC, BM", y = V1), size = .9, alpha = .4) + 
  ggtitle('Percent of platelet') +
  xlab("") +
  ylab("")

qcplots <- ggarrange(p1, p2, p3, p4, p5, p6, nrow = 1)
ggsave("qc_violins_filt_2.png", qcplots, device = 'png', dpi = "retina", width = 21, height = 12, bg = "white")


feats <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb", "percent_plat")

VlnPlot(apf, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

# RNA vs antibody filtered ----

chk_f <- data.table(colnames(apf@assays$RNA), apf$nCount_RNA, colnames(apf@assays$ab), apf$nCount_ab, apf@active.ident)

rna_ab_f <- ggplot(chk_f[V4 < 10000]) +
  geom_jitter(aes(x = V2, y = V4, color = V5),
              size = 1.5, alpha = .4,
              position = position_jitter(height = .05)) +
  scale_color_manual(values = c("#ffa23e", "green")) +
  xlab("RNA") +
  ylab("antibody tag") +
  theme(
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  ggtitle('RNA read count vs antibody tag count')

ggsave("RNA_vs_ab_filt.png", rna_ab_f, device = 'png', dpi = "retina", width = 16, height = 12, bg = "white")



# Data normalization and clustering ----

apf <- NormalizeData(apf, normalization.method='RC', scale.factor=1e4)

apf <- FindVariableFeatures(apf, nfeatures = 10000)
length(apf@assays$RNA@var.features)

apf <- ScaleData(apf)
apf <- RunPCA(apf, verbose = T, npcs = 20)
apf <- RunUMAP(apf, dims = 1:20, verbose = T)
apf <- RunTSNE(apf)


DimPlot(apf, group.by = "orig.ident", pt.size = 2)
DimPlot(apf, group.by = "orig.ident", pt.size = 2, reduction = "pca")
DimPlot(apf, group.by = "orig.ident", pt.size = 2, reduction = "tsne")

p4f <- subset(apf, idents = "pat4")
p4f <- NormalizeData(p4f, normalization.method='RC', scale.factor=1e4)
p4f <- FindVariableFeatures(p4f, nfeatures = 10000)
length(p4f@assays$RNA@var.features)

p4f <- ScaleData(p4f)
p4f <- RunPCA(p4f, verbose = T, npcs = 30)
p4f <- RunUMAP(p4f, dims = 1:10, verbose = T, n.neighbors = 10)

DimPlot(p4f, pt.size = 2)
DimPlot(p4f, pt.size = 2, reduction = "pca")


# Subset most highly expressed
C <- p4f@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
below5 <- names(which(apply(C, 1, median) <= .05))

p4f <- subset(p4f, features = below5)

p4f <- FindVariableFeatures(p4f, nfeatures = 10000)
p4f <- ScaleData(p4f)
p4f <- RunPCA(p4f, verbose = T, npcs = 30)
p4f <- RunUMAP(p4f, dims = 1:10, verbose = T, n.neighbors = 10)
DimPlot(p4f, pt.size = 2)
DimPlot(p4f, pt.size = 2, reduction = "pca")

# most expressed gene
d <- p4f@assays$RNA@data

most <- order(apply(d, 1, median), decreasing = T)[1:20]
rownames(d)[most]



out <- p4f@assays$RNA@scale.data
out.n <- rownames(out)
# We create the .csv file
write.csv(out, "fisk1.csv")
write.csv(out.n, "fisk1names.csv")



# remove outliers from the PCA plot

outl <- p4f@reductions$pca@cell.embeddings
outl <- as.data.table(outl, keep.rownames = TRUE)

outl[order(PC_1), .(rn, PC_1, PC_2)][1:50]

# remove outliers with PCA_1 lower than -100
keep <- outl[PC_1 > -100, rn]
p4fp <- subset(p4f, cells = keep)

p4fp <- FindVariableFeatures(p4fp, nfeatures = 10000)
feat <- rownames(p4fp)
p4fp <- ScaleData(p4fp, features = feat)
p4fp <- RunPCA(p4fp, verbose = T, npcs = 30, features = VariableFeatures(object = p4fp))
p4fp <- RunUMAP(p4fp, dims = 1:10, verbose = T, n.neighbors = 4)
DimPlot(p4fp, pt.size = 2)
DimPlot(p4fp, pt.size = 2, reduction = "pca")

FeaturePlot(p4fp, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 2, reduction = 'pca')

top10 <- head(VariableFeatures(p4fp), 10)
plot1 <- VariableFeaturePlot(p4fp)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Explore the dimensionality - which PCs are good
p4fp <- JackStraw(p4fp, num.replicate = 100)
p4fp <- ScoreJackStraw(p4fp, dims = 1:20)
JackStrawPlot(p4fp, dims = 1:15)
ElbowPlot(p4fp)

p4fp <- FindNeighbors(p4fp, dims = 1:10)
p4fp <- FindClusters(p4fp, resolution = 0.5)

# Inspect different features on a dim reduc plot
FeaturePlot(p4fp, features = c("percent_plat", "percent_hb", "percent_mt", "percent_ribo", "nCount_RNA", "nFeature_RNA"), pt.size = 2, ncol = 3)

save.image(file = "fiskesund_hsc_10x.RData")

p4fp <- ProjectUMAP(p4fp )
p4fp <- FindClusters(p4fp, reduction.type='pca', dims.use=1:9, resolution=c(0.6,0.7,0.8,0.9,1,0.5), print.output=0, save.SNN=TRUE)



