
setwd("E:/project_data_ssd/hpk-1_dataset")

library(Biobase)
library(openxlsx)
library(gProfileR)
library(ggplot2)
library(VennDiagram)
library(gplots)


DE_hpk1.all <- read.xlsx("DE summary hpk1_US_vs_N2_US.xlsx", sheet = 1)
# replace NA p-values with 1
DE_hpk1.all$DESeq2_fdrpval[is.na(DE_hpk1.all$DESeq2_fdrpval)] <- 1
DE_hpk1 <- DE_hpk1.all[ (abs(DE_hpk1.all$DESeq2_log2FC) >= 1) & (DE_hpk1.all$DESeq2_fdrpval < 0.05)  , ]
DE_hpk1 <- DE_hpk1[ , c("GeneSymbol", "DESeq2_log2FC", "DESeq2_fdrpval", "count_range_HPK-1_US", "count_range_N2_US",
                        "avg_TPM.HPK-1_US", "avg_TPM.N2_US") ]

wb_path <- "E:/public_data_annotation/wormbase/WS277"


# read in VST - TODO: HAVE TO SAVE OUT FROM RDATA OBJECT


#gprofiler_res <- gprofiler(DE_hpk1$GeneSymbol, organism = "celegans")
#write.table(DE_hpk1$GeneSymbol, file = "hpk1_DEGs.txt", quote = FALSE, sep = "\t", row.names = FALSE)


# we're also interested in comparing to the Berber microarray dataset
berber_path <- "F:/project_data/MBP_wormlab_analysis_all_11JAN2018/analysis/hpk_study/GEO DATA/berber"
# original paper DE gene list:
#  note these gene names are molecular names
berber_deg_original <- rbind(read.xlsx(file.path(berber_path, "paper_gene_list.xlsx"), sheet = 1)[,c(1,3,4)],
      read.xlsx(file.path(berber_path, "paper_gene_list.xlsx"), sheet = 2))
# my reanalysis of this dataset:
berber_deg_reanalysis <- read.xlsx(file.path(berber_path,"output", "DE results - berber 2016 p01 (nofdr).xlsx"), sheet = 1)


# ------Back to some basics----
#
#

# volcano plot for DE genes with FDR 0.05 and 1 LFC cutoff

volc <- ggplot(data = as.data.frame(DE_hpk1.all), 
               aes(x = DESeq2_log2FC, 
                   y = -log10(DESeq2_fdrpval),
                   colour = as.factor((abs(DESeq2_log2FC) >= 1) & (DESeq2_fdrpval < 0.05)))) +
  geom_point(alpha=0.4, size = 1.75) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot- hpk-1 US vs N2 US") +
  xlab("log2 foldchange")+
  ylab("-log10(adjusted pvalue)")+
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text( size = 14)) +
  geom_hline(yintercept = -log10(0.05), colour = "red", linetype = "dashed")+
  geom_vline(xintercept = 1, colour = "red", linetype = "dotdash")+
  geom_vline(xintercept = -1, colour = "red", linetype = "dotdash")
ggsave(plot = volc, filename = "DE volcano plot for hpk-1 US vs N2 US.pdf", device = "pdf", width = 8, height = 8)
rm(volc)

# DE counts
sum(DE_hpk1$DESeq2_log2FC > 0)
sum(DE_hpk1$DESeq2_log2FC < 0)

#  overlap with basic functional lists of interest
degs_annot <- DE_hpk1

# THIS IS MY OLD LIST
lifespan_genes <- read.xlsx(file.path("E:/project_data_ssd/myc_project/rna-seq_spelunking",
                                      "lifespan_genes", "aggregated lifespan genes (ABC v1).xlsx"))
# NEW GEROGENE LIST
lifespan_genes_v2 <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/gerogenes", 
                    "Gerogene table with orthologs (RNAi and var combined)- Wormbase WS282.xlsx"))

# and some new additions from recent papers that haven't been incorporated into databases
lifespan_genes_additions <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/gerogenes", 
                                         "gerogenes_our_manual_additions.xlsx"))

protNet <- read.table(file.path("F:/project_data/MBP_wormlab_analysis_all_11JAN2018/analysis/hpk_study/docs from andy/proteostasis network",
                                "proteostasis network full from paper.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)[ , c(-1,-2) ]

wTF <- read.xlsx(file.path("E:/public_data_annotation", "wTF-3.0.xlsx"))

degs_annot <- merge(degs_annot, lifespan_genes, by.x = "GeneSymbol", by.y = "gene", all.x = TRUE)
degs_annot <- merge(degs_annot, protNet, by.x = "GeneSymbol", by.y = "Gene.name", all.x = TRUE)
degs_annot$is_TF <- degs_annot$GeneSymbol %in% wTF$Public.name

# how many of the proteostasis network genes have been renamed?
# the vast majority are still found by gene name
sum(!(protNet$Gene.name %in% DE_hpk1.all$GeneSymbol))

# venn with berber
pdf("venn- samuelson hpk-1 vs berber DEGs.pdf", width = 10, height = 8)
# all DE
plot.new()
grid.draw(venn.diagram(
  x=
    list("Samuelson RNA-Seq-\nAll DE" = DE_hpk1$GeneSymbol,
         "Berber Arrays-\nAll DE" = berber_deg_reanalysis$GeneSymbol),
  filename = NULL,
  fill = c("cornflowerblue", "tomato"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1.5,
  cat.col = c("darkblue",  "tomato3"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.05)
)

# upregulated
plot.new()
grid.draw(venn.diagram(
  x=
    list("Samuelson RNA-Seq-\nupregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "Berber Arrays-\nupregulated" = berber_deg_reanalysis$GeneSymbol[berber_deg_reanalysis$logFC > 0]),
  filename = NULL,
  fill = c("cornflowerblue", "tomato"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1.5,
  cat.col = c("darkblue",  "tomato3"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.05)
)

graphics.off()


# functional enrichment for shared DE genes in same direction
berberShared <- c(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],  berber_deg_reanalysis$GeneSymbol[berber_deg_reanalysis$logFC > 0]),
                  intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],  berber_deg_reanalysis$GeneSymbol[berber_deg_reanalysis$logFC < 0]))

berberShared[berberShared %in% lifespan_genes$gene]


gp <- gprofiler(query = list("rnaseq_berber_shared" = berberShared),
                             organism = "celegans", evcodes = TRUE,
                             src_filter = c("GO", "KEGG", "REAC"), correction_method = "gSCS",
                             max_set_size = 1000, max_p_value = 0.05)


gp <- gprofiler(query = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                organism = "celegans", evcodes = TRUE,
                src_filter = c("GO", "KEGG", "REAC"), correction_method = "gSCS",
                max_set_size = 1000, max_p_value = 0.05)
write.table(gp, file = "gprofiler- RNAseq upregulated genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Dillin lab hsf-1 dataset
dillin_hsf1 <- read.xlsx(file.path("F:/project_data/MBP_wormlab_analysis_all_11JAN2018/analysis/hpk_study/GEO DATA/dillin_lab/output",
                                   "DE results - sy441 vs N2.xlsx"))
sum(dillin_hsf1$logFC > 0)
sum(dillin_hsf1$logFC < 0)
hsf1Shared <- c(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],  dillin_hsf1$GeneSymbol[dillin_hsf1$logFC > 0]),
                  intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],  dillin_hsf1$GeneSymbol[dillin_hsf1$logFC < 0]))


# comparison plot
library(MASS)
datMerge <- merge(DE_hpk1.all[,c(1:2)], dillin_hsf1[,c(2,4)], by.x = 1, by.y = 1)
png(file.path("corrplot- dillin hsf1 dataset compare.png"), width=800, height=600)
plot(datMerge$DESeq2_log2FC, datMerge$logFC, xlab="hpk-1 RNA-Seq logfc", ylab = "hsf-1 sy441 atrray logfc")
abline(rlm(datMerge$logFC ~ datMerge$DESeq2_log2FC), col = "red")

cortest <- cor.test(x = datMerge$DESeq2_log2FC, y = datMerge$logFC, alternative="greater", method = "pearson")
legend(x = "topleft", legend = c(paste("Pearson corr:", signif(cortest$estimate,5)),
                                 paste("Pearson corr p-value:", signif(cortest$p.value,5))), bty="n")
rm(cortest)
graphics.off()

# heatmap for neuropeptides?



# ------Neuron Expression--------------------
#
#



# look at expression information from WormBase
#  The expression ontology is of course hierarchical!
source("OBOparse.r")
obo_text <- readLines(file.path(wb_path, "anatomy_ontology.WS277.obo"))
# note: this prints a ton of stuff to stdout as it parses
obo_parsed <- readOboFile(obo_text, cols = c("name", "id")) 
rm(obo_text)
gc()
gene_anatomy <- read.table(file.path(wb_path, "anatomy_association.WS277.wb"), sep = "\t",
                           quote = "", stringsAsFactors = FALSE, comment.char = "!")
table(gene_anatomy$V13) # all celegans
table(gene_anatomy$V4)


# broadest neuron term in pheno ontology is WBbt:0003679
# probably have to collapse the neuron_related ontology terms to get all neuron-expressed genes?
#   or use the API, which I did once before.
gene_anatomy.neuron_expressed <- gene_anatomy[ gene_anatomy$V5 == "WBbt:0003679" , ]
length(unique(gene_anatomy.neuron_expressed$V3))
table(gene_anatomy.neuron_expressed$V4)
gene_anatomy.non_neuron <- gene_anatomy[ !(gene_anatomy$V3 %in% gene_anatomy.neuron_expressed$V3) , ]
length(unique(gene_anatomy.non_neuron$V3))

unique(gene_anatomy[-which(gene_anatomy$V5 == "WBbt:0003679") , "V3" ]) # genes with expression in tissues outside neurons
onlyNeuron <- setdiff(unique(gene_anatomy.neuron_expressed$V3),
        unique(gene_anatomy[-which(gene_anatomy$V5 == "WBbt:0003679") , "V3" ]))


intersect(DE_hpk1$GeneSymbol, onlyNeuron)
intersect(DE_hpk1$GeneSymbol, gene_anatomy.neuron_expressed$V3)


#
# ---> DUE TO MEMORY ISSUES- WILL DEAL WITH THESE ONE AT A TIME! <----
# 


# # Load data from CenGen
# # so these are cell x feature, 52412 cells x 20842 features total for neurons
# # and 65450 cells x 21393 features for all other cell types
# neurons <- readRDS("081519_L4_neuron_cds.rds")
# allCells <- readRDS("081519_L4_all_cells_cds.rds")
# 
# # what other relevant info is here?
# head(neurons$total_counts_Mito) # presumably they already filtered out cells with high mito counts but maybe not?
# summary(neurons$total_counts) # some of these are pretty darn low for total and expressed
# head(neurons$num_genes_expressed) # this seems to be based just on counts > 0
# table(neurons$Neuron.type) # neuron cell identity
# table(neurons$ident) # neuron type in terms of neurotransmitters
# table(neurons$Experiment) # related to neuron ident- label for flow used to pull out the type of neuron
# table(neurons$outlier) # some of these are "true" for outlier- remove? what is this based on?
# head(neurons$total_counts_endogenous) # this is total counts - mito counts
# 
# table(neurons$Tissue.type)
# #table(allCells$Tissue.type)
# 
# # ?
# neurons$ident

# want to try to summarize this into what cell types have expresion of which genes
# this information somehow is NOT available in "bulk" format from CenGen
# we'll call a count greater than 1 as expressed??
# to deal with memory limitations, we'll preallocate a matrix, and then replace values, rather than
#   reassigning to the variable
# binning will be done per-sample based on quantiles; as the data is of course quite sparse this will be done with excluding
#   the 0-counts. 0 = no counts, 1 = possible expression, 2 = low/med expression, 3 = expression
# it seems like there might be some cells with absurdly low expression across the board, and maybe these should be filtered out?
#  (there's one for example with a max count of 2 across all genes) it would probably be good to summarize that
#  if necessary, since clearly the single cell dataset object is able to handle a very large matrix
#   without memory issues, we could duplicate that object and reassign the expression values with the binned values

# gc()
# expressed_neurons <- mat.or.vec(dim(exprs(neurons))[1], dim(exprs(neurons))[2])
# 
# # this is going to take a while.
# # might need to try without sapply if there are still memory issues
# for(i in 1:ncol(expressed_neurons)){
#   qs = quantile(exprs(neurons)[,i][exprs(neurons)[,i] > 0], probs = c(0.5, 0.9))
#   expressed_neurons[,i] <- sapply(exprs(neurons[,1]), function(x){
#     if(x == 0){0}
#     else if((x > 0) && (x <= qs[1])){1}
#     else if((x > qs[1]) && (x <= qs[2])){2}
#     else if(x > qs[2]){3}
#     })
#   rm(qs)
#   gc()
# }
# rm(i)
# gc()
# save.image()
# class(expressed_neurons)
# typeof(expressed_neurons)
# object.size(expressed_neurons)
# mode(expressed_neurons) <- "integer"
# object.size(expressed_neurons)
# gc()

#save(expressed_neurons, file = "cengen081519_L4_neuron_expressed_summary_binned_ABC.RData")
# load(file = "cengen081519_L4_neuron_expressed_summary_binned_ABC.RData") # object is expressed_neurons


# "all cells" set 

# expressed_allCells <- mat.or.vec(dim(exprs(allCells))[1], dim(exprs(allCells))[2])
# mode(expressed_allCells) <- "integer"
# 
# # this is (also) going to take a while.
# for(i in 1:ncol(expressed_allCells)){
#   qs = quantile(exprs(allCells)[,i][exprs(allCells)[,i] > 0], probs = c(0.5, 0.9))
#   expressed_allCells[,i] <- sapply(exprs(allCells[,1]), function(x){
#     if(x == 0){0}
#     else if((x > 0) && (x <= qs[1])){1}
#     else if((x > qs[1]) && (x <= qs[2])){2}
#     else if(x > qs[2]){3}
#   })
#   rm(qs)
# }
# rm(i)
# 
# save(expressed_allCells, file = "cengen081519_L4_allcells_expressed_summary_binned_ABC.RData")
# load(file = "cengen081519_L4_allcells_expressed_summary_binned_ABC.RData") # object is expressed_allCells
# 
# # NEED TO NAME THE DIMENSIONS!
# 
# 
# gc()

#  summarize further by tissue/cell type to build gene sets and reduce the dimensionality
# SHOULD ALSO FILTER OUT CELLS that had a very low number of total genes, and the remaining ones
#  that had a pretty high number of mito read counts- it seemed like there might still be some in here.


#
# Re-re-revisiting this, Sep 2022
#


# DE summary N2_S_vs_N2_US.xlsx
# DE summary hpk1_S_vs_N2_S.xlsx
# DE summary hpk1_S_vs_hpk1_US.xlsx

DE.all <- list()


DE.all$N2_S_vs_N2_US <- read.xlsx("DE summary N2_S_vs_N2_US.xlsx", sheet = 1)
DE.all$hpk1_S_vs_hpk_US <- read.xlsx("DE summary hpk1_S_vs_hpk1_US.xlsx", sheet = 1)
DE.all$hpk1_S_vs_N2_S <- read.xlsx("DE summary hpk1_S_vs_N2_S.xlsx", sheet = 1)


# replace NA p-values with 1
DE.all <- lapply(DE.all, function(x){
  x$DESeq2_fdrpval[is.na(x$DESeq2_fdrpval)] <- 1
  x <- x[ , grepl("GeneSymbol|DESeq2_log2FC|DESeq2_fdrpval|count_range_*|avg_TPM", colnames(x)) ]
  return(x)
})

lapply(DE.all, nrow)

DE.degs <- lapply(DE.all, function(x){
  return(x[ (abs(x$DESeq2_log2FC) >= 1) & (x$DESeq2_fdrpval < 0.05)  , ])
})

DE.all$hpk1_US_vs_N2_US <- DE_hpk1.all
DE.degs$hpk1_US_vs_N2_US <- DE_hpk1

lapply(DE.all, nrow)
lapply(DE.degs, nrow)

# actually let's get those WBGeneIDs in there after all...
gtf <- rtracklayer::import.gff(file.path("E:/public_data_annotation/ensembl/celegans/82", "Caenorhabditis_elegans.WBcel235.82.gtf"), 
           format="gtf", genome="WBcel235",feature.type="gene")
gtf <- data.frame(gtf)

DE_hpk1 <- cbind("WBGeneID" = gtf$gene_id[match(DE_hpk1$GeneSymbol, gtf$gene_name)], DE_hpk1)
DE_hpk1.all <- cbind("WBGeneID" = gtf$gene_id[match(DE_hpk1.all$GeneSymbol, gtf$gene_name)], DE_hpk1.all)

#
# OK! So, it seems somehow more genes change expression substantially in hpk-1 with stress than N2.
#  What of the N2 ones also change expression in hpk-1 though?
#

pdf("venn- N2 stress responsive vs hpk-1 null stress responsive.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("N2 stress responsive, down" = with(DE.degs$N2_S_vs_N2_US, GeneSymbol[DESeq2_log2FC < 0]),
         "N2 stress responsive, up" = with(DE.degs$N2_S_vs_N2_US, GeneSymbol[DESeq2_log2FC > 0]),
         "hpk-1 null stress responsive, down" = with(DE.degs$hpk1_S_vs_hpk_US, GeneSymbol[DESeq2_log2FC < 0]),
         "hpk-1 null stress responsive, up" = with(DE.degs$hpk1_S_vs_hpk_US, GeneSymbol[DESeq2_log2FC > 0])),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold", "grey66"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4", "grey33"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()


pdf("venn- N2 vs Hpk-1 stress vs baseline.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 vs N2, baseline, down" = with(DE.degs$hpk1_US_vs_N2_US, GeneSymbol[DESeq2_log2FC < 0]),
         "hpk-1 vs N2, baseline, up" = with(DE.degs$hpk1_US_vs_N2_US, GeneSymbol[DESeq2_log2FC > 0]),
         "hpk-1 vs N2, stress, down" = with(DE.degs$hpk1_S_vs_N2_S, GeneSymbol[DESeq2_log2FC < 0]),
         "hpk-1 vs N2, stress, up" = with(DE.degs$hpk1_S_vs_N2_S, GeneSymbol[DESeq2_log2FC > 0])),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold", "grey66"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4", "grey33"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()

#
# overlap with mxl-2 DE results??
#
de_mxl2 <- read.xlsx(file.path("C:/Users/Adam Cornwell/Box/adam/DR myc manuscript reboot/code/output/DE/MXL-2_EV_vs_N2_EV",
                               "DE_results_2022_MXL-2_EV_vs_N2_EV.xlsx"), sheet = 3)

pdf("venn- N2 vs Hpk-1 baseline vs mxl-2 DE.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 vs N2, baseline, down" = with(DE.degs$hpk1_US_vs_N2_US, GeneSymbol[DESeq2_log2FC < 0]),
         "hpk-1 vs N2, baseline, up" = with(DE.degs$hpk1_US_vs_N2_US, GeneSymbol[DESeq2_log2FC > 0]),
         "mxl-2 vs N2, down" = with(de_mxl2, GeneSymbol[log2FoldChange < 0]),
         "mxl-2 vs N2, up" = with(de_mxl2, GeneSymbol[log2FoldChange > 0])),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold", "grey66"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4", "grey33"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()


#
# load TPM expression for these samples,
#  going to want that for some plots
#
tmpFiles <- list.files(path = "F:/project_data/MBP_wormlab_analysis_all_11JAN2018/analysis/hpk_study/rsem_gene", pattern = "_US_", full.names = TRUE)
tpm <- lapply(tmpFiles, function(x){
  read.table(x, sep = "\t", header = TRUE)[,c("gene_id", "TPM")]
})
names(tpm) <- gsub("F:/project_data/MBP_wormlab_analysis_all_11JAN2018/analysis/hpk_study/rsem_gene/clt_|_rsem.out.genes.results", "", tmpFiles)
# gene ID vectors identical? They should be.
identical(tpm$`HPK-1_US_64`$gene_id, tpm$N2_US_63$gene_id)

tmpGeneVec <- tpm$`HPK-1_US_64`$gene_id # any is fine
tpm <- sapply(names(tpm), function(x){tpm[[x]][,2]})
rownames(tpm) <- tmpGeneVec
rm(tmpFiles, tmpGeneVec)

# subset to relevant genes
tpm <- merge(tpm, wb_gene_info[,c("WBGeneID", "gene_name")], by.x = 0, by.y = 1)
# a few gene names have been updated...
tpm <- tpm[ tpm$gene_name %in% c(DE.all$N2_S_vs_N2_US$GeneSymbol, "pals-1", "pals-32", "exc-14")  , ]
row.names(tpm) <- tpm[,1]

tpm_wbid <- as.matrix(tpm[ , -c(1, ncol(tpm)) ])
tpm_sym <- tpm
rownames(tpm_sym) <- tpm_sym[,ncol(tpm_sym)]
tpm_sym <- as.matrix(tpm_sym[ , -c(1, ncol(tpm_sym)) ])
rm(tpm)

# build a simple sample metadata matrix
sampleMetadata <- data.frame(rbind(
  c(SampleName = "N2_US_63", SampleGroup = "N2_US", SampleColor = "#6AD885"),
  c(SampleName = "N2_US_67", SampleGroup = "N2_US", SampleColor = "#6AD885"),
  c(SampleName = "N2_US_71", SampleGroup = "N2_US", SampleColor = "#6AD885"),
  c(SampleName = "HPK-1_US_64", SampleGroup = "HPK1_US", SampleColor = "#DC8E36"),
  c(SampleName = "HPK-1_US_68", SampleGroup = "HPK1_US", SampleColor = "#DC8E36"),
  c(SampleName = "HPK-1_US_72", SampleGroup = "HPK1_US", SampleColor = "#DC8E36")
))


source("E:/project_data_ssd/heatmap.3.5.R") # read in heatmap.3 function
library(viridis)
#
# load genage lifespan annotation that I recurated a bit
#

genage_simple <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/genage", "genage_edited_gerogene_categories.xlsx"))
colnames(genage_simple)[1] <- "GeneSymbol"


DE_hpk1_gerogenes <- merge(DE_hpk1, genage_simple, by = "GeneSymbol")
# A few of these gene names have been updated
# C31B8.4 is now pals-32
# F15D3.8 is now pals-1
# K11D12.9 is now exc-14
DE_hpk1_gerogenes$GeneSymbol[DE_hpk1_gerogenes$GeneSymbol == "C31B8.4" ] <- "pals-32"
DE_hpk1_gerogenes$GeneSymbol[DE_hpk1_gerogenes$GeneSymbol == "F15D3.8" ] <- "pals-1"
DE_hpk1_gerogenes$GeneSymbol[DE_hpk1_gerogenes$GeneSymbol == "K11D12.9" ] <- "exc-14"

# going to try to make a heatmap showing expression across samples for gerogenes,
#  along with row-side annotation for DE direction and gerogene category.
heatmapData <- tpm_sym[ DE_hpk1_gerogenes$GeneSymbol , ]
identical(rownames(heatmapData), DE_hpk1_gerogenes$GeneSymbol) # should be true
# need as shifted log2
heatmapData <- log2(heatmapData + 1)

tmpColData <- as.matrix(sampleMetadata$SampleColor)

tmpRowData <- cbind(
  # gerogene color vector- assumes only two conditions
  gerogene_category = ifelse(DE_hpk1_gerogenes$genage_gerogene_cats == "progeric_gene", "plum4", "wheat"),
  # DE direction color vector
  FC_direction = ifelse(DE_hpk1_gerogenes$DESeq2_log2FC > 0, "red3", "royalblue3")
)

pdf("hpk-1 and N2 basal expression for DE gerogenes.pdf", width = 10, height = 17)
heatmap.3(x = heatmapData, 
          ColSideColors = tmpColData, 
          RowSideColors = t(tmpRowData),
          dendrogram = "row",
          col = viridis(option = "plasma", n = 64),
          cexCol = 1.2, cexRow = 1, margins = c(12,10), keysize = 0.7, KeyValueName = "log2 (TPM + 1)")
graphics.off()

rm(heatmapData, tmpColData, tmpRowData)

#
# LOAD SUMMARIZED CENGEN DATA
#

cengen_sum_level2 <- read.csv(file.path("E:/public_data_annotation/c_elegans_other/cengen_summarized", "021821_medium_threshold2.csv"))
rownames(cengen_sum_level2) <- cengen_sum_level2$gene_name
cengen_sum_level2 <- cengen_sum_level2[ , -c(1,2,3) ]
cengen_sum_level2 <- as.matrix(cengen_sum_level2)

hist(cengen_sum_level2["hpk-1",], main = "hpk-1 expression level across neurons in Cengen SC", breaks = 20)
sort(cengen_sum_level2["hpk-1",], decreasing = TRUE)[1:10]

# serotonin- tph-1 
# gabaergic- unc-47
# glutamatergic- eat-4
# dopaminergic- cat-2
# cholinergic- unc-17
# thermosensory- ttx-3

# filter step one- how many neurons have hpk-1 expressed?

cengen_sum_level2.hpk1 <- cengen_sum_level2[  , cengen_sum_level2["hpk-1",] > 0 ]
ncol(cengen_sum_level2)
ncol(cengen_sum_level2.hpk1)
setdiff(colnames(cengen_sum_level2), colnames(cengen_sum_level2.hpk1))

# which of the DE genes have expression in any neurons?

DE.degs_neuron_expressed <- lapply(DE.degs, function(x){
  
  return(x[x$GeneSymbol %in% rownames(cengen_sum_level2.hpk1) , ])
  
})

lapply(DE.degs, nrow)
lapply(DE.degs_neuron_expressed, nrow)


dim(cengen_sum_level2[ , cengen_sum_level2["tph-1",] > 0 ])

dim(cengen_sum_level2[ , cengen_sum_level2["unc-47",] > 0 ])

dim(cengen_sum_level2[ , cengen_sum_level2["eat-4",] > 0 ])

dim(cengen_sum_level2[ , cengen_sum_level2["cat-2",] > 0 ])

dim(cengen_sum_level2[ , cengen_sum_level2["unc-17",] > 0 ])

hpk1_neuron_type_expr <- Reduce(f = rbind, x = lapply(c("tph-1", "unc-47", "eat-4", "cat-2", "unc-17", "ttx-3"), function(x){
  data.frame(hpk1_expr = cengen_sum_level2[ "hpk-1" , cengen_sum_level2[x,] > 0 ], marker = x)
}))
# add "everything else" 
hpk1_neuron_type_expr = rbind(hpk1_neuron_type_expr, data.frame("hpk1_expr" = cengen_sum_level2["hpk-1", !(colnames(cengen_sum_level2) %in% rownames(hpk1_neuron_type_expr))],
           "marker" = "all other neurons"))

# concat neuron type info (er concat by replacement, because there aren't many categories)
hpk1_neuron_type_expr$marker <- gsub("tph-1", "serotonergic (tph-1)", hpk1_neuron_type_expr$marker)
hpk1_neuron_type_expr$marker <- gsub("unc-47", "GABAergic (unc-47)", hpk1_neuron_type_expr$marker)
hpk1_neuron_type_expr$marker <- gsub("eat-4", "glutamatergic (eat-4)", hpk1_neuron_type_expr$marker)
hpk1_neuron_type_expr$marker <- gsub("cat-2", "dopaminergic (cat-2)", hpk1_neuron_type_expr$marker)
hpk1_neuron_type_expr$marker <- gsub("unc-17", "cholinergic (unc-17)", hpk1_neuron_type_expr$marker)
hpk1_neuron_type_expr$marker <- gsub("ttx-3", "thermosensory (ttx-3)", hpk1_neuron_type_expr$marker)

ggplot(data = hpk1_neuron_type_expr, aes(x = hpk1_expr)) +
  geom_histogram(bins = 20) +
  facet_wrap(vars(marker)) +
  theme_bw(base_size = 20) +
  xlab("hpk-1 expression") + ylab("# of neurons")

#
# As defined by these markers, some neurons fall into multiple categories
#  Try to make an UpSet plot to look at the size of these overlaps
#  This should be a simple case, there are only five groups.
#
library(UpSetR)

markerExprList <- list("serotonergic(tph-1)" = colnames(cengen_sum_level2[ , cengen_sum_level2["tph-1",] > 100 ]),
                       "GABAergic(unc-47)" = colnames(cengen_sum_level2[ , cengen_sum_level2["unc-47",] > 100 ]),
                       "glutamatergic(eat-4)" = colnames(cengen_sum_level2[ , cengen_sum_level2["eat-4",] > 100 ]),
                       "dopaminergic(cat-2)" = colnames(cengen_sum_level2[ , cengen_sum_level2["cat-2",] > 100 ]),
                       "cholinergic(unc-17)" = colnames(cengen_sum_level2[ , cengen_sum_level2["unc-17",] > 100 ]))

upset(fromList(markerExprList), order.by = "freq", text.scale = 1.5)


#
# General functional enrichment (running again to update)
#
detach("package:gProfileR", unload = TRUE)
library(gprofiler2)

# note- user-defined gmt files can also be uploaded now
# with upload_GMT_file(gmtfile = "")

hpk1_de_gprofiled <- gost(query = list("hpk-1 null up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                  "hpk-1 null down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]),
     organism = "celegans", evcodes = TRUE, domain_scope = "custom", custom_bg = DE_hpk1.all$GeneSymbol,
     sources = c("GO", "KEGG", "REAC"))

View(hpk1_de_gprofiled$result)

gostplot(gostres = hpk1_de_gprofiled)

# these are useful, allows highlighting certain terms of interest
publish_gostplot()
publish_gosttable()

# they also provide a suggestion for getting GProfiler results into a form that
#  can be loaded for visualization with EnrichmentMap in Cytoscape-
#  see "reating a Generic Enrichment Map (GEM) file for EnrichmentMap" section in 
#  the vignette https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html


#
# Maybe use a treemap?
#
library(treemap)

# Helper function to deal with text wrapping in treemaps
# Also modified to transform p-values to -log10
wrapify <- function(myDat = NULL, pvalCol = "Adjusted.P.value", term="Term"){
  myDat[,term] <- sapply(myDat[,term], function(x){
    paste0(strwrap(gsub("_", " ", x), width = 22, 
                   prefix = "\n", initial = ""), collapse = "")})
  myDat$pTransform <- -log10(as.numeric(myDat[,pvalCol]))
  return(myDat)
}


# ---- gprofiler results ----
unique(hpk1_de_gprofiled$result$query)
gp_res <- list(hpk1_up = hpk1_de_gprofiled$result[hpk1_de_gprofiled$result$query == "hpk-1 null up" , ], 
               hpk1_down = hpk1_de_gprofiled$result[hpk1_de_gprofiled$result$query == "hpk-1 null down" , ])
gp_res <- lapply(gp_res, function(x){
  x$source <- factor(x$source, levels = c("GO:BP","GO:MF","GO:CC","KEGG"));
  return(x);
})

pdf("Treemap- significant gprofiler results.pdf",
    width = 12, height = 8)
for(currRes in names(gp_res)){
  treemap(dtf = wrapify(gp_res[[currRes]], "p_value", "term_name"),
          index = c("source", "term_name"),
          vSize = "pTransform",
          vColor = "source",
          type = "categorical", 
          force.print.labels = FALSE,
          aspRatio = 1.5,
          fontsize.labels = c(0, 8),
          position.legend = "right", 
          fontsize.legend = 6,
          drop.unused.levels = FALSE,
          title = paste0("Significant gprofileR results\n", currRes))
}
graphics.off()
rm(currRes)

#
# Which significant terms have genes that are enriched or uniquely expressed in neurons?
#

kaletsky_unique_neurons <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/murphylab_tissue_expression",
                                               "table s9- tissue unique genes.xlsx"), sheet = "neurons")[,"Gene"]

kaletsky_enriched_neurons <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/murphylab_tissue_expression",
                                               "table S8- tissue enriched genes.xlsx"), sheet = "neurons")[,"Public.Name"]

length(kaletsky_unique_neurons)
length(kaletsky_enriched_neurons)
length(union(kaletsky_enriched_neurons, kaletsky_unique_neurons))

neuron_unique_enriched <- union(kaletsky_enriched_neurons, kaletsky_unique_neurons)

# what fraction of the DE genes overall overlap with the neuron enriched/specific genes?
length(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0], neuron_unique_enriched))
length(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0])

length(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0], neuron_unique_enriched))
length(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0])
intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0], neuron_unique_enriched)

hpk1_up_neuron_assoc <- gp_res$hpk1_up[ (sapply(gp_res$hpk1_up$intersection, function(gene_isect){
  sum(strsplit(gene_isect, split = ",")[[1]] %in% neuron_unique_enriched)
}) / gp_res$hpk1_up$intersection_size) > 0.5 , ]

pdf("Treemap- significant gprofiler results- neuron associated.pdf",
    width = 12, height = 8)
treemap(dtf = wrapify(hpk1_up_neuron_assoc, "p_value", "term_name"),
        index = c("source", "term_name"),
        vSize = "pTransform",
        vColor = "source",
        type = "categorical", 
        force.print.labels = FALSE,
        aspRatio = 1.5,
        fontsize.labels = c(0, 8),
        position.legend = "right", 
        fontsize.legend = 6,
        drop.unused.levels = FALSE,
        title = paste0("Significant gprofileR results\n", "upregulated, neuron-associated"))
graphics.off()

#
# How many unique genes represented in the enrichment overlap?
#
enrich_overlap_unique <- unique(as.vector(unlist(sapply(hpk1_up_neuron_assoc$intersection, function(gene_isect){
  as.vector(strsplit(gene_isect, split = ",")[[1]])
}))))
length(enrich_overlap_unique)
sum(enrich_overlap_unique %in% neuron_unique_enriched)

# export for enrichmentmap
gem <- hpk1_up_neuron_assoc[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
write.table(gem, file = "hpk-1 up neuron-associated enrichment - gProfiler_gem.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(hpk1_up_neuron_assoc[,c("term_id", "source")], file = "hpk-1 up neuron-associated enrichment addtl info.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


#
# OK great, we've got a reasonable looking term map which I've grouped by some broad themes
#  Now to build a table of the associated genes, with the foldchange in our dataset, and some additional annotation
# A gene may appear in multiple groups, but we only want it to show up once per group. The gene entry will also include the
#   terms in the group that it's included in.
# Annotation to include:
#   - lifespan association
#   - human ortholog
#   - protein domains
#   - nervous system gene category, per Oliver Hobert's wormbook tables
#   - associated enriched GO terms
#   - PGP?
#   - possible binding by TFs of interest?
#

ortholist <- read.xlsx(file.path("E:/public_data_annotation/homology/celegans", "ortholist_master.xlsx"))
# collapse into a separated string for multiple human symbols per WBGeneID
ortholist_collapsed <- Reduce(x = lapply(unique(ortholist$WormBase.ID), function(geneID){
  data.frame(Human_GeneSymbol = paste(na.omit(ortholist$HGNC.Symbol[ortholist$WormBase.ID == geneID]),
                                      collapse = ","),
             Human_EnsemblID = paste(na.omit(ortholist$Ensembl.ID[ortholist$WormBase.ID == geneID]),
                                     collapse = ","))
}), f = rbind)
rownames(ortholist_collapsed) <- unique(ortholist$WormBase.ID)
# merge in c elegans "gene names"
ortholist_collapsed <- merge(ortholist_collapsed, wb_gene_info[,c("WBGeneID", "gene_name")], by.x = 0, by.y = "WBGeneID", all.x = TRUE)

# PGP
pgpList <- read.xlsx(file.path("C:/Users/Adam Cornwell/Box/adam/DR myc manuscript reboot/omics_data_sources/genesets",
                               "PGP List and Addresses.xlsx"))
pgpList <- pgpList$WBID[!is.na(pgpList$WBID)]
pgpList <- merge(pgpList, wb_gene_info[,c("WBGeneID", "gene_name")], by.x = 1, by.y = "WBGeneID")

# protein domains. Why in the world did this file get a csv extension when it's tab-delimited??
annot_prot_domains <- read.table(file.path(paths$wormbase, "c_elegans.PRJNA13758.WS284.protein_domains.csv"),
                                 sep = "\n")

# how many split records per line? Do all have at least three? YES - so every line is at least:
# WBGeneID -> gene_name -> proteinID -- and everything after is some kind of domain annotation.
table(sapply(1:nrow(annot_prot_domains), function(x){length(unlist(strsplit(annot_prot_domains[x,], "\t")))}))
# Find all the unique entries that include some kind of description
prot_domain_desc <- na.omit(unlist(sapply(1:nrow(annot_prot_domains), function(x){
  splitted = unlist(strsplit(annot_prot_domains[x,], "\t"))
  if(any(grepl(" ", splitted))){
    return(unique(splitted[grepl(" ", splitted)]))
  }
})))
# split up by source
prot_domain_desc_interpro <- unique(na.omit(sapply(1:length(prot_domain_desc), function(x){
  splitted = unlist(strsplit(prot_domain_desc[x], ":"))
  if(splitted[1] == "INTERPRO"){return(prot_domain_desc[x])}else{return(NA)}
})))
# now build the gene associations
# interpro
protdomain_interpro_long <- Reduce(lapply(1:nrow(annot_prot_domains), function(x){
  splitted = unlist(strsplit(annot_prot_domains[x,], "\t"))
  gene_id = splitted[1]
  protein_id = splitted[3]
  if(length(splitted) > 3){
    if(any(prot_domain_desc_interpro %in% splitted[4:length(splitted)])){
      Reduce(lapply(prot_domain_desc_interpro[prot_domain_desc_interpro %in% splitted[4:length(splitted)]], function(z){
        c(gene_id, protein_id, z)
      }), f = rbind)
    }
  }
}), f = rbind)
protdomain_interpro_long <- data.frame(protdomain_interpro_long)
colnames(protdomain_interpro_long) <- c("gene_id", "protein_id", "domain_desc")
rownames(protdomain_interpro_long) <- NULL


library(ontologyIndex)
obo_parsed <- get_OBO(file.path("E:/public_data_annotation/wormbase/WS284", "anatomy_ontology.WS284.obo"))

gene_anatomy <- read.table(file.path("E:/public_data_annotation/wormbase/WS284", "anatomy_association.WS284.wb"), sep = "\t",
                           quote = "", stringsAsFactors = FALSE, comment.char = "!")[ , 2:5]
colnames(gene_anatomy) <- c("WBGeneID", "gene_name", "certainty", "anatomy_ontology")

# try to get an id :: term index from the ontology
anatomy_ontology_index <- data.frame(id = obo_parsed$id, name = obo_parsed$name)
anatomy_ontology_index_neuro <- anatomy_ontology_index[ anatomy_ontology_index$name %in% 
                                                          c("interneuron", "GABAergic neuron", "motor neuron", "sensory neuron",
                                                            "dopaminergic neuron", "glutamatergic neuron", "serotonergic neuron",
                                                            "cholinergic neuron") , ]
#merge(gene_anatomy[gene_anatomy$gene_name == "tph-1" , ], anatomy_ontology_index, by.x = "anatomy_ontology", by.y = "id")
# ha tph-1 isn't directly associated with the term for serotonergic neuron.... what.

grouped_go <- read.xlsx("Grouped GO terms corresponding to hpk-1 up and neuron-associated enrichment map.xlsx")


#
#
#

hpk1_up_neuronenrich_summary <- Reduce(f = rbind, lapply(unique(grouped_go$broad_group), function(curr_group){
  # get all unique genes across terms in the group
  group_genes = unique(unlist(strsplit(hpk1_up_neuron_assoc$intersection[hpk1_up_neuron_assoc$term_id %in% 
                                                                           grouped_go$term[grouped_go$broad_group == curr_group]],
                                       split = ",")))
  Reduce(f = rbind, x = lapply(group_genes, function(currGene){
    # get WBGeneID
    currWBGeneID = NA
    if(currGene %in% wb_gene_info$gene_name){
      currWBGeneID = wb_gene_info$WBGeneID[wb_gene_info$gene_name == currGene]
    }else if(currGene %in% wb_gene_info$molecular_name){
      # This means the gene name should also be updated...
      currWBGeneID = wb_gene_info$WBGeneID[which(wb_gene_info$molecular_name == currGene)]
    }else{
      currWBGeneID = NA
    }
    
    # get differential expression info
    currFC = with(DE_hpk1, DESeq2_log2FC[GeneSymbol == currGene])
    currPval = with(DE_hpk1, DESeq2_fdrpval[GeneSymbol == currGene])
    # is it neuron-specific or neuron-enriched?
    currNeuronExprClass = ifelse(currGene %in% kaletsky_unique_neurons, "neuron-specific",
                                 ifelse(currGene %in% kaletsky_enriched_neurons, "neuron-enriched",
                                        ""))
    # anatomy ontology association
    currNeuroAnat <- unique(merge(gene_anatomy[gene_anatomy$WBGeneID == currWBGeneID , ],
                           anatomy_ontology_index_neuro, by.x = "anatomy_ontology", by.y = "id")[,"name"])
    currNeuroAnat <- ifelse(length(currNeuroAnat) == 0, "", paste(currNeuroAnat, sep = "; ", collapse = "; "))
    # get human homologs- may be none
    currHomologs = ortholist_collapsed$Human_GeneSymbol[ortholist_collapsed$gene_name == currGene]
    # get lifespan category- may be none
    currLifespan = genage_simple$genage_gerogene_cats[genage_simple$GeneSymbol == currGene]
    # get if it's in the PGP- pretty unlikely
    currPGP = currGene %in% pgpList$gene_name
    # is a TF? what domain?
    currTF = wTF$DBD[wTF$WBGeneID == currWBGeneID]
    currTF = ifelse(length(currTF) == 0, "", ifelse(is.na(currTF), "", paste0("yes, family ", currTF)))
    # get any associated protein domains- could be a few results, need to combine
    currProtDomain = paste(protdomain_interpro_long$domain_desc[protdomain_interpro_long$gene_id %in% currWBGeneID],
                           sep = "; ", collapse = "; ")
    # # Associated GO terms (names) from this "group" ...  this is so ugly I'm sorry
    # currGoTerms = paste(hpk1_up_neuron_assoc$term_name[ hpk1_up_neuron_assoc$term_id %in%
    #                                                       grouped_go$term[unlist(
    #                                                         lapply(
    #                                                           strsplit(
    #                                                             hpk1_up_neuron_assoc$intersection[hpk1_up_neuron_assoc$term_id %in% 
    #                                                                                                 grouped_go$term[grouped_go$broad_group == curr_group]],
    #                                                             split = ","), function(geneVec){currGene %in% geneVec}))]], sep = " ||| ", collapse = " ||| ")
    return(
      data.frame(
        "broad_group" = curr_group,
        "gene_name" = currGene,
        "WBGeneID" =  ifelse(length(currWBGeneID) == 0, "", currWBGeneID), # THESE ARE NAMES THAT NEED TO BE CHANGED
        "Foldchange_log2" = currFC,
        "adj_pval" = currPval,
        "Neuron_expr_class" = currNeuronExprClass,
        "WB_anatomy_association_neuron_type" = currNeuroAnat,
        "human_homologs" = ifelse(length(currHomologs) == 0, "", currHomologs),
        "gerogene_category" = ifelse(length(currLifespan) == 0, "", currLifespan),
        "is_in_pgp" = currPGP,
        "is_TF" = currTF,
        "protein_domains" = ifelse(length(currProtDomain) == 0, "", currProtDomain)#,
        #"in_GO_terms_in_group" = currGoTerms
      )
    )
    
  }))
}))

# collapse group for genes that appear in multiple groups
hpk1_up_neuronenrich_summary_collapse <- 
  Reduce(f = rbind,
         x = lapply(unique(hpk1_up_neuronenrich_summary$gene_name),
                    function(currGene){
                      
                      if(nrow(hpk1_up_neuronenrich_summary[ hpk1_up_neuronenrich_summary$gene_name == currGene , ]) == 1){
                        return(hpk1_up_neuronenrich_summary[ hpk1_up_neuronenrich_summary$gene_name == currGene , ])
                      }else{
                        currDat = hpk1_up_neuronenrich_summary[ hpk1_up_neuronenrich_summary$gene_name == currGene , ]
                        retVal = currDat[1,] # pick one, only the first column is unique
                        retVal[1,1] = paste(currDat[,1], sep = " ; ", collapse = " ; ")
                        return(retVal)
                      }
                    }))


write.xlsx(hpk1_up_neuronenrich_summary_collapse, "hpk-1 upregulated gene enrichment neuron-associated term summary.xlsx", asTable = TRUE)




#
# TF enrichment
#

paths <- list()
paths$seq_data <- "C:/Users/Adam Cornwell/Box/adam/DR myc manuscript reboot/omics_data_sources"
paths$wormbase <- "E:/public_data_annotation/wormbase/WS284"
library(qusage)

wb_gene_info <- read.csv(file.path(paths$wormbase, "c_elegans.PRJNA13758.WS284.geneIDs.txt"), header = FALSE)[ , 2:6 ]
colnames(wb_gene_info) <- c("WBGeneID", "gene_name", "molecular_name", "status", "gene_biotype")
# keep live genes
wb_gene_info <- wb_gene_info[ wb_gene_info$status == "Live" , ]
# if gene name is NA, and molecular name is not NA, assign the molecular name as the gene name
wb_gene_info[ wb_gene_info == "" ] <- NA
wb_gene_info$gene_name[is.na(wb_gene_info$gene_name)] <- wb_gene_info$molecular_name[is.na(wb_gene_info$gene_name)]
any(is.na(wb_gene_info$gene_name))
any(duplicated(wb_gene_info$WBGeneID))
any(duplicated(wb_gene_info$gene_name))

# all predicted binding interactions
sets_tfs_all_ce <- qusage::read.gmt(file.path(paths$seq_data, "genesets",
                                      "motif-based binding predictions- all C elegans hits (no homology filter).gmt"))
sets_tfs_all_ce_long <- Reduce(x = lapply(names(sets_tfs_all_ce), function(x){
  cbind("motif" = x, "gene_name" = sets_tfs_all_ce[[x]])}), f = rbind)
sets_tfs_all_ce_long <- merge(sets_tfs_all_ce_long, wb_gene_info[ , c("WBGeneID", "gene_name") ], by.x = "gene_name")

# predicted binding interactions after filtering for homology
sets_tfs_species_filt <- qusage::read.gmt(file.path(paths$seq_data, "genesets",
                                            "motif-based binding predictions- matches in at least two other species.gmt"))
sets_tfs_species_filt_long <- Reduce(x = lapply(names(sets_tfs_species_filt), function(x){
  cbind("motif" = x, "gene_name" = sets_tfs_species_filt[[x]])}), f = rbind)
sets_tfs_species_filt_long <- merge(sets_tfs_species_filt_long, wb_gene_info[ , c("WBGeneID", "gene_name") ], by.x = "gene_name")

# add TF names to set names
metadata_tf_cisbp <- read.table(file.path("E:/project_data_ssd/myc_project/reboot_2022/TF binding prediction/",
                                          "Metadata for motif predictions- cisbp.txt"),
                                sep = "\t", quote = "", header = TRUE)
metadata_tf_cisbp$Motif_ID <- paste0("cisbp_", metadata_tf_cisbp$Motif_ID)
metadata_tf_jaspar <- read.table(file.path("E:/project_data_ssd/myc_project/reboot_2022/TF binding prediction/",
                                           "Metadata for motif predictions- jaspar.txt"),
                                 sep = "\t", quote = "", header = TRUE)
metadata_tf_jaspar$jasparID <- paste0("JASPAR_", metadata_tf_jaspar$jasparID)

sets_tfs_species_filt_long = merge(sets_tfs_species_filt_long, metadata_tf_cisbp, by.x = "motif", by.y = "Motif_ID", all.x = TRUE)
sets_tfs_species_filt_long = merge(sets_tfs_species_filt_long, metadata_tf_jaspar, by.x = "motif", by.y = "jasparID", all.x = TRUE)


sets_tfs_species_filt_long$motif_tf <- apply(sets_tfs_species_filt_long, MARGIN = 1, function(row){
  # motif/tf name is DBID or geneSymbol, pasted with motif, OR motif alone
  motifName = ifelse(is.na(row["DBID"]) && is.na(row["geneSymbol"]),
                     row["motif"],
                     paste(row["motif"], 
                           ifelse(is.na(row["DBID"]), 
                                  row["geneSymbol"], 
                                  wb_gene_info$gene_name[wb_gene_info$WBGeneID == row["DBID"]]), sep = "_"))
  return(motifName)
})

sets_tfs_species_filt_long$geneSymbol[is.na(sets_tfs_species_filt_long$geneSymbol)] <- sapply(
  sets_tfs_species_filt_long$motif_tf[is.na(sets_tfs_species_filt_long$geneSymbol)], function(x){
    ifelse(grepl("cisbp", x), gsub("^cisbp_[A-Z0-9]+_2.00_", "", x), NA)
  })

#  so yes- the "geneSymbol" column is the TF gene symbol, not the predicted binding gene
# then there are some additional ones we can set
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_unc130_core"] <- "unc-130"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "OH2004_CHE-1"] <- "che-1"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "OH2011_UNC-3"] <- "unc-3"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "OH2018_HLH-4"] <- "hlh-4"
# there are also a bunch of complexes, in which case we need to use a delimiter to include both genes
#  and then also check for the expression of both complex members
# I will *not* assign any gene symbols for the ebox motifs I manually defined, as there are too many possible
#  factors that bind eboxes.
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "OH2007_CEH-10_TTX-3"] <- "ceh-10_&&_ttx-3"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_CND-1"] <- "hlh-2_&&_cnd-1"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-10"] <- "hlh-2_&&_hlh-10"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-14"] <- "hlh-2_&&_hlh-14"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-15"] <- "hlh-2_&&_hlh-15"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-19"] <- "hlh-2_&&_hlh-19"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-3"] <- "hlh-2_&&_hlh-3"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-4"] <- "hlh-2_&&_hlh-4"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-8"] <- "hlh-2_&&_hlh-8"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_HLH-2_LIN-32"] <- "hlh-2_&&_lin-32"
sets_tfs_species_filt_long$geneSymbol[sets_tfs_species_filt_long$motif_tf == "uniprobe_grove_sw_MXL-1_MDL-1"] <- "mxl-1_&&_mdl-1"

#
# Ok, now we'll do the same thing but in the "non-species-filtered" object
#  probably could have done this one first and then just subset it but oh well
#

# add TF names to set names
sets_tfs_all_ce_long = merge(sets_tfs_all_ce_long, metadata_tf_cisbp, by.x = "motif", by.y = "Motif_ID", all.x = TRUE)
sets_tfs_all_ce_long = merge(sets_tfs_all_ce_long, metadata_tf_jaspar, by.x = "motif", by.y = "jasparID", all.x = TRUE)


sets_tfs_all_ce_long$motif_tf <- apply(sets_tfs_all_ce_long, MARGIN = 1, function(row){
  # motif/tf name is DBID or geneSymbol, pasted with motif, OR motif alone
  motifName = ifelse(is.na(row["DBID"]) && is.na(row["geneSymbol"]),
                     row["motif"],
                     paste(row["motif"], 
                           ifelse(is.na(row["DBID"]), 
                                  row["geneSymbol"], 
                                  wb_gene_info$gene_name[wb_gene_info$WBGeneID == row["DBID"]]), sep = "_"))
  return(motifName)
})

sets_tfs_all_ce_long$geneSymbol[is.na(sets_tfs_all_ce_long$geneSymbol)] <- sapply(
  sets_tfs_all_ce_long$motif_tf[is.na(sets_tfs_all_ce_long$geneSymbol)], function(x){
    ifelse(grepl("cisbp", x), gsub("^cisbp_[A-Z0-9]+_2.00_", "", x), NA)
  })

#  so yes- the "geneSymbol" column is the TF gene symbol, not the predicted binding gene
# then there are some additional ones we can set
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_unc130_core"] <- "unc-130"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "OH2004_CHE-1"] <- "che-1"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "OH2011_UNC-3"] <- "unc-3"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "OH2018_HLH-4"] <- "hlh-4"
# there are also a bunch of complexes, in which case we need to use a delimiter to include both genes
#  and then also check for the expression of both complex members
# I will *not* assign any gene symbols for the ebox motifs I manually defined, as there are too many possible
#  factors that bind eboxes.
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "OH2007_CEH-10_TTX-3"] <- "ceh-10_&&_ttx-3"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_CND-1"] <- "hlh-2_&&_cnd-1"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-10"] <- "hlh-2_&&_hlh-10"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-14"] <- "hlh-2_&&_hlh-14"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-15"] <- "hlh-2_&&_hlh-15"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-19"] <- "hlh-2_&&_hlh-19"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-3"] <- "hlh-2_&&_hlh-3"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-4"] <- "hlh-2_&&_hlh-4"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_HLH-8"] <- "hlh-2_&&_hlh-8"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_HLH-2_LIN-32"] <- "hlh-2_&&_lin-32"
sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == "uniprobe_grove_sw_MXL-1_MDL-1"] <- "mxl-1_&&_mdl-1"


# how many of these TFs have expression in neurons?
sets_tfs_species_filt_long_neuro <- sets_tfs_species_filt_long[ sets_tfs_species_filt_long$geneSymbol %in% rownames(cengen_sum_level2.hpk1) , ]
nrow(sets_tfs_species_filt_long)
nrow(sets_tfs_species_filt_long_neuro)
unique(sets_tfs_species_filt_long_neuro$geneSymbol)

library(goseq)
library(GenomicFeatures)
library(GenomicRanges)
transcriptDb <- makeTxDbFromGFF("E:/public_data_annotation/ensembl/celegans/106/Caenorhabditis_elegans.WBcel235.106.gtf", format="gtf")
txsByGene=transcriptsBy(transcriptDb,"gene")
lengthData=median(width(txsByGene))

lengthData <- merge(wb_gene_info[ wb_gene_info$gene_name %in%  DE_hpk1.all$GeneSymbol,  ], lengthData, by.x = "WBGeneID", by.y = 0)
colnames(lengthData)[ncol(lengthData)] <- "length"

currDE <- DE.degs_neuron_expressed$hpk1_US_vs_N2_US
currDE_up <- currDE$GeneSymbol[ currDE$DESeq2_log2FC > 0]
currDE_down <- currDE$GeneSymbol[ currDE$DESeq2_log2FC < 0]

# merge with length info
lengthData_currDE <- lengthData
lengthData_currDE$de_up <- lengthData_currDE$gene_name %in% currDE_up
lengthData_currDE$de_down <- lengthData_currDE$gene_name %in% currDE_down

# --DE upreg--
pwf=nullp(lengthData_currDE$de_up, bias.data = lengthData_currDE$length)
rownames(pwf) <- lengthData_currDE$WBGeneID
goEnrich <- goseq(pwf, gene2cat = sets_tfs_species_filt_long_neuro[, c("WBGeneID", "motif_tf") ])
goEnrich[,"over_represented_pvalue"] <- p.adjust(goEnrich$over_represented_pvalue, method="BH")
goEnrich[,"under_represented_pvalue"] <- p.adjust(goEnrich$under_represented_pvalue, method="BH")

goEnrich.sig.up <- goEnrich[ goEnrich$over_represented_pvalue < 0.05 , ]
    
goEnrich.sig.up$intersect_WBGeneID <- unlist(lapply(goEnrich.sig.up$category, function(x){
  paste(intersect(sets_tfs_species_filt_long_neuro[ sets_tfs_species_filt_long_neuro[, "motif_tf"] == x , "WBGeneID"], # genes in cat
                  lengthData_currDE$gene_name[lengthData_currDE$de_up]), sep = ", ", collapse = ", ") # genes sig in query set
}))

# --DE downreg--
pwf=nullp(lengthData_currDE$de_down, bias.data = lengthData_currDE$length)
rownames(pwf) <- lengthData_currDE$WBGeneID
goEnrich <- goseq(pwf, gene2cat = sets_tfs_species_filt_long_neuro[, c("WBGeneID", "motif_tf") ])
goEnrich[,"over_represented_pvalue"] <- p.adjust(goEnrich$over_represented_pvalue, method="BH")
goEnrich[,"under_represented_pvalue"] <- p.adjust(goEnrich$under_represented_pvalue, method="BH")

goEnrich.sig.down <- goEnrich[ goEnrich$over_represented_pvalue < 0.05 , ]

goEnrich.sig.down$intersect_WBGeneID <- unlist(lapply(goEnrich.sig.down$category, function(x){
  paste(intersect(sets_tfs_species_filt_long_neuro[ sets_tfs_species_filt_long_neuro[, "motif_tf"] == x , "WBGeneID"], # genes in cat
                  lengthData_currDE$gene_name[lengthData_currDE$de_down]), sep = ", ", collapse = ", ") # genes sig in query set
}))


# connections to longevity, stress response, proteostatic network
#  serotonergic and gabaeregic neuron <function>,
# hpk-1 interacting TFs
# 



#
# CONNECTOME
#
# read in the broad metadata- this is part of what what used to define styles in cytoscape
wormWiring_meta <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wormwiring connectome", "H_Attributes, Master.xlsx"))

table(wormWiring_meta$Neuron.Type)
table(wormWiring_meta$Neuron.Group)
# are the names unique??
sum(duplicated(wormWiring_meta$name))
sort(table(wormWiring_meta$name), decreasing = TRUE)
# some of these are just completely duplicate rows?
wormWiring_meta <- unique(wormWiring_meta)
sum(duplicated(wormWiring_meta$name)) # seems to have been the case.

#View(as.matrix(colnames(cengen_sum_level2)))

# now... we need to map cell identity between this and CenGen- which includes dealing with neuron pairs that were summarized in CenGen
# a bunch of these won't actually map since they're meta-nodes and not actually neurons.
# Any neuron that ends in an L or R that has a corresponding terminal L AND R entry
#   should correspond to something in cengen with the L/R stripped. Otherwise it's probably standalone, or not a neuron.

# we'll just iterate over the names in the connectome, identify the corresponding cengen entry, or NA if there doesn't see to be one.
#  Then I suppose we can go through the gaps manually to see if anything was missed.

# # THIS WAS ALREADY DONE AND THE RESULTS REVISED MANUALLY
# pairedNeuronNames <- sapply(wormWiring_meta$name, function(currName){
#   if(substr(currName, nchar(currName), nchar(currName)) == "L"){
#     # look for entry with R
#     if(paste0(substr(currName, 1, nchar(currName)-1), "R") %in% wormWiring_meta$name){
#       # return curr name as part of a pair
#       return(currName)
#     }
#   }
#   # this is kind of lazy, but now do the opposite.
#   else if(substr(currName, nchar(currName), nchar(currName)) == "R"){
#     # look for entry with L
#     if(paste0(substr(currName, 1, nchar(currName)-1), "L") %in% wormWiring_meta$name){
#       # return curr name as part of a pair
#       return(currName)
#     }
#   } else{return(NA)}
# })
# 
# pairedNeuronNames <- unlist(pairedNeuronNames)
# pairedNeuronNames <- pairedNeuronNames[!is.na(pairedNeuronNames)]
# names(pairedNeuronNames) <- NULL
# 
# wormWiring_meta$isPaired <- sapply(wormWiring_meta$name, function(x){x %in% pairedNeuronNames})
# 
# wormWiring_meta$mergeName <- sapply(1:nrow(wormWiring_meta), function(i){
#   if(wormWiring_meta$isPaired[i]){
#     return(substr(wormWiring_meta$name[i], 1, nchar(wormWiring_meta$name[i])-1))
#   }else{
#     return(wormWiring_meta$name[i])
#   }
# })

# I'm just going to fix these manually and read back in
#write.xlsx(wormWiring_meta, file.path("E:/public_data_annotation/c_elegans_other/wormwiring connectome", "name_mapping_to_cengen.xlsx"))
wormWiring_meta2 <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wormwiring connectome", "name_mapping_to_cengen.xlsx"), sheet = 1)

cengenToWireMap <- merge(colnames(cengen_sum_level2), wormWiring_meta2,
      by.x = 1, by.y = "CengenName", all.x = TRUE)

cengenToWireMap <- cengenToWireMap[ !is.na(cengenToWireMap$ConnectomeName) , ]
colnames(cengenToWireMap)[1] <- "CengenName"
table(cengenToWireMap$NeuronType)

#
# Ok, now we'll see if we can export something that we can incorporate into the cytoscape connectome map
# There are 473 nodes in that network, which is... fewer than their master annotation excel file
#  so I've just exported the node table from cytoscape, since ultimately we need things to match up to that

connectome_nodetable <- read.csv(file.path("E:/public_data_annotation/c_elegans_other/wormwiring connectome",
                                           "Hermaphrodite_map_node_table_export.csv"))

connectome_nodetable_merge <- merge(connectome_nodetable, cengenToWireMap[,c("CengenName", "ConnectomeName")],
                                    by.x = "name", by.y = "ConnectomeName", all.x = TRUE)

write.xlsx(connectome_nodetable_merge,
           file.path("E:/public_data_annotation/c_elegans_other/wormwiring connectome", "name_mapping_to_cengen_v2.xlsx"))


# merge into hpk-1 expression object, and plot

hpk1_neuron_type_expr <- merge(hpk1_neuron_type_expr, connectome_nodetable_merge[,c("CengenName", "Neuron.Type")],
                               by.x = 0, by.y = "CengenName", all.x = TRUE)

ggplot(data = hpk1_neuron_type_expr[ (hpk1_neuron_type_expr$Neuron.Type != "Endorgan") & !is.na(hpk1_neuron_type_expr$Neuron.Type)  , ],
       aes(x = hpk1_expr)) +
  geom_histogram(bins = 20) +
  facet_wrap(vars(Neuron.Type)) +
  theme_bw(base_size = 20) +
  xlab("hpk-1 expression") + ylab("# of neurons")


# go the other direction to add the hpk-1 expression info to the nodetable
connectome_nodetable_hpk1 <- merge(connectome_nodetable_merge, cengen_sum_level2["hpk-1",],
                                   by.x = "CengenName", by.y = 0, all.x = TRUE)
colnames(connectome_nodetable_hpk1)[ncol(connectome_nodetable_hpk1)] <- "expr_hpk1"

write.csv(connectome_nodetable_hpk1, file.path("E:/public_data_annotation/c_elegans_other/wormwiring connectome",
                                               "nodetable_hpk1expr_merge.csv"), row.names = FALSE)

summary(cengen_sum_level2["hpk-1",])


#
# Oliver Hobert has a great compendium in wormbook on the gene classes that are functionally important for neurons
#  but these tables are only in html... so... we'll scrape it...
#

library(rvest)

htmlDoc <- read_html("http://www.wormbook.org/chapters/www_neuronalgenome/neuronalgenome.html")

divClass <- html_attr(html_elements(htmlDoc, "div"), "class")
divClass[is.na(divClass)] <- ""
tableNodes <- html_elements(htmlDoc, "div")[ divClass == "table"]
tableNodes <- tableNodes[-1] # don't care about table 1

tables <- html_table(tableNodes)

# can we pull the table names out?
tableNames <- html_text2(html_elements(htmlDoc, "h3"))
tableNames <- tableNames[33:length(tableNames)]

names(tables) <- tableNames

# write each of these to an excel file?
lapply(names(tables), function(currTable){
  write.xlsx(tables[[currTable]], file = file.path("E:/project_data_ssd/hpk-1_dataset/hobert_neuron_genes",
                                                   paste0(gsub(":", "-", currTable), ".xlsx")))
})

write.xlsx(tables[[11]], file = file.path("E:/project_data_ssd/hpk-1_dataset/hobert_neuron_genes",
                                                 paste0(gsub(":|/", "-", names(tables)[11]), ".xlsx")))

write.xlsx(tables[[16]], file = file.path("E:/project_data_ssd/hpk-1_dataset/hobert_neuron_genes",
                                          paste0(gsub(":|/", "-", names(tables)[16]), ".xlsx")))

#
# Doing a bit of manual editing and then will read back in...
#
# most of these tables have a bit of a different set of headers- we'll have to define which are most relevant here!
# okay, got it! I think.

hobert_table_idx <- read.xlsx(file.path("E:/project_data_ssd/hpk-1_dataset", "Hobert neuron genes table inventory.xlsx"))
hobert_table_filenames <- list.files("E:/project_data_ssd/hpk-1_dataset/hobert_neuron_genes", pattern = "*.xlsx")
hobert_table_filenames <- cbind(hobert_table_filenames,
                                table_shortname = sapply(hobert_table_filenames, function(x){strsplit(x, "-")[[1]][1]}))
row.names(hobert_table_filenames) <- NULL
hobert_table_idx <- merge(hobert_table_idx, hobert_table_filenames)
rm(hobert_table_filenames)

# ok, I want to assemble a single giant list from all of these
# Starting with three columns- the table description,
#   the category (if applicable, some don't have subcategories), the gene name

# Ack, my original version of this was assigning the same major and sub categories to every row in the table!
#  The subcategory column may be NA, in which case all the results for subcategory for that column will be NA,
#   but otherwise, it should just return whatever the subcategory column is.

neuron_gene_classes <- Reduce(f = rbind, x = lapply(1:nrow(hobert_table_idx), function(idx){
  currTable = read.xlsx(file.path("E:/project_data_ssd/hpk-1_dataset/hobert_neuron_genes",
                                  hobert_table_idx$hobert_table_filenames[idx]))
  currTable = data.frame(currTable[,na.omit(c(hobert_table_idx$category_column[idx], hobert_table_idx$gene_column[idx]))])
  return(data.frame(
    broad_class = rep(hobert_table_idx$table_description[idx], nrow(currTable)),
    subcategory = if(is.na(hobert_table_idx$category_column[idx])){
      rep("N/A", nrow(currTable))
    }else{
      currTable[,1]
    },
    gene = if(ncol(currTable) == 1){currTable[,1]}else{currTable[,2]}
  ))
}))




# write out this master table
write.xlsx(neuron_gene_classes, file = "Oliver Hobert wormbook neuron transcriptome tables aggregated and simplified (FIXED).xlsx")


neuron_gene_classes_hpk1DE <- merge(neuron_gene_classes, DE_hpk1[ , c("GeneSymbol", "DESeq2_log2FC", "count_range_HPK-1_US","count_range_N2_US")],
                                    by.x = "gene", by.y = "GeneSymbol")

#
# Aging atlas info
#
aging_atlas_cluster_shared_up <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/supplemental materials",
                                                     "Table S8 aging gene changes- shared across clusters (ABC edit).xlsx"),
                                           sheet = "common_aging_up")
aging_atlas_cluster_shared_down <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/supplemental materials",
                                                     "Table S8 aging gene changes- shared across clusters (ABC edit).xlsx"),
                                           sheet = "common_aging_down")

neuron_gene_classes_hpk1DE$aging_atlas_common <- sapply(neuron_gene_classes_hpk1DE$gene, function(currGene){
  if(currGene %in% aging_atlas_cluster_shared_up$gene){return("up_with_age")
  }else if(currGene %in% aging_atlas_cluster_shared_down$gene){return("down_with_age")
      }else{return(NA)}
})


#
#  What's the overlap with the neuron-associated GO term list I already assembled?
# 

length(hpk1_up_neuronenrich_summary_collapse$gene_name)
length(neuron_gene_classes_hpk1DE$gene)
length(intersect(neuron_gene_classes_hpk1DE$gene, hpk1_up_neuronenrich_summary_collapse$gene_name))
length(union(neuron_gene_classes_hpk1DE$gene, hpk1_up_neuronenrich_summary_collapse$gene_name))

#
# Aging atlas neuron type mapping-
#  I already mapped the cengen cell clusters to the connectome, for additional neuron annotation on top of what
#   what was already available from Cengen.
#  We also need to do this for the aging atlas cell clusters. 
#   They provide a table, supplemental table 3 ("media-4.xlsx") where they provide the correspondence
#     between their clusters, and the highest-correlation match in the cengen L4 dataset. So we can start with that.
#
#  I previously had un-grouped some of the cengen cluster names to map the ambiguous ones to the connectome neurons
#   and I'll need to do that again here, so I'll have to modify their table first (done).
#  To reduce duplication, I'll actually just take the cengen and atlas names for now, and not merge in the connectome names
#
#  TEN cell clusters in Cengen don't have a mappable equivalent or class in the aging atlas dataset
#
atlas_cengen_neuron_map <- read.xlsx("aging_atlas_to_cengen_cell_map (adapted from table s3 in atlas preprint).xlsx", colNames = TRUE)
atlas_cengen_neuron_map <- unique(merge(cengenToWireMap[,c("CengenName", "NeuronType")], atlas_cengen_neuron_map[,1:2],
                                 by.x = "CengenName", by.y = "cengen_name"))

# we can try to use marker gene expression to annotate neurotransmitter class
#  and see if there's a consistent result across datasets and age
# serotonin- tph-1 
# gabaergic- unc-47
# glutamatergic- eat-4
# dopaminergic- cat-2
# cholinergic- unc-17
# thermosensory- ttx-1 (note NOT ttx-3)

atlas_cengen_neuron_map$cengen_signal_class <- sapply(1:nrow(atlas_cengen_neuron_map), function(i){
  
  currCengen = atlas_cengen_neuron_map$CengenName[i]
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-1")
  
  return(paste(names(typeMarkers)[which(cengen_sum_level2[typeMarkers,currCengen] > 75)], collapse = ", "))
  
})

hist(as.vector(aging_atlas_cluster_mean_counts$D1[c("serotonergic" = "tph-1",
                                                    "GABAergic" = "unc-47",
                                                    "glutamatergic" = "eat-4",
                                                    "dopaminergic" = "cat-2",
                                                    "cholinergic" = "unc-17",
                                                    "thermosensory" = "ttx-3") , ]), breaks = seq(0, 5, by = 0.1))
hist(as.vector(aging_atlas_cluster_mean_counts$D3[c("serotonergic" = "tph-1",
                                                    "GABAergic" = "unc-47",
                                                    "glutamatergic" = "eat-4",
                                                    "dopaminergic" = "cat-2",
                                                    "cholinergic" = "unc-17",
                                                    "thermosensory" = "ttx-3") , ]), breaks = seq(0, 6, by = 0.1))

# atlas, D1
atlas_cengen_neuron_map$atlas_signal_class_d1 <- sapply(1:nrow(atlas_cengen_neuron_map), function(i){
  currAtlas = atlas_cengen_neuron_map$atlas_cluster[i]
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D1[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
})

# atlas, D3
atlas_cengen_neuron_map$atlas_signal_class_d3 <- sapply(1:nrow(atlas_cengen_neuron_map), function(i){
  currAtlas = atlas_cengen_neuron_map$atlas_cluster[i]
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D3[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
})

# atlas, D5
atlas_cengen_neuron_map$atlas_signal_class_d5 <- sapply(1:nrow(atlas_cengen_neuron_map), function(i){
  currAtlas = atlas_cengen_neuron_map$atlas_cluster[i]
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D5[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
})

# atlas, D8
atlas_cengen_neuron_map$atlas_signal_class_d8 <- sapply(1:nrow(atlas_cengen_neuron_map), function(i){
  currAtlas = atlas_cengen_neuron_map$atlas_cluster[i]
  typeMarkers = c("serotonin" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D8[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
})

#
# For neuron clusters that are NOT able to be mapped to a cengen cluster,
#  use the aging atlas expression to make a similar call
#
tmp_atlas_unique_neurons <- setdiff(na.omit(aging_atlas_clusters_all$shortID[aging_atlas_clusters_all$`Tissue.type.(higher.degree)` == "neuron"]),
        atlas_cengen_neuron_map$atlas_cluster)

atlas_unique_neuron_signal_class = data.frame(atlas_D1 = sapply(tmp_atlas_unique_neurons, function(currAtlas){
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D1[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
}),
atlas_D3 = sapply(tmp_atlas_unique_neurons, function(currAtlas){
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D3[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
}),
atlas_D5 = sapply(tmp_atlas_unique_neurons, function(currAtlas){
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D5[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
}),
atlas_D8 = sapply(tmp_atlas_unique_neurons, function(currAtlas){
  typeMarkers = c("serotonergic" = "tph-1",
                  "GABAergic" = "unc-47",
                  "glutamatergic" = "eat-4",
                  "dopaminergic" = "cat-2",
                  "cholinergic" = "unc-17",
                  "thermosensory" = "ttx-3")
  return(paste(names(typeMarkers)[which(aging_atlas_cluster_mean_counts$D8[typeMarkers,currAtlas] > 0.5)], collapse = ", "))
}))

rm(tmp_atlas_unique_neurons)
# make consensus call
atlas_unique_neuron_signal_class$consensus <- unlist(apply(atlas_unique_neuron_signal_class, MARGIN = 1, function(x){
  y = unique(x)
  y = y[y != ""]
  if(length(y) > 1){
    return("")
  }else if(length(y) == 0){
    return("")
  }else{
    return(y)
  }
}))


#
# Now that we have the rest of the aging atlas info, we can include more detail on which genes change with age
#
# first, read in table I made which has the clusters with DE results (not all of them did)
aging_atlas_deg_clusters <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                "Mean expr for young and old adapted from DEtable by cluster.xlsx"), sheet = 1)
# read in the cluster-level DE gene results. All of these are OLD VS YOUNG
aging_atlas_degs <- lapply(aging_atlas_deg_clusters$shortID, function(x){
  read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas","Copy of deg_table.xlsx"), sheet = x)
})

# convert this to long format. Maybe useful for plots later, but also useful for not having to page
#  through a list of dataframes to pull out the FCs for every gene of interest.
names(aging_atlas_degs) <- aging_atlas_deg_clusters$shortID
# all these should have length > 0
summary(unlist(lapply(aging_atlas_degs, nrow)))
# I've previously confirmed that column names in this file are the same in every sheet

aging_atlas_degs_long <- Reduce(f = rbind, x = lapply(names(aging_atlas_degs), function(clustName){
  return(data.frame(
    cluster_shortID = clustName,
    gene_name = aging_atlas_degs[[clustName]][,1],
    logFoldChange = aging_atlas_degs[[clustName]][,3],
    padj = aging_atlas_degs[[clustName]][,5]
  ))
})
)

aging_atlas_degs_long.unfilt <- Reduce(f = rbind, x = lapply(names(aging_atlas_degs), function(clustName){
  return(data.frame(
    cluster_shortID = clustName,
    gene_name = aging_atlas_degs[[clustName]][,1],
    logFoldChange = aging_atlas_degs[[clustName]][,3],
    padj = aging_atlas_degs[[clustName]][,5]
  ))
})
)

# this is not just genes with significant p-values-
#  FILTER
nrow(aging_atlas_degs_long)
aging_atlas_degs_long <- aging_atlas_degs_long[ aging_atlas_degs_long$padj < 0.05 , ]
nrow(aging_atlas_degs_long)

plot(density(aging_atlas_degs_long$logFoldChange))
plot(density(abs(aging_atlas_degs_long$logFoldChange)))
sum(abs(aging_atlas_degs_long$logFoldChange) > 1)
sum(abs(aging_atlas_degs_long$logFoldChange) > 0.5) # will go with this
nrow(aging_atlas_degs_long)
aging_atlas_degs_long <- aging_atlas_degs_long[ abs(aging_atlas_degs_long$logFoldChange) >= 0.5 , ]
nrow(aging_atlas_degs_long)

# merge in tissue annotation for convenience
aging_atlas_degs_long <- merge(aging_atlas_degs_long, aging_atlas_deg_clusters[,c("shortID","Clusters.explanatory.name")], 
                               by.x = "cluster_shortID", by.y = "shortID", all.x = TRUE)
nrow(aging_atlas_degs_long)

# how many genes in this are on my hobert intersect list? 
length(neuron_gene_classes_hpk1DE$gene)
sum(neuron_gene_classes_hpk1DE$gene %in% aging_atlas_degs_long$gene_name)

# how many hpk-1 DE genes are DE in the aging atlast but NOT in the neuron transcriptome list?
sum(DE_hpk1$GeneSymbol %in% aging_atlas_degs_long$gene_name)
sum(DE_hpk1$GeneSymbol %in% setdiff(aging_atlas_degs_long$gene_name, neuron_gene_classes_hpk1DE$gene))
# ok actually a lot.
# maybe we'll look at tissue aging DE association and make an upset plot that way

# but first, finish the neuron focus component
# for each gene in neuron transcriptome hpk-1 DE overlap list,
#  in what tissues do they change with age? need a separate down and up column probably
neuron_gene_classes_hpk1DE <- cbind(neuron_gene_classes_hpk1DE, t(sapply(neuron_gene_classes_hpk1DE$gene, function(currGene){
  if(currGene %in% aging_atlas_degs_long$gene_name){
    
    with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
         if((sum(logFoldChange > 0) > 0) & (sum(logFoldChange < 0) > 0)){ # gene is both up and down in different places
           return(data.frame(up = paste(Clusters.explanatory.name[logFoldChange > 0], sep = "; ", collapse = "; "),
                             down = paste(Clusters.explanatory.name[logFoldChange < 0], sep = "; ", collapse = "; ")))
         }else if(sum(logFoldChange > 0) > 0){ # gene is just upregulated
          return(data.frame(up = paste(Clusters.explanatory.name[logFoldChange > 0], sep = "; ", collapse = "; "), down = NA))
         }else if(sum(logFoldChange < 0) > 0){ # gene is just downregulated
           return(data.frame(up = NA, down = paste(Clusters.explanatory.name[logFoldChange < 0], sep = "; ", collapse = "; ")))
         }else{stop("something unexpected occurred")}
    )
  }else{
    return(data.frame(up = NA, down = NA))
  }
})))

neuron_gene_classes_hpk1DE$num_clusters_aging_up <- sapply(neuron_gene_classes_hpk1DE$gene, function(currGene){
  if(currGene %in% aging_atlas_degs_long$gene_name){
    with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
        sum(logFoldChange > 0))
  }else{ return(0)}
})

neuron_gene_classes_hpk1DE$num_clusters_aging_down <- sapply(neuron_gene_classes_hpk1DE$gene, function(currGene){
  if(currGene %in% aging_atlas_degs_long$gene_name){
    with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
         sum(logFoldChange < 0))
  }else{ return(0)}
})


#
# Also include co-expression with hpk-1 across cells in this
#  Looking at what genes are co-expressed with hpk-1 across clusters within a timepoint gives us an idea
#   of what genes have an spatial expression pattern that strongly resembles that of hpk-1
#  This does not necessarily get us to the influence hpk-1 has on expression.
#  On the other hand, correlation across age might be more interesting for inference of causal or regulatory
#   relationshops, since we know that hpk-1 expression levels *do* change with age.
# Would we want to look at that *within* each cluster? across all clusters?
#   Could be different in each tissue, maybe want to start with within-cluster.
#
#
# Not including D15 for now
#  have to skip the first couple rows since I had manually included additional headers
aging_atlas_cluster_mean_counts <- list(D1 = read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                       "Mean counts for all genes per cluster by day.xlsx"), sheet = "D1",
                                                       startRow = 4, check.names = FALSE, rowNames = TRUE),
                                        D3 = read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                                 "Mean counts for all genes per cluster by day.xlsx"), sheet = "D3",
                                                       startRow = 4, check.names = FALSE, rowNames = TRUE),
                                        D5 = read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                                 "Mean counts for all genes per cluster by day.xlsx"), sheet = "D5",
                                                       startRow = 4, check.names = FALSE, rowNames = TRUE),
                                        D8 = read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                                 "Mean counts for all genes per cluster by day.xlsx"), sheet = "D8",
                                                       startRow = 4, check.names = FALSE, rowNames = TRUE),
                                        D11 = read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                                 "Mean counts for all genes per cluster by day.xlsx"), sheet = "D11",
                                                        startRow = 4, check.names = FALSE, rowNames = TRUE)
                                        )

head(aging_atlas_cluster_mean_counts$D1)
class(aging_atlas_cluster_mean_counts$D1)
mode(aging_atlas_cluster_mean_counts$D1)
aging_atlas_cluster_mean_counts <- lapply(aging_atlas_cluster_mean_counts, function(x){
  x = as.matrix(x)
  mode(x) = "numeric"
  return(x)
})
class(aging_atlas_cluster_mean_counts$D1)
mode(aging_atlas_cluster_mean_counts$D1)


#
# alright, so now try to look for genes with correlation to hpk-1 expression across clusters
#   We'll start simple and use Pearson. We can do this for each day separately.
#
hpk1_corr_pearson_d1 <- Reduce(f = rbind, x = lapply(rownames(aging_atlas_cluster_mean_counts$D1), function(geneName){
    corRes = cor.test(x = aging_atlas_cluster_mean_counts$D1["hpk-1",],
                      y = aging_atlas_cluster_mean_counts$D1[geneName,],
                      method = "pearson")
    return(data.frame(gene_name = geneName, time = "D1", cor = corRes$estimate, p.value = corRes$p.value))
}))

hpk1_corr_pearson_d3 <- Reduce(f = rbind, x = lapply(rownames(aging_atlas_cluster_mean_counts$D3), function(geneName){
  corRes = cor.test(x = aging_atlas_cluster_mean_counts$D3["hpk-1",],
                    y = aging_atlas_cluster_mean_counts$D3[geneName,],
                    method = "pearson")
  return(data.frame(gene_name = geneName, time = "D3", cor = corRes$estimate, p.value = corRes$p.value))
}))

hpk1_corr_pearson_d5 <- Reduce(f = rbind, x = lapply(rownames(aging_atlas_cluster_mean_counts$D5), function(geneName){
  corRes = cor.test(x = aging_atlas_cluster_mean_counts$D5["hpk-1",],
                    y = aging_atlas_cluster_mean_counts$D5[geneName,],
                    method = "pearson")
  return(data.frame(gene_name = geneName, time = "D5", cor = corRes$estimate, p.value = corRes$p.value))
}))

hpk1_corr_pearson_d8 <- Reduce(f = rbind, x = lapply(rownames(aging_atlas_cluster_mean_counts$D8), function(geneName){
  corRes = cor.test(x = aging_atlas_cluster_mean_counts$D8["hpk-1",],
                    y = aging_atlas_cluster_mean_counts$D8[geneName,],
                    method = "pearson")
  return(data.frame(gene_name = geneName, time = "D8", cor = corRes$estimate, p.value = corRes$p.value))
}))

hpk1_corr_pearson_d11 <- Reduce(f = rbind, x = lapply(rownames(aging_atlas_cluster_mean_counts$D11), function(geneName){
  corRes = cor.test(x = aging_atlas_cluster_mean_counts$D11["hpk-1",],
                    y = aging_atlas_cluster_mean_counts$D11[geneName,],
                    method = "pearson")
  return(data.frame(gene_name = geneName, time = "D11", cor = corRes$estimate, p.value = corRes$p.value))
}))

# double check row order is the same
identical(hpk1_corr_pearson_d1$gene_name, hpk1_corr_pearson_d11$gene_name) # should be true!

# adjust all the p-values for multiple testing.
#  p-values >= 0.01 get set to NA
hpk1_corr_pearson <- list(D1 = hpk1_corr_pearson_d1,
                          D3 = hpk1_corr_pearson_d3,
                          D5 = hpk1_corr_pearson_d5,
                          D8 = hpk1_corr_pearson_d8,
                          D11 = hpk1_corr_pearson_d11)

hpk1_corr_pearson <- lapply(hpk1_corr_pearson, function(dat){
  dat$p.value = p.adjust(dat$p.value, method = "fdr")
  dat$cor[dat$p.value >= 0.01] = NA
  return(dat)
})

# combine cor across timepoints into one matrix
hpk1_corrMat_pearson <- Reduce(f = cbind, x =lapply(hpk1_corr_pearson, function(dat){
  dat$cor
}))
colnames(hpk1_corrMat_pearson) <- names(hpk1_corr_pearson)
rownames(hpk1_corrMat_pearson) <- hpk1_corr_pearson$D1$gene_name # they should all be the same gene order

hpk1_corrMat_pearson <- cbind(hpk1_corrMat_pearson,
                              "meanCorr" = rowMeans(hpk1_corrMat_pearson, na.rm = TRUE),
                              "numCorrDays" = rowSums(!is.na(hpk1_corrMat_pearson)))
#
# combine with the neuron info
# 
neuron_gene_classes_hpk1DE <- merge(neuron_gene_classes_hpk1DE, hpk1_corrMat_pearson[,c("meanCorr", "numCorrDays")],
                                    by.x = "gene", by.y = 0, all.x = TRUE)

# write this out for now
write.xlsx(x = neuron_gene_classes_hpk1DE, asTable = TRUE,
           file = "hpk-1 DE genes with neuron-associated function, with aging atlas DE info and hpk1 corr results across tissues v2.1.xlsx")

#
# Quick and dirty permutation test to determine if the overlap with the neuron gene compendium and the hpk-1 DE genes
#  is much more than we'd expect by chance- and it is.
#
set.seed(12345)
# sample from the genes expressed in the actual rna-seq experiment.
tmp_permTest <- sapply(1:100000, function(i){sum(sample(DE_hpk1.all$GeneSymbol, size = 2201) %in% neuron_gene_classes$gene)})
summary(tmp_permTest)
plot(density(tmp_permTest ))

#
# Heatmap for correlation across cell clusters
#
heatmapData <- Reduce(f = cbind, x =lapply(hpk1_corr_pearson, function(dat){
  dat$cor
}))
colnames(heatmapData) <- names(hpk1_corr_pearson)
rownames(heatmapData) <- hpk1_corr_pearson$D1$gene_name
table(rowSums(!is.na(heatmapData)))
# keep results with at least 3 non-NA results
heatmapData <- heatmapData[ rowSums(!is.na(heatmapData)) >= 3  , ]
any(heatmapData < 0, na.rm = TRUE) # all values are positive, so we can use a different color scale?
sum(abs(rowMeans(heatmapData, na.rm = TRUE)) > 0.4) 
heatmapData <- heatmapData[ abs(rowMeans(heatmapData, na.rm = TRUE)) > 0.4 , ]
# set NA to 0
heatmapData[is.na(heatmapData)] <- 0


# we don't have any annotation to add here for columns since it's across clusters

# are any of these genes differentially expressed in our dataset?
rownames(heatmapData)[rownames(heatmapData) %in% DE_hpk1$GeneSymbol]

pdf(file = "Correlation heatmap with hpk-1 across cell clusters for age x genes.pdf", width = 10, height = 14)
heatmap.3(heatmapData,
          col = viridis(option = "magma", n = 128, direction = 1), #scale = "row",
          cexCol = 3, cexRow = 1, margins = c(10,7), keysize = 0.8,
          Colv = FALSE, dendrogram = "row",
          KeyValueName = "Expr corr across cell clusters")
graphics.off()

rm(heatmapData)


#
# ok now we look at correlation with hpk-1 across time, within clusters.
#
# for this we have a file with 211 sheets
aging_atlas_clusters_all <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                                                "Mean counts for all genes per cluster by day.xlsx"), sheet = 1)

aging_atlas_within_tissue <- lapply(aging_atlas_clusters_all$shortID, function(clusterName){
  read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/calico_aging_atlas/wrangled",
                      "Calico_single_cell_summarized-mean_counts_by_day_per_cluster.xlsx"), sheet = clusterName, rowNames = TRUE)
})
names(aging_atlas_within_tissue) <- aging_atlas_clusters_all$shortID

# pearson correlation between hpk-1 and all other genes, within each cell cluster
#  and then we'll try to filter, since it's unlikely that there will be informative correlations
#  in every cluster. We expect informative correlations across time where hpk-1 actually changes
#  expression with age.
# There are many more cases here where hpk-1 has no expression in a cluster, will produce NA correlations.

# this takes a while, could be parallelized.

# should only look in clusters where hpk-1 is actually expressed!
hpk1_expressed_clusters <- data.frame(Reduce(f = rbind, lapply(aging_atlas_within_tissue, function(x){x["hpk-1",]})))
hpk1_expressed_clusters$cluster_name <- names(aging_atlas_within_tissue)
hpk1_expressed_clusters$num_zero <- rowSums(hpk1_expressed_clusters[ , c("d1", "d3", "d5", "d8", "d11", "d15") ] == 0, na.rm = TRUE)
hpk1_expressed_clusters$num_na <- rowSums(is.na(hpk1_expressed_clusters[ , c("d1", "d3", "d5", "d8", "d11", "d15") ]), na.rm = TRUE)
hpk1_expressed_clusters$total_zero_na <- hpk1_expressed_clusters$num_zero + hpk1_expressed_clusters$num_na


# next time add a status print to this... and parallelize....

aging_atlas_within_tissue.cors <- lapply(
  aging_atlas_within_tissue[hpk1_expressed_clusters$cluster_name[hpk1_expressed_clusters$total_zero_na < 3]], function(clusterDat){
  
  Reduce(f = rbind, x = lapply(rownames(clusterDat), function(geneName){
    corRes = cor.test(x = as.numeric(clusterDat["hpk-1",]),
                      y = as.numeric(clusterDat[geneName,]),
                      method = "pearson")
    return(data.frame(gene_name = geneName, cor = corRes$estimate, p.value = corRes$p.value))
  }))
  
})
names(aging_atlas_within_tissue.cors)

aging_atlas_within_tissue.cors <- lapply(aging_atlas_within_tissue.cors, function(x){
  x$adj.p.value = p.adjust(x$p.value, method = "fdr")
  return(x)
})

# due to the small number of datapoints (six ages, with most missing day 15)
#  most of the p-values are not particularly small, and are not significant after multiple testing
#  so we probably want to use a low p-value cutoff on the uncorrected p-values

# include correlation across tissues as another possible dimension for each gene?
#  could even do PCA then, hypothetically.

# combine results across tissues
aging_atlas_cors <- Reduce(f = cbind, lapply(aging_atlas_within_tissue.cors, function(x){
  x$cor[x$p.value >= 0.01] = NA
  return(x$cor)
}))
colnames(aging_atlas_cors) <- names(aging_atlas_within_tissue.cors)
identical(aging_atlas_within_tissue.cors$`0_0`$gene_name, aging_atlas_within_tissue.cors$`143_0`$gene_name) # sanity check
rownames(aging_atlas_cors) <- aging_atlas_within_tissue.cors$`0_0`$gene_name # all should be the same

# remove any rows with NA in every column
table(rowSums(is.na(aging_atlas_cors)))
dim(aging_atlas_cors)
aging_atlas_cors <- aging_atlas_cors[ rowSums(is.na(aging_atlas_cors)) < ncol(aging_atlas_cors) , ]
dim(aging_atlas_cors)
# any columns that only have 1 result? (would be hpk-1 self-correlation)
aging_atlas_cors <- aging_atlas_cors[ , colSums(!is.na(aging_atlas_cors)) > 1 ]
dim(aging_atlas_cors)

# with original cutoff- 
#  11_1 (vulva)- 1157 genes positively correlated with hpk-1 across age?
#  others with 50 or more genes (less to more0: 6_3, 28_0, 84_0, 5_1, 4_0, 17_0, 142_0, 35_4, 4_1, 30_0, 38_0, 113_0, 101_0
# but this is too sparse to be visually informative, so we can try to slim it down more

# filter rows
aging_atlas_cors <- aging_atlas_cors[ rowSums(!is.na(aging_atlas_cors)) > 10 , ]
dim(aging_atlas_cors)

# filter columns
aging_atlas_cors <- aging_atlas_cors[ , colSums(!is.na(aging_atlas_cors)) > 4 ]
dim(aging_atlas_cors) # much more reasonable

#
# Plot heatmap of correlations within cell clusters across aging
#

# indicate neurons vs non-neurons for expression clusters at least
tissue_type_vec <- aging_atlas_clusters_all[,c("Tissue.type.(higher.degree)", "shortID")]
tissue_type_vec$`Tissue.type.(higher.degree)`[tissue_type_vec$`Tissue.type.(higher.degree)` == "various"] <- "various_or_unknown"
tissue_type_vec$`Tissue.type.(higher.degree)`[is.na(tissue_type_vec$`Tissue.type.(higher.degree)`)] <-  "various_or_unknown"

tissue_type_color <- setNames(c(colorBlindness::paletteMartin, "#4AAA8A", "#D7442E", "#A9998F"),
                              unique(tissue_type_vec$`Tissue.type.(higher.degree)`))
swatch2(tissue_type_color)
tissue_type_vec <- merge(tissue_type_color, tissue_type_vec, by.x = 0, by.y = "Tissue.type.(higher.degree)")
colnames(tissue_type_vec) <- c("tissue_type", "type_color", "shortID")
tissue_type_vec <- tissue_type_vec[ order(tissue_type_vec$shortID, decreasing = FALSE) , ]

# need to eliminate NA values to allow clustering
# otherwise need to cluster separately and then provide the hclust object
heatmapData <- aging_atlas_cors
heatmapData[is.na(heatmapData)] <- 0

tissue_type_color_vec <- tissue_type_vec$type_color[match(colnames(heatmapData), tissue_type_vec$shortID)]

# are any of these genes differentially expressed in our dataset?
rownames(heatmapData)[rownames(heatmapData) %in% DE_hpk1$GeneSymbol]

pdf(file = "Correlation heatmap with hpk-1 across age for cell clusters x genes.pdf", width = 10, height = 14)
heatmap.3(heatmapData,
          ColSideColors = as.matrix(tissue_type_color_vec),
          col = viridis(option = "cividis", n = 128, direction = 1), #scale = "row",
          cexCol = 0.75, cexRow = 1.1, margins = c(2,7), keysize = 0.8,
          labCol = NA,
          KeyValueName = "Expr corr across age")
graphics.off()

rm(heatmapData, aging_atlas_cors, tissue_type_color, tissue_type_color_vec)

#
# What if we go the other direction- what do the DE genes in our dataset look like in terms of
#   correlations with hpk-1 expression across age, perhaps with a slightly less conservative filter?
#
aging_atlas_cors <- Reduce(f = cbind, lapply(aging_atlas_within_tissue.cors, function(x){
  x$cor[x$p.value >= 0.01] = NA
  return(x$cor)
}))
colnames(aging_atlas_cors) <- names(aging_atlas_within_tissue.cors)
identical(aging_atlas_within_tissue.cors$`0_0`$gene_name, aging_atlas_within_tissue.cors$`143_0`$gene_name) # sanity check
rownames(aging_atlas_cors) <- aging_atlas_within_tissue.cors$`0_0`$gene_name # all should be the same

# remove any rows with NA in every column
table(rowSums(is.na(aging_atlas_cors)))
dim(aging_atlas_cors)
aging_atlas_cors <- aging_atlas_cors[ rowSums(is.na(aging_atlas_cors)) < ncol(aging_atlas_cors) , ]
dim(aging_atlas_cors)
# any columns that only have 1 result? (would be hpk-1 self-correlation)
aging_atlas_cors <- aging_atlas_cors[ , colSums(!is.na(aging_atlas_cors)) > 1 ]
dim(aging_atlas_cors)

# how many of the DE genes are still here? 
aging_atlas_cors <- aging_atlas_cors[ rownames(aging_atlas_cors) %in% DE_hpk1$GeneSymbol , ]
dim(aging_atlas_cors)

# filter rows
aging_atlas_cors <- aging_atlas_cors[ rowSums(!is.na(aging_atlas_cors)) > 4 , ]
dim(aging_atlas_cors)

# filter columns
aging_atlas_cors <- aging_atlas_cors[ , colSums(!is.na(aging_atlas_cors)) > 0 ]
dim(aging_atlas_cors) # much more reasonable

tissue_type_vec <- aging_atlas_clusters_all[,c("Tissue.type.(higher.degree)", "shortID")]
tissue_type_vec$`Tissue.type.(higher.degree)`[tissue_type_vec$`Tissue.type.(higher.degree)` == "various"] <- "various_or_unknown"
tissue_type_vec$`Tissue.type.(higher.degree)`[is.na(tissue_type_vec$`Tissue.type.(higher.degree)`)] <-  "various_or_unknown"

tissue_type_color <- setNames(c(colorBlindness::paletteMartin, "#4AAA8A", "#D7442E", "#A9998F"),
                              unique(tissue_type_vec$`Tissue.type.(higher.degree)`))
swatch2(tissue_type_color)
tissue_type_vec <- merge(tissue_type_color, tissue_type_vec, by.x = 0, by.y = "Tissue.type.(higher.degree)")
colnames(tissue_type_vec) <- c("tissue_type", "type_color", "shortID")
tissue_type_vec <- tissue_type_vec[ order(tissue_type_vec$shortID, decreasing = FALSE) , ]

# need to eliminate NA values to allow clustering
# otherwise need to cluster separately and then provide the hclust object
heatmapData <- aging_atlas_cors
heatmapData[is.na(heatmapData)] <- 0

tissue_type_color_vec <- tissue_type_vec$type_color[match(colnames(heatmapData), tissue_type_vec$shortID)]

# sort rows by # non-0 results high to low, and don't cluster
heatmapData <- heatmapData[ order(rowSums(heatmapData != 0), decreasing = TRUE) , ]

# add indication of foldchange in our dataset
geneCols <- setNames(DE_hpk1$DESeq2_log2FC, DE_hpk1$GeneSymbol)
geneCols <- geneCols[rownames(heatmapData)]
geneCols[geneCols > 0] <- "darkred"
geneCols[geneCols < 0] <- "darkblue"

pdf(file = "Correlation heatmap with hpk-1 based on DE genes in hpk-1 null across age for cell clusters x genes.pdf",
    width = 10, height = 14)
heatmap.3(heatmapData,
          ColSideColors = as.matrix(tissue_type_color_vec),
          RowSideColors = t(as.matrix(geneCols)),
          col = viridis(option = "cividis", n = 128, direction = 1), #scale = "row",
          cexCol = 0.75, cexRow = 1.2, margins = c(2,10), keysize = 0.8,
          labCol = NA, Rowv = FALSE, dendrogram = "col",
          KeyValueName = "Expr corr across age")
graphics.off()

rm(heatmapData, aging_atlas_cors, tissue_type_color, tissue_type_color_vec, geneCols)

#
# hm, there's a somewhat interesting pattern of correlations across cell clusters
#  negative correlation in a few, positive correlation in a few
#  and it's not specific to tissue type
# --> Is there any specificity to this gene set with that?
#  What if we just pick a random set of the same size?
# 1803 genes were in the overlap with the correlation set

# RANDOM GENE SET CONTROL
aging_atlas_cors <- Reduce(f = cbind, lapply(aging_atlas_within_tissue.cors, function(x){
  x$cor[x$p.value >= 0.01] = NA
  return(x$cor)
}))
colnames(aging_atlas_cors) <- names(aging_atlas_within_tissue.cors)
identical(aging_atlas_within_tissue.cors$`0_0`$gene_name, aging_atlas_within_tissue.cors$`143_0`$gene_name) # sanity check
rownames(aging_atlas_cors) <- aging_atlas_within_tissue.cors$`0_0`$gene_name # all should be the same

# remove any rows with NA in every column
table(rowSums(is.na(aging_atlas_cors)))
dim(aging_atlas_cors)
aging_atlas_cors <- aging_atlas_cors[ rowSums(is.na(aging_atlas_cors)) < ncol(aging_atlas_cors) , ]
dim(aging_atlas_cors)
# any columns that only have 1 result? (would be hpk-1 self-correlation)
aging_atlas_cors <- aging_atlas_cors[ , colSums(!is.na(aging_atlas_cors)) > 1 ]
dim(aging_atlas_cors)

# how many of the DE genes are still here? 
aging_atlas_cors <- aging_atlas_cors[ sample(rownames(aging_atlas_cors), size = 1803) , ]
dim(aging_atlas_cors)

# filter rows
aging_atlas_cors <- aging_atlas_cors[ rowSums(!is.na(aging_atlas_cors)) > 4 , ]
dim(aging_atlas_cors)

# filter columns
aging_atlas_cors <- aging_atlas_cors[ , colSums(!is.na(aging_atlas_cors)) > 0 ]
dim(aging_atlas_cors) 


tissue_type_vec <- aging_atlas_clusters_all[,c("Tissue.type.(higher.degree)", "shortID")]
tissue_type_vec$`Tissue.type.(higher.degree)`[tissue_type_vec$`Tissue.type.(higher.degree)` == "various"] <- "various_or_unknown"
tissue_type_vec$`Tissue.type.(higher.degree)`[is.na(tissue_type_vec$`Tissue.type.(higher.degree)`)] <-  "various_or_unknown"

tissue_type_color <- setNames(c(colorBlindness::paletteMartin, "#4AAA8A", "#D7442E", "#A9998F"),
                              unique(tissue_type_vec$`Tissue.type.(higher.degree)`))
swatch2(tissue_type_color)
tissue_type_vec <- merge(tissue_type_color, tissue_type_vec, by.x = 0, by.y = "Tissue.type.(higher.degree)")
colnames(tissue_type_vec) <- c("tissue_type", "type_color", "shortID")
tissue_type_vec <- tissue_type_vec[ order(tissue_type_vec$shortID, decreasing = FALSE) , ]

# need to eliminate NA values to allow clustering
# otherwise need to cluster separately and then provide the hclust object
heatmapData <- aging_atlas_cors
heatmapData[is.na(heatmapData)] <- 0

tissue_type_color_vec <- tissue_type_vec$type_color[match(colnames(heatmapData), tissue_type_vec$shortID)]

# sort rows by # non-0 results high to low, and don't cluster
heatmapData <- heatmapData[ order(rowSums(heatmapData != 0), decreasing = TRUE) , ]


pdf(file = "Correlation heatmap with hpk-1 RANDOM GENE SET CONTROL.pdf",
    width = 10, height = 14)
heatmap.3(heatmapData,
          ColSideColors = as.matrix(tissue_type_color_vec),
          #RowSideColors = t(as.matrix(geneCols)),
          col = viridis(option = "cividis", n = 128, direction = 1), #scale = "row",
          cexCol = 0.75, cexRow = 0.75, margins = c(2,10), keysize = 0.8,
          labCol = NA, Rowv = FALSE, dendrogram = "col",
          KeyValueName = "Expr corr across age")
graphics.off()

rm(heatmapData, aging_atlas_cors, tissue_type_color, tissue_type_color_vec)


#
# Are there genes that are both spatially and temporally correlated with hpk-1 expression?
#
# already filtered based on adj p < 0.01
hpk1_corrMat_spatial <- Reduce(f = cbind, x =lapply(hpk1_corr_pearson, function(dat){
  dat$cor
}))
colnames(hpk1_corrMat_spatial) <- names(hpk1_corr_pearson)
rownames(hpk1_corrMat_spatial) <- hpk1_corr_pearson$D1$gene_name

# keep results with at least two non-NA results this time
hpk1_corrMat_spatial <- hpk1_corrMat_spatial[ rowSums(!is.na(hpk1_corrMat_spatial)) >= 2  , ]
any(hpk1_corrMat_spatial < 0, na.rm = TRUE) # all values are positive
sum(abs(rowMeans(hpk1_corrMat_spatial, na.rm = TRUE)) > 0.35) 
hpk1_corrMat_spatial <- hpk1_corrMat_spatial[ abs(rowMeans(hpk1_corrMat_spatial, na.rm = TRUE)) > 0.35, ]
# set NA to 0
hpk1_corrMat_spatial[is.na(hpk1_corrMat_spatial)] <- 0

# across aging
hpk1_corrMat_aging <- Reduce(f = cbind, lapply(aging_atlas_within_tissue.cors, function(x){
  x$cor[x$p.value >= 0.0001] = NA
  return(x$cor)
}))
colnames(hpk1_corrMat_aging) <- names(aging_atlas_within_tissue.cors)
rownames(hpk1_corrMat_aging) <- aging_atlas_within_tissue.cors$`0_0`$gene_name # all should be the same

# remove any rows with NA in every column
hpk1_corrMat_aging <- hpk1_corrMat_aging[ rowSums(is.na(hpk1_corrMat_aging)) < ncol(hpk1_corrMat_aging) , ]
# any columns that only have 1 result? (would be hpk-1 self-correlation)
hpk1_corrMat_aging <- hpk1_corrMat_aging[ , colSums(!is.na(hpk1_corrMat_aging)) > 1 ]
dim(hpk1_corrMat_aging)


heatmapData1 <- hpk1_corrMat_spatial[intersect(rownames(hpk1_corrMat_aging), rownames(hpk1_corrMat_spatial)) , ]
heatmapData2 <- hpk1_corrMat_aging[intersect(rownames(hpk1_corrMat_aging), rownames(hpk1_corrMat_spatial)) , ]

# remove columns that are all NA for this subset
heatmapData1 <- heatmapData1[ , colSums(!is.na(heatmapData1)) > 1 ]
heatmapData2 <- heatmapData2[ , colSums(!is.na(heatmapData2)) > 1 ]
dim(heatmapData1)
dim(heatmapData2)

# tissue annotation for across aging
tissue_type_vec <- aging_atlas_clusters_all[,c("Tissue.type.(higher.degree)", "shortID")]
tissue_type_vec$`Tissue.type.(higher.degree)`[tissue_type_vec$`Tissue.type.(higher.degree)` == "various"] <- "various_or_unknown"
tissue_type_vec$`Tissue.type.(higher.degree)`[is.na(tissue_type_vec$`Tissue.type.(higher.degree)`)] <-  "various_or_unknown"

tissue_type_color <- setNames(c(colorBlindness::paletteMartin, "#4AAA8A", "#D7442E", "#A9998F"),
                              unique(tissue_type_vec$`Tissue.type.(higher.degree)`))
tissue_type_vec <- merge(tissue_type_color, tissue_type_vec, by.x = 0, by.y = "Tissue.type.(higher.degree)")
colnames(tissue_type_vec) <- c("tissue_type", "type_color", "shortID")
tissue_type_vec <- tissue_type_vec[ order(tissue_type_vec$shortID, decreasing = FALSE) , ]
tissue_type_color_vec <- tissue_type_vec$type_color[match(colnames(heatmapData2), tissue_type_vec$shortID)]

# set NA values to 0
heatmapData1[is.na(heatmapData1)] <- 0
heatmapData2[is.na(heatmapData2)] <- 0

pdf(file = "Correlation heatmap with hpk-1- intersect of spatial and age for spatial.pdf",
    width = 10, height = 14)
heatmap.3(heatmapData1,
          col = viridis(option = "magma", n = 128, direction = 1), #scale = "row",
          cexCol = 1, cexRow = 1.5, margins = c(8,10), keysize = 0.8,
           Rowv = FALSE, dendrogram = "col",
          KeyValueName = "Expr corr across cell cluters")
graphics.off()

pdf(file = "Correlation heatmap with hpk-1- intersect of spatial and age for age.pdf",
    width = 10, height = 14)
heatmap.3(heatmapData2,
          ColSideColors = as.matrix(tissue_type_color_vec),
          col = viridis(option = "cividis", n = 128, direction = 1), #scale = "row",
          cexCol = 0.75, cexRow = 1.5, margins = c(10,10), keysize = 0.8,
          labCol = NA, Rowv = FALSE, dendrogram = "col",
          KeyValueName = "Expr corr across age")
graphics.off()



#
# Looking for correspondence among known aging-associated TFs
#  --> first of all, across aging, within neurons
# Our starting list of primary interest:
# "hpk-1", "hlh-30", "daf-16", "skn-1", "hsf-1", "nhr-49", "pha-4", "nhr-62", "mml-1", "mxl-2", "mdl-1", "mxl-1", ("mxl-3"), ("daf-19")
# 
# Aging atlas paper, figure 7A- TFs "most up-regulated with age": (NOTE: SOME OVERLAP HERE)
goi_atlas_age_up_tfs <- c("nhr-25", "eyg-1", "ztf-14", "hlh-30", "daf-16", "skn-1", "gei-3", "daf-12", "elt-7", "fkh-9", "fkh-7", "dpy-27",
                          "efl-2", "nhr-128", "ztf-29", "nhr-79", "sem-2", "nhr-23", "lsl-1", "egl-5")
# Aging atlas paper, figure 7B- TFs "most down-regulated with age":
goi_atlas_age_down_tfs <- c("fkh-2", "ces-1", "che-1", "ceh-53", "ceh-28", "ceh-36", "xbp-1", "dsc-1", "pha-2", "ces-2", "vab-7", "lin-32",
                            "odd-2", "dmd-4", "ceh-34", "ceh-40", "unc-4", "ceh-22", "ceh-6", "ceh-33")
#
# Andy saw mention of these TFs on the aging atlas site and thought it was from an aging-associated analysis, but it's
#  not, they're from figure 2 in the paper which was based on estimating TF activity based on variability between
#  the TF expression and expression of predicted target genes *after regressing out age as a factor*.
goi_atlas_tf_pred_high_activity <- c("crh-1", "dpy-27", "mef-2", "jun-1", "gei-3", "daf-12", "skn-1", "sptf-3",
                                     "fkn-7", "hlh-30", "egl-27,", "daf-16", "nhr-34", "atf-7", "nhr-142", "nhr-46",
                                     "nhr-66", "fkn-9", "daf-19", "lin-1", "atf-2", "nhr-84", "nhr-86",
                                     "mdl-1,", "syd-9", "nhr-49", "dmd-6", "nhr-3", "mxl-3", "xbp-1", "egl-13", "ztf-3", "mnm-2",
                                     "lin-11", "Y44A6D.3", "dsc-1", "ceh-20", "dmd-4", "hlh-8", "ceh-9", "nhr-67", "cey-1",
                                     "atf-5", "vab-7", "pag-3", "mab-23")
#  We're interested in looking at all of these in conjunction with hpk-1
#
# Note that the TF activity correlations/predictions in Figure 2D in their paper are across cells of each cluster,
#  not considering differences through age (cells of each age were integrated for the analysis)

ggVennDiagram::ggVennDiagram(list("Aging atlas\nTFs down with age" = goi_atlas_age_down_tfs,
                                  "Aging atlas\nTFs up with age" = goi_atlas_age_up_tfs,
                                  "Aging atlas\nTFs with high predicted\nactivity within clusters" = goi_atlas_tf_pred_high_activity))

#
# I have also now mapped the aging atlas neurons to the cengen neurons, and the associated annotations
#  so we can use all of that information here.
# I also have predicted neurotransmitter associations, but since a neuron can be associated with more than one, 
#  we need to have a color indication for each (glutamatergic, serotonergic, GABAergic, dopaminergic, cholinergic)
#

heatmapData <- Reduce(f = cbind, lapply(aging_atlas_within_tissue.cors, function(x){
  x$cor[x$p.value >= 0.1] = NA
  return(x$cor)
}))
colnames(heatmapData) <- names(aging_atlas_within_tissue.cors)
rownames(heatmapData) <- aging_atlas_within_tissue.cors$`0_0`$gene_name # all should be the same
heatmapData[is.na(heatmapData)] <- 0

genesOfInterest <- c("hpk-1", "hlh-30", "daf-16", "skn-1", "hsf-1", "nhr-49", "pha-4", "nhr-62", "mml-1", "mxl-2", "mdl-1", "mxl-1")

atlas_cengen_neuron_map.sub <- atlas_cengen_neuron_map[ atlas_cengen_neuron_map$atlas_cluster %in% colnames(heatmapData)  , ]

# just neurons and genes of interest
heatmapData <- heatmapData[ genesOfInterest , atlas_cengen_neuron_map.sub$atlas_cluster ]

tmp_type = atlas_cengen_neuron_map.sub$NeuronType
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]

swatch2(c("Interneuron" = as.vector(colorBlindness::paletteMartin[2]),
        "Sensory" = as.vector(colorBlindness::paletteMartin[4]),
        "Endorgan" = as.vector(colorBlindness::paletteMartin[6]),
        "Motorneuron" = as.vector(colorBlindness::paletteMartin[8])))

tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

pdf(file = "Correlation heatmap with hpk-1 across age for aging TFs on neuron cell clusters x genes.pdf", width = 10, height = 12)
heatmap.3(heatmapData,
          ColSideColors = neuron_colors,
          col = viridis(option = "cividis", n = 128, direction = 1), #scale = "row",
          cexCol = 0.75, cexRow = 1.5, margins = c(3,7), keysize = 0.8,
          labCol = NA, main = "Correlations across age with hpk-1, Pearson, adj p < 0.1", 
          KeyValueName = "Corr across age")
graphics.off()

rm(heatmapData, atlas_cengen_neuron_map.sub, tmp_cholinergic,
   tmp_dopaminergic, tmp_GABAergic, tmp_glutamatergic, tmp_serotonergic, tmp_type, genesOfInterest)

# 6_2 shows positive corr in a few cases, 114_0 negative in a few case, plot some examples

# colors- younger 3 black, older 3 red.
plot(x = as.numeric(aging_atlas_within_tissue$`6_2`["hpk-1" , ]) ,
     y = as.numeric(aging_atlas_within_tissue$`6_2`["skn-1" , ]),
     xlab = "hpk-1", ylab = "skn-1", main = "Cluster 6_2", pch = 3, col = c(rep("black", 3), rep("red", 3)))

plot(x = as.numeric(aging_atlas_within_tissue$`6_2`["hpk-1" , ]) ,
     y = as.numeric(aging_atlas_within_tissue$`6_2`["daf-16" , ]),
     xlab = "hpk-1", ylab = "daf-16", main = "Cluster 6_2", pch = 3, col = c(rep("black", 3), rep("red", 3)))

plot(x = as.numeric(aging_atlas_within_tissue$`114_0`["hpk-1" , ]) ,
     y = as.numeric(aging_atlas_within_tissue$`114_0`["daf-16" , ]),
     xlab = "hpk-1", ylab = "daf-16", main = "Cluster 114_0", pch = 3, col = c(rep("black", 3), rep("red", 3)))

plot(x = as.numeric(aging_atlas_within_tissue$`114_0`["hpk-1" , ]) ,
     y = as.numeric(aging_atlas_within_tissue$`114_0`["hsf-1" , ]),
     xlab = "hpk-1", ylab = "hsf-1", main = "Cluster 114_0", pch = 3, col = c(rep("black", 3), rep("red", 3)))


#
# Do the same with their differential expression results
#
# heatmapData_list <- aging_atlas_degs[names(aging_atlas_degs) %in% atlas_cengen_neuron_map$atlas_cluster]

#
# WITH ADDITIONAL SIGNIFICANCE FILTERING
# 
#heatmapData_long <- aging_atlas_degs_long.unfilt[ aging_atlas_degs_long.unfilt$cluster_shortID %in% atlas_cengen_neuron_map$atlas_cluster , ]
heatmapData_long <- aging_atlas_degs_long[ aging_atlas_degs_long$cluster_shortID %in% atlas_cengen_neuron_map$atlas_cluster , ]
heatmapData <- t(reshape2::acast(heatmapData_long[,1:3], cluster_shortID ~ gene_name, drop = FALSE, fill = 0))
# andy actually wants to see the neuron identity, so we'll actually have to combine 
#   neuron names for the corresponding clusters, since they don't all have 1:1 equivalence
atlas_cengen_neuron_map.sub <- unique(atlas_cengen_neuron_map[ atlas_cengen_neuron_map$atlas_cluster %in% colnames(heatmapData)  , 1:4 ])
atlas_cengen_neuron_map.sub_agg <- Reduce(f = rbind, lapply(unique(atlas_cengen_neuron_map.sub$atlas_cluster), function(clustID){
  if(sum(atlas_cengen_neuron_map.sub$atlas_cluster == clustID) == 1){
    return(atlas_cengen_neuron_map.sub[ atlas_cengen_neuron_map.sub$atlas_cluster == clustID , ]) # return existing
  }else{ # multiple cells for cluster, need to aggregate
    return(with(atlas_cengen_neuron_map.sub[ atlas_cengen_neuron_map.sub$atlas_cluster == clustID , ],
         data.frame(CengenName = paste(unique(CengenName), sep = ", ", collapse = ", "),
                    NeuronType = paste(unique(NeuronType), sep = ", ", collapse = ", "),
                    atlas_cluster = clustID,
                    cengen_signal_class = paste(unique(cengen_signal_class), sep = ", ", collapse = ", "))
         ))
  }
}))

genesOfInterest <- c("hpk-1", "hlh-30", "daf-16", "skn-1", "hsf-1", "nhr-49", "pha-4", "nhr-62", "mml-1", "mxl-2", "mdl-1", "mxl-1")

genesOfInterest <- genesOfInterest[genesOfInterest %in% heatmapData_long$gene_name]

heatmapData <- heatmapData[ genesOfInterest  , atlas_cengen_neuron_map.sub_agg$atlas_cluster ]

tmp_type = atlas_cengen_neuron_map.sub_agg$NeuronType
# now we're aggregating neurons, so some clusters don't represent just one type
tmp_type[!(tmp_type %in% c("Interneuron", "Sensory", "Endorgan", "Motorneuron")) ] = "white"
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]


tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

pdf(file = "Aging atlas DE across age for aging TFs in neurons (Adam version- with filtering) v1.1.pdf", width = 16, height = 12)
heatmap.3(heatmapData,
          ColSideColors = neuron_colors,
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 0.75, cexRow = 1.5, margins = c(8,7), keysize = 0.8,
          labCol = atlas_cengen_neuron_map.sub_agg$CengenName, # neuron names in cengen
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5, neurons only",
          KeyValueName = "LogFoldChange")
graphics.off()

rm(heatmapData, atlas_cengen_neuron_map.sub, tmp_cholinergic, atlas_cengen_neuron_map.sub_agg,
   tmp_dopaminergic, tmp_GABAergic, tmp_glutamatergic, tmp_serotonergic, tmp_type, heatmapData_long)

#
#  WITHOUT ADDITIONAL SIGNIFICANCE FILTERING!!
#
heatmapData_long <- aging_atlas_degs_long.unfilt[ aging_atlas_degs_long.unfilt$cluster_shortID %in% atlas_cengen_neuron_map$atlas_cluster , ]
heatmapData <- t(reshape2::acast(heatmapData_long[,1:3], cluster_shortID ~ gene_name, drop = FALSE, fill = 0))
atlas_cengen_neuron_map.sub <- unique(atlas_cengen_neuron_map[ atlas_cengen_neuron_map$atlas_cluster %in% colnames(heatmapData)  , 2:4 ])

genesOfInterest <- c("hpk-1", "hlh-30", "daf-16", "skn-1", "hsf-1", "nhr-49", "pha-4", "nhr-62", "mml-1", "mxl-2", "mdl-1", "mxl-1")

genesOfInterest <- genesOfInterest[genesOfInterest %in% heatmapData_long$gene_name]

heatmapData <- heatmapData[ genesOfInterest  , atlas_cengen_neuron_map.sub$atlas_cluster ]

tmp_type = atlas_cengen_neuron_map.sub$NeuronType
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]

tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub$cengen_signal_class),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

pdf(file = "Aging atlas DE across age for aging TFs in neurons (Adam version- no filtering).pdf", width = 10, height = 12)
heatmap.3(heatmapData,
          ColSideColors = neuron_colors,
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 0.75, cexRow = 1.5, margins = c(3,7), keysize = 0.8,
          labCol = NA, main = "Aging atlas DE, no filter- as provided, neurons only",
          KeyValueName = "LogFoldChange")
graphics.off()

rm(heatmapData, atlas_cengen_neuron_map.sub, tmp_cholinergic,
   tmp_dopaminergic, tmp_GABAergic, tmp_glutamatergic, tmp_serotonergic, tmp_type)
#
# ---- aging DE plots for Other genes of interest: TFs most up and down regulated with age ----
#
heatmapData_long <- aging_atlas_degs_long[ aging_atlas_degs_long$cluster_shortID %in% atlas_cengen_neuron_map$atlas_cluster , ]
heatmapData <- t(reshape2::acast(heatmapData_long[,1:3], cluster_shortID ~ gene_name, drop = FALSE, fill = 0))
# andy actually wants to see the neuron identity, so we'll actually have to combine 
#   neuron names for the corresponding clusters, since they don't all have 1:1 equivalence
atlas_cengen_neuron_map.sub <- unique(atlas_cengen_neuron_map[ atlas_cengen_neuron_map$atlas_cluster %in% colnames(heatmapData)  , 1:4 ])
atlas_cengen_neuron_map.sub_agg <- Reduce(f = rbind, lapply(unique(atlas_cengen_neuron_map.sub$atlas_cluster), function(clustID){
  if(sum(atlas_cengen_neuron_map.sub$atlas_cluster == clustID) == 1){
    return(atlas_cengen_neuron_map.sub[ atlas_cengen_neuron_map.sub$atlas_cluster == clustID , ]) # return existing
  }else{ # multiple cells for cluster, need to aggregate
    return(with(atlas_cengen_neuron_map.sub[ atlas_cengen_neuron_map.sub$atlas_cluster == clustID , ],
                data.frame(CengenName = paste(unique(CengenName), sep = ", ", collapse = ", "),
                           NeuronType = paste(unique(NeuronType), sep = ", ", collapse = ", "),
                           atlas_cluster = clustID,
                           cengen_signal_class = paste(unique(cengen_signal_class), sep = ", ", collapse = ", "))
    ))
  }
}))

genesOfInterest <- c(goi_atlas_age_up_tfs, goi_atlas_age_down_tfs, "hpk-1")
genesOfInterest <- genesOfInterest[genesOfInterest %in% heatmapData_long$gene_name]

heatmapData <- heatmapData[ genesOfInterest  , atlas_cengen_neuron_map.sub_agg$atlas_cluster ]

tmp_type = atlas_cengen_neuron_map.sub_agg$NeuronType
# now we're aggregating neurons, so some clusters don't represent just one type
tmp_type[!(tmp_type %in% c("Interneuron", "Sensory", "Endorgan", "Motorneuron")) ] = "white"
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]


tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

pdf(file = "Aging atlas DE across age for tfs most up or down with age (fig 7A+B) adam filtered.pdf", width = 16, height = 12)
heatmap.3(heatmapData,
          ColSideColors = neuron_colors,
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 0.75, cexRow = 1.5, margins = c(8,10), keysize = 0.8,
          labCol = atlas_cengen_neuron_map.sub_agg$CengenName, # neuron names in cengen
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5, neurons only",
          KeyValueName = "LogFoldChange")
graphics.off()

rm(heatmapData, atlas_cengen_neuron_map.sub, tmp_cholinergic, atlas_cengen_neuron_map.sub_agg,
   tmp_dopaminergic, tmp_GABAergic, tmp_glutamatergic, tmp_serotonergic, tmp_type, heatmapData_long)

#
# ----  aging DE plots for Other genes of interest: TFs with highest predicted activity within cell clusters ----
#
heatmapData_long <- aging_atlas_degs_long[ aging_atlas_degs_long$cluster_shortID %in% atlas_cengen_neuron_map$atlas_cluster , ]
heatmapData <- t(reshape2::acast(heatmapData_long[,1:3], cluster_shortID ~ gene_name, drop = FALSE, fill = 0))
# andy actually wants to see the neuron identity, so we'll actually have to combine 
#   neuron names for the corresponding clusters, since they don't all have 1:1 equivalence
atlas_cengen_neuron_map.sub <- unique(atlas_cengen_neuron_map[ atlas_cengen_neuron_map$atlas_cluster %in% colnames(heatmapData)  , 1:4 ])
atlas_cengen_neuron_map.sub_agg <- Reduce(f = rbind, lapply(unique(atlas_cengen_neuron_map.sub$atlas_cluster), function(clustID){
  if(sum(atlas_cengen_neuron_map.sub$atlas_cluster == clustID) == 1){
    return(atlas_cengen_neuron_map.sub[ atlas_cengen_neuron_map.sub$atlas_cluster == clustID , ]) # return existing
  }else{ # multiple cells for cluster, need to aggregate
    return(with(atlas_cengen_neuron_map.sub[ atlas_cengen_neuron_map.sub$atlas_cluster == clustID , ],
                data.frame(CengenName = paste(unique(CengenName), sep = ", ", collapse = ", "),
                           NeuronType = paste(unique(NeuronType), sep = ", ", collapse = ", "),
                           atlas_cluster = clustID,
                           cengen_signal_class = paste(unique(cengen_signal_class), sep = ", ", collapse = ", "))
    ))
  }
}))

genesOfInterest <- c(goi_atlas_tf_pred_high_activity, "hpk-1")

genesOfInterest <- genesOfInterest[genesOfInterest %in% heatmapData_long$gene_name]

heatmapData <- heatmapData[ genesOfInterest  , atlas_cengen_neuron_map.sub_agg$atlas_cluster ]

tmp_type = atlas_cengen_neuron_map.sub_agg$NeuronType
# now we're aggregating neurons, so some clusters don't represent just one type
tmp_type[!(tmp_type %in% c("Interneuron", "Sensory", "Endorgan", "Motorneuron")) ] = "white"
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]


tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub_agg$cengen_signal_class),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

pdf(file = "Aging atlas DE across age for TFs with high predicted activity (from preprint fig 2).pdf", width = 16, height = 12)
heatmap.3(heatmapData,
          ColSideColors = neuron_colors,
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 0.75, cexRow = 1.5, margins = c(8,7), keysize = 0.8,
          labCol = atlas_cengen_neuron_map.sub_agg$CengenName, # neuron names in cengen
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5, neurons only",
          KeyValueName = "LogFoldChange")
graphics.off()

rm(heatmapData, atlas_cengen_neuron_map.sub, tmp_cholinergic, atlas_cengen_neuron_map.sub_agg,
   tmp_dopaminergic, tmp_GABAergic, tmp_glutamatergic, tmp_serotonergic, tmp_type, heatmapData_long)


#
# Another idea Andy had was to try to demonstrate if hpk-1 is distinct as a kinase in its
#  being transcriptionally regulated as a function of age. So let's see how many other
#  kinases are differentially expressed with aging....
# I got this list from Zaru et al, which was an update of the Gerard Wormbook article from 2005
#

elegans_kinases <- read.xlsx(file.path("E:/project_data_ssd/hpk-1_dataset/kinases",
                                       "Zaru et al elegans kinome supplement S10.xlsx"), sheet = 1)
elegans_kinases <- elegans_kinases[,c("Group", "Family", "Gene.sequence.name")]
elegans_kinases <- merge(elegans_kinases, wb_gene_info[,c("gene_name", "molecular_name")],
                         by.x = "Gene.sequence.name", by.y = "molecular_name", all.x = TRUE)
nrow(elegans_kinases)
elegans_kinases <- elegans_kinases[ !is.na(elegans_kinases$gene_name) , ]
nrow(elegans_kinases)

# how many genes per kinase group?
table(elegans_kinases$Group)

# so... Andy wanted a heatmap with just all of them.
#  the least we can do is also indicate the class 
aging_kinase_degs <- aging_atlas_degs_long[ aging_atlas_degs_long$gene_name %in% elegans_kinases$gene_name , ]
any(aging_kinase_degs$padj > 0.05)
sum((aging_kinase_degs$logFoldChange > 0 ) & (aging_kinase_degs$logFoldChange < 0.5 ))
# transform to a gene x cluster matrix of foldchanges, since these are all significant results
aging_kinase_degs_mat <- mat.or.vec(nr = length(unique(aging_kinase_degs$gene_name)), nc = length(unique(aging_kinase_degs$cluster_shortID)))
rownames(aging_kinase_degs_mat) <- unique(aging_kinase_degs$gene_name)
colnames(aging_kinase_degs_mat) <- unique(aging_kinase_degs$cluster_shortID)

for(i in rownames(aging_kinase_degs_mat)){
  for(j in colnames(aging_kinase_degs_mat)){
    if(any(with(aging_kinase_degs, (cluster_shortID == j) & (gene_name == i) ))){
    aging_kinase_degs_mat[i,j] <- with(aging_kinase_degs, logFoldChange[(cluster_shortID == j) & (gene_name == i)] )
    }else{
      aging_kinase_degs_mat[i,j] <- 0 # I'd feel better with NA but for now this is ok
    }
  }
}
rm(i, j)

sum(aging_kinase_degs_mat["hpk-1" , ] != 0)

library(viridis)
library(hues)

swatch2 <- function(x){
  par(mai = c(0.2, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 
              0.2, 0.4))
  barplot(rep(1, length(x)), col = rev(x), space = 0.1, axes = FALSE, 
          names.arg = rev(names(x)), cex.names = 0.8, horiz = T, las = 1)
  return(invisible(NULL))
}

# create color vector for kinase classes
kinase_class_color <- iwanthue(length(unique(elegans_kinases$Group)))
names(kinase_class_color) <- unique(elegans_kinases$Group)
swatch2(kinase_class_color)

kinase_class_color_vec <- kinase_class_color[match(elegans_kinases$Group[match(rownames(aging_kinase_degs_mat), elegans_kinases$gene_name)],
      names(kinase_class_color))]

ncol(aging_kinase_degs_mat)
aging_atlas_clusters_all$`Tissue.type.(higher.degree)`[aging_atlas_clusters_all$shortID %in% colnames(aging_kinase_degs_mat)]

# indicate broad cell type for cell cluster
tissue_type_vec <- aging_atlas_clusters_all[,c("Tissue.type.(higher.degree)", "shortID")]
tissue_type_vec$`Tissue.type.(higher.degree)`[tissue_type_vec$`Tissue.type.(higher.degree)` == "various"] <- "various_or_unknown"
tissue_type_vec$`Tissue.type.(higher.degree)`[is.na(tissue_type_vec$`Tissue.type.(higher.degree)`)] <-  "various_or_unknown"

tissue_type_color <- setNames(c(colorBlindness::paletteMartin, "#4AAA8A", "#D7442E", "#A9998F"),
                              unique(tissue_type_vec$`Tissue.type.(higher.degree)`))
#swatch2(tissue_type_color)
tissue_type_vec <- merge(tissue_type_color, tissue_type_vec, by.x = 0, by.y = "Tissue.type.(higher.degree)")
colnames(tissue_type_vec) <- c("tissue_type", "type_color", "shortID")

setdiff(colnames(aging_kinase_degs_mat), tissue_type_vec$shortID) # SHOULD BE NO ITEMS IN DIFF

# re-order so that all the neurons are together, move all neurons to one end
#  KEEP tissue_type_vec WE CAN USE IT AGAIN and avoid redundancy
tissue_type_vec <- tissue_type_vec[ order(tissue_type_vec$tissue_type, decreasing = FALSE) , ]
tissue_type_vec <- rbind(tissue_type_vec[ tissue_type_vec$tissue_type != "neuron"  , ], tissue_type_vec[ tissue_type_vec$tissue_type == "neuron"  , ])

setdiff(colnames(aging_kinase_degs_mat), tissue_type_vec$shortID) # SHOULD BE NO ITEMS IN DIFF

svg(file = "Elegans kinome and aging differential expression in cell clusters v3.svg", width = 10, height = 25)
heatmap.3(aging_kinase_degs_mat[,tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(aging_kinase_degs_mat)]],
          RowSideColors = t(as.matrix(kinase_class_color_vec)),
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(aging_kinase_degs_mat)]),
          dendrogram = "row", Colv = FALSE,
          col = colorRampPalette(c("blue4", "white", "red4")), #scale = "row",
          cexCol = 0.75, cexRow = 0.75, margins = c(2,7), keysize = 0.8,
          labCol = NA,
           KeyValueName = "Aging atlas old vs yound LFC")

graphics.off()

pdf(file = "Elegans kinome and aging differential expression in cell clusters COLOR KEY.pdf", width = 2, height = 3)
swatch2(kinase_class_color)
swatch2(with(unique(tissue_type_vec[,c("type_color", "tissue_type")]), setNames(type_color, tissue_type)))
graphics.off()


#
# For Andy's reference- check DE of TOR components (along with hpk-1):
#   rsks-1, aak-1, smg-1, let-363, mlst-8, raga-1, ragc-1, rheb-1, nprl-2, nprl-3, rict-1, sinh-1
#
tor_related_genes_and_hpk1 <- c("rsks-1", "aak-1", "smg-1", "let-363", "mlst-8", "raga-1", "ragc-1",
                       "rheb-1", "nprl-2", "nprl-3", "rict-1", "sinh-1", "hpk-1")

aging_tor_degs <- aging_atlas_degs_long[ aging_atlas_degs_long$gene_name %in% tor_related_genes_and_hpk1 , ]
# transform to a gene x cluster matrix of foldchanges, since these are all significant results
aging_tor_degs_mat <- mat.or.vec(nr = length(unique(aging_tor_degs$gene_name)), nc = length(unique(aging_tor_degs$cluster_shortID)))
rownames(aging_tor_degs_mat) <- unique(aging_tor_degs$gene_name)
colnames(aging_tor_degs_mat) <- unique(aging_tor_degs$cluster_shortID)

for(i in rownames(aging_tor_degs_mat)){
  for(j in colnames(aging_tor_degs_mat)){
    if(any(with(aging_tor_degs, (cluster_shortID == j) & (gene_name == i) ))){
      aging_tor_degs_mat[i,j] <- with(aging_tor_degs, logFoldChange[(cluster_shortID == j) & (gene_name == i)] )
    }else{
      aging_tor_degs_mat[i,j] <- 0 # I'd feel better with NA but for now this is ok
    }
  }
}
rm(i, j)


setdiff(colnames(aging_tor_degs_mat), tissue_type_vec$shortID) # SHOULD BE NO ITEMS IN DIFF

pdf(file = "Elegans TOR-related genes aging differential expression in cell clusters.pdf", width = 10, height = 12)
heatmap.3(aging_tor_degs_mat[,tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(aging_tor_degs_mat)]],
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(aging_tor_degs_mat)]),
          dendrogram = "row", Colv = FALSE,
          col = colorRampPalette(c("blue4", "white", "red4")), #scale = "row",
          cexCol = 0.75, cexRow = 0.75, margins = c(2,7), keysize = 0.8,
          labCol = NA,
          colsep = seq(1, ncol(aging_tor_degs_mat), by = 1),
          rowsep = seq(1, nrow(aging_tor_degs_mat), by = 1),
          sepcolor = "grey90",
          KeyValueName = "Aging atlas old vs yound LFC")

graphics.off()



# for reference- how many cases of up/down regulation for each?
aging_kinase_degs_occurrence <- data.frame("num_up_cluster" = rowSums(aging_kinase_degs_mat > 0),
                                           "num_down_cluster" = rowSums(aging_kinase_degs_mat < 0))
write.xlsx(aging_kinase_degs_occurrence, file = "Elegans kinome aging DE occurrence table.xlsx", rowNames = TRUE)

# add more info to this: WBGeneID, longevity associations
aging_kinase_degs_occurrence_2 <- merge(aging_kinase_degs_occurrence, wb_gene_info[,1:2], by.x = 0, by.y = 2,  all.x = TRUE)
# !!!!!!!!!!!!!!!!!!!!!!
# TODO- USE NEWER WORMBASE LIFESPAN GENE LIST
# TODO- ADD new "manually" curated additional lifespan genes
# !!!!!!!!!!!!!!!!!!!!!!
aging_kinase_degs_occurrence_2 <- merge(aging_kinase_degs_occurrence_2, lifespan_genes, by.x = 1, by.y = 1, all.x = TRUE)
aging_kinase_degs_occurrence_2$in_pgp <- aging_kinase_degs_occurrence_2$WBGeneID %in% pgpList$x
aging_kinase_degs_occurrence_2 <- merge(aging_kinase_degs_occurrence_2, genage_simple, by.x = 1, by.y = 1, all.x = TRUE)
colnames(aging_kinase_degs_occurrence_2)[c(1,5)] <- c("gene_name","wormbase_lifespan_annotation")
write.xlsx(aging_kinase_degs_occurrence_2, file = "Elegans kinome aging DE occurrence table (v2 with lifespan info).xlsx", rowNames = TRUE)

#
# summarize the aging-associated gene info more broadly, since we'll use it again
#

# added additional gerogenes from recent papers
#  and for now just a true/false indication of wormbase lifespan phenotype as of ws282

aging_genes_merged <- merge(wb_gene_info[,1:2], lifespan_genes, by.x = 2, by.y = 1, all.x = TRUE)
head(aging_genes_merged)
aging_genes_merged <- merge(aging_genes_merged, genage_simple, by.x = "gene_name", by.y = 1, all.x = TRUE)
head(aging_genes_merged)
aging_genes_merged$inPGP <- aging_genes_merged$WBGeneID %in% pgpList$x
aging_genes_merged$lifespanPheno_WB.WS282 <- aging_genes_merged$WBGeneID %in% lifespan_genes_v2$WBGeneID
aging_genes_merged$new_lifespan_genes <- aging_genes_merged$WBGeneID %in% lifespan_genes_additions$WBGeneID
head(aging_genes_merged)

#
# This is any gene:
#  1) associated with aging phenotypes in wormbase
#  2) in genage
#  3) in the PGP
#
aging_genes_merged$any_aging_pheno <- (!is.na(aging_genes_merged$wb_annotation) |
                                         !is.na(aging_genes_merged$genage_gerogene_cats) |
                                         aging_genes_merged$inPGP |
                                         aging_genes_merged$lifespanPheno_WB.WS282 |
                                         aging_genes_merged$new_lifespan_genes)
sum(aging_genes_merged$any_aging_pheno)

tmp_aging_genes <- setNames(aging_genes_merged$any_aging_pheno, aging_genes_merged$gene_name)
tmp_aging_genes <- tmp_aging_genes[rownames(aging_kinase_degs_mat)]
identical(names(tmp_aging_genes), rownames(aging_kinase_degs_mat))

svg(file = "Elegans kinome and aging differential expression in cell clusters v3.4.svg", width = 10, height = 25)
heatmap.3(aging_kinase_degs_mat[,tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(aging_kinase_degs_mat)]],
          RowSideColors = rbind("kinase_class" = kinase_class_color_vec,
                                "aging_implicated" = ifelse(tmp_aging_genes, "gold4", "black")),
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(aging_kinase_degs_mat)]),
          dendrogram = "row", Colv = FALSE,
          col = colorRampPalette(c("blue4", "white", "red4")), #scale = "row",
          cexCol = 0.75, cexRow = 0.75, margins = c(10,7), keysize = 0.8,
          labCol = NA,
          KeyValueName = "Aging atlas old vs yound LFC")

graphics.off()

rm(tmp_aging_genes)

table(rowSums(aging_kinase_degs_mat != 0))
#
# Filtering kinase result for main text figure
#  DE in ten or more clusters (Since Andy must have counted them and wrote it that way)
#
dim(aging_kinase_degs_mat)

# remove kinases DE in less than 10 clusters
aging_kinase_degs_mat.filt <- aging_kinase_degs_mat[ rowSums(aging_kinase_degs_mat > 0) >= 10 , ]
# remove empty columns
aging_kinase_degs_mat.filt <- aging_kinase_degs_mat.filt[ , colSums(aging_kinase_degs_mat.filt == 0) < nrow(aging_kinase_degs_mat.filt) ]
dim(aging_kinase_degs_mat.filt)

# order by average across rows
#aging_kinase_degs_mat.filt <- aging_kinase_degs_mat.filt[ order(rowMeans(aging_kinase_degs_mat.filt), decreasing = TRUE) , ]
# order by NUMBER OF UPREGULATED occurrences
aging_kinase_degs_mat.filt <- aging_kinase_degs_mat.filt[ order(rowSums(aging_kinase_degs_mat.filt > 0), decreasing = TRUE) , ]

tmp_aging_genes <- setNames(aging_genes_merged$any_aging_pheno, aging_genes_merged$gene_name)
tmp_aging_genes <- tmp_aging_genes[rownames(aging_kinase_degs_mat.filt)]
identical(names(tmp_aging_genes), rownames(aging_kinase_degs_mat.filt))

kinase_class_color_vec_filt <- kinase_class_color[match(elegans_kinases$Group[match(rownames(aging_kinase_degs_mat.filt), elegans_kinases$gene_name)],
                                                   names(kinase_class_color))]

svg(file = "Elegans kinases DE with age in GTE 10 clusters .svg", width = 12, height = 8)
heatmap.3(aging_kinase_degs_mat.filt[,tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(aging_kinase_degs_mat.filt)]],
          RowSideColors = rbind("kinase_class" = kinase_class_color_vec_filt,
                                "aging_implicated" = ifelse(tmp_aging_genes, "gold4", "black")),
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(aging_kinase_degs_mat.filt)]),
          dendrogram = "none", Colv = FALSE, Rowv = FALSE,
          col = colorRampPalette(c("blue4", "white", "red4")), #scale = "row",
          cexCol = 0.75, cexRow = 1.4, margins = c(10,7), keysize = 0.8,
          labCol = NA,
          rowsep = seq(1, nrow(aging_kinase_degs_mat.filt)),
          #colsep = seq(1, ncol(aging_kinase_degs_mat.filt)),
          sepcolor = "grey90",sepwidth = c(0.001, 0.001),
          KeyValueName = "Aging atlas old vs yound LFC")

graphics.off()

rm(tmp_aging_genes)

#
# Kinase table update 3
#
# add kinase family info

intersect(aging_kinase_degs_occurrence_2$gene_name, paper_SF06_2$GeneSymbol)
aging_kinase_degs_occurrence_3 <- merge(aging_kinase_degs_occurrence_2, elegans_kinases[,2:4], by.x = "gene_name", by.y = "gene_name", all.x = TRUE)
aging_kinase_degs_occurrence_3 <- merge(aging_kinase_degs_occurrence_3, paper_SF06_2[,2:3], by.x = "gene_name", by.y = "GeneSymbol", all.x = TRUE)
aging_kinase_degs_occurrence_3 <- merge(aging_kinase_degs_occurrence_3, aging_genes_merged[,c(1,6)],
                                        by.x = "gene_name", by.y = "gene_name", all.x = TRUE)
write.xlsx(aging_kinase_degs_occurrence_3, file = file.path("E:/project_data_ssd/hpk-1_dataset",
                                                            "Elegans kinome aging DE occurrence table (v3 with lifespan info).xlsx"), asTable = TRUE)

## !! TODO !!
#
# Kinases- how about genes that are expressed in the same cells as hpk-1, and not necessarily differentially expressed?
#  Anything for known PPIs with hipk1 in STRING or similar?
# Andy kept saying "co-expression" and there aren't that many kinases that show similar trends in the same cell with aging
#  as hpk-1- but we also are stating that hpk-1 is pretty unusual for being regulated at the level of expression-
#  so we'd expect that other kinases could participate in a signal cascade and just be expressed rather than differentially expressed
#
# So this would be across ages
# Should also just check if I happened to have pulled any kinases out of my correlation analysis previously
#
elegans_kinases



#
# How about phosphatases?
#
#
phosphatase_table <- read.xlsx(file.path("E:/project_data_ssd/hpk-1_dataset/phosphatases",
                                         "Uniprot-C elegans protein phosphatases 10Nov2022.xlsx"))
# have to map gene IDs
# try to use the molecular names
phosphatase_table <- phosphatase_table[ !is.na(phosphatase_table$`Gene.Names.(ORF)`) ,]
phosphatase_table$`Gene.Names.(ORF)` <- gsub("CELE_", "", unlist(lapply(strsplit(phosphatase_table$`Gene.Names.(ORF)`, " "), function(x){x[1]})))
phosphatase_table$`Gene.Names.(ORF)`[phosphatase_table$`Gene.Names.(ORF)` == "C09D8.1/C09D8.2"] <- "C09D8.1" # C09D8.2 is an ncRNA
# some of these are duplicated- representing different isoforms!
# I feel like some things might be missing here, but it's most of them.
length(unique(phosphatase_table$`Gene.Names.(ORF)`))
sort(table(phosphatase_table$`Gene.Names.(ORF)`))
# for ones that have multiple entries- 
#  if there's an entry that's been reviewed, keep it
#  otherwise, keep the longest isoform
# the whole sapply just picks the indices to keep
phosphatase_table_uq <- phosphatase_table[ sapply(1:length(unique(phosphatase_table$`Gene.Names.(ORF)`)), function(i){
  currGene = unique(phosphatase_table$`Gene.Names.(ORF)`)[i]
  if(sum(phosphatase_table$`Gene.Names.(ORF)` == currGene) == 1){
    return(which(phosphatase_table$`Gene.Names.(ORF)` == currGene))
  }else if(any(phosphatase_table$Reviewed[phosphatase_table$`Gene.Names.(ORF)` == currGene] == "reviewed")){
    return(which((phosphatase_table$Reviewed == "reviewed") & (phosphatase_table$`Gene.Names.(ORF)` == currGene))[1])
  }else{
    return(
      which(phosphatase_table$`Gene.Names.(ORF)` == currGene)[
        order(phosphatase_table$Length[phosphatase_table$`Gene.Names.(ORF)` == currGene], decreasing = TRUE)][1]
    )
  }
}) , ]

phosphatase_table_uq <- merge(phosphatase_table_uq, wb_gene_info, by.x = "Gene.Names.(ORF)", by.y = "molecular_name", all.x = TRUE)
# select columns to keep
phosphatase_table_uq <- phosphatase_table_uq[ , c("WBGeneID", "gene_name", "gene_biotype", "Protein.names", "Entry", "Reviewed") ]

sum(phosphatase_table_uq$gene_name %in% aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno])

#
# NOTE- SINCE THE CONTEXT HERE IS REGARDING hpk-1 WE WILL ALSO INCLUDE IT
# SO THIS IS NOT STRICTLY JUST PHOSPHATASES IN THIS PARTICULAR OBJECT
aging_phosphatase_degs <- aging_atlas_degs_long[ aging_atlas_degs_long$gene_name %in% c(phosphatase_table_uq$gene_name, "hpk-1") , ]

# how many phosphatases have any DE with age in any cluster?
length(unique(aging_phosphatase_degs$gene_name)) # only 54 of 177 phosphatases

any(aging_phosphatase_degs$padj > 0.05)
sum((aging_phosphatase_degs$logFoldChange > 0 ) & (aging_phosphatase_degs$logFoldChange < 0.5 ))
# transform to a gene x cluster matrix of foldchanges, since these are all significant results
aging_phosphatase_degs_mat <- mat.or.vec(nr = length(unique(aging_phosphatase_degs$gene_name)),
                                         nc = length(unique(aging_phosphatase_degs$cluster_shortID)))
rownames(aging_phosphatase_degs_mat) <- unique(aging_phosphatase_degs$gene_name)
colnames(aging_phosphatase_degs_mat) <- unique(aging_phosphatase_degs$cluster_shortID)

for(i in rownames(aging_phosphatase_degs_mat)){
  for(j in colnames(aging_phosphatase_degs_mat)){
    if(any(with(aging_phosphatase_degs, (cluster_shortID == j) & (gene_name == i) ))){
      aging_phosphatase_degs_mat[i,j] <- with(aging_phosphatase_degs, logFoldChange[(cluster_shortID == j) & (gene_name == i)] )
    }else{
      aging_phosphatase_degs_mat[i,j] <- 0 # I'd feel better with NA but for now this is ok
    }
  }
}
rm(i, j)


# ORDER ROWS BY AVERAGE VALUE (LFC)
# and don't cluster the heatmap
aging_phosphatase_degs_mat <- aging_phosphatase_degs_mat[ order(rowMeans(aging_phosphatase_degs_mat), decreasing = TRUE) , ]

tmp_aging_genes <- setNames(aging_genes_merged$any_aging_pheno, aging_genes_merged$gene_name)
tmp_aging_genes <- tmp_aging_genes[rownames(aging_phosphatase_degs_mat)]
identical(names(tmp_aging_genes), rownames(aging_phosphatase_degs_mat)) # should be true

svg(file = "Elegans phosphatome (with hpk1) and aging differential expression in cell clusters v1.svg", width = 11, height = 16)
heatmap.3(aging_phosphatase_degs_mat[,tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(aging_phosphatase_degs_mat)]],
          RowSideColors =t(as.matrix(ifelse(tmp_aging_genes, "gold4", "black"))),
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(aging_phosphatase_degs_mat)]),
          dendrogram = "none", Colv = FALSE, Rowv = FALSE,
          col = colorRampPalette(c("blue4", "white", "red4")), #scale = "row",
          cexCol = 0.75, cexRow = 1, margins = c(10,7), keysize = 0.8,
          labCol = NA,
          rowsep = seq(1, nrow(aging_phosphatase_degs_mat)),
          colsep = seq(1, ncol(aging_phosphatase_degs_mat)),
          sepcolor = "grey90",sepwidth = c(0.001, 0.001),
          KeyValueName = "Aging atlas old vs yound LFC")

graphics.off()

rm(tmp_aging_genes)

#
# Now restrict to only clusters with hpk-1 DE
#
aging_phosphatase_degs_mat.2 <- aging_phosphatase_degs_mat[ , aging_phosphatase_degs_mat["hpk-1",] != 0 ]
# remove rows without any DE results in these clusters
aging_phosphatase_degs_mat.2 <- aging_phosphatase_degs_mat.2[ rowSums(aging_phosphatase_degs_mat.2 > 0) > 0 , ]

# sort again
aging_phosphatase_degs_mat.2 <- aging_phosphatase_degs_mat.2[ order(rowMeans(aging_phosphatase_degs_mat.2), decreasing = TRUE) , ]

tmp_aging_genes <- setNames(aging_genes_merged$any_aging_pheno, aging_genes_merged$gene_name)
tmp_aging_genes <- tmp_aging_genes[rownames(aging_phosphatase_degs_mat.2)]
identical(names(tmp_aging_genes), rownames(aging_phosphatase_degs_mat.2)) # should be true



svg(file = "Elegans phosphatome (with hpk1) and aging differential expression ONLY IN HPK-1 DE CLUSTERS v1.svg", width = 14, height = 12)
heatmap.3(aging_phosphatase_degs_mat.2[,tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(aging_phosphatase_degs_mat.2)]],
          RowSideColors =t(as.matrix(ifelse(tmp_aging_genes, "gold4", "black"))),
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(aging_phosphatase_degs_mat.2)]),
          dendrogram = "none", Colv = FALSE, Rowv = FALSE,
          col = colorRampPalette(c("blue4", "white", "red4")), #scale = "row",
          cexCol = 1, cexRow = 1.4, margins = c(20,10), keysize = 0.8,
          labCol =  with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(aging_phosphatase_degs_mat.2[ , tissue_type_vec$shortID[tissue_type_vec$shortID %in% 
                                                                                     colnames(aging_phosphatase_degs_mat.2)]]), shortID)]),
          rowsep = seq(1, nrow(aging_phosphatase_degs_mat.2)),
          colsep = seq(1, ncol(aging_phosphatase_degs_mat.2)),
          sepcolor = "grey90",sepwidth = c(0.001, 0.001),
          KeyValueName = "Aging atlas old vs yound LFC")

graphics.off()

rm(tmp_aging_genes)




# What's next?
#  - TF enrichment, based on neuron-associated genes DE in hpk-1
#  - co-expression with hpk-1 in Cengen? (TFs, as well as genes that are DE and enriched for a TF)
#     - This is done for aging atlas
#  - heatmap across cells, labeled by classes/categories, for autophagy and other genes of interest.

#
# look at all (ALL) TFs across aging in nervous system (and elsewhere) ONLY in cell clusters with significant change in hpk-1 (out of 48)
#   hpk-1 expression trends in non-neuronal cells with age?
#   ** predicted target gene expression with age?
#

hpk1_age_de_clusters <- aging_atlas_degs_long[aging_atlas_degs_long$gene_name == "hpk-1" , ]
hpk1_age_de_clusters <- merge(hpk1_age_de_clusters, aging_atlas_deg_clusters[,c(2,5)],
                              by.x = "cluster_shortID", by.y = "shortID", all.x = TRUE)
sort(table(hpk1_age_de_clusters$`Tissue.type.(higher.degree)`))

#  Question numero ichi: what TFs are differentially expressed in the same clusters as hpk-1 across normal aging?
#    Just going to go through wTF in atlas_degs_long
tfs_de_in_hpk1_clusters <- aging_atlas_degs_long[ with(aging_atlas_degs_long, 
                                                       (gene_name %in% c(wTF$Public.name, "hpk-1")) & # ADD hpk-1
                                                         (cluster_shortID %in% hpk1_age_de_clusters$cluster_shortID)) ,    ]
sum((tfs_de_in_hpk1_clusters$logFoldChange > 0 ) & (tfs_de_in_hpk1_clusters$logFoldChange < 0.5 )) # should have already been filtered...

length(unique(tfs_de_in_hpk1_clusters$cluster_shortID))
length(unique(tfs_de_in_hpk1_clusters$gene_name))

#  Turn the "long" data object into a matrix
# ah yes, I was indeed using an unfiltered version of the DE results before

tfs_de_in_hpk1_clusters_fc_mat <- mat.or.vec(nr = length(unique(tfs_de_in_hpk1_clusters$gene_name)),
                                             nc = length(unique(tfs_de_in_hpk1_clusters$cluster_shortID)))
rownames(tfs_de_in_hpk1_clusters_fc_mat) <- unique(tfs_de_in_hpk1_clusters$gene_name)
colnames(tfs_de_in_hpk1_clusters_fc_mat) <- unique(tfs_de_in_hpk1_clusters$cluster_shortID)

for(i in 1:nrow(tfs_de_in_hpk1_clusters_fc_mat)){
  for(j in 1:ncol(tfs_de_in_hpk1_clusters_fc_mat)){
    rv = tfs_de_in_hpk1_clusters$logFoldChange[
      (tfs_de_in_hpk1_clusters$gene_name == rownames(tfs_de_in_hpk1_clusters_fc_mat)[i]) &
        (tfs_de_in_hpk1_clusters$cluster_shortID == colnames(tfs_de_in_hpk1_clusters_fc_mat)[j])]
    if(length(rv) == 0){
      tfs_de_in_hpk1_clusters_fc_mat[i,j] <- NA
    }else{
      tfs_de_in_hpk1_clusters_fc_mat[i,j] <- rv 
    }
  }
}
rm(i,j, rv)
hist(tfs_de_in_hpk1_clusters_fc_mat)
dim(tfs_de_in_hpk1_clusters_fc_mat )

#
# Summarize results for each resulting TF, number up, number down, and split by neuron or non-neuron
#  ---------> also include aging/lifespan associations
#
tfs_de_in_hpk1_clust_summary <- Reduce(f = rbind, lapply(rownames(tfs_de_in_hpk1_clusters_fc_mat), function(geneName){
  
  up_clust = na.omit(colnames(tfs_de_in_hpk1_clusters_fc_mat)[tfs_de_in_hpk1_clusters_fc_mat[geneName, ] > 0])
  down_clust = na.omit(colnames(tfs_de_in_hpk1_clusters_fc_mat)[tfs_de_in_hpk1_clusters_fc_mat[geneName, ] < 0])
  
  up_clust_types = sort(table(aging_atlas_deg_clusters$`Tissue.type.(higher.degree)`[aging_atlas_deg_clusters$shortID %in% up_clust]),
                        decreasing = TRUE)
  up_clust_types <- paste(apply(cbind(names(up_clust_types), up_clust_types), MARGIN = 1, function(x){
    paste(x, collapse = ":")
    }),
    collapse = ", ")
  
  down_clust_types = sort(table(aging_atlas_deg_clusters$`Tissue.type.(higher.degree)`[aging_atlas_deg_clusters$shortID %in% down_clust]),
                        decreasing = TRUE)
  down_clust_types <- paste(apply(cbind(names(down_clust_types), down_clust_types), MARGIN = 1, function(x){
    paste(x, collapse = ":")
  }),
  collapse = ", ")
  
  return(data.frame(
    "gene_name" = geneName,
    "WBGeneID" = ifelse(length(wb_gene_info$WBGeneID[wb_gene_info$gene_name == geneName]) == 0, NA,
                        wb_gene_info$WBGeneID[wb_gene_info$gene_name == geneName]),
    "total_clust" = length(up_clust) + length(down_clust),
    "num_clust_up" = length(up_clust),
    "num_clust_down" = length(down_clust),
    "up_clust_types" = up_clust_types,
    "down_clust_types" = down_clust_types,
    "implicated_in_aging" = geneName %in% aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno]))
}))

write.xlsx(x = tfs_de_in_hpk1_clust_summary, file = "TFs DE in hpk-1 DE clusters summary table v2.1.xlsx", asTable = TRUE)


#
# Heatmap for TFs DE in old vs young in same clusters as hpk-1
#
heatmapData <- as.matrix(tfs_de_in_hpk1_clusters_fc_mat)
mode(heatmapData) <- "numeric"
heatmapData[is.na(heatmapData)] <- 0

heatmapData <- heatmapData[ , tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData)]]

# order by number of non-zero rows
# heatmapData <- heatmapData[ order(rowSums(heatmapData == 0), decreasing = FALSE) , ]
# order by overall mean of row FC
heatmapData <- heatmapData[ order(rowMeans(heatmapData), decreasing = TRUE) , ]

tmp_aging_genes <- setNames(aging_genes_merged$any_aging_pheno, aging_genes_merged$gene_name)
setdiff(rownames(heatmapData), names(tmp_aging_genes))
# A BUNCH OF GENES HAVE CHANGED NAMES so we'll fix that manually for now
# atf-5 changed to atf-4
# C34D10.2 is unk-1
# T18D3.7 is tsct-1
# C52E12.1 is znf-598
# C01F6.9 is znf-706
# F47E1.3 is prdm-14
# F43G9.12 is mog-7
rownames(heatmapData)[rownames(heatmapData) == "atf-5"] <- "atf-4"
rownames(heatmapData)[rownames(heatmapData) == "C34D10.2"] <- "unk-1"
rownames(heatmapData)[rownames(heatmapData) == "T18D3.7"] <- "tsct-1"
rownames(heatmapData)[rownames(heatmapData) == "C52E12.1"] <- "znf-598"
rownames(heatmapData)[rownames(heatmapData) == "C01F6.9"] <- "znf-706"
rownames(heatmapData)[rownames(heatmapData) == "F47E1.3"] <- "prdm-14"
rownames(heatmapData)[rownames(heatmapData) == "F43G9.12"] <- "mog-7"
setdiff(rownames(heatmapData), names(tmp_aging_genes))
tmp_aging_genes <- tmp_aging_genes[rownames(heatmapData)]
identical(names(tmp_aging_genes), rownames(heatmapData)) # SHOULD BE TRUE

# CHECK ORDER OF COLUMNS TO MAKE SURE IT MATCHES COLOR ANNOTATION
identical(colnames(heatmapData),
          tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData) ])


svg(file = "Aging atlas DE across age for TFs in same clusters as hpk-1 v2.2.svg", width = 12, height = 32)
heatmap.2(heatmapData, trace = 'none', density.info = "none",
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(heatmapData)]),
          RowSideColors = t(as.matrix(ifelse(tmp_aging_genes, "gold4", "black"))),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(32), #scale = "row",
          cexCol = 1.1, cexRow = 0.85, margins = c(1,8), keysize = 0.8,
          Colv = FALSE, Rowv = FALSE, dendrogram = "none",
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          labCol = NA,
          colsep = seq(1, ncol(heatmapData)),
          rowsep = seq(1, nrow(heatmapData)),
          sepcol = "grey80", sepwidth = c(0.01, 0.01), offsetRow = 0,
          key.title = "old vs young", key.xlab = "log2 FC")
graphics.off()


#
# OK now we'll subset this into neurons and non-neurons
#  and filter out ones that don't have any change in each case
# 

heatmapData_neurons <- heatmapData[ , colnames(heatmapData) %in% tissue_type_vec$shortID[tissue_type_vec$tissue_type == "neuron"] ]
heatmapData_notneurons <- heatmapData[ , colnames(heatmapData) %in% tissue_type_vec$shortID[tissue_type_vec$tissue_type != "neuron"] ]

# remove rows that are empty for each
heatmapData_neurons <- heatmapData_neurons[ rowSums(heatmapData_neurons == 0) != ncol(heatmapData_neurons) , ]
heatmapData_notneurons <- heatmapData_notneurons[ rowSums(heatmapData_notneurons == 0) != ncol(heatmapData_notneurons) , ]

dim(heatmapData)
dim(heatmapData_neurons)
dim(heatmapData_notneurons)

# order by overall mean of row FC
heatmapData_neurons <- heatmapData_neurons[ order(rowMeans(heatmapData_neurons), decreasing = TRUE) , ]
heatmapData_notneurons <- heatmapData_notneurons[ order(rowMeans(heatmapData_notneurons), decreasing = TRUE) , ]

# Additional filtered object- keep rows that are non-zero in three or more clusters in either neurons OR non-neurons
#  This strategy means that a TF has three total DE clusters split between neurons and non-neurons won't show up
heatmapData_neurons_filt2 <- heatmapData_neurons[ rowSums(heatmapData_neurons == 0) < (ncol(heatmapData_neurons)-2) , ]
heatmapData_notneurons_filt2 <- heatmapData_notneurons[ rowSums(heatmapData_notneurons == 0) < (ncol(heatmapData_notneurons)-2) , ]

dim(heatmapData_neurons)
dim(heatmapData_neurons_filt2)

dim(heatmapData_notneurons)
dim(heatmapData_notneurons_filt2)


# now we can color neurons by neuron type
atlas_cengen_neuron_map.sub <- merge(colnames(heatmapData_neurons), atlas_cengen_neuron_map, by.x = 1, by.y = "atlas_cluster", all.x = TRUE)
colnames(atlas_cengen_neuron_map.sub)[1] <- "atlas_cluster"
atlas_cengen_neuron_map.sub <- merge(atlas_cengen_neuron_map.sub, setNames(atlas_unique_neuron_signal_class[,"consensus"],
                                                                           rownames(atlas_unique_neuron_signal_class)),
                                     by.x = "atlas_cluster", by.y = 0, all.x = TRUE)
colnames(atlas_cengen_neuron_map.sub)[ncol(atlas_cengen_neuron_map.sub)] <- "atlas_unique_signal_cons"
atlas_cengen_neuron_map.sub$NeuronType[is.na(atlas_cengen_neuron_map.sub$NeuronType)] <- aging_atlas_clusters_all$Clusters.explanatory.name[match( 
  atlas_cengen_neuron_map.sub$atlas_cluster[is.na(atlas_cengen_neuron_map.sub$NeuronType)],
  aging_atlas_clusters_all$shortID)]

atlas_cengen_neuron_map.sub$NeuronType[grepl("*inter neuron*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*inter neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*interneurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motor neuron*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motor neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motorneurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"

atlas_cengen_neuron_map.sub <- atlas_cengen_neuron_map.sub[ match(colnames(heatmapData_neurons),
                                                                  atlas_cengen_neuron_map.sub$atlas_cluster) , ] 
identical(colnames(heatmapData_neurons), atlas_cengen_neuron_map.sub$atlas_cluster) # should be true

tmp_type = atlas_cengen_neuron_map.sub$NeuronType
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]

swatch2(c("Interneuron" = as.vector(colorBlindness::paletteMartin[2]),
          "Sensory" = as.vector(colorBlindness::paletteMartin[4]),
          "Endorgan" = as.vector(colorBlindness::paletteMartin[6]),
          "Motorneuron" = as.vector(colorBlindness::paletteMartin[8])))

# UPDATE- for cells where there was no call in cengen expression, look for a call in the aging atlas expression


tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub$cengen_signal_class) | 
                             grepl("glutamatergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                            grepl("serotonergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                         grepl("GABAergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                            grepl("dopaminergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                           grepl("cholinergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

# USE CELL NAMES AS COLUMN LABELS (Clusters.explanatoryname in aging_atlas_clusters_all)

# so now we also want to sort the columns by neuron type.

# NEURONS ONLY
svg(file = "Aging atlas DE across age for TFs in same NEURON clusters as hpk-1 v2.1.svg", width = 12, height = 28)
heatmap.3(heatmapData_neurons[ , order(neuron_colors[,1])],
          ColSideColors = as.matrix(neuron_colors[order(neuron_colors[,1]) , ]),
          RowSideColors = t(as.matrix(ifelse(tmp_aging_genes[rownames(heatmapData_neurons)], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_neurons[ , order(neuron_colors[,1])]), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 1.1, cexRow = 0.7, margins = c(8,10), keysize = 0.8,
          colsep = seq(1, ncol(heatmapData_neurons)),
          rowsep = seq(1, nrow(heatmapData_neurons)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", 
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          KeyValueName = "LogFoldChange")
graphics.off()

# NEURONS ONLY - ADDITIONAL FILTER
svg(file = "Aging atlas DE across age for TFs in same NEURON clusters as hpk-1 (ADDITIONAL GENE FILT) v2.svg", width = 12, height = 15)
heatmap.2(heatmapData_neurons_filt2[ , order(neuron_colors[,1])],
          density.info = "none", trace = "none",
          #ColSideColors = as.matrix(neuron_colors[order(neuron_colors[,1]) , ]),
          RowSideColors = t(as.matrix(ifelse(tmp_aging_genes[rownames(heatmapData_neurons_filt2)], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_neurons_filt2[ , order(neuron_colors[,1])]), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 1.1, cexRow = 1, margins = c(8,10), keysize = 0.8,
          colsep = seq(1, ncol(heatmapData_neurons_filt2)),
          rowsep = seq(1, nrow(heatmapData_neurons_filt2)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", 
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          key.xlab = "LogFoldChange")
graphics.off()


# NON-NEURONS ONLY
# CHECK ORDER OF COLUMNS TO MAKE SURE IT MATCHES COLOR ANNOTATION
identical(colnames(heatmapData_notneurons),
          tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData_notneurons) ]) 
svg(file = "Aging atlas DE across age for TFs in same NON-NEURON clusters as hpk-1 v2.1.svg", width = 12, height = 18)
heatmap.3(heatmapData_notneurons,
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(heatmapData_notneurons)]),
          RowSideColors =  t(as.matrix(ifelse(tmp_aging_genes[rownames(heatmapData_notneurons)], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_notneurons), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 1.1, cexRow = 0.7, margins = c(8,10), keysize = 0.8,
          Colv = FALSE, Rowv = FALSE, dendrogram = "none",
          colsep = seq(1, ncol(heatmapData_notneurons)),
          rowsep = seq(1, nrow(heatmapData_notneurons)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          KeyValueName = "LogFoldChange")
graphics.off()

# NON-NEURONS ONLY - additional filter
identical(colnames(heatmapData_notneurons_filt2),
          tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData_notneurons_filt2) ]) 
svg(file = "Aging atlas DE across age for TFs in same NON-NEURON clusters as hpk-1 (ADDITIONAL GENE FILT) v1.1.svg", width = 12, height = 12)
heatmap.3(heatmapData_notneurons_filt2,
          ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(heatmapData_notneurons_filt2)]),
          RowSideColors =  t(as.matrix(ifelse(tmp_aging_genes[rownames(heatmapData_notneurons_filt2)], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_notneurons_filt2), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 1.1, cexRow = 1.1, margins = c(8,10), keysize = 0.8,
          Colv = FALSE, Rowv = FALSE, dendrogram = "none",
          colsep = seq(1, ncol(heatmapData_notneurons_filt2)),
          rowsep = seq(1, nrow(heatmapData_notneurons_filt2)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          KeyValueName = "LogFoldChange")
graphics.off()
# for fixing grid
svg(file = "Aging atlas DE across age for TFs in same NON-NEURON clusters as hpk-1 (ADDITIONAL GENE FILT)(USE FOR GRID).svg", width = 12, height = 12)
heatmap.2(heatmapData_notneurons_filt2,
          trace = "none", density.info = "none", 
          #ColSideColors = as.matrix(tissue_type_vec$type_color[tissue_type_vec$shortID %in% colnames(heatmapData_notneurons_filt2)]),
          RowSideColors =  t(as.matrix(ifelse(tmp_aging_genes[rownames(heatmapData_notneurons_filt2)], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_notneurons_filt2), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 1.1, cexRow = 1.1, margins = c(8,10), keysize = 0.8,
          Colv = FALSE, Rowv = FALSE, dendrogram = "none",
          colsep = seq(1, ncol(heatmapData_notneurons_filt2)),
          rowsep = seq(1, nrow(heatmapData_notneurons_filt2)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          key.xlab = "LogFoldChange")
graphics.off()


# Venn diagram for TFs in neurons vs TFs in non-neurons?

venn.diagram(list("TFs DE with hpk-1 in neurons" = setdiff(rownames(heatmapData_neurons), "hpk-1"),
                  "TFs DE with hpk-1 in other tissues" = setdiff(rownames(heatmapData_notneurons), "hpk-1")),
             filename = "Venn for TFs DE in hpk-1 in neurons vs non-neurons.svg",
             imagetype = "svg", width = 10, height = 8, units = "in")

venn.diagram(list("TFs DE with hpk-1 in neurons" = setdiff(rownames(heatmapData_neurons_filt2), "hpk-1"),
                  "TFs DE with hpk-1 in other tissues" = setdiff(rownames(heatmapData_notneurons_filt2), "hpk-1")),
             filename = "Venn for TFs DE in hpk-1 in neurons vs non-neurons after extra filt.svg",
             imagetype = "svg", width = 10, height = 8, units = "in")

# bonus- intersection with aging-associated genes
venn.diagram(list("TFs DE with hpk-1 in neurons" = setdiff(rownames(heatmapData_neurons), "hpk-1"),
                  "TFs DE with hpk-1 in other tissues" = setdiff(rownames(heatmapData_notneurons), "hpk-1"),
                  "Aging-associated TFs" = wTF$Public.name[wTF$WBGeneID %in% aging_genes_merged$WBGeneID[aging_genes_merged$any_aging_pheno]]),
             filename = "Venn for TFs DE in hpk-1 in neurons vs non-neurons vs aging TFs.svg",
             imagetype = "svg", width = 10, height = 8, units = "in")

venn.diagram(list("TFs DE with hpk-1 in neurons" = setdiff(rownames(heatmapData_neurons_filt2), "hpk-1"),
                  "TFs DE with hpk-1 in other tissues" = setdiff(rownames(heatmapData_notneurons_filt2), "hpk-1"),
                  "Aging-associated TFs" = wTF$Public.name[wTF$WBGeneID %in% aging_genes_merged$WBGeneID[aging_genes_merged$any_aging_pheno]]),
             filename = "Venn for TFs DE in hpk-1 in neurons vs non-neurons vs aging TFs after extra filt.svg",
             imagetype = "svg", width = 10, height = 8, units = "in")



rm(heatmapData, heatmapData_neurons, heatmapData_neurons_filt2, heatmapData_notneurons, heatmapData_notneurons_filt2,
    tmp_cholinergic, atlas_cengen_neuron_map.sub, neuron_colors,
   tmp_dopaminergic, tmp_GABAergic, tmp_glutamatergic, tmp_serotonergic, tmp_type)


# 
# Even further filtering for neurons- just those DE in 9 or more clusters
#
heatmapData_neurons_filt3 <- heatmapData_neurons_filt2[ rowSums(heatmapData_neurons_filt2 != 0) >= 9 , ]

# NEURONS ONLY - ADDITIONAL FILTER
svg(file = "Aging atlas DE across age for TFs in same NEURON clusters as hpk-1 (only top results) v1.svg", width = 12, height = 8)
heatmap.3(heatmapData_neurons_filt3[ , order(neuron_colors[,1])],
          ColSideColors = as.matrix(neuron_colors[order(neuron_colors[,1]) , ]),
          RowSideColors = t(as.matrix(ifelse(tmp_aging_genes[rownames(heatmapData_neurons_filt3)], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_neurons_filt3[ , order(neuron_colors[,1])]), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 1.1, cexRow = 1, margins = c(8,10), keysize = 0.8,
          colsep = seq(1, ncol(heatmapData_neurons_filt3)),
          rowsep = seq(1, nrow(heatmapData_neurons_filt3)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", 
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5, TFs DE in 9+ clusts",
          KeyValueName = "LogFoldChange")
graphics.off()

#
# Treemap for the functions associated with this smaller subset of TFs (heatmapData_neurons_filt3)
#  I loaded the library earlier in the file/workspace, so I won't load it again here
# I was not assuming that we'd have very much to work from for enrichment, but 
#  we might be able to just show all associated terms and pathways
# up and down separately???
library(GO.db)


#
#  !!! These objects are actually loaded later in the file!
#
# sets_kegg_all_long
# sets_reactome_long
# elegansGO

# add ontology branch to GO object
# really should do this in bulk rather than one at a time
#  unfortunately some of these terms need an update of the GO DB object.
# This takes a REALLY long time time run doing lookups individually, 
#  I'm just going to do it for the tiny set of GO terms that actually overlap!

# setdiff(unique(elegansGO$GO), keys(GO.db))
# 
# elegansGO$onto_branch <- sapply(elegansGO$GO, function(x){
#   tryCatch( # return NA if not found
#     {as.vector(select(GO.db, keys = x, columns = "ONTOLOGY", keytype = "GOID")[,"ONTOLOGY"])},
#            error = function(e){NA}
#   )
# })

# DO NOT INCLUDE hpk-1 here
neuron_aging_tf_sub3 <- merge(data.frame(gene_name = rownames(heatmapData_neurons_filt3),
                                         aging_de_mean = rowMeans(heatmapData_neurons_filt3, na.rm = TRUE)), 
      wb_gene_info[,c("WBGeneID", "gene_name")], all.x = TRUE)
# remove hpk-1 since it is not a TF!
neuron_aging_tf_sub3 <- neuron_aging_tf_sub3[ !(neuron_aging_tf_sub3$gene_name == "hpk-1") , ]
#--------------------------------------------
# !!! WITH UPREGULATED ON AVERAGE FIRST !!! 
#--------------------------------------------
tmp_isect_kegg <- sets_kegg_all_long[ sets_kegg_all_long$gene_id %in% neuron_aging_tf_sub3$WBGeneID[neuron_aging_tf_sub3$aging_de_mean > 0] , ]
tmp_isect_reac <- sets_reactome_long[ sets_reactome_long$gene_id %in% neuron_aging_tf_sub3$WBGeneID[neuron_aging_tf_sub3$aging_de_mean > 0] , ]
tmp_isect_go <- elegansGO[(elegansGO$GeneID %in% neuron_aging_tf_sub3$WBGeneID[neuron_aging_tf_sub3$aging_de_mean > 0]) , ]
tmp_isect_go$onto_branch <- sapply(tmp_isect_go$GO, function(x){
  tryCatch( # return NA if not found
    {as.vector(select(GO.db, keys = x, columns = "ONTOLOGY", keytype = "GOID")[,"ONTOLOGY"])},
    error = function(e){NA}
  )
})
sum(is.na(tmp_isect_go$onto_branch)) # only 1, remove it
tmp_isect_go <- tmp_isect_go[ !is.na(tmp_isect_go$onto_branch) , ]
sum(is.na(tmp_isect_go$onto_branch))
# add the GO definitions
tmp_isect_go$term_def <- sapply(tmp_isect_go$GO, function(x){
  tryCatch( # return NA if not found
    {as.vector(select(GO.db, keys = x, columns = "TERM", keytype = "GOID")[,"TERM"])},
    error = function(e){NA}
  )
})

# GO needs more work to get unique associations than the others
tmp_isect_go <- tmp_isect_go[ , -4 ]
tmp_isect_go <- unique(tmp_isect_go)
tmp_isect_go.bp <- tmp_isect_go[ tmp_isect_go$onto_branch == "BP" , ]
tmp_isect_go.mf <- tmp_isect_go[ tmp_isect_go$onto_branch == "MF" , ]
tmp_isect_go.cc <- tmp_isect_go[ tmp_isect_go$onto_branch == "CC" , ]

# Andy suggested removing CC

neuron_aging_tf_sub3_up_func_tabl <- rbind(
  data.frame(DataSource = "GO:BP", table(tmp_isect_go.bp$term_def)),
  data.frame(DataSource = "GO:MF", table(tmp_isect_go.mf$term_def)),
  #data.frame(DataSource = "GO:CC", table(tmp_isect_go.cc$term_def)),
  data.frame(DataSource = "KEGG", table(tmp_isect_kegg$pathway)),
  data.frame(DataSource = "Reactome", table(tmp_isect_reac$pathway))
)
rm(tmp_isect_go, tmp_isect_go.bp, tmp_isect_go.cc, tmp_isect_go.mf, tmp_isect_kegg, tmp_isect_reac)

#--------------------------------------------
# !!! NOW WITH DOWNREGULATED ON AVERAGE !!! 
#--------------------------------------------
tmp_isect_kegg <- sets_kegg_all_long[ sets_kegg_all_long$gene_id %in% neuron_aging_tf_sub3$WBGeneID[neuron_aging_tf_sub3$aging_de_mean < 0] , ]
tmp_isect_reac <- sets_reactome_long[ sets_reactome_long$gene_id %in% neuron_aging_tf_sub3$WBGeneID[neuron_aging_tf_sub3$aging_de_mean < 0] , ]
tmp_isect_go <- elegansGO[(elegansGO$GeneID %in% neuron_aging_tf_sub3$WBGeneID[neuron_aging_tf_sub3$aging_de_mean < 0]) , ]
tmp_isect_go$onto_branch <- sapply(tmp_isect_go$GO, function(x){
  tryCatch( # return NA if not found
    {as.vector(select(GO.db, keys = x, columns = "ONTOLOGY", keytype = "GOID")[,"ONTOLOGY"])},
    error = function(e){NA}
  )
})
sum(is.na(tmp_isect_go$onto_branch)) # only 1, remove it
tmp_isect_go <- tmp_isect_go[ !is.na(tmp_isect_go$onto_branch) , ]
sum(is.na(tmp_isect_go$onto_branch))
# add the GO definitions
tmp_isect_go$term_def <- sapply(tmp_isect_go$GO, function(x){
  tryCatch( # return NA if not found
    {as.vector(select(GO.db, keys = x, columns = "TERM", keytype = "GOID")[,"TERM"])},
    error = function(e){NA}
  )
})

# GO needs more work to get unique associations than the others
tmp_isect_go <- tmp_isect_go[ , -4 ]
tmp_isect_go <- unique(tmp_isect_go)
tmp_isect_go.bp <- tmp_isect_go[ tmp_isect_go$onto_branch == "BP" , ]
tmp_isect_go.mf <- tmp_isect_go[ tmp_isect_go$onto_branch == "MF" , ]
tmp_isect_go.cc <- tmp_isect_go[ tmp_isect_go$onto_branch == "CC" , ]

# Andy suggested removing CC
neuron_aging_tf_sub3_down_func_tabl <- rbind(
  data.frame(DataSource = "GO:BP", table(tmp_isect_go.bp$term_def)),
  data.frame(DataSource = "GO:MF", table(tmp_isect_go.mf$term_def)),
  #data.frame(DataSource = "GO:CC", table(tmp_isect_go.cc$term_def)),
  data.frame(DataSource = "KEGG", table(tmp_isect_kegg$pathway)),
  data.frame(DataSource = "Reactome", table(tmp_isect_reac$pathway))
)

rm(tmp_isect_go, tmp_isect_go.bp, tmp_isect_go.cc, tmp_isect_go.mf, tmp_isect_kegg, tmp_isect_reac)

colnames(neuron_aging_tf_sub3_up_func_tabl) <- c("DataSource", "term", "count")
colnames(neuron_aging_tf_sub3_down_func_tabl) <- c("DataSource", "term", "count")
#
# Ok so it's probably not going to be great, but we'll use the terms with more than one gene associated
#  and use the # associated to determine the box size.
#


aging_tf_func_tabl <- list(up = neuron_aging_tf_sub3_up_func_tabl, 
                           down = neuron_aging_tf_sub3_down_func_tabl)
# also filter
aging_tf_func_tabl <- lapply(aging_tf_func_tabl, function(x){
  x$DataSource <- factor(x$DataSource, levels = c("GO:BP","GO:MF","KEGG", "Reactome"))
  x <- x[ x$count > 1 , ]
  return(x)
})

pdf("Treemap- functions associated with top TFs DE with aging in hpk-1 DE clusters.pdf",
    width = 12, height = 8)
for(currRes in names(aging_tf_func_tabl)){
  treemap(dtf = wrapify(aging_tf_func_tabl[[currRes]], "count", "term"),
          index = c("DataSource", "term"),
          vSize = "count",
          vColor = "DataSource",
          type = "categorical", 
          force.print.labels = FALSE,
          aspRatio = 1.5,
          fontsize.labels = c(0, 16),
          position.legend = "right", 
          fontsize.legend = 6,
          drop.unused.levels = FALSE,
          title = paste0("Functional terms associated with TFs that are on average \n", currRes))
}
graphics.off()
rm(currRes)

#
#  FOR OUR RNA-SEQ DATA-
# Since I hadn't plotted anything for the GOseq enrichment results
# based on all hpk-1 null DE genes- doing that here.
# Not including CC here
#

hpk1_de_all_enriched <- list(
  up = rbind(data.frame(DataSource = "GO:BP", hpk1_up_go[hpk1_up_go$ontology == "BP" ,c("term", "over_represented_pvalue")]),
             data.frame(DataSource = "GO:MF", hpk1_up_go[hpk1_up_go$ontology == "MF" ,c("term", "over_represented_pvalue")]),
             with(hpk1_up_kegg[ , c("category", "over_represented_pvalue") ],
                  data.frame(DataSource = "KEGG", term = category, over_represented_pvalue )),
             with(hpk1_up_reac[ , c("category", "over_represented_pvalue") ],
                  data.frame(DataSource = "Reactome", term = category, over_represented_pvalue ))),
  
  down = rbind(data.frame(DataSource = "GO:BP", hpk1_down_go[hpk1_down_go$ontology == "BP" ,c("term", "over_represented_pvalue")]),
               data.frame(DataSource = "GO:MF", hpk1_down_go[hpk1_down_go$ontology == "MF" ,c("term", "over_represented_pvalue")]),
               with(hpk1_down_kegg[ , c("category", "over_represented_pvalue") ],
                    data.frame(DataSource = "KEGG", term = category, over_represented_pvalue )),
               with(hpk1_down_reac[ , c("category", "over_represented_pvalue") ],
                    data.frame(DataSource = "Reactome", term = category, over_represented_pvalue )))
)

plot(density(hpk1_de_all_enriched$up$over_represented_pvalue))
plot(density(hpk1_de_all_enriched$down$over_represented_pvalue))

# also filter for p-value < 0.01
hpk1_de_all_enriched <- lapply(hpk1_de_all_enriched, function(x){
  x$DataSource <- factor(x$DataSource, levels = c("GO:BP","GO:MF","KEGG", "Reactome"))
  x = x[x$over_represented_pvalue < 0.01 , ]
  return(x)
})

pdf("Treemap- enriched functional annotations for all hpk1 DE genes.pdf",
    width = 12, height = 8)
currRes="up"
  treemap(dtf = wrapify(hpk1_de_all_enriched[[currRes]], "over_represented_pvalue", "term"),
          index = c("DataSource", "term"),
          vSize = "pTransform",
          vColor = "DataSource",
          type = "categorical", 
          force.print.labels = FALSE,
          aspRatio = 1,
          fontsize.labels = c(0, 16),
          position.legend = "right", 
          fontsize.legend = 6,
          drop.unused.levels = FALSE,
          title = paste0("Functional enrichment for all hpk-1 DE genes: \n", currRes))

currRes="down"
treemap(dtf = wrapify(hpk1_de_all_enriched[[currRes]], "over_represented_pvalue", "term"),
        index = c("DataSource", "term"),
        vSize = "pTransform",
        vColor = "DataSource",
        type = "categorical", 
        force.print.labels = FALSE,
        aspRatio = 0.2,
        fontsize.labels = c(0, 16),
        position.legend = "right", 
        fontsize.legend = 6,
        drop.unused.levels = FALSE,
        title = paste0("Functional enrichment for all hpk-1 DE genes: \n", currRes))
  
graphics.off()
rm(currRes)

#
# --- Also output the functional enrichment as barplots ---

table(hpk1_de_all_enriched$up$DataSource)
library(viridis)
library(ggpubr)

# what's a reasonable p-value range to use for the scale?
summary(-1*log10(c(hpk1_up_go$over_represented_pvalue, hpk1_up_kegg$over_represented_pvalue, hpk1_up_reac$over_represented_pvalue)))

# UP AND DOWN SEPARATELY

# -- up --
dat <- hpk1_up_go # UP
dat <- data.frame(dat)

dat.bp <- dat[dat$ontology == "BP" , ]
dat.bp <- dat.bp[ order(dat.bp$over_represented_pvalue, decreasing = TRUE) , ] # note - reverse order
dat.bp$term <- factor(dat.bp$term, levels = dat.bp$term, ordered = TRUE)
dat.bp$ratio = dat.bp$numDEInCat/dat.bp$numInCat

dat.mf <- dat[dat$ontology == "MF" , ]
dat.mf <- dat.mf[ order(dat.mf$over_represented_pvalue, decreasing = TRUE) , ] # note - reverse order
dat.mf$term <- factor(dat.mf$term, levels = dat.mf$term, ordered = TRUE)
dat.mf$ratio = dat.mf$numDEInCat/dat.mf$numInCat
rm(dat)

dat <- rbind(hpk1_up_kegg, hpk1_up_reac) # COMBINE PATHWAY DATASETS
dat.path <- dat[ order(dat$over_represented_pvalue, decreasing = TRUE) , ] # note - reverse order
dat.path$category <- factor(dat.path$category, levels = dat.path$category, ordered = TRUE)
dat.path$ratio = dat.path$numDEInCat/dat.path$numInCat
rm(dat)

gp.up_bp <- ggplot(data = dat.bp, aes(y = term, x = ratio, fill = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 18) +
  ylab("GO Term") +
  xlab("Ratio (# DE / # in term)") +
  labs(fill = "-log10(p-value)") +
  ggtitle("GO Biological Process terms enriched in \nhpk-1 vs N2 DE genes (upregulated)") +
  scale_fill_gradientn(colors = turbo(6), limits = c(1.3, 75), breaks = c(10, 30, 50, 70), values = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  geom_text(aes(x = 0, label = numDEInCat), color = "white",  hjust = 0, nudge_x = 0.01 ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))

gp.up_mf <- ggplot(data = dat.mf, aes(y = term, x = ratio, fill = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 18) +
  ylab("GO Term") +
  xlab("Ratio (# DE / # in term)") +
  labs(fill = "-log10(p-value)") +
  ggtitle("GO Molecular Function terms enriched in \nhpk-1 vs N2 DE genes (upregulated)") +
  scale_fill_gradientn(colors = turbo(6), limits = c(1.3, 75), breaks = c(10, 30, 50, 70), values = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  geom_text(aes(x = 0, label = numDEInCat), color = "white",  hjust = 0, nudge_x = 0.01 ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))

gp.up_path <- ggplot(data = dat.path, aes(y = category, x = ratio, fill = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 18) +
  ylab("Pathway") +
  xlab("Ratio (# DE / # in pathway)") +
  labs(fill = "-log10(p-value)") +
  ggtitle("KEGG and Reactome pathways enriched in \nhpk-1 vs N2 DE genes (upregulated)") +
  scale_fill_gradientn(colors = turbo(6), limits = c(1.3, 75), breaks = c(10, 30, 50, 70), values = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  geom_text(aes(x = 0, label = numDEInCat), color = "white",  hjust = 0, nudge_x = 0.01 ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))

# ggsave(plot = gp.up, filename = file.path(analysis_sets$main$comps$`EAT-2_EV_vs_N2_EV`$outputPathDE,
#                                           "GO BP enrichment- DE upregulated hpk-1 vs N2.pdf"), device = "pdf", width = 12, height = 10, units = "in")
# rm(dat, dat.bp)

# plot these in one view
ggsave(plot = ggarrange(gp.up_bp, gp.up_mf, gp.up_path, ncol = 1, align = "v", heights = c(1.05, 0.35, 0.6), common.legend = TRUE),
       filename = "Enrichment barplot- hpk-1 vs N2 all DE upregulated.pdf", device = "pdf", width = 12, height = 20, units = "in")

# -- down --
dat <- hpk1_down_go
dat <- data.frame(dat)

dat.bp <- dat[dat$ontology == "BP" , ]
dat.bp <- dat.bp[ order(dat.bp$over_represented_pvalue, decreasing = TRUE) , ] # note - reverse order
dat.bp$term <- factor(dat.bp$term, levels = dat.bp$term, ordered = TRUE)
dat.bp$ratio = dat.bp$numDEInCat/dat.bp$numInCat

dat.mf <- dat[dat$ontology == "MF" , ]
dat.mf <- dat.mf[ order(dat.mf$over_represented_pvalue, decreasing = TRUE) , ] # note - reverse order
dat.mf$term <- factor(dat.mf$term, levels = dat.mf$term, ordered = TRUE)
dat.mf$ratio = dat.mf$numDEInCat/dat.mf$numInCat
rm(dat)

dat <- rbind(hpk1_down_kegg, hpk1_down_reac) # COMBINE PATHWAY DATASETS
dat.path <- dat[ order(dat$over_represented_pvalue, decreasing = TRUE) , ] # note - reverse order
dat.path$category <- factor(dat.path$category, levels = dat.path$category, ordered = TRUE)
dat.path$ratio = dat.path$numDEInCat/dat.path$numInCat
rm(dat)

gp.down_bp <- ggplot(data = dat.bp, aes(y = term, x = ratio, fill = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 18) +
  ylab("GO Term") +
  xlab("Ratio (# DE / # in term)") +
  labs(fill = "-log10(p-value)") +
  ggtitle("GO Biological Process terms enriched in \nhpk-1 vs N2 DE genes (downregulated)") +
  scale_fill_gradientn(colors = turbo(6), limits = c(1.3, 75), breaks = c(10, 30, 50, 70), values = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  geom_text(aes(x = 0, label = numDEInCat), color = "white",  hjust = 0, nudge_x = 0.01 ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))

gp.down_mf <- ggplot(data = dat.mf, aes(y = term, x = ratio, fill = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 18) +
  ylab("GO Term") +
  xlab("Ratio (# DE / # in term)") +
  labs(fill = "-log10(p-value)") +
  ggtitle("GO Molecular Function terms enriched in \nhpk-1 vs N2 DE genes (downregulated)") +
  scale_fill_gradientn(colors = turbo(6), limits = c(1.3, 75), breaks = c(10, 30, 50, 70), values = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  geom_text(aes(x = 0, label = numDEInCat), color = "white",  hjust = 0, nudge_x = 0.01 ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))

gp.down_path <- ggplot(data = dat.path, aes(y = category, x = ratio, fill = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 18) +
  ylab("Pathway") +
  xlab("Ratio (# DE / # in pathway)") +
  labs(fill = "-log10(p-value)") +
  ggtitle("KEGG and Reactome pathways enriched in \nhpk-1 vs N2 DE genes (downregulated)") +
  scale_fill_gradientn(colors = turbo(6), limits = c(1.3, 75), breaks = c(10, 30, 50, 70), values = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  geom_text(aes(x = 0, label = numDEInCat), color = "white",  hjust = 0, nudge_x = 0.01 ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))
# ggsave(plot = gp.down, filename = file.path(analysis_sets$main$comps$`EAT-2_EV_vs_N2_EV`$outputPathDE,
#                                             "GO BP enrichment- DE downregulated hpk-1 vs N2.pdf"), device = "pdf", width = 12, height = 10, units = "in")
# rm(dat, dat.bp)
ggsave(plot = ggarrange(gp.down_bp, gp.down_mf, gp.down_path, ncol = 1, align = "v", heights = c(0.3, 0.3, 0.4),
                        common.legend = TRUE),
       filename = "Enrichment barplot- hpk-1 vs N2 all DE downregulated.pdf", device = "pdf", width = 12, height = 8, units = "in")

# ggsave(plot = ggarrange(gp.up_bp, gp.up_mf, gp.up_path,
#                         gp.down_bp, gp.down_mf, gp.down_path,
#                         ncol = 1, align = "v", heights = c(1, 0.3, 0.3, 0.1, 0.1, 0.2), common.legend = TRUE),
#        filename = "Enrichment barplot- hpk-1 vs N2 all DE up and down regulated.pdf", device = "pdf", width = 12, height = 24, units = "in")


#
# For the aging-associated TFs DE in neurons-
# how many *cells* have shared expression of them, based on cengen?
#

# convert this to a list of cells, where each entry has a vector of which of these TFs are expressed in that cell
tmp_cengen_tf_mat <- cengen_sum_level2[rownames(heatmapData_neurons_filt3)[rownames(heatmapData_neurons_filt3) %in% aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno]] , ]

# cengen_TOI_list <- lapply(colnames(tmp_cengen_tf_mat), function(cellName){
#   rownames(tmp_cengen_tf_mat)[tmp_cengen_tf_mat[,cellName] > 0]
# })
# names(cengen_TOI_list) <- colnames(tmp_cengen_tf_mat)

# genewise?
cengen_TOI_list_t <- lapply(rownames(tmp_cengen_tf_mat), function(tfName){
  colnames(tmp_cengen_tf_mat)[tmp_cengen_tf_mat[tfName,] > 0]
})
names(cengen_TOI_list_t) <- rownames(tmp_cengen_tf_mat)

library(UpSetR)

pdf("Num neurons with expression of top aging TFs that are DE with hpk-1.pdf", width = 10, height = 6)
upset(fromList(cengen_TOI_list_t), nsets = 100, order.by = "freq", text.scale = 1.6, set_size.show = TRUE, set_size.scale_max  = 150,
      point.size = 3)
graphics.off()


#
# For TFs which are expressed by do not necessarily change expression levels with age in the same cells as hpk-1 is expressed in:
#   * pull out the predicted target genes, where such information is available
#   * Is there an association between the predicted target genes and age-associated changes?
# The idea being that we can identify putative modulation of TF activity by hpk-1 even when TF expression isn't changing
#

# I already have the TF information loaded.
#  We need to get the corresponding TF gene info, and identify the cells in which they're expressed with hpk-1 at the same time
#   that expression could be at any time point.

# start with all predicted binding interactions in C elegans
#   see how this goes first.
# !!! there are a bunch of motifs that were defined based on complexes of TFs using PBMs, so we need to check if BOTH of
#      those genes are expressed. I've delimited them as "gene-1+&&_gene-2"
# don't have to do anything for the cases where the TF gene symbol is NA

# maybe we can make a motif x tissue matrix to indicate association (needs to be for each time point?)
# this also needs to be just the subset of cell clusters that have DE data

hpk1_expressed_clusters <- lapply(aging_atlas_cluster_mean_counts, function(x){
  as.vector(na.omit(colnames(x)[x["hpk-1",] >= 0.2])) # We'll use some kind of expression threshold, though lenient...
})
# nevermind, we're just using the 48 clusters that hpk-1 is DE with age in :)
#hpk1_age_de_clusters
tf_hpk1_coex_mat <- mat.or.vec(nr = length(unique(sets_tfs_all_ce_long$motif_tf[!is.na(sets_tfs_all_ce_long$geneSymbol)])),
                               nc = length(hpk1_age_de_clusters$cluster_shortID))
rownames(tf_hpk1_coex_mat) <-unique(sets_tfs_all_ce_long$motif_tf[!is.na(sets_tfs_all_ce_long$geneSymbol)])
colnames(tf_hpk1_coex_mat) <- hpk1_age_de_clusters$cluster_shortID

for(motif in rownames(tf_hpk1_coex_mat)){
  for(cluster in colnames(tf_hpk1_coex_mat)){
    curr_tf = unique(sets_tfs_all_ce_long$geneSymbol[sets_tfs_all_ce_long$motif_tf == motif])
    curr_tf = strsplit(curr_tf, split = "_&&_", fixed = TRUE)[[1]]
    ### EDGE CASE- atf-5 name changed to atf-4 recently, so substitute...
    if(curr_tf == "atf-4"){curr_tf = "atf-5"}
    if(!(curr_tf %in% rownames(aging_atlas_cluster_mean_counts$D1))){
      message(paste0(curr_tf, " not in expression set")) # many developmentally expressed TFs
    }else{
      # determine how many timepoints have expression of the relevant TF genes-
      #  if there are two (because they're in a complex) then both need to be expressed
      if(length(curr_tf) == 1){
        tf_hpk1_coex_mat[motif, cluster] <- sum(unlist(lapply(aging_atlas_cluster_mean_counts, function(x){
          # We'll use some kind of expression threshold, though lenient
          as.vector(na.omit(x[curr_tf,cluster] >= 0.2)) 
        })))
      }else if(length(curr_tf) == 2){
        tf_hpk1_coex_mat[motif, cluster] <- length(intersect(unlist(lapply(aging_atlas_cluster_mean_counts, function(x){
          as.vector(na.omit(x[curr_tf[1],cluster] >= 0.2)) 
        })),
        unlist(lapply(aging_atlas_cluster_mean_counts, function(x){
          as.vector(na.omit(x[curr_tf[2],cluster] >= 0.2)) 
        }))))
      }else{
        stop("unexpected condition")
      }
    }
  }
}


table(as.vector(tf_hpk1_coex_mat))
# if we consider only the ones expressed at 4 or 5 time points
# then we have 956 combinations of cell cluster and binding profile to look at in the aging DE data
# which TFs have the most overlap in cluster expression?
sort(rowSums(tf_hpk1_coex_mat >= 4), decreasing = TRUE)
# and which clusters have the most overlapping TFs?
#  The top ones could just be particularly transcriptionally active cells expressing more TFs than other cells
sort(colSums(tf_hpk1_coex_mat >= 4), decreasing = TRUE)

# a number of the "top coexpressed" TFs are also DE with age, hlh-30 and daf-16 for example.

#
# Ok so now we try to come up with some way to evaluate if there's enrichment of the predicted target genes
#  in the aging DE set....
#  easiest would to just look at mean foldchange, or actually look at fisher's exact test for the overlap
#  between the predicted binding set and the aging-DE up or down regulated genes
#  I also considered possibly doing something like how the GSEA score is calculated, by ranking the aging DE genes
# 

# I'll try using fgsea first
#  note that it seems that in the fgsea docs, when they say that the gene names need
#   to be the same between the gene sets to test and the ranked list,
#   they mean that they do not do identifier re-mapping, and so alias or 
#   other types of gene ID mis-matches will be treated as missing.
library(fgsea)

# need to make gene sets of DE genes by cluster that change with age
# These are the ranked gene input lists. Since we're actually using the up/down info, we don't
#  need to separate them by foldchange direction.

# should already be filtered
any(aging_atlas_degs_long$padj > 0.05)
any(abs(aging_atlas_degs_long$logFoldChange) < 0.5)
# vector of ranked LFCs per cluster
rankedAgeClustDegs <- lapply(unique(aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$cluster_shortID %in% 
                                                      colnames(tf_hpk1_coex_mat)[colSums(tf_hpk1_coex_mat >= 4) > 0]]),
       function(clustName){
         datSub = aging_atlas_degs_long[aging_atlas_degs_long$cluster_shortID == clustName,]
         datSub = datSub[ order(datSub$logFoldChange, decreasing = TRUE) , ]
         return(with(datSub,setNames(as.vector(logFoldChange), gene_name)))
       })
names(rankedAgeClustDegs) <- unique(aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$cluster_shortID %in% 
                                                                            colnames(tf_hpk1_coex_mat)[colSums(tf_hpk1_coex_mat >= 4) > 0]])
unlist(lapply(rankedAgeClustDegs, length))

# These will then be tested against the TF binding prediction genes for the given TF
# We're really only interested in running one test per TF and cluster
#   so if we need to consider the p-values, we can run multiple testing correction on them all after

# going to end up with a two-dimensional structure of lists
#  outer list: clusters
#  inner list: motifs
# Since each cluster may have a different set of relevant motifs
# Then we can summarize this on the level of clusters (etc) after
# The adjusted p-values in the original results will not be relevant,
#  since it's being run as an n=1 for all

# took about 4 minutes to run, for TFs with >= 4 timepoints with coexpression
#  (parallelization is supported for fgsea, but only when multiple genesets are specified per ranked list,
#    otherwise we'd have to parallelize across the lists)
gsea_clust_tf_aging <- lapply(names(rankedAgeClustDegs), function(clustName){
  
  rv = lapply(rownames(tf_hpk1_coex_mat)[tf_hpk1_coex_mat[,clustName] >= 4], function(motifName){
    
    predBindGenes = list(sets_tfs_all_ce[[unique(sets_tfs_all_ce_long$motif[sets_tfs_all_ce_long$motif_tf == motifName])]])
    names(predBindGenes) = motifName
    return(fgsea(pathways=predBindGenes, stats=rankedAgeClustDegs[[clustName]], eps = 0)) # eps = 0 for more accurate pvals
  })
  names(rv) = rownames(tf_hpk1_coex_mat)[tf_hpk1_coex_mat[,clustName] >= 4]
  return(rv)
})
names(gsea_clust_tf_aging) <- names(rankedAgeClustDegs)

# example of making a GSEA-style score plot 
#plotEnrichment(tmpPredGenes, rankedAgeClustDegs[[tmpName]])

# summarize all results
# all these results are actually just datatables
#  and aren't custom classes or anything so this is actually trivial
# NOTE THAT THESE ARE WITH *ALL* C ELEGANS PREDICTED TARGETS AND NOT THE HOMOLOGY FILTERED SET
gsea_clust_tf_aging.summary = Reduce(f = rbind, lapply(names(gsea_clust_tf_aging), function(clust_name){
  cbind(cluster = clust_name, data.frame(Reduce(f = rbind, gsea_clust_tf_aging[[clust_name]])))
}))
nrow(gsea_clust_tf_aging.summary)

# so then we can calculate the actual adjusted p-values based on all the tests that were run
gsea_clust_tf_aging.summary$padj <- p.adjust(gsea_clust_tf_aging.summary$pval, method = "fdr")

sum(gsea_clust_tf_aging.summary$padj < 0.05)
sum(gsea_clust_tf_aging.summary$pval < 0.01)

unique(gsea_clust_tf_aging.summary$cluster[ gsea_clust_tf_aging.summary$pval < 0.01])
unique(gsea_clust_tf_aging.summary$pathway[ gsea_clust_tf_aging.summary$pval < 0.01])

unique(gsea_clust_tf_aging.summary$cluster[ gsea_clust_tf_aging.summary$padj < 0.1])
unique(gsea_clust_tf_aging.summary$pathway[ gsea_clust_tf_aging.summary$padj < 0.1])

# add some more cluster info for convenience
gsea_clust_tf_aging.summary$cluster_info <- apply(aging_atlas_deg_clusters[match(gsea_clust_tf_aging.summary$cluster,
                                                                                 aging_atlas_deg_clusters$shortID) , c(3,4,5) ],
                                                  MARGIN = 1, function(x){
                                                    paste(unique(x), sep = "_", collapse = "_")
                                                  })

# plot TF by cluster
#  first create a subset of "significant-ish" results
gsea_clust_tf_aging.sum_sig <- gsea_clust_tf_aging.summary[ (gsea_clust_tf_aging.summary$pval < 0.01) & 
                                                              (gsea_clust_tf_aging.summary$size >= 10), ] 

ggplot(gsea_clust_tf_aging.sum_sig, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=size)) +
  coord_flip() +
  labs(x="TF Motif", y="Normalized Enrichment Score", title="") + 
  facet_wrap(vars(cluster_info)) + 
  theme_bw()

ggplot(gsea_clust_tf_aging.sum_sig, aes(cluster_info, NES)) +
  geom_col(aes(fill=size)) +
  coord_flip() +
  labs(x="Cell Cluster", y="Normalized Enrichment Score", title="") + 
  facet_wrap(vars(pathway), nrow = 2) + 
  theme_bw()

# Just make it a heatmap...
heatmapData <- mat.or.vec(nr = length(unique(gsea_clust_tf_aging.sum_sig$pathway)),
                          nc = length(unique(gsea_clust_tf_aging.sum_sig$cluster)))
rownames(heatmapData) = unique(gsea_clust_tf_aging.sum_sig$pathway)
colnames(heatmapData) = unique(gsea_clust_tf_aging.sum_sig$cluster_info)

for(i in rownames(heatmapData)){
  for(j in colnames(heatmapData)){
    rv = gsea_clust_tf_aging.sum_sig$NES[(gsea_clust_tf_aging.sum_sig$pathway == i) & 
                                           (gsea_clust_tf_aging.sum_sig$cluster_info == j)]
    heatmapData[i, j] <- ifelse(length(rv) == 0, 0, rv)
  }
}
rm(rv)

sum(as.vector(heatmapData != 0))

pdf(file = "GSEA score for TF targets and aging DEGs for TFs in hpk-1 clusters.pdf", width = 14, height = 12)
heatmap.3(t(heatmapData),
          col = colorRampPalette(c("purple4", "white", "gold4"))(64), #scale = "row",
          cexCol = 1, cexRow = 0.8, margins = c(14,30), keysize = 0.8,
          main = "",
          KeyValueName = "GSEA NES",
          rowsep = seq(1, ncol(heatmapData), by = 1), 
          colsep = seq(1, nrow(heatmapData), by = 1),
          sepcolor = "grey90", sepwidth = 0.001)
graphics.off()
rm(heatmapData)

# make a couple of example enrichment plots for significant results
plotEnrichment(sets_tfs_all_ce$cisbp_M00709_2.00, rankedAgeClustDegs[["81_0"]]) # daf-19 - negative example
plotEnrichment(sets_tfs_all_ce$cisbp_M01903_2.00, rankedAgeClustDegs[["124_0"]]) # cey-1 in I4 interneurons - positive example
plotEnrichment(sets_tfs_all_ce$cisbp_M00709_2.00, rankedAgeClustDegs[["124_0"]]) # no enrichment example


# how many of the "significant" TFs here overlap with the aging DE genes?
# note again that atf-4 is known at atf-4 in their dataset


tfs_de_in_hpk1_clust_summary[ tfs_de_in_hpk1_clust_summary$gene_name %in%  
                                c("lim-7", "fos-1", "atf-5", "xbp-1", "crh-1", "daf-16", "let-381", "fkh-7", "cey-1",
                                "syd-9", "gei-3", "egl-27", "klu-2", "atf-7", "daf-19", "nhr-88", "skn-1") , ]

setdiff( c("lim-7", "fos-1", "atf-5", "xbp-1", "crh-1", "daf-16", "let-381", "fkh-7", "cey-1",
           "syd-9", "gei-3", "egl-27", "klu-2", "atf-7", "daf-19", "nhr-88", "skn-1"),
         tfs_de_in_hpk1_clust_summary$gene_name)

#
# TFs expressed in hpk-1-DE clusters and aging target gene DE- Version 2
#  * Don't filter by number of timepoints with expression of the TF
#  * Use Fisher's Exact test instead of GSEA (so don't consider up/down direction)
# So we won't necessarily need ranked genes for this
#

# including all clusters with any TF that co-expresses with hpk-1
#  built the relevant lists of DE genes for each cluster
ageClustDegs <- lapply(unique(aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$cluster_shortID %in% 
                                                                            colnames(tf_hpk1_coex_mat)[colSums(tf_hpk1_coex_mat >= 1) > 0]]),
                             function(clustName){
                               datSub = aging_atlas_degs_long[aging_atlas_degs_long$cluster_shortID == clustName,]
                               datSub = datSub[ order(datSub$logFoldChange, decreasing = TRUE) , ]
                               return(with(datSub,setNames(as.vector(logFoldChange), gene_name)))
                             })
names(ageClustDegs) <- unique(aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$cluster_shortID %in% 
                                                                            colnames(tf_hpk1_coex_mat)[colSums(tf_hpk1_coex_mat >= 1) > 0]])
unlist(lapply(rankedAgeClustDegs, length))
unlist(lapply(ageClustDegs, length)) # only adds one cluster

# These will then be tested against the TF binding prediction genes for the given TF
# We're really only interested in running one test per TF and cluster
#   so if we need to consider the p-values, we can run multiple testing correction on them all after
# going to end up with a two-dimensional structure of lists
#  outer list: clusters
#  inner list: motifs
# Since each cluster may have a different set of relevant motifs
# Then we can summarize this on the level of clusters (etc) after
# The adjusted p-values in the original results will not be relevant,
#  since it's being run as an n=1 for all
#
# hypergeometric enrichment
#  "universe": # of genes with detectable expression at any time point in the cluster
#  Set A: DE genes in cluster
#  Set B: Genes with predicted binding for the given TF motif (that are expressed at any time point in the cluster)
#  specified as formulated for phyper with lower.tail as false: hitInSample-1 (q), hitInPop (m), failInPop (n), sampleSize (k)
#  m = the # of genes of Set A 
#  n = length(universe) - length(set A)
#  k = length(set B)
#  q = length(intersect(Set A, Set B)) - 1

hypergeo_clust_tf_aging <- lapply(names(ageClustDegs), function(clustName){
  
  # determine set of genes with expression in the cluster at any timepoint
  # (rownames are the same for all matrices in the list)
  clust_expressed_genes <- rownames(aging_atlas_cluster_mean_counts$D1)[(aging_atlas_cluster_mean_counts$D1[,clustName] > 0) |
    (aging_atlas_cluster_mean_counts$D3[,clustName] > 0) |
    (aging_atlas_cluster_mean_counts$D5[,clustName] > 0) |
    (aging_atlas_cluster_mean_counts$D8[,clustName] > 0) |
    (aging_atlas_cluster_mean_counts$D11[,clustName] > 0)]
  
  # Cluster age DE genes
  currDegs = names(ageClustDegs[[clustName]])
  # presumably these are all expressed in the cluster at some point, but we'll make sure
  message(paste0("total curr degs in cluster: ", length(currDegs)))
  currDegs = currDegs[currDegs %in% clust_expressed_genes]
  message(paste0("Curr degs with expression in cluster (should be same): ", length(currDegs)))
  
  rv = lapply(rownames(tf_hpk1_coex_mat)[tf_hpk1_coex_mat[,clustName] >= 1], function(motifName){ # any TF with coexpression
    
    # get all the genes with binding predictions for the current TF
    predBindGenes = sets_tfs_all_ce[[unique(sets_tfs_all_ce_long$motif[sets_tfs_all_ce_long$motif_tf == motifName])]]
    predBindGenes2 = predBindGenes[predBindGenes %in% clust_expressed_genes]
    message(paste0("# genes with pred binding for ", motifName, ": ", length(predBindGenes),
                   " After considering expression in cluster: ", length(predBindGenes2)))
    
    m = length(currDegs)
    n = length(setdiff(clust_expressed_genes, currDegs))
    k = length(predBindGenes2)
    q = length(intersect(currDegs, predBindGenes2))

    return(phyper(q-1, m, n, k, lower.tail = FALSE))
  })
  names(rv) = rownames(tf_hpk1_coex_mat)[tf_hpk1_coex_mat[,clustName] >= 1]
  return(rv)
})
names(hypergeo_clust_tf_aging) <- names(ageClustDegs)

# example
hypergeo_clust_tf_aging$`100_0`$`cisbp_M00342_2.00_vab-3`


# summarize all results
# NOTE THAT THESE ARE WITH *ALL* C ELEGANS PREDICTED TARGETS AND NOT THE HOMOLOGY FILTERED SET
hypergeo_clust_tf_aging.summary = Reduce(f = rbind, lapply(names(hypergeo_clust_tf_aging), function(clust_name){
  data.frame(cluster = clust_name,
        motif = names(hypergeo_clust_tf_aging[[clust_name]]),
        pval = as.numeric(unlist(hypergeo_clust_tf_aging[[clust_name]])))
}))
nrow(hypergeo_clust_tf_aging.summary)

# so then we can calculate the actual adjusted p-values based on all the tests that were run
hypergeo_clust_tf_aging.summary$padj <- p.adjust(hypergeo_clust_tf_aging.summary$pval, method = "fdr")

# how many uncorrected vs corrected pvalues < 0.05?
sum(hypergeo_clust_tf_aging.summary$pval < 0.05)
sum(hypergeo_clust_tf_aging.summary$padj < 0.05)

hypergeo_clust_tf_aging.summary$cluster_info <- apply(aging_atlas_deg_clusters[match(hypergeo_clust_tf_aging.summary$cluster,
                                                                                 aging_atlas_deg_clusters$shortID) , c(3,4,5) ],
                                                  MARGIN = 1, function(x){
                                                    paste(unique(x), sep = "_", collapse = "_")
                                                  })
#  First, FDR p-value < 0.05
hypergeo_clust_tf_aging.sum_sig <- hypergeo_clust_tf_aging.summary[ hypergeo_clust_tf_aging.summary$padj < 0.05 , ]

heatmapData <- mat.or.vec(nr = length(unique(hypergeo_clust_tf_aging.sum_sig$motif)),
                          nc = length(unique(hypergeo_clust_tf_aging.sum_sig$cluster)))
rownames(heatmapData) = unique(hypergeo_clust_tf_aging.sum_sig$motif)
colnames(heatmapData) = unique(hypergeo_clust_tf_aging.sum_sig$cluster_info)

for(i in rownames(heatmapData)){
  for(j in colnames(heatmapData)){
    rv = hypergeo_clust_tf_aging.sum_sig$padj[(hypergeo_clust_tf_aging.sum_sig$motif == i) & 
                                           (hypergeo_clust_tf_aging.sum_sig$cluster_info == j)]
    heatmapData[i, j] <- ifelse(length(rv) == 0, NA, rv) # we don't really want to have to assume that everything non-significant is 1
  }
}
rm(rv)

sum(as.vector(heatmapData != 0), na.rm = TRUE)

heatmapData <- -1*log10(heatmapData)
# Sort by # non-NA per row
heatmapData <- heatmapData[order(rowSums(!is.na(heatmapData)), decreasing = TRUE) , ]
heatmapData <- heatmapData[ , order(colSums(!is.na(heatmapData)), decreasing = FALSE) ]

# add some annotation- 
#  indicate if a TF was also differentially expressed with age (since this is per-cluster-
#       indicate # of clusters among the ones shown??)
#  indicate a TF aging association?
#  indicate tissue types for columns, and then drop the last concatenated description in the names
uq_motifs <- unique(sets_tfs_all_ce_long[,c("geneSymbol", "motif_tf")])
uq_motifs <- uq_motifs[ uq_motifs$motif_tf %in% rownames(heatmapData) , ]
uq_motifs <- uq_motifs[ match(rownames(heatmapData), uq_motifs$motif_tf) ,  ]
motif_tf_hmap_de_counts <- sapply(uq_motifs$geneSymbol, function(gene){
  gene = unlist(strsplit(gene, "_&&_"))
  if(length(gene) == 1){
    sum(gene == aging_atlas_degs_long$gene_name[ aging_atlas_degs_long$cluster_shortID %in%
                                       hypergeo_clust_tf_aging.sum_sig$cluster ])
  }else if(length(gene) == 2){
    length(unique(aging_atlas_degs_long$cluster_shortID[ 
      (aging_atlas_degs_long$cluster_shortID %in% hypergeo_clust_tf_aging.sum_sig$cluster) & 
        (aging_atlas_degs_long$gene_name %in% gene)]))
  }
})

# from https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

tmp = unique(data.frame(map2color(motif_tf_hmap_de_counts, viridis::mako(33)),motif_tf_hmap_de_counts))
tmp = tmp[ order(as.numeric(tmp[,2]), decreasing = FALSE) , ]

pdf(file = "color key for num clusters de for TFs in tf targets and aging degs.pdf", width = 2, height = 3)
swatch2(setNames(tmp[,1], tmp[,2]))
graphics.off()
rm(tmp)

tmp = unique(cbind(tissue_type_vec$type_color, tissue_type_vec$tissue_type))
swatch2(setNames(tmp[,1], tmp[,2]))

# ok and now the tissue indication like we've been adding
tissue_type_vec$type_color
temp_tissue_type <- tissue_type_vec[ match(hypergeo_clust_tf_aging.sum_sig$cluster[match(colnames(heatmapData),
                                              hypergeo_clust_tf_aging.sum_sig$cluster_info)], tissue_type_vec$shortID) , ]

svg(file = "Hypergeometric enrichment for TF targets and aging DEGs for TFs in hpk-1 clusters (FDR 05).svg", width = 18, height = 18)
heatmap.3(heatmapData,
          col = viridis::viridis(64), #scale = "row",
          cexCol = 1, cexRow = 0.8, margins = c(23,29), keysize = 0.8,
          main = "",
          Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          RowSideColors = t(as.matrix(map2color(motif_tf_hmap_de_counts, viridis::mako(33)))),
          ColSideColors = as.matrix(temp_tissue_type$type_color),
          labCol = aging_atlas_clusters_all[ match(hypergeo_clust_tf_aging.sum_sig$cluster[match(colnames(heatmapData),
                                                                                                 hypergeo_clust_tf_aging.sum_sig$cluster_info)], 
                                                   aging_atlas_clusters_all$shortID), ]$Clusters.explanatory.name,
          KeyValueName = "-log10(FDR p-value)",
          rowsep = seq(1, nrow(heatmapData), by = 1), 
          colsep = seq(1, ncol(heatmapData), by = 1),
          sepcolor = "grey90", sepwidth = c(0.01, 0.01))
graphics.off()
rm(heatmapData)


#
# Andy asked about if there are any motifs enriched for binding predictions to the
#  TFs that are DE with age in the same clusters as hpk-1... which I didn't do yet.
#  start simple- TFs upregulated and downregulated (separately) with age in any 
#  cluster that hpk-1 is DE in.
#

# subset to TFs that are upregulated
ageClustDegs_tf_up <- lapply(ageClustDegs, function(clustDEGs){
  return(names(clustDEGs)[(clustDEGs > 0) & (names(clustDEGs) %in% wTF$Public.name)])
})
sort(table(unlist(ageClustDegs_tf_up)))
length(unique(unlist(ageClustDegs_tf_up)))

ageClustDegs_tf_down <- lapply(ageClustDegs, function(clustDEGs){
  return(names(clustDEGs)[(clustDEGs < 0) & (names(clustDEGs) %in% wTF$Public.name)])
})
sort(table(unlist(ageClustDegs_tf_down)))
length(unique(unlist(ageClustDegs_tf_down)))

# ok, now run enrichment vs the TF motifs for the unique genes across these
# hypergeometric enrichment
#  "universe": # of genes with detectable expression at any time point
#  Set A: DE genes upregulated (or downregulated) across any of the clusters that are TFs
#  Set B: Genes with predicted binding for the given TF motif 
#  specified as formulated for phyper with lower.tail as false: hitInSample-1 (q), hitInPop (m), failInPop (n), sampleSize (k)
#  m = the # of genes of Set A 
#  n = length(universe) - length(set A)
#  k = length(set B)
#  q = length(intersect(Set A, Set B)) - 1

expressed_genes <- unique(c(rownames(aging_atlas_cluster_mean_counts$D1)[rowSums(aging_atlas_cluster_mean_counts$D1 > 0.5,
                                                     na.rm = TRUE) > 0],
         rownames(aging_atlas_cluster_mean_counts$D3)[rowSums(aging_atlas_cluster_mean_counts$D3 > 0.5,
                                                              na.rm = TRUE) > 0],
         rownames(aging_atlas_cluster_mean_counts$D5)[rowSums(aging_atlas_cluster_mean_counts$D5 > 0.5,
                                                              na.rm = TRUE) > 0],
         rownames(aging_atlas_cluster_mean_counts$D8)[rowSums(aging_atlas_cluster_mean_counts$D8 > 0.5,
                                                              na.rm = TRUE) > 0],
         rownames(aging_atlas_cluster_mean_counts$D11)[rowSums(aging_atlas_cluster_mean_counts$D11 > 0.5,
                                                              na.rm = TRUE) > 0]
))


DE_TF_motifEnrich <- Reduce(f = rbind, x = lapply(names(sets_tfs_all_ce), function(motifName){
  retVal = list()
  
  ### UP
  m = length(unique(unlist(ageClustDegs_tf_up)))
  n = length(setdiff(expressed_genes, unique(unlist(ageClustDegs_tf_up))))
  k = length(sets_tfs_all_ce[[motifName]])
  q = length(intersect(unique(unlist(ageClustDegs_tf_up)), sets_tfs_all_ce[[motifName]]))
  
  retVal$up = phyper(q-1, m, n, k, lower.tail = FALSE)
  rm(m,n,k,q)
  
  ### DOWN
  m = length(unique(unlist(ageClustDegs_tf_down)))
  n = length(setdiff(expressed_genes, unique(unlist(ageClustDegs_tf_down))))
  k = length(sets_tfs_all_ce[[motifName]])
  q = length(intersect(unique(unlist(ageClustDegs_tf_down)), sets_tfs_all_ce[[motifName]]))
  
  retVal$down = phyper(q-1, m, n, k, lower.tail = FALSE)
  rm(m,n,k,q)
  
  return(c("motif" = motifName, "upTFs_pval" = retVal$up, "downTFs_pval" = retVal$down))
  
}))

DE_TF_motifEnrich <- merge(DE_TF_motifEnrich, unique(sets_tfs_species_filt_long[, c("motif", "motif_tf")]),
      by = "motif")

DE_TF_motifEnrich.sig <- DE_TF_motifEnrich[ (DE_TF_motifEnrich$upTFs_pval < 0.05) |
                                              (DE_TF_motifEnrich$downTFs_pval < 0.05), ]

#
# Actually, what he wanted was just the known aging-associated TFs that
#  that are regulated in any neurons, plus hpk-1
#

expressed_genes_neurons <- unique(unlist(lapply(aging_atlas_cluster_mean_counts, function(datMat){
  rownames(datMat)[
    rowSums(datMat[,na.omit(aging_atlas_clusters_all$shortID[
      aging_atlas_clusters_all$`Tissue.type.(higher.degree)` == "neuron"])] > 0.5,
      na.rm = TRUE) > 0]
})))

age_tfs_neuron_up <- unique(c(with(aging_atlas_degs_long, gene_name[(gene_name %in% lifespan_genes_v2$GeneSymbol) & 
                                                                      (gene_name %in% wTF$Public.name) &
                                                                      (cluster_shortID %in%
                                                                         na.omit(aging_atlas_clusters_all$shortID[
                                                                           aging_atlas_clusters_all$`Tissue.type.(higher.degree)` == "neuron"])) &
                                                                      (logFoldChange > 0) ]), "hpk-1")) # andy asked for hpk-1 to be in there

age_tfs_neuron_down <- unique(with(aging_atlas_degs_long, gene_name[(gene_name %in% lifespan_genes_v2$GeneSymbol) & 
                                                                      (gene_name %in% wTF$Public.name) &
                                                                      (cluster_shortID %in%
                                                                         na.omit(aging_atlas_clusters_all$shortID[
                                                                           aging_atlas_clusters_all$`Tissue.type.(higher.degree)` == "neuron"])) &
                                                                      (logFoldChange < 0) ]))
length(age_tfs_neuron_up)
length(age_tfs_neuron_down)


DE_TF_motifEnrich_2 <- Reduce(f = rbind, x = lapply(names(sets_tfs_all_ce), function(motifName){
  retVal = list()
  
  ### UP
  m = length(age_tfs_neuron_up)
  n = length(setdiff(expressed_genes_neurons, age_tfs_neuron_up))
  k = length(sets_tfs_all_ce[[motifName]])
  q = length(intersect(age_tfs_neuron_up, sets_tfs_all_ce[[motifName]]))
  
  retVal$up = phyper(q-1, m, n, k, lower.tail = FALSE)
  rm(m,n,k,q)
  
  ### DOWN
  m = length(age_tfs_neuron_down)
  n = length(setdiff(expressed_genes_neurons, age_tfs_neuron_down))
  k = length(sets_tfs_all_ce[[motifName]])
  q = length(intersect(age_tfs_neuron_down, sets_tfs_all_ce[[motifName]]))
  
  retVal$down = phyper(q-1, m, n, k, lower.tail = FALSE)
  rm(m,n,k,q)
  
  ### UP AND DOWN
  
  m = length(unique(c(age_tfs_neuron_up, age_tfs_neuron_down)))
  n = length(setdiff(expressed_genes_neurons, unique(c(age_tfs_neuron_up, age_tfs_neuron_down))))
  k = length(sets_tfs_all_ce[[motifName]])
  q = length(intersect(unique(c(age_tfs_neuron_up, age_tfs_neuron_down)), sets_tfs_all_ce[[motifName]]))
  
  retVal$combined = phyper(q-1, m, n, k, lower.tail = FALSE)
  rm(m,n,k,q)
  
  return(c("motif" = motifName, "upTFs_pval" = retVal$up,
           "downTFs_pval" = retVal$down,
           "combinedTFs_pval" = retVal$combined))
  
}))

DE_TF_motifEnrich_2 <- merge(DE_TF_motifEnrich_2, unique(sets_tfs_species_filt_long[, c("motif", "motif_tf")]),
                           by = "motif")

DE_TF_motifEnrich_2.sig <- DE_TF_motifEnrich_2[ (DE_TF_motifEnrich_2$upTFs_pval < 0.05) |
                                              (DE_TF_motifEnrich_2$downTFs_pval < 0.05) |
                                                (DE_TF_motifEnrich_2$combinedTFs_pval < 0.05), ]





#
# At what time point does hpk-1 expression start to increase?
# Is the change linear?
# Consistent across cell clusters?
#
#... parallel coordinate plot?
# 
library(GGally)

# clusters of interest
coi <- aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$gene_name == "hpk-1"]

hpk1_expr_by_cluster <- lapply(coi, function(clust){
  aging_atlas_within_tissue[[clust]][ "hpk-1" , ]
})
hpk1_expr_by_cluster <- Reduce(hpk1_expr_by_cluster, f = rbind)
rownames(hpk1_expr_by_cluster) <- coi

# convert to long format
hpk1_expr_by_cluster_long <- reshape2::melt(as.matrix(hpk1_expr_by_cluster))
colnames(hpk1_expr_by_cluster_long) <- c("cluster", "time", "mean_count")

# exclude day 15- too many missing values
hpk1_expr_by_cluster_df <- data.frame(hpk1_expr_by_cluster)
#hpk1_expr_by_cluster_df <- hpk1_expr_by_cluster_df[,-ncol(hpk1_expr_by_cluster_df)]
hpk1_expr_by_cluster_df$cluster = rownames(hpk1_expr_by_cluster_df)

hpk1_expr_by_cluster_df <- merge(hpk1_expr_by_cluster_df, aging_atlas_clusters_all[,c("shortID", "Tissue.type.(higher.degree)")],
                                 by.x = "cluster", by.y = "shortID", all.x = TRUE)
colnames(hpk1_expr_by_cluster_df)[ncol(hpk1_expr_by_cluster_df)] <- "tissue_type"

ggsave(plot = ggparcoord(hpk1_expr_by_cluster_df, showPoints = TRUE, groupColumn = "tissue_type",
           columns = 2:6, scale = "globalminmax", boxplot = TRUE, alphaLines = 0.5, missing = "median") +
  facet_wrap(~tissue_type) + ggtitle("hpk-1 mean counts across clusters with hpk-1 DE in old vs young") +
    scale_color_manual(name = "tissue_type",
                       values = setNames(tissue_type_vec[tissue_type_vec$tissue_type %in% hpk1_expr_by_cluster_df$tissue_type,2],
                                         tissue_type_vec[tissue_type_vec$tissue_type %in% hpk1_expr_by_cluster_df$tissue_type,1])) +
    theme_gray(),
  filename = "hpk-1 expression across DE clusters, parallel coordinate plot.pdf", device = "pdf", width = 10, height = 8)

# order by p-value?
hpk1_expr_by_cluster <- hpk1_expr_by_cluster[ aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$gene_name == "hpk-1"][
  order(aging_atlas_degs_long$padj[aging_atlas_degs_long$gene_name == "hpk-1"], decreasing = FALSE)] , ]
# order by FC?
hpk1_expr_by_cluster <- hpk1_expr_by_cluster[ aging_atlas_degs_long$cluster_shortID[aging_atlas_degs_long$gene_name == "hpk-1"][
  order(aging_atlas_degs_long$logFoldChange[aging_atlas_degs_long$gene_name == "hpk-1"], decreasing = TRUE)] , ]


pdf(file = "hpk-1 expression across time in hpk-1 DE clusters.pdf", width = 12, height = 10)
heatmap.3(hpk1_expr_by_cluster, Colv = FALSE, Rowv = FALSE, dendrogram = "none", 
          col = viridis::viridis(64),
          cexCol = 1, cexRow = 0.8, margins = c(8,25), keysize = 0.8,
          colsep = 3,
          main = "",
          RowSideColors = rbind("log2FC" = map2color(with(aging_atlas_degs_long[aging_atlas_degs_long$gene_name == "hpk-1" , ], 
                                                     logFoldChange[match(rownames(hpk1_expr_by_cluster), cluster_shortID)]), viridis::plasma(33),
                                                limits = c(0, 1.5)),
                            "tissue_type" = tissue_type_vec$type_color[match(rownames(hpk1_expr_by_cluster), tissue_type_vec$shortID)]),
          labRow = aging_atlas_deg_clusters$Curated.annotations[match(rownames(hpk1_expr_by_cluster), aging_atlas_deg_clusters$shortID)])
graphics.off()


tmp = unique(data.frame(map2color(with(aging_atlas_degs_long[aging_atlas_degs_long$gene_name == "hpk-1" , ], 
                                       logFoldChange[match(rownames(hpk1_expr_by_cluster), cluster_shortID)]),
                                  viridis::plasma(33),
                                  limits = c(0, 1.5)) ,
             with(aging_atlas_degs_long[aging_atlas_degs_long$gene_name == "hpk-1" , ], 
                  logFoldChange[match(rownames(hpk1_expr_by_cluster), cluster_shortID)])))
tmp = tmp[ order(as.numeric(tmp[,2]), decreasing = TRUE) , ]
tmp[,2] <- signif(tmp[,2], 2)
tmp <- unique(tmp)

pdf(file = "color key for hpk-1 DE log2FC sidebar in heatmap.pdf", width = 2, height = 3)
swatch2(setNames(tmp[,1], tmp[,2]))
graphics.off()
rm(tmp)


#
# ok so now break that up into neuron functional subtypes-
#  motor neurons, interneurons, sensory neurons
#
# start by taking the existing object, subset it to neurons, then add the neuron subtype info

hpk1_expr_by_cluster_df_neuron <- hpk1_expr_by_cluster_df[ hpk1_expr_by_cluster_df$tissue_type == "neuron" , ]

hpk1_expr_by_cluster_df_neuron <- merge(hpk1_expr_by_cluster_df_neuron, atlas_cengen_neuron_map[,c("atlas_cluster", "NeuronType")],
                                        by.x = "cluster", by.y = "atlas_cluster", all.x = TRUE)

# some of these aren't annotated in cengen, add annotation manually
hpk1_expr_by_cluster_df_neuron$cluster[is.na(hpk1_expr_by_cluster_df_neuron$NeuronType)]
# 113_0 = Interneuron
# 126_0 = Motorneuron
# 128_0 = Interneuron
# 138_0 = Interneuron
# 29_0 = Interneuron
# 58_1 = Interneuron
# 68_0 = Motorneuron
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "113_0"] <- "Interneuron"
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "126_0"] <- "Motorneuron"
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "128_0"] <- "Interneuron"
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "138_0"] <- "Interneuron"
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "29_0"] <- "Interneuron"
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "58_1"] <- "Interneuron"
hpk1_expr_by_cluster_df_neuron$NeuronType[hpk1_expr_by_cluster_df_neuron$cluster == "68_0"] <- "Motorneuron"


table(hpk1_expr_by_cluster_df_neuron$NeuronType) # looks ok

ggsave(plot = ggparcoord(hpk1_expr_by_cluster_df_neuron, showPoints = TRUE, groupColumn = "NeuronType",
                         columns = 2:6, scale = "globalminmax", boxplot = TRUE, alphaLines = 0.3, missing = "median") +
         facet_wrap(~NeuronType) + ggtitle("hpk-1 mean counts across Neuron subtypes with hpk-1 DE in old vs young")  +
         theme_gray(),
       filename = "hpk-1 expression across DE neuron subtype clusters, parallel coordinate plot.pdf", device = "pdf", width = 10, height = 8)


hpk1_expr_by_cluster_df_neuron$young <- rowMeans(hpk1_expr_by_cluster_df_neuron[,c("d1", "d3", "d5")], na.rm = TRUE)
hpk1_expr_by_cluster_df_neuron$old <- rowMeans(hpk1_expr_by_cluster_df_neuron[,c("d8", "d11", "d15")], na.rm = TRUE)

# average young and old
ggsave(plot =
ggparcoord(hpk1_expr_by_cluster_df_neuron, showPoints = TRUE, groupColumn = "NeuronType",
           columns = 10:11, scale = "globalminmax", boxplot = TRUE, alphaLines = 0, missing = "median") +
  facet_wrap(~NeuronType) + ggtitle("hpk-1 mean counts across Neuron subtypes with hpk-1 DE in old vs young")  +
  theme_gray(),
filename = "hpk-1 expression age AVERAGE across DE neuron subtype clusters, parallel coordinate plot.pdf", device = "pdf", width = 10, height = 8)

# split into age class



#  spatial correlation in nervous system only vs non-neuronal (instead of correlations across all clusters)
#  phosphatases that change with age, any association with hpk-1? cross-reference Zack's list

#
# Back to basics- pathway enrichment in all hpk-1 DE genes
#  We've been focusing mostly on neurons, but we do need to show more than just neuronal functions to set the stage
#
#

sets_kegg_all <- qusage::read.gmt(file.path(paths$seq_data, "genesets", "KEGG C elegans gene sets- all pathways and genes- 2022-06-22.gmt"))
sets_kegg_all_long <-  data.frame(Reduce(x = lapply(names(sets_kegg_all), function(x){
  cbind("pathway" = x, "gene_id" = sets_kegg_all[[x]])}), f = rbind))

# +++ REACTOME - pathways as gene sets - WBGENEID +++
# I updated this set 07Sept2022
sets_reactome <- qusage::read.gmt(file.path("E:/public_data_annotation/reactome", "NCBI2Reactome_07Sept2022.gmt"))
# check to make sure all pathways have a sane number of members
# just keep ones with at least 10 members
#table(unlist(lapply(sets_reactome, length)))
#names(sets_reactome)[unlist(lapply(sets_reactome, length)) < 10]
sets_reactome <- sets_reactome[unlist(lapply(sets_reactome, length)) >= 10]

sets_reactome_long <-  data.frame(Reduce(x = lapply(names(sets_reactome), function(x){
  cbind("pathway" = x, "gene_id" = sets_reactome[[x]])}), f = rbind))

# +++ Gene ontology +++
elegansGO <- read.table(file.path("E:/public_data_annotation/wormbase/WS284",
                                  "c_elegans.PRJNA13758.WS284.gene_association.wb"),
                        comment.char="!", sep = "\t", stringsAsFactors = FALSE, quote="")[,c(2,3,5,7)]
colnames(elegansGO) <- c("GeneID", "GeneSymbol", "GO", "Evidence")

# assumes lengthdata is already defined
goSeqEnrich_singleSet <- function(geneVec, setLongObj, setSelCols){

  # wb <- createWorkbook()
  # gene sig vector is 1 for the genes in the geneVec and 0 otherwise
  sig = numeric(length(lengthData$length))
  names(sig) = lengthData$WBGeneID
  sig[names(sig) %in% geneVec] = 1
  pwf=nullp(sig, bias.data = lengthData$length)
  enrichRes <- goseq(pwf = pwf, #method = "hypergeometric",
                     gene2cat = setLongObj[, setSelCols ])  
  enrichRes[,"over_represented_pvalue"] <- p.adjust(enrichRes$over_represented_pvalue, method="BH")
  enrichRes[,"under_represented_pvalue"] <- p.adjust(enrichRes$under_represented_pvalue, method="BH")
  enrichRes.sig.all <- enrichRes[ enrichRes$over_represented_pvalue < 0.05 , ]
  if(nrow(enrichRes.sig.all) > 0){
    
    enrichRes.sig.all$intersect_WBGeneID <- unlist(lapply(enrichRes.sig.all$category, function(x){
      paste(intersect(setLongObj[ setLongObj[, setSelCols[2]] == x , setSelCols[1]], # genes in cat
                      names(sig)[sig == 1]), sep = ", ", collapse = ", ") # genes sig in query set
    }))
    # for convenience, also get the gene symbols
    enrichRes.sig.all$intersect_genename <- unlist(lapply(enrichRes.sig.all$category, function(x){
      paste(wb_gene_info$gene_name[wb_gene_info$WBGeneID %in% intersect(setLongObj[ setLongObj[, setSelCols[2]] == x , setSelCols[1]], # genes in cat
                                                     names(sig)[sig == 1])], sep = ", ", collapse = ", ") # genes sig in query set
    }))
  }
  # addWorksheet(wb, sheetName = paste("enriched"))
  # writeDataTable(wb, sheet = 1, x = as.data.frame(enrichRes.sig.all), rowNames = FALSE)
  # 
  # saveWorkbook(wb, file.path(paths$output, "set_enrichment", paste0(filePrefix, ".xlsx")),
  #              overwrite = TRUE)
  return(enrichRes.sig.all)
}



hpk1_up_kegg = goSeqEnrich_singleSet(geneVec = DE_hpk1$WBGeneID[DE_hpk1$DESeq2_log2FC > 0],
                             setLongObj = sets_kegg_all_long, setSelCols = c("gene_id", "pathway"))

hpk1_down_kegg = goSeqEnrich_singleSet(geneVec = DE_hpk1$WBGeneID[DE_hpk1$DESeq2_log2FC < 0],
                                     setLongObj = sets_kegg_all_long, setSelCols = c("gene_id", "pathway"))

hpk1_up_reac = goSeqEnrich_singleSet(geneVec = DE_hpk1$WBGeneID[DE_hpk1$DESeq2_log2FC > 0],
                                     setLongObj = sets_reactome_long, setSelCols = c("gene_id", "pathway"))

hpk1_down_reac = goSeqEnrich_singleSet(geneVec = DE_hpk1$WBGeneID[DE_hpk1$DESeq2_log2FC < 0],
                                       setLongObj = sets_reactome_long, setSelCols = c("gene_id", "pathway"))

hpk1_up_go = goSeqEnrich_singleSet(geneVec = DE_hpk1$WBGeneID[DE_hpk1$DESeq2_log2FC > 0],
                                     setLongObj = elegansGO, setSelCols = c("GeneID", "GO"))

hpk1_down_go = goSeqEnrich_singleSet(geneVec = DE_hpk1$WBGeneID[DE_hpk1$DESeq2_log2FC < 0],
                                       setLongObj = elegansGO, setSelCols = c("GeneID", "GO"))


write.xlsx(x = rbind(hpk1_up_kegg,hpk1_up_reac), file = "hpk-1 DE up all enrichment in KEGG and Reactome (goseq).xlsx")

write.xlsx(x = rbind(hpk1_down_kegg,hpk1_down_reac), file = "hpk-1 DE down all enrichment in KEGG and Reactome (goseq).xlsx")


write.xlsx(x = hpk1_up_go, file = "hpk-1 DE up all enrichment in GO (goseq).xlsx")
write.xlsx(x = hpk1_down_go, file = "hpk-1 DE down all enrichment in GO (goseq).xlsx")



# TF enrichment for neuron associated genes DE in hpk-1
#  going to use the union of the two neuron lists I have- upregulated only- there are a few down genes on the second list
#  but only 10, so not enough to see any significant enrichment

hpk1_de_up_neuron_list_union <- union(neuron_gene_classes_hpk1DE$gene[neuron_gene_classes_hpk1DE$DESeq2_log2FC > 0],
                                      hpk1_up_neuronenrich_summary_collapse$gene_name)
length(hpk1_de_up_neuron_list_union)

hpk1_de_up_neuron_list_union <- merge(hpk1_de_up_neuron_list_union, wb_gene_info[, c("WBGeneID", "gene_name")],
                                      by.x = 1, by.y = 2)
colnames(hpk1_de_up_neuron_list_union)[1] <- "gene_name"

#
# !!! ASSUMING I ALREADY HAVE WHAT I NEED LOADED FROM EARLIER
#
head(lengthData)
nrow(lengthData)

# merge with length info
lengthData_currDE <- lengthData
lengthData_currDE$GOI <- lengthData$WBGeneID %in% hpk1_de_up_neuron_list_union$WBGeneID
sum(lengthData_currDE$GOI) # genes of interest


pwf=nullp(lengthData_currDE$GOI, bias.data = lengthData_currDE$length)
rownames(pwf) <- lengthData_currDE$WBGeneID
tfEnrich_neuronDE_hpk1 <- goseq(pwf, gene2cat = sets_tfs_species_filt_long[, c("WBGeneID", "motif_tf") ])
tfEnrich_neuronDE_hpk1[,"over_represented_pvalue"] <- p.adjust(tfEnrich_neuronDE_hpk1$over_represented_pvalue, method="BH")
tfEnrich_neuronDE_hpk1[,"under_represented_pvalue"] <- p.adjust(tfEnrich_neuronDE_hpk1$under_represented_pvalue, method="BH")

tfEnrich_neuronDE_hpk1 <- tfEnrich_neuronDE_hpk1[ tfEnrich_neuronDE_hpk1$over_represented_pvalue < 0.05 , ]

tfEnrich_neuronDE_hpk1$intersect_WBGeneID <- unlist(lapply(tfEnrich_neuronDE_hpk1$category, function(x){
  paste(intersect(sets_tfs_species_filt_long[ sets_tfs_species_filt_long[, "motif_tf"] == x , "gene_name"], # genes in cat
                  lengthData_currDE$gene_name[lengthData_currDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
}))


tfEnrich_neuronDE_hpk1$category
tfEnrich_neuronDE_sigTFs <- c("scrt-1", "klf-3", "klf-1", "sptf-3",
                              "cey-1", "unc-3", "daf-19", "svh-5", "daf-12", "nhr-48", "nhr-8", "ceh-22")

tfEnrich_neuronDE_sigTFs[tfEnrich_neuronDE_sigTFs %in% neuron_unique_enriched] # only daf-19 is enriched or unique to neurons
tfEnrich_neuronDE_sigTFs[tfEnrich_neuronDE_sigTFs %in% rownames(cengen_sum_level2)] # 9 of 13 have neuron expression

DE_hpk1[ DE_hpk1$GeneSymbol %in% tfEnrich_neuronDE_sigTFs , ]

tfEnrich_neuronDE_sigTFs %in% aging_atlas_cluster_shared_down
tfEnrich_neuronDE_sigTFs %in% aging_atlas_cluster_shared_up
genage_simple[genage_simple$GeneSymbol %in% tfEnrich_neuronDE_sigTFs , ]

write.xlsx(tfEnrich_neuronDE_hpk1, file = "tfEnrich_neuronDE_hpk1.xlsx")



heatmapData <- cengen_sum_level2[ c(tfEnrich_neuronDE_sigTFs[ tfEnrich_neuronDE_sigTFs %in% rownames(cengen_sum_level2) ], "hpk-1") , ]

tmpColData <- connectome_nodetable_merge$Neuron.Type[ match(colnames(cengen_sum_level2), connectome_nodetable_merge$CengenName)]
tmpColData[tmpColData == "Sensory"] <- "sienna3"
tmpColData[tmpColData == "Interneuron"] <- "turquoise3"
tmpColData[tmpColData == "Motorneuron"] <- "ivory2"
tmpColData[tmpColData == "Endorgan"] <- "seashell4"
tmpColData[is.na(tmpColData)] <- "seashell4"
tmpColData <- as.matrix(tmpColData)
colnames(tmpColData) <- "neuron_type"

pdf("Heatmap- cengen- TFs enriched for neuron-associated hpk-1 DE genes.pdf", width = 24, height = 8)
heatmap.3(x = heatmapData, 
          ColSideColors = tmpColData, 
          # RowSideColors = t(tmpRowData),
          dendrogram = "both",
          #scale = "row",
          col = viridis(option = "turbo", n = 64), #scale = "column",
          cexCol = 1, cexRow = 1.2, margins = c(9,5), keysize = 0.7, KeyValueName = "Cengen expression level")
graphics.off()





#
# Aging transcriptional profiles and hpk-1 null DE genes
#
#

#
# Wang et al 2022 EMBO
# 
#  It's  unclear if the foldchanges in the file are log transformed or not.
#  Actually the fold-changes are fold-ratios, so downregulated genes have foldchanges between 0 and 1,
#   also indicating these are not log transformed values. So the "FC < 0.67" for down-regulated is also 1.5.
#  That would also mean that their FC threshold on a log2 scale was just 0.58.
#  The methods don't mention p-value correction, or a specific approach that was used for DE analysis.
#   They also may have attempted to change the distribution from negative binomial to normal??
#     It mentions DE analysis was performed after transforming data with Box-Cox, which is intended to
#       induce "normality". So that's not a typical analysis workflow by any means.
#       The Box-Cox transform was applied by some people a few years ago- there's a paper on the application-
#       but that's not the best accepted workflow now by any means.
#   Will go with it for now, and possibly just filter this more stringently.
#
wang_aging <- list()
wang_aging$neuron <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wang_2022_bulk_aging",
                    "embj2021109633-sup-0003-datasetev1.xlsx"), sheet = "Neuron")
wang_aging$intestine <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wang_2022_bulk_aging",
                                         "embj2021109633-sup-0003-datasetev1.xlsx"), sheet = "Intestine")
wang_aging$bwm <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wang_2022_bulk_aging",
                                         "embj2021109633-sup-0003-datasetev1.xlsx"), sheet = "Body Wall Muscle")
wang_aging$hypodermis <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wang_2022_bulk_aging",
                                         "embj2021109633-sup-0003-datasetev1.xlsx"), sheet = "Hypodermis")
wang_aging$coelomocyte <- read.xlsx(file.path("E:/public_data_annotation/c_elegans_other/wang_2022_bulk_aging",
                                         "embj2021109633-sup-0003-datasetev1.xlsx"), sheet = "Coelomocyte")

# first, some annotation cleanup... some of the genes here are "dead"
# looks like they used both the sequence and the annotation from WB235 which is... old
wang_aging <- lapply(wang_aging, function(x){
  x = x[ x$Wormbase.ID %in% wb_gene_info$WBGeneID , ]
})

# check foldchanges 
plot(density(wang_aging$neuron$foldchange))
plot(density(wang_aging$intestine$foldchange))
plot(density(wang_aging$bwm$foldchange))
plot(density(wang_aging$hypodermis$foldchange))
plot(density(wang_aging$coelomocyte$foldchange))

# check p-values
plot(density(wang_aging$neuron$pvalue))
plot(density(wang_aging$intestine$pvalue))
plot(density(wang_aging$bwm$pvalue))
plot(density(wang_aging$hypodermis$pvalue))
plot(density(wang_aging$coelomocyte$pvalue))

# how do the p-values look after correction?
plot(density(p.adjust(wang_aging$neuron$pvalue, method = "fdr")))
sum(wang_aging$neuron$pvalue < 0.05)
sum(p.adjust(wang_aging$neuron$pvalue, method = "fdr") < 0.05)

# transform fold-ratio to fold change
wang_aging <- lapply(wang_aging, function(x){
  fc = x$foldchange
  fc[fc < 1] = -1*(1/fc[fc < 1])
  x$foldchange = fc
  return(x)
})

plot(density(wang_aging$neuron$foldchange))
plot(density(wang_aging$intestine$foldchange))
plot(density(wang_aging$bwm$foldchange))
plot(density(wang_aging$hypodermis$foldchange))
plot(density(wang_aging$coelomocyte$foldchange))
# yeah, those are not log-transformed foldchanges. Anyway, moving on.

# "high threshold" differential expression
#  criteria I've been using for our own data
wang_aging_highThreshDE <- lapply(wang_aging, function(x){
  x <- x[ p.adjust(x$pvalue, method = "fdr") < 0.05 , ]
  x <- x[ abs(x$foldchange) >= 2 , ]
  return(x)
})

# "low threshold" differential expression-
#  same criteria they used
wang_aging_lowThreshDE <- lapply(wang_aging, function(x){
    x <- x[ x$pvalue < 0.05 , ]
    x <- x[ abs(x$foldchange) > 1.5 , ]
    return(x)
})



# number of DE genes
cbind(names(wang_aging), Reduce(f = rbind, x = lapply(wang_aging_highThreshDE, function(x){
  data.frame("upregulated" = sum(x$foldchange > 0),
             "downregulated" = sum(x$foldchange < 0),
             "total" = nrow(x))
})))

cbind(names(wang_aging), Reduce(f = rbind, x = lapply(wang_aging_lowThreshDE, function(x){
  data.frame("upregulated" = sum(x$foldchange > 0),
             "downregulated" = sum(x$foldchange < 0),
             "total" = nrow(x))
})))

### ok so now overlap with hpk-1 DE
#  note that these foldchanges are D8 vs D1, so upregulated = aging-associated, presumably.

tissue_res_ht <- Reduce(f = cbind, x = lapply(names(wang_aging), function(currTissue){
  sapply(DE_hpk1$WBGeneID, function(currGene){
    if(currGene %in% wang_aging_highThreshDE[[currTissue]]$Wormbase.ID[wang_aging_highThreshDE[[currTissue]]$foldchange < 0]){
      return("down_D8_HT")
    }else if(currGene %in% wang_aging_highThreshDE[[currTissue]]$Wormbase.ID[wang_aging_highThreshDE[[currTissue]]$foldchange > 0]){
      return("up_D8_HT")
    }else{
      NA
    }
  })
}))
colnames(tissue_res_ht) <- paste0(names(wang_aging), "_D8_vs_D1")

tissue_res_lt <- Reduce(f = cbind, x = lapply(names(wang_aging), function(currTissue){
  sapply(DE_hpk1$WBGeneID, function(currGene){
    if(currGene %in% wang_aging_lowThreshDE[[currTissue]]$Wormbase.ID[wang_aging_lowThreshDE[[currTissue]]$foldchange < 0]){
      return("down_D8_LT")
    }else if(currGene %in% wang_aging_lowThreshDE[[currTissue]]$Wormbase.ID[wang_aging_lowThreshDE[[currTissue]]$foldchange > 0]){
      return("up_D8_LT")
    }else{
      NA
    }
  })
}))
colnames(tissue_res_lt) <- paste0(names(wang_aging), "_D8_vs_D1")

identical(DE_hpk1$WBGeneID, rownames(tissue_res_lt)) # should be true

DE_hpk1_wang_aging <- cbind(DE_hpk1, tissue_res_lt)

write.xlsx(x = DE_hpk1_wang_aging, "hpk-1 US vs N2 US DE and wang2022 D8 vs D1 results.xlsx")

# TF prediction on those concordant upregulated ones?

# merge with length info
lengthData_currDE <- lengthData
lengthData_currDE$GOI <- lengthData$WBGeneID %in% 
  DE_hpk1_wang_aging$WBGeneID[(DE_hpk1_wang_aging$DESeq2_log2FC > 0) & grepl("up", DE_hpk1_wang_aging$neuron_D8_vs_D1)]
sum(lengthData_currDE$GOI) # genes of interest


pwf=nullp(lengthData_currDE$GOI, bias.data = lengthData_currDE$length)
rownames(pwf) <- lengthData_currDE$WBGeneID
wang_neuron_up_olap <- goseq(pwf, gene2cat = sets_tfs_species_filt_long[, c("WBGeneID", "motif_tf") ])
wang_neuron_up_olap[,"over_represented_pvalue"] <- p.adjust(wang_neuron_up_olap$over_represented_pvalue, method="BH")
wang_neuron_up_olap[,"under_represented_pvalue"] <- p.adjust(wang_neuron_up_olap$under_represented_pvalue, method="BH")

wang_neuron_up_olap <- wang_neuron_up_olap[ wang_neuron_up_olap$over_represented_pvalue < 0.05 , ]

wang_neuron_up_olap$intersect_WBGeneID <- unlist(lapply(wang_neuron_up_olap$category, function(x){
  paste(intersect(sets_tfs_species_filt_long[ sets_tfs_species_filt_long[, "motif_tf"] == x , "gene_name"], # genes in cat
                  lengthData_currDE$gene_name[lengthData_currDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
}))


write.xlsx(x = wang_neuron_up_olap, file = "tf_enrich_hpk1US_up_Wang_neuron_aging_up_overlap.xlsx")


# what do we get with the ones that are up in our dataset and down in D8 neurons?
#  Is there a different set of TFs, or are we effectively enriching for TFs that are associated with
#  genes that are strongly enriched in neurons?
lengthData_currDE <- lengthData
lengthData_currDE$GOI <- lengthData$WBGeneID %in% 
  DE_hpk1_wang_aging$WBGeneID[(DE_hpk1_wang_aging$DESeq2_log2FC > 0) & grepl("down", DE_hpk1_wang_aging$neuron_D8_vs_D1)]
sum(lengthData_currDE$GOI) # genes of interest


pwf=nullp(lengthData_currDE$GOI, bias.data = lengthData_currDE$length)
rownames(pwf) <- lengthData_currDE$WBGeneID
wang_neuron_up_down_olap <- goseq(pwf, gene2cat = sets_tfs_species_filt_long[, c("WBGeneID", "motif_tf") ])
wang_neuron_up_down_olap[,"over_represented_pvalue"] <- p.adjust(wang_neuron_up_down_olap$over_represented_pvalue, method="BH")
wang_neuron_up_down_olap[,"under_represented_pvalue"] <- p.adjust(wang_neuron_up_down_olap$under_represented_pvalue, method="BH")

wang_neuron_up_down_olap <- wang_neuron_up_down_olap[ wang_neuron_up_down_olap$over_represented_pvalue < 0.05 , ]

wang_neuron_up_down_olap$intersect_WBGeneID <- unlist(lapply(wang_neuron_up_down_olap$category, function(x){
  paste(intersect(sets_tfs_species_filt_long[ sets_tfs_species_filt_long[, "motif_tf"] == x , "gene_name"], # genes in cat
                  lengthData_currDE$gene_name[lengthData_currDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
}))


#
# See if maybe GeneOverlap can help present the intersections between hpk-1 null and the tissue-specific aging data from wang et al
#  --> note that I had looked at two different thresholds for their data- "low", the same they used for the paper (non-adj p < 0.05, FC > 1.5), 
#         and "high", the same that I've been using for our projects (p adj < 0.05, |FC| >= 2 )
#
library(GeneOverlap)

GOM_hpk1_wang_lt <- newGOM(gsetA = list("hpk-1 null D2 up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                                     "hpk-1 null D2 down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]),
                        gsetB = list("neurons up" = with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange > 0]),
                                     "neurons down" = with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange < 0]),
                                     "intestine up" = with(wang_aging_lowThreshDE$intestine, Public.Name[foldchange > 0]),
                                     "intestine down" = with(wang_aging_lowThreshDE$intestine, Public.Name[foldchange < 0]),
                                     "BWM up" = with(wang_aging_lowThreshDE$bwm, Public.Name[foldchange > 0]),
                                     "BWM down" = with(wang_aging_lowThreshDE$bwm, Public.Name[foldchange < 0]),
                                     "hyp up" = with(wang_aging_lowThreshDE$hypodermis, Public.Name[foldchange > 0]),
                                     "hyp down" = with(wang_aging_lowThreshDE$hypodermis, Public.Name[foldchange < 0]),
                                     "coel up" = with(wang_aging_lowThreshDE$coelomocyte, Public.Name[foldchange > 0]),
                                     "coel down" = with(wang_aging_lowThreshDE$coelomocyte, Public.Name[foldchange < 0])),
                        genome.size =   length(unique(c(wang_aging$neuron$Public.Name, wang_aging$intestine$Public.Name,
                                                        wang_aging$bwm$Public.Name, wang_aging$hypodermis$Public.Name, 
                                                        wang_aging$coelomocyte$Public.Name,
                                                        DE_hpk1.all$GeneSymbol))))
pdf("GeneOverlapMatrix- Wang et al D8 vs D1 tissue DE LT and hpk-1 null DE.pdf", width = 20, height = 8)
drawHeatmap(GOM_hpk1_wang_lt, what = "Jaccard", adj.p = TRUE)
graphics.off()

#
# How large are the most significant overlaps?
#

# hpk-1 upregulated and wang neurons D8 vs D1
intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
          with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange < 0])) # 416 genes (out of 1898)

intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
          with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange > 0])) # 168 genes (out of 1821)

# hpk-1 downregulated and wang neurons D8 vs D1
intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],
          with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange < 0])) # 17 genes (out of 1898)

intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],
          with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange > 0])) # 55 genes (out of 1821)

# hpk-1 down and wang intestine down
intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],
          with(wang_aging_lowThreshDE$intestine, Public.Name[foldchange < 0])) # 35 genes (out of 1358)


GOM_hpk1_wang_ht <- newGOM(gsetA = list("hpk-1 null D2 up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                                        "hpk-1 null D2 down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]),
                           gsetB = list("neurons up" = with(wang_aging_highThreshDE$neuron, Public.Name[foldchange > 0]),
                                        "neurons down" = with(wang_aging_highThreshDE$neuron, Public.Name[foldchange < 0]),
                                        "intestine up" = with(wang_aging_highThreshDE$intestine, Public.Name[foldchange > 0]),
                                        "intestine down" = with(wang_aging_highThreshDE$intestine, Public.Name[foldchange < 0]),
                                        "BWM up" = with(wang_aging_highThreshDE$bwm, Public.Name[foldchange > 0]),
                                        "BWM down" = with(wang_aging_highThreshDE$bwm, Public.Name[foldchange < 0]),
                                        "hyp up" = with(wang_aging_highThreshDE$hypodermis, Public.Name[foldchange > 0]),
                                        "hyp down" = with(wang_aging_highThreshDE$hypodermis, Public.Name[foldchange < 0]),
                                        "coel up" = with(wang_aging_highThreshDE$coelomocyte, Public.Name[foldchange > 0]),
                                        "coel down" = with(wang_aging_highThreshDE$coelomocyte, Public.Name[foldchange < 0])),
                           genome.size =   length(unique(c(wang_aging$neuron$Public.Name, wang_aging$intestine$Public.Name,
                                                           wang_aging$bwm$Public.Name, wang_aging$hypodermis$Public.Name, 
                                                           wang_aging$coelomocyte$Public.Name,
                                                           DE_hpk1.all$GeneSymbol))))
pdf("GeneOverlapMatrix- Wang et al D8 vs D1 tissue DE HT and hpk-1 null DE.pdf", width = 20, height = 8)
drawHeatmap(GOM_hpk1_wang_ht, what = "Jaccard", adj.p = TRUE)
graphics.off()


#
# hpk-1 and aging-associated genes
#

table(aging_genes_merged$wb_annotation)
table(aging_genes_merged$genage_gerogene_cats)

GOM_hpk1_gerogenes <- newGOM(gsetA = list("hpk-1 null D2 up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                                        "hpk-1 null D2 down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]),
                           gsetB = list("Wormbase- extended lifespan" = 
                                          as.vector(na.omit(aging_genes_merged$gene_name[aging_genes_merged$wb_annotation == "extended lifespan"])),
                                        "Wormbase- shortened lifespan" = 
                                          as.vector(na.omit(aging_genes_merged$gene_name[aging_genes_merged$wb_annotation == "shortened lifespan"])),
                                        "GenAge- longevity gene" =
                                          as.vector(na.omit(aging_genes_merged$gene_name[aging_genes_merged$genage_gerogene_cats == "longevity_gene"])),
                                        "GenAge- progeric gene" =
                                          as.vector(na.omit(aging_genes_merged$gene_name[aging_genes_merged$genage_gerogene_cats == "progeric_gene"])),
                                        "any lifespan phenotype from all sources" = 
                                          aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno]),
                           genome.size =   20000)

# NOTHING SIGNIFICANT
pdf("GeneOverlapMatrix-hpk-1 null DE and gerogenes.pdf", width = 20, height = 20)
drawHeatmap(GOM_hpk1_gerogenes, what = "Jaccard", adj.p = FALSE, cutoff = 0.5)
graphics.off()

pdf("venn- hpk-1 null and lifespan genes.pdf", width = 10, height = 8)
grid.draw(venn.diagram(
  x=
    list("hpk-1 null up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "hpk-1 null down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],
         "Genes with any lifespan phenotype" = aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno]),
  filename = NULL,
  fill = c("tomato", "cornflowerblue", "cornsilk"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1.5,
  cat.col = c("tomato3", "darkblue", "cornsilk3"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.05)
)
graphics.off()


#
# While we're here, since it might be useful for the DR paper:
# * TF enrichment in the tissue-defined aging genes, since their approach was different
# * Write out the significant results after I filtered and such
#
#  This will be, as usual, split into up and down regulated sets

agingTissue_TFenrich <- lapply(wang_aging_lowThreshDE, function(currDat){
  
  currDat.up = currDat[ currDat$foldchange > 0 , ]
  currDat.down = currDat[ currDat$foldchange < 0 , ]
  
  lengthData_upDE <- lengthData # Assuming lengthData was defined earlier
  lengthData_downDE <- lengthData
  
  lengthData_upDE$GOI <- lengthData$WBGeneID %in% currDat.up$Wormbase.ID
  lengthData_downDE$GOI <- lengthData$WBGeneID %in% currDat.down$Wormbase.ID
  
  # up
  pwf=with(data = lengthData_upDE, nullp(GOI, bias.data = length))
  rownames(pwf) <- lengthData_upDE$WBGeneID
  up_enrich <- goseq(pwf, gene2cat = sets_tfs_species_filt_long[, c("WBGeneID", "motif_tf") ])
  up_enrich[,"over_represented_pvalue"] <- p.adjust(up_enrich$over_represented_pvalue, method="BH")
  up_enrich[,"under_represented_pvalue"] <- p.adjust(up_enrich$under_represented_pvalue, method="BH")
  
  up_enrich <- up_enrich[ up_enrich$over_represented_pvalue < 0.05 , ]
  
  up_enrich$intersect_genes <- unlist(lapply(up_enrich$category, function(x){
    paste(intersect(sets_tfs_species_filt_long$gene_name[ sets_tfs_species_filt_long$motif_tf == x ], # genes in cat
                    lengthData_upDE$gene_name[lengthData_upDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
  }))
  
  # down
  pwf=with(data = lengthData_downDE, nullp(GOI, bias.data = length))
  rownames(pwf) <- lengthData_downDE$WBGeneID
  down_enrich <- goseq(pwf, gene2cat = sets_tfs_species_filt_long[, c("WBGeneID", "motif_tf") ])
  down_enrich[,"over_represented_pvalue"] <- p.adjust(down_enrich$over_represented_pvalue, method="BH")
  down_enrich[,"under_represented_pvalue"] <- p.adjust(down_enrich$under_represented_pvalue, method="BH")
  
  down_enrich <- down_enrich[ down_enrich$over_represented_pvalue < 0.05 , ]
  
  down_enrich$intersect_genes <- unlist(lapply(down_enrich$category, function(x){
    paste(intersect(sets_tfs_species_filt_long$gene_name[ sets_tfs_species_filt_long$motif_tf == x ], # genes in cat
                    lengthData_downDE$gene_name[lengthData_downDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
  }))
  
  if((nrow(up_enrich) == 0) && (nrow(down_enrich) == 0)){
    return(NULL)
  }else if((nrow(up_enrich) == 0)){
    return(data.frame("Direction" = "downregulated", down_enrich))
  }else if((nrow(down_enrich) == 0)){
    return(data.frame("Direction" = "upregulated", up_enrich))
  }else{
    return(rbind(data.frame("Direction" = "upregulated", up_enrich),
                 data.frame("Direction" = "downregulated", down_enrich)))
  }

})

#
# TF enrichment, on the set without species-based filtering
#
agingTissue_TFenrich_allCeTgts <- lapply(wang_aging_lowThreshDE, function(currDat){
  
  currDat.up = currDat[ currDat$foldchange > 0 , ]
  currDat.down = currDat[ currDat$foldchange < 0 , ]
  
  lengthData_upDE <- lengthData # Assuming lengthData was defined earlier
  lengthData_downDE <- lengthData
  
  lengthData_upDE$GOI <- lengthData$WBGeneID %in% currDat.up$Wormbase.ID
  lengthData_downDE$GOI <- lengthData$WBGeneID %in% currDat.down$Wormbase.ID
  
  # up
  pwf=with(data = lengthData_upDE, nullp(GOI, bias.data = length))
  rownames(pwf) <- lengthData_upDE$WBGeneID
  up_enrich <- goseq(pwf, gene2cat = sets_tfs_all_ce_long[, c("WBGeneID", "motif") ])
  up_enrich[,"over_represented_pvalue"] <- p.adjust(up_enrich$over_represented_pvalue, method="BH")
  up_enrich[,"under_represented_pvalue"] <- p.adjust(up_enrich$under_represented_pvalue, method="BH")
  
  up_enrich <- up_enrich[ up_enrich$over_represented_pvalue < 0.05 , ]
  
  up_enrich$intersect_genes <- unlist(lapply(up_enrich$category, function(x){
    paste(intersect(sets_tfs_species_filt_long$gene_name[ sets_tfs_species_filt_long$motif == x ], # genes in cat
                    lengthData_upDE$gene_name[lengthData_upDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
  }))
  
  # down
  pwf=with(data = lengthData_downDE, nullp(GOI, bias.data = length))
  rownames(pwf) <- lengthData_downDE$WBGeneID
  down_enrich <- goseq(pwf, gene2cat = sets_tfs_all_ce_long[, c("WBGeneID", "motif") ])
  down_enrich[,"over_represented_pvalue"] <- p.adjust(down_enrich$over_represented_pvalue, method="BH")
  down_enrich[,"under_represented_pvalue"] <- p.adjust(down_enrich$under_represented_pvalue, method="BH")
  
  down_enrich <- down_enrich[ down_enrich$over_represented_pvalue < 0.05 , ]
  
  down_enrich$intersect_genes <- unlist(lapply(down_enrich$category, function(x){
    paste(intersect(sets_tfs_species_filt_long$gene_name[ sets_tfs_species_filt_long$motif == x ], # genes in cat
                    lengthData_downDE$gene_name[lengthData_downDE$GOI]), sep = ", ", collapse = ", ") # genes sig in query set
  }))
  
  if((nrow(up_enrich) == 0) && (nrow(down_enrich) == 0)){
    return(NULL)
  }else if((nrow(up_enrich) == 0)){
    return(data.frame("Direction" = "downregulated", down_enrich))
  }else if((nrow(down_enrich) == 0)){
    return(data.frame("Direction" = "upregulated", up_enrich))
  }else{
    return(rbind(data.frame("Direction" = "upregulated", up_enrich),
                 data.frame("Direction" = "downregulated", down_enrich)))
  }
})

#
# looks ok, now write out the tissue aging/enrichment results
#  There's a lot of replication here because openxlsx uses a lot of pointers
#   so there's some strange behaviour with copying objects.
#
wb <- createWorkbook()
addWorksheet(wb, "neuron")
addWorksheet(wb, "intestine")
addWorksheet(wb, "BWM")
addWorksheet(wb, "hypodermis")
addWorksheet(wb, "coelomocyte")

wb1 <- createWorkbook()
addWorksheet(wb1, "neuron")
addWorksheet(wb1, "intestine")
addWorksheet(wb1, "BWM")
addWorksheet(wb1, "hypodermis")
addWorksheet(wb1, "coelomocyte")

wb2 <- createWorkbook()
addWorksheet(wb2, "neuron")
addWorksheet(wb2, "intestine")
addWorksheet(wb2, "BWM")
addWorksheet(wb2, "hypodermis")
addWorksheet(wb2, "coelomocyte")


# filtered/updated gene names for DE aging tissue genes
writeDataTable(wb, sheet = "neuron", x = wang_aging_lowThreshDE$neuron)
writeDataTable(wb, sheet = "intestine", x = wang_aging_lowThreshDE$intestine)
writeDataTable(wb, sheet = "BWM", x = wang_aging_lowThreshDE$bwm)
writeDataTable(wb, sheet = "hypodermis", x = wang_aging_lowThreshDE$hypodermis)
writeDataTable(wb, sheet = "coelomocyte", x = wang_aging_lowThreshDE$coelomocyte)

# TF enrichment, species-filtered
writeDataTable(wb1, sheet = "neuron", x = agingTissue_TFenrich$neuron)
writeDataTable(wb1, sheet = "intestine", x = agingTissue_TFenrich$intestine)
writeDataTable(wb1, sheet = "BWM", x = agingTissue_TFenrich$bwm)
writeDataTable(wb1, sheet = "hypodermis", x = agingTissue_TFenrich$hypodermis)
writeDataTable(wb1, sheet = "coelomocyte", x = agingTissue_TFenrich$coelomocyte)

# TF enrichment, all Ce predicted targets
writeDataTable(wb2, sheet = "neuron", x = agingTissue_TFenrich_allCeTgts$neuron)
writeDataTable(wb2, sheet = "intestine", x = agingTissue_TFenrich_allCeTgts$intestine)
writeDataTable(wb2, sheet = "BWM", x = agingTissue_TFenrich_allCeTgts$bwm)
writeDataTable(wb2, sheet = "hypodermis", x = agingTissue_TFenrich_allCeTgts$hypodermis)
writeDataTable(wb2, sheet = "coelomocyte", x = agingTissue_TFenrich_allCeTgts$coelomocyte)

# save files
saveWorkbook(wb, file = "Wang 2022 tissue aging DE genes filtered.xlsx")
saveWorkbook(wb1, file = "Wang 2022 tissue aging Adam TF enrichment (species filtered targets).xlsx")
saveWorkbook(wb2, file = "Wang 2022 tissue aging Adam TF enrichment (all Ce targets).xlsx")

rm(wb, wb1, wb2)


#
# Now we'll do something similar with the aging atlas dataset:
#  What is the overlap in aging-associated DEGs in the hpk-1 dataset?
# Andy suggested leaving out things that have mixed aging changes (both up down across different tissues)
#  but we won't throw those out immediately.
aging_atlas_deg_clusters



lapply(DE_hpk1$GeneSymbol, function(gene_name){
  
  aging_atlas_degs
  
  data.frame()
  
})

aging_atlas_degs


aging_atlas_hpk1_DE_overlap_all <-  t(sapply(DE_hpk1$GeneSymbol, function(currGene){
  if(currGene %in% aging_atlas_degs_long$gene_name){
    
    with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
         if((sum(logFoldChange > 0) > 0) & (sum(logFoldChange < 0) > 0)){ # gene is both up and down in different places
           return(data.frame(up = paste(Clusters.explanatory.name[logFoldChange > 0], sep = "; ", collapse = "; "),
                             down = paste(Clusters.explanatory.name[logFoldChange < 0], sep = "; ", collapse = "; ")))
         }else if(sum(logFoldChange > 0) > 0){ # gene is just upregulated
           return(data.frame(up = paste(Clusters.explanatory.name[logFoldChange > 0], sep = "; ", collapse = "; "), down = NA))
         }else if(sum(logFoldChange < 0) > 0){ # gene is just downregulated
           return(data.frame(up = NA, down = paste(Clusters.explanatory.name[logFoldChange < 0], sep = "; ", collapse = "; ")))
         }else{stop("something unexpected occurred")}
    )
  }else{
    return(data.frame(up = NA, down = NA))
  }
}))

aging_atlas_hpk1_DE_overlap_all <- merge(aging_atlas_hpk1_DE_overlap_all, DE_hpk1[,c("GeneSymbol", "DESeq2_log2FC")],
                                         by.x = 0, by.y = "GeneSymbol")
colnames(aging_atlas_hpk1_DE_overlap_all)[1] <- "gene"

aging_atlas_hpk1_DE_overlap_all$num_clusters_aging_up <- sapply(aging_atlas_hpk1_DE_overlap_all$gene, function(currGene){
  if(currGene %in% aging_atlas_degs_long$gene_name){
    with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
         sum(logFoldChange > 0))
  }else{ return(0)}
})

aging_atlas_hpk1_DE_overlap_all$num_clusters_aging_down <- sapply(aging_atlas_hpk1_DE_overlap_all$gene, function(currGene){
  if(currGene %in% aging_atlas_degs_long$gene_name){
    with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
         sum(logFoldChange < 0))
  }else{ return(0)}
})

# number of cases with upregulation in hpk-1 null
with(aging_atlas_hpk1_DE_overlap_all[ aging_atlas_hpk1_DE_overlap_all$DESeq2_log2FC > 0 , ], c("all_up" = sum((num_clusters_aging_up > 0)&
                                                                                                                (num_clusters_aging_down) == 0),
                                                                                               "all_down" = sum((num_clusters_aging_up == 0)&
                                                                                                                  (num_clusters_aging_down) > 0),
                                                                                               "mixed" = sum((num_clusters_aging_up > 0)&
                                                                                                               (num_clusters_aging_down) > 0),
                                                                                               "no_aging_de" = sum((num_clusters_aging_up == 0)&
                                                                                                                  (num_clusters_aging_down) == 0)))
with(aging_atlas_hpk1_DE_overlap_all[ aging_atlas_hpk1_DE_overlap_all$DESeq2_log2FC < 0 , ], c("all_up" = sum((num_clusters_aging_up > 0)&
                                                                                                                (num_clusters_aging_down) == 0),
                                                                                               "all_down" = sum((num_clusters_aging_up == 0)&
                                                                                                                  (num_clusters_aging_down) > 0),
                                                                                               "mixed" = sum((num_clusters_aging_up > 0)&
                                                                                                               (num_clusters_aging_down) > 0),
                                                                                               "no_aging_de" = sum((num_clusters_aging_up == 0)&
                                                                                                                     (num_clusters_aging_down) == 0)))
# for reference, what's the size of each of these categories for all the aging atlas genes?
aging_atlas_DE_up_down_count_all <- data.frame(
  num_clusters_aging_up = sapply(as.character(unique(unlist(lapply(aging_atlas_degs, function(x){x[,1]})))), 
                                 function(currGene){
                                   if(currGene %in% aging_atlas_degs_long$gene_name){
                                     with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
                                          sum(logFoldChange > 0))
                                   }else{ return(0)}
                                 }),
  num_clusters_aging_down = sapply(as.character(unique(unlist(lapply(aging_atlas_degs, function(x){x[,1]})))), 
                                   function(currGene){
                                     if(currGene %in% aging_atlas_degs_long$gene_name){
                                       with(data = aging_atlas_degs_long[aging_atlas_degs_long$gene_name == currGene , ],
                                            sum(logFoldChange < 0))
                                     }else{ return(0)}
                                   })
  )


aging_atlas_DE_up_down_count_all$count_sum_cat <- sapply(1:nrow(aging_atlas_DE_up_down_count_all), function(i){
  x = aging_atlas_DE_up_down_count_all
  if((x[i,1] > 0) &
     (x[i,2] == 0)){
    return("all_up")
  }else if((x[i,1] == 0) &
           (x[i,2] > 0)){
    return("all_down")
  }else if((x[i,1] > 0) &
           (x[i,2] > 0)){
    return("mixed")
  }else if((x[i,1] == 0) &
           (x[i,2] == 0)){
    return("no_aging_DE")
  }
})

table(aging_atlas_DE_up_down_count_all$count_sum_cat)


#
# GeneOverlap with hpk-1 DE genes and the aging atlas
#
GOM_hpk1_atlas <- newGOM(gsetA = list("hpk-1 null D2 up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                                        "hpk-1 null D2 down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]),
                           gsetB = 
                           list("only upregulated" = 
                                  rownames(aging_atlas_DE_up_down_count_all)[aging_atlas_DE_up_down_count_all$count_sum_cat == "all_up"],
                                "only downregulated" = 
                                  rownames(aging_atlas_DE_up_down_count_all)[aging_atlas_DE_up_down_count_all$count_sum_cat == "all_down"],
                                "mixed up down" = 
                                  rownames(aging_atlas_DE_up_down_count_all)[aging_atlas_DE_up_down_count_all$count_sum_cat == "mixed"],
                                "no aging DE" = 
                                  rownames(aging_atlas_DE_up_down_count_all)[aging_atlas_DE_up_down_count_all$count_sum_cat == "no_aging_DE"]),
                           genome.size =   length(unique(c(unlist(lapply(aging_atlas_cluster_mean_counts, function(x){rownames(x)})),
                                                           DE_hpk1.all$GeneSymbol))))
pdf("GeneOverlapMatrix- aging atlas cluster direction summary and hpk-1 null DE.pdf", width = 20, height = 8)
drawHeatmap(GOM_hpk1_atlas, what = "Jaccard", adj.p = TRUE)
graphics.off()

# let's try to do this more specifically, by tissue type, as specified by "Tissue.type.(higher.degree)"
#  neuron, muscle, epithelium, glia, hypodermis, intestine, gland, valve, mesoderm, seam, various | coelomocyte (these are both coelomocytes!)
table(aging_atlas_clusters_all$`Tissue.type.(higher.degree)`)

aging_atlas_degs_long

aging_atlas_degs_long[aging_atlas_degs_long$cluster_shortID %in%
                        aging_atlas_clusters_all$shortID[aging_atlas_clusters_all$`Tissue.type.(higher.degree)` == "mesoderm"] , ]

#
# without specifying a minimum # of DE clusters...
atlas_tissue_type_DEbroad <- lapply(
  c("neuron", "muscle", "epithelium", "glia",
    "hypodermis", "intestine", "gland", "valve",
    "mesoderm", "seam", "various|coelomocyte"), function(tissue){
      # get clusters for current tissue
      currClusts = aging_atlas_clusters_all$shortID[grepl(tissue, aging_atlas_clusters_all$`Tissue.type.(higher.degree)`)]
      # check if any of them are in the DE set
      if(any(currClusts %in% unique(aging_atlas_degs_long$cluster_shortID))){
        # for each unique DE gene across the relevant clusters, determine if it's up, down, or mixed across the clusters
        if(length(currClusts) == 1){
          # can only be all up or all down
          return(setNames(ifelse(aging_atlas_degs_long$logFoldChange[ aging_atlas_degs_long$cluster_shortID == currClusts] > 0,
                                 "all_up", "all_down"), aging_atlas_degs_long$gene_name[ aging_atlas_degs_long$cluster_shortID == currClusts]))
        }else if(length(currClusts) > 1){
          sapply(unique(aging_atlas_degs_long$gene_name[aging_atlas_degs_long$cluster_shortID %in% currClusts]),
                 function(gene){
                   currFCs = aging_atlas_degs_long$logFoldChange[ (aging_atlas_degs_long$cluster_shortID %in% currClusts) &
                                                                    (aging_atlas_degs_long$gene_name %in% gene)]
                   if(sum(currFCs > 0) == length(currFCs)){
                     return("all_up")
                   }else if(sum(currFCs < 0) == length(currFCs)){
                     return("all_down")
                   }else{
                     return("mixed")
                   }
                 })
        }
      }
    })
names(atlas_tissue_type_DEbroad) <- c("neuron", "muscle", "epithelium", "glia",
                                      "hypodermis", "intestine", "gland", "valve",
                                      "mesoderm", "seam", "various|coelomocyte")


GOM_hpk1_atlas_tissueSum <- newGOM(gsetB = list("hpk-1 null D2 up" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
                                      "hpk-1 null D2 down" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]),
                         gsetA = 
                           list("neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"],
                                "neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"],
                                "neuron mixed" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "mixed"],
                                
                                "muscle upregulated" = names(atlas_tissue_type_DEbroad$muscle)[atlas_tissue_type_DEbroad$muscle == "all_up"],
                                "muscle downregulated" = names(atlas_tissue_type_DEbroad$muscle)[atlas_tissue_type_DEbroad$muscle == "all_down"],
                                "muscle mixed" = names(atlas_tissue_type_DEbroad$muscle)[atlas_tissue_type_DEbroad$muscle == "mixed"],
                                
                                "epithelium upregulated" = names(atlas_tissue_type_DEbroad$epithelium)[atlas_tissue_type_DEbroad$epithelium == "all_up"],
                                "epithelium downregulated" = names(atlas_tissue_type_DEbroad$epithelium)[atlas_tissue_type_DEbroad$epithelium == "all_down"],
                                "epithelium mixed" = names(atlas_tissue_type_DEbroad$epithelium)[atlas_tissue_type_DEbroad$epithelium == "mixed"],
                                
                                "glia upregulated" = names(atlas_tissue_type_DEbroad$glia)[atlas_tissue_type_DEbroad$glia == "all_up"],
                                "glia downregulated" = names(atlas_tissue_type_DEbroad$glia)[atlas_tissue_type_DEbroad$glia == "all_down"],
                                "glia mixed" = names(atlas_tissue_type_DEbroad$glia)[atlas_tissue_type_DEbroad$glia == "mixed"],
                                
                                "hypodermis upregulated" = names(atlas_tissue_type_DEbroad$hypodermis)[atlas_tissue_type_DEbroad$hypodermis == "all_up"],
                                "hypodermis downregulated" = names(atlas_tissue_type_DEbroad$hypodermis)[atlas_tissue_type_DEbroad$hypodermis == "all_down"],
                                "hypodermis mixed" = names(atlas_tissue_type_DEbroad$hypodermis)[atlas_tissue_type_DEbroad$hypodermis == "mixed"],
                                
                                "intestine upregulated" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "all_up"],
                                "intestine downregulated" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "all_down"],
                                "intestine mixed" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "mixed"],
                                
                                "gland upregulated" = names(atlas_tissue_type_DEbroad$gland)[atlas_tissue_type_DEbroad$gland == "all_up"],
                                "gland downregulated" = names(atlas_tissue_type_DEbroad$gland)[atlas_tissue_type_DEbroad$gland == "all_down"],
                                "gland mixed" = names(atlas_tissue_type_DEbroad$gland)[atlas_tissue_type_DEbroad$gland == "mixed"],
                                
                                "valve upregulated" = names(atlas_tissue_type_DEbroad$valve)[atlas_tissue_type_DEbroad$valve == "all_up"],
                                "valve downregulated" = names(atlas_tissue_type_DEbroad$valve)[atlas_tissue_type_DEbroad$valve == "all_down"],
                                "valve mixed" = names(atlas_tissue_type_DEbroad$valve)[atlas_tissue_type_DEbroad$valve == "mixed"],
                                
                                "mesoderm upregulated" = names(atlas_tissue_type_DEbroad$mesoderm)[atlas_tissue_type_DEbroad$mesoderm == "all_up"],
                                "mesoderm downregulated" = names(atlas_tissue_type_DEbroad$mesoderm)[atlas_tissue_type_DEbroad$mesoderm == "all_down"],
                                "mesoderm mixed" = names(atlas_tissue_type_DEbroad$mesoderm)[atlas_tissue_type_DEbroad$mesoderm == "mixed"],
                                
                                "seam upregulated" = names(atlas_tissue_type_DEbroad$seam)[atlas_tissue_type_DEbroad$seam == "all_up"],
                                "seam downregulated" = names(atlas_tissue_type_DEbroad$seam)[atlas_tissue_type_DEbroad$seam == "all_down"],
                                "seam mixed" = names(atlas_tissue_type_DEbroad$seam)[atlas_tissue_type_DEbroad$seam == "mixed"],
                                
                                "coelomocyte upregulated" = names(atlas_tissue_type_DEbroad$`various|coelomocyte`)[atlas_tissue_type_DEbroad$`various|coelomocyte` == "all_up"],
                                "coelomocyte downregulated" = names(atlas_tissue_type_DEbroad$`various|coelomocyte`)[atlas_tissue_type_DEbroad$`various|coelomocyte` == "all_down"],
                                "coelomocyte mixed" = names(atlas_tissue_type_DEbroad$`various|coelomocyte`)[atlas_tissue_type_DEbroad$`various|coelomocyte` == "mixed"]),
                         genome.size =   length(unique(c(unlist(lapply(aging_atlas_cluster_mean_counts, function(x){rownames(x)})),
                                                         DE_hpk1.all$GeneSymbol))))


pdf("GeneOverlapMatrix- aging atlas tissue type cluster direction summary and hpk-1 null DE.pdf", width = 20, height = 20)
drawHeatmap(GOM_hpk1_atlas_tissueSum, what = "Jaccard", adj.p = TRUE, cutoff = 0.05)
graphics.off()

#
# get intersect sizes
#

# hpk-1 up
unlist(lapply(list("neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"],
"neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"],
"neuron mixed" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "mixed"],
"intestine downregulated" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "all_down"]),
function(x){length(intersect(x, DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0]))}))

# hpk-1 down
unlist(lapply(list("neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"],
                   "neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"],
                   "neuron mixed" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "mixed"],
                   "intestine downregulated" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "all_down"]),
              function(x){length(intersect(x, DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0]))}))

#
# how well do the up/down genes overlap between the two aging datasets? for neurons?
#
# wang et al neurons up
unlist(lapply(list("neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"],
                   "neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"],
                   "neuron mixed" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "mixed"]),
              function(x){length(intersect(x, with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange > 0])))}))

# wang et al neurons down
unlist(lapply(list("neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"],
                   "neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"],
                   "neuron mixed" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "mixed"]),
              function(x){length(intersect(x, with(wang_aging_lowThreshDE$neuron, Public.Name[foldchange < 0])))}))

#
# overlap matrix- wang tissue-specific vs aging atlas single-cell
#

GOM_wang_atlas <- newGOM(gsetB = list("neurons up" = with(wang_aging_highThreshDE$neuron, Public.Name[foldchange > 0]),
                    "neurons down" = with(wang_aging_highThreshDE$neuron, Public.Name[foldchange < 0]),
                    "intestine up" = with(wang_aging_highThreshDE$intestine, Public.Name[foldchange > 0]),
                    "intestine down" = with(wang_aging_highThreshDE$intestine, Public.Name[foldchange < 0]),
                    "BWM up" = with(wang_aging_highThreshDE$bwm, Public.Name[foldchange > 0]),
                    "BWM down" = with(wang_aging_highThreshDE$bwm, Public.Name[foldchange < 0]),
                    "hyp up" = with(wang_aging_highThreshDE$hypodermis, Public.Name[foldchange > 0]),
                    "hyp down" = with(wang_aging_highThreshDE$hypodermis, Public.Name[foldchange < 0]),
                    "coel up" = with(wang_aging_highThreshDE$coelomocyte, Public.Name[foldchange > 0]),
                    "coel down" = with(wang_aging_highThreshDE$coelomocyte, Public.Name[foldchange < 0])),
       gsetA = 
         list("neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"],
              "neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"],
              "neuron mixed" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "mixed"],
              
              "muscle upregulated" = names(atlas_tissue_type_DEbroad$muscle)[atlas_tissue_type_DEbroad$muscle == "all_up"],
              "muscle downregulated" = names(atlas_tissue_type_DEbroad$muscle)[atlas_tissue_type_DEbroad$muscle == "all_down"],
              "muscle mixed" = names(atlas_tissue_type_DEbroad$muscle)[atlas_tissue_type_DEbroad$muscle == "mixed"],
              
              "epithelium upregulated" = names(atlas_tissue_type_DEbroad$epithelium)[atlas_tissue_type_DEbroad$epithelium == "all_up"],
              "epithelium downregulated" = names(atlas_tissue_type_DEbroad$epithelium)[atlas_tissue_type_DEbroad$epithelium == "all_down"],
              "epithelium mixed" = names(atlas_tissue_type_DEbroad$epithelium)[atlas_tissue_type_DEbroad$epithelium == "mixed"],
              
              "glia upregulated" = names(atlas_tissue_type_DEbroad$glia)[atlas_tissue_type_DEbroad$glia == "all_up"],
              "glia downregulated" = names(atlas_tissue_type_DEbroad$glia)[atlas_tissue_type_DEbroad$glia == "all_down"],
              "glia mixed" = names(atlas_tissue_type_DEbroad$glia)[atlas_tissue_type_DEbroad$glia == "mixed"],
              
              "hypodermis upregulated" = names(atlas_tissue_type_DEbroad$hypodermis)[atlas_tissue_type_DEbroad$hypodermis == "all_up"],
              "hypodermis downregulated" = names(atlas_tissue_type_DEbroad$hypodermis)[atlas_tissue_type_DEbroad$hypodermis == "all_down"],
              "hypodermis mixed" = names(atlas_tissue_type_DEbroad$hypodermis)[atlas_tissue_type_DEbroad$hypodermis == "mixed"],
              
              "intestine upregulated" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "all_up"],
              "intestine downregulated" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "all_down"],
              "intestine mixed" = names(atlas_tissue_type_DEbroad$intestine)[atlas_tissue_type_DEbroad$intestine == "mixed"],
              
              "gland upregulated" = names(atlas_tissue_type_DEbroad$gland)[atlas_tissue_type_DEbroad$gland == "all_up"],
              "gland downregulated" = names(atlas_tissue_type_DEbroad$gland)[atlas_tissue_type_DEbroad$gland == "all_down"],
              "gland mixed" = names(atlas_tissue_type_DEbroad$gland)[atlas_tissue_type_DEbroad$gland == "mixed"],
              
              "valve upregulated" = names(atlas_tissue_type_DEbroad$valve)[atlas_tissue_type_DEbroad$valve == "all_up"],
              "valve downregulated" = names(atlas_tissue_type_DEbroad$valve)[atlas_tissue_type_DEbroad$valve == "all_down"],
              "valve mixed" = names(atlas_tissue_type_DEbroad$valve)[atlas_tissue_type_DEbroad$valve == "mixed"],
              
              "mesoderm upregulated" = names(atlas_tissue_type_DEbroad$mesoderm)[atlas_tissue_type_DEbroad$mesoderm == "all_up"],
              "mesoderm downregulated" = names(atlas_tissue_type_DEbroad$mesoderm)[atlas_tissue_type_DEbroad$mesoderm == "all_down"],
              "mesoderm mixed" = names(atlas_tissue_type_DEbroad$mesoderm)[atlas_tissue_type_DEbroad$mesoderm == "mixed"],
              
              "seam upregulated" = names(atlas_tissue_type_DEbroad$seam)[atlas_tissue_type_DEbroad$seam == "all_up"],
              "seam downregulated" = names(atlas_tissue_type_DEbroad$seam)[atlas_tissue_type_DEbroad$seam == "all_down"],
              "seam mixed" = names(atlas_tissue_type_DEbroad$seam)[atlas_tissue_type_DEbroad$seam == "mixed"],
              
              "coelomocyte upregulated" = names(atlas_tissue_type_DEbroad$`various|coelomocyte`)[atlas_tissue_type_DEbroad$`various|coelomocyte` == "all_up"],
              "coelomocyte downregulated" = names(atlas_tissue_type_DEbroad$`various|coelomocyte`)[atlas_tissue_type_DEbroad$`various|coelomocyte` == "all_down"],
              "coelomocyte mixed" = names(atlas_tissue_type_DEbroad$`various|coelomocyte`)[atlas_tissue_type_DEbroad$`various|coelomocyte` == "mixed"]),
       genome.size =   length(unique(c(unlist(lapply(aging_atlas_cluster_mean_counts, function(x){rownames(x)})),
                                       unlist(lapply(wang_aging, function(x){x$Public.Name}))))))

pdf("GeneOverlapMatrix- aging atlas tissue type cluster direction summary and wang 2022 tissue specific aging.pdf", width = 20, height = 20)
drawHeatmap(GOM_wang_atlas, what = "Jaccard", adj.p = TRUE, cutoff = 0.05)
graphics.off()


#
# try to get the specific numbers for all of these overlap, to be added to the plot
# Andy suggested this as a thesis figure, could be included in the intro or discussion on limitations or future directions?
#


#
# ok so not considering the wang 2022 dataset for now-
# For the hpk-1 DE up genes:
#  For the up ones that are also up in neurons with age- what are the functions? Is this associated with progeria?
#  For the up ones that are actually down in neurons with age- is this part of a compensatory response? What are these genes?

# also need to look at how broad these changes are


# --- response matches direction with aging ---
hpk1_atlas_neurons_trend <- list()

# down
hpk1_atlas_neurons_trend$both_down <- intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0] ,
          names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"])

# up
hpk1_atlas_neurons_trend$both_up <- intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0] ,
          names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"])

# --- response in inverted direction from aging ---

# down
hpk1_atlas_neurons_trend$hpk1_down_atlas_up <- intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0] ,
          names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"])

# up
hpk1_atlas_neurons_trend$hpk1_up_atlas_down <- intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0] ,
          names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"])

unlist(lapply(hpk1_atlas_neurons_trend, length))

# overlap with lifespan genes?

intersect(hpk1_atlas_neurons_trend$both_up, aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno])
table(aging_genes_merged$wb_annotation[ aging_genes_merged$gene_name %in% hpk1_atlas_neurons_trend$both_up ])
table(aging_genes_merged$genage_gerogene_cats[ aging_genes_merged$gene_name %in% hpk1_atlas_neurons_trend$both_up ])
table(aging_genes_merged$inPGP[ aging_genes_merged$gene_name %in% hpk1_atlas_neurons_trend$both_up ])

intersect(hpk1_atlas_neurons_trend$both_down, aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno])

intersect(hpk1_atlas_neurons_trend$hpk1_down_atlas_up, aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno])

intersect(hpk1_atlas_neurons_trend$hpk1_up_atlas_down, aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno])
table(aging_genes_merged$wb_annotation[ aging_genes_merged$gene_name %in% hpk1_atlas_neurons_trend$hpk1_up_atlas_down ])
table(aging_genes_merged$genage_gerogene_cats[ aging_genes_merged$gene_name %in% hpk1_atlas_neurons_trend$hpk1_up_atlas_down ])
table(aging_genes_merged$inPGP[ aging_genes_merged$gene_name %in% hpk1_atlas_neurons_trend$hpk1_up_atlas_down ])

# ENRICHMENT - across all the lists
hpk1_atlas_neurons_enrichKegg <- lapply(hpk1_atlas_neurons_trend, function(geneVec){
  goSeqEnrich_singleSet(geneVec = wb_gene_info$WBGeneID[wb_gene_info$gene_name %in% geneVec],
                        setLongObj = sets_kegg_all_long, setSelCols = c("gene_id", "pathway"))
})
unlist(lapply(hpk1_atlas_neurons_enrichKegg, nrow))

# Reactome
hpk1_atlas_neurons_enrichReac <- lapply(hpk1_atlas_neurons_trend, function(geneVec){
  goSeqEnrich_singleSet(geneVec = wb_gene_info$WBGeneID[wb_gene_info$gene_name %in% geneVec],
                        setLongObj = sets_reactome_long, setSelCols = c("gene_id", "pathway"))
})
unlist(lapply(hpk1_atlas_neurons_enrichReac, nrow))

# GO enrichment
hpk1_atlas_neurons_enrichGO <- lapply(hpk1_atlas_neurons_trend, function(geneVec){
  goSeqEnrich_singleSet(geneVec = wb_gene_info$WBGeneID[wb_gene_info$gene_name %in% geneVec],
                        setLongObj = elegansGO, setSelCols = c("GeneID", "GO"))
})
unlist(lapply(hpk1_atlas_neurons_enrichGO, nrow))

# TF enrichment 
hpk1_atlas_neurons_enrichTFs <- lapply(hpk1_atlas_neurons_trend, function(geneVec){
  goSeqEnrich_singleSet(geneVec = wb_gene_info$WBGeneID[wb_gene_info$gene_name %in% geneVec],
                        setLongObj = sets_tfs_species_filt_long, setSelCols = c("WBGeneID", "motif_tf"))
})
unlist(lapply(hpk1_atlas_neurons_enrichTFs, nrow))


# are those TFs usually expressed in neurons? from hpk1_up_atlas_down- 
# Seven are NHRs that have the same motif, and two different motifs for blmp-1 and unc-30
tmpTFs <- c("ceh-22", "blmp-1", "eor-1", "ceh-27", "tab-1", "klf-1", "nhr-110", "nhr-117", "nhr-125", "nhr-143", "nhr-177", "nhr-180",
  "nhr-55", "dmd-4", "ets-7", "F19F10.1", "ceh-24", "unc-30")

intersect(tmpTFs, aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno])

intersect(tmpTFs, union(kaletsky_enriched_neurons, kaletsky_unique_neurons))

View(aging_atlas_degs_long[ aging_atlas_degs_long$gene_name %in% tmpTFs  , ])
# 10 of these are differentially expressed with age, almost all downregulated
#  yet are upregulated in hpk-1 loss
DE_hpk1[ DE_hpk1$GeneSymbol %in% tmpTFs , ]

table(atlas_tissue_type_DEbroad$neuron[wTF$Public.name])
atlas_tissue_type_DEbroad$neuron[tmpTFs]

#
# we'll try to use heatmaps to show how broadly these genes are differentially expressed in
#  neurons with aging, to see if hpk-1 "up" genes that show similar or different aging patterns
#   have different activity in different neurons, etc. To make this comparable, we'll keep
#   the neurons in the same order between the plots.

atlasDE_UpUp <- aging_atlas_degs_long[ aging_atlas_degs_long$gene_name %in% hpk1_atlas_neurons_trend$both_up , ]

atlasDE_UpDown <- aging_atlas_degs_long[ aging_atlas_degs_long$gene_name %in% hpk1_atlas_neurons_trend$hpk1_up_atlas_down , ]

atlasDE_UpUp_fc_mat <- mat.or.vec(nr = length(unique(atlasDE_UpUp$gene_name)),
                                             nc = length(unique(atlasDE_UpUp$cluster_shortID)))

atlasDE_UpDown_fc_mat <- mat.or.vec(nr = length(unique(atlasDE_UpDown$gene_name)),
                                  nc = length(unique(atlasDE_UpDown$cluster_shortID)))

rownames(atlasDE_UpUp_fc_mat) <- unique(atlasDE_UpUp$gene_name)
colnames(atlasDE_UpUp_fc_mat) <- unique(atlasDE_UpUp$cluster_shortID)

rownames(atlasDE_UpDown_fc_mat) <- unique(atlasDE_UpDown$gene_name)
colnames(atlasDE_UpDown_fc_mat) <- unique(atlasDE_UpDown$cluster_shortID)

for(i in 1:nrow(atlasDE_UpUp_fc_mat)){
  for(j in 1:ncol(atlasDE_UpUp_fc_mat)){
    rv = atlasDE_UpUp$logFoldChange[
      (atlasDE_UpUp$gene_name == rownames(atlasDE_UpUp_fc_mat)[i]) &
        (atlasDE_UpUp$cluster_shortID == colnames(atlasDE_UpUp_fc_mat)[j])]
    if(length(rv) == 0){
      atlasDE_UpUp_fc_mat[i,j] <- NA
    }else{
      atlasDE_UpUp_fc_mat[i,j] <- rv 
    }
  }
}
rm(i,j, rv)

for(i in 1:nrow(atlasDE_UpDown_fc_mat)){
  for(j in 1:ncol(atlasDE_UpDown_fc_mat)){
    rv = atlasDE_UpDown$logFoldChange[
      (atlasDE_UpDown$gene_name == rownames(atlasDE_UpDown_fc_mat)[i]) &
        (atlasDE_UpDown$cluster_shortID == colnames(atlasDE_UpDown_fc_mat)[j])]
    if(length(rv) == 0){
      atlasDE_UpDown_fc_mat[i,j] <- NA
    }else{
      atlasDE_UpDown_fc_mat[i,j] <- rv 
    }
  }
}
rm(i,j, rv)


heatmapData_UpUp <- as.matrix(atlasDE_UpUp_fc_mat)
mode(heatmapData_UpUp) <- "numeric"
heatmapData_UpUp[is.na(heatmapData_UpUp)] <- 0

heatmapData_UpDown <- as.matrix(atlasDE_UpDown_fc_mat)
mode(heatmapData_UpDown) <- "numeric"
heatmapData_UpDown[is.na(heatmapData_UpDown)] <- 0

heatmapData_UpUp <- heatmapData_UpUp[ , tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData_UpUp)]]
heatmapData_UpDown <- heatmapData_UpDown[ , tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData_UpDown)]]

# CHECK ORDER OF COLUMNS TO MAKE SURE IT MATCHES COLOR ANNOTATION
identical(colnames(heatmapData_UpUp),
          tissue_type_vec$shortID[tissue_type_vec$shortID %in% colnames(heatmapData_UpUp) ])


heatmapData_UpUp_neurons <- heatmapData_UpUp[ , colnames(heatmapData_UpUp) %in%
                                                tissue_type_vec$shortID[tissue_type_vec$tissue_type == "neuron"] ]
heatmapData_UpDown_neurons <- heatmapData_UpDown[ , colnames(heatmapData_UpDown) %in%
                                                tissue_type_vec$shortID[tissue_type_vec$tissue_type == "neuron"] ]

# FILTER to include at least three non-zero rows
heatmapData_UpUp_neurons <- heatmapData_UpUp_neurons[ rowSums(heatmapData_UpUp_neurons == 0) < (ncol(heatmapData_UpUp_neurons)-4) , ]
heatmapData_UpDown_neurons <- heatmapData_UpDown_neurons[ rowSums(heatmapData_UpDown_neurons == 0) < (ncol(heatmapData_UpDown_neurons)-4) , ]

# remove empty columns
heatmapData_UpUp_neurons <- heatmapData_UpUp_neurons[ , colSums(heatmapData_UpUp_neurons == 0) != nrow(heatmapData_UpUp_neurons) ]
heatmapData_UpDown_neurons <- heatmapData_UpDown_neurons[ , colSums(heatmapData_UpDown_neurons == 0) != nrow(heatmapData_UpDown_neurons) ]

dim(heatmapData_UpUp_neurons)
dim(heatmapData_UpDown_neurons)

# order by number of non-zero rows
# heatmapData <- heatmapData[ order(rowSums(heatmapData == 0), decreasing = FALSE) , ]
# order by overall mean of row FC
heatmapData_UpUp_neurons <- heatmapData_UpUp_neurons[ order(rowMeans(heatmapData_UpUp_neurons), decreasing = TRUE) , ]
heatmapData_UpDown_neurons <- heatmapData_UpDown_neurons[ order(rowMeans(heatmapData_UpDown_neurons), decreasing = FALSE) , ]

# ---- finishing with the "up/up" set ----

# now we can color neurons by neuron type
atlas_cengen_neuron_map.sub <- merge(colnames(heatmapData_UpUp_neurons), atlas_cengen_neuron_map, by.x = 1, by.y = "atlas_cluster", all.x = TRUE)
colnames(atlas_cengen_neuron_map.sub)[1] <- "atlas_cluster"
atlas_cengen_neuron_map.sub <- merge(atlas_cengen_neuron_map.sub, setNames(atlas_unique_neuron_signal_class[,"consensus"],
                                                                           rownames(atlas_unique_neuron_signal_class)),
                                     by.x = "atlas_cluster", by.y = 0, all.x = TRUE)
colnames(atlas_cengen_neuron_map.sub)[ncol(atlas_cengen_neuron_map.sub)] <- "atlas_unique_signal_cons"
atlas_cengen_neuron_map.sub$NeuronType[is.na(atlas_cengen_neuron_map.sub$NeuronType)] <- aging_atlas_clusters_all$Clusters.explanatory.name[match( 
  atlas_cengen_neuron_map.sub$atlas_cluster[is.na(atlas_cengen_neuron_map.sub$NeuronType)],
  aging_atlas_clusters_all$shortID)]

atlas_cengen_neuron_map.sub$NeuronType[grepl("*inter neuron*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*inter neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*interneurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motor neuron*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motor neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motorneurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*sensory neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Sensory"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*PHC*|*PVM*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Sensory"
atlas_cengen_neuron_map.sub$NeuronType[!(atlas_cengen_neuron_map.sub$NeuronType %in% c("Interneuron", 
                                                                                       "Sensory", "Endorgan", 
                                                                                      "Motorneuron"))] <- "Other Neuron"

atlas_cengen_neuron_map.sub <- atlas_cengen_neuron_map.sub[ match(colnames(heatmapData_UpUp_neurons),
                                                                  atlas_cengen_neuron_map.sub$atlas_cluster) , ] 
identical(colnames(heatmapData_UpUp_neurons), atlas_cengen_neuron_map.sub$atlas_cluster) # should be true

tmp_type = atlas_cengen_neuron_map.sub$NeuronType
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]
tmp_type[tmp_type == "Other Neuron"] = colorBlindness::paletteMartin[10]

pdf(file = "Color key- Neuron functional type 2.pdf", width = 2, height = 4)
swatch2(c("Interneuron" = as.vector(colorBlindness::paletteMartin[2]),
          "Sensory" = as.vector(colorBlindness::paletteMartin[4]),
          "Endorgan" = as.vector(colorBlindness::paletteMartin[6]),
          "Motorneuron" = as.vector(colorBlindness::paletteMartin[8]),
          "Other neuron" = as.vector(colorBlindness::paletteMartin[10])))
graphics.off()

# UPDATE- for cells where there was no call in cengen expression, look for a call in the aging atlas expression
tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub$cengen_signal_class) | 
                             grepl("glutamatergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                            grepl("serotonergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                         grepl("GABAergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                            grepl("dopaminergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                           grepl("cholinergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

# USE CELL NAMES AS COLUMN LABELS (Clusters.explanatoryname in aging_atlas_clusters_all)

# so now we also want to sort the columns by neuron type.

svg(file = "Aging atlas DE across age for genes up in hpk-1 and only upregulated with age in neurons.svg", width = 18, height = 16)
heatmap.3(heatmapData_UpUp_neurons[ , order(neuron_colors[,1])],
          ColSideColors = as.matrix(neuron_colors[order(neuron_colors[,1]) , ]),
          RowSideColors = t(as.matrix(ifelse(rownames(heatmapData_UpUp_neurons) %in% 
                                               aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_UpUp_neurons[ , order(neuron_colors[,1])]), shortID)]),
          col = colorRampPalette(c("white", "darkred"))(64), #scale = "row",
          cexCol = 0.9, cexRow = 1.1, margins = c(12,8), keysize = 0.8, symkey = FALSE,
          colsep = seq(1, ncol(heatmapData_UpUp_neurons)),
          rowsep = seq(1, nrow(heatmapData_UpUp_neurons)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", 
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          KeyValueName = "LogFoldChange")
graphics.off()

# ---- finishing with the "up/down" set ----

# now we can color neurons by neuron type
atlas_cengen_neuron_map.sub <- merge(colnames(heatmapData_UpDown_neurons), atlas_cengen_neuron_map, by.x = 1, by.y = "atlas_cluster", all.x = TRUE)
colnames(atlas_cengen_neuron_map.sub)[1] <- "atlas_cluster"
atlas_cengen_neuron_map.sub <- merge(atlas_cengen_neuron_map.sub, setNames(atlas_unique_neuron_signal_class[,"consensus"],
                                                                           rownames(atlas_unique_neuron_signal_class)),
                                     by.x = "atlas_cluster", by.y = 0, all.x = TRUE)
colnames(atlas_cengen_neuron_map.sub)[ncol(atlas_cengen_neuron_map.sub)] <- "atlas_unique_signal_cons"
atlas_cengen_neuron_map.sub$NeuronType[is.na(atlas_cengen_neuron_map.sub$NeuronType)] <- aging_atlas_clusters_all$Clusters.explanatory.name[match( 
  atlas_cengen_neuron_map.sub$atlas_cluster[is.na(atlas_cengen_neuron_map.sub$NeuronType)],
  aging_atlas_clusters_all$shortID)]

atlas_cengen_neuron_map.sub$NeuronType[grepl("*inter neuron*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*inter neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*interneurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Interneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motor neuron*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motor neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*motorneurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Motorneuron"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*sensory neurons*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Sensory"
atlas_cengen_neuron_map.sub$NeuronType[grepl("*PHC*|*PVM*", atlas_cengen_neuron_map.sub$NeuronType)] <- "Sensory"
atlas_cengen_neuron_map.sub$NeuronType[!(atlas_cengen_neuron_map.sub$NeuronType %in% c("Interneuron", 
                                                                                       "Sensory", "Endorgan", 
                                                                                       "Motorneuron"))] <- "Other Neuron"

atlas_cengen_neuron_map.sub <- atlas_cengen_neuron_map.sub[ match(colnames(heatmapData_UpDown_neurons),
                                                                  atlas_cengen_neuron_map.sub$atlas_cluster) , ] 
identical(colnames(heatmapData_UpDown_neurons), atlas_cengen_neuron_map.sub$atlas_cluster) # should be true

tmp_type = atlas_cengen_neuron_map.sub$NeuronType
tmp_type[tmp_type == "Interneuron"] = colorBlindness::paletteMartin[2]
tmp_type[tmp_type == "Sensory"] = colorBlindness::paletteMartin[4]
tmp_type[tmp_type == "Endorgan"] = colorBlindness::paletteMartin[6]
tmp_type[tmp_type == "Motorneuron"] = colorBlindness::paletteMartin[8]
tmp_type[tmp_type == "Other Neuron"] = colorBlindness::paletteMartin[10]

swatch2(c("Interneuron" = as.vector(colorBlindness::paletteMartin[2]),
          "Sensory" = as.vector(colorBlindness::paletteMartin[4]),
          "Endorgan" = as.vector(colorBlindness::paletteMartin[6]),
          "Motorneuron" = as.vector(colorBlindness::paletteMartin[8]),
          "Other neuron" = as.vector(colorBlindness::paletteMartin[10])))

# UPDATE- for cells where there was no call in cengen expression, look for a call in the aging atlas expression
tmp_glutamatergic = ifelse(grepl("glutamatergic", atlas_cengen_neuron_map.sub$cengen_signal_class) | 
                             grepl("glutamatergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                           "black", "white")
tmp_serotonergic = ifelse(grepl("serotonergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                            grepl("serotonergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                          "black", "white")
tmp_GABAergic = ifelse(grepl("GABAergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                         grepl("GABAergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                       "black", "white")
tmp_dopaminergic = ifelse(grepl("dopaminergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                            grepl("dopaminergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                          "black", "white")
tmp_cholinergic = ifelse(grepl("cholinergic", atlas_cengen_neuron_map.sub$cengen_signal_class) |
                           grepl("cholinergic", atlas_cengen_neuron_map.sub$atlas_unique_signal_cons),
                         "black", "white")

neuron_colors <- cbind("neuron_type" = tmp_type,
                       "is_glutamatergic" = tmp_glutamatergic,
                       "is_serotonergic" = tmp_serotonergic,
                       "is_GABAergic" = tmp_GABAergic,
                       "is_dopaminergic" = tmp_dopaminergic,
                       "is_cholinergic" = tmp_cholinergic)

# USE CELL NAMES AS COLUMN LABELS (Clusters.explanatoryname in aging_atlas_clusters_all)

# so now we also want to sort the columns by neuron type.

svg(file = "Aging atlas DE across age for genes up in hpk-1 and only downregulated with age in neurons.svg", width = 20, height = 18)
heatmap.3(heatmapData_UpDown_neurons[ , order(neuron_colors[,1])],
          ColSideColors = as.matrix(neuron_colors[order(neuron_colors[,1]) , ]),
          RowSideColors = t(as.matrix(ifelse(rownames(heatmapData_UpDown_neurons) %in% 
                               aging_genes_merged$gene_name[aging_genes_merged$any_aging_pheno], "gold4", "black"))),
          labCol = with(aging_atlas_clusters_all, Curated.annotations[
            match(colnames(heatmapData_UpDown_neurons[ , order(neuron_colors[,1])]), shortID)]),
          col = colorRampPalette(c("darkblue", "white", "darkred"))(64), #scale = "row",
          cexCol = 0.9, cexRow = 0.9, margins = c(12,8), keysize = 0.8, symkey = FALSE,
          colsep = seq(1, ncol(heatmapData_UpDown_neurons)),
          rowsep = seq(1, nrow(heatmapData_UpDown_neurons)),
          sepcol = "grey80", sepwidth = c(0.0001, 0.0001),
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", 
          main = "Aging atlas DE, adj p-value < 0.05 and LFC > 0.5",
          KeyValueName = "LogFoldChange")
graphics.off()


pdf("venn- hpk-1 DE up vs atlas neuron all up.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 upregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "atlas neuron upregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_up"]),
  filename = NULL,
  fill = c("tomato", "gold"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("tomato3", "gold4"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)
graphics.off()

pdf("venn- hpk-1 DE up vs atlas neuron all down.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 upregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "atlas neuron downregulated" = names(atlas_tissue_type_DEbroad$neuron)[atlas_tissue_type_DEbroad$neuron == "all_down"]),
  filename = NULL,
  fill = c("tomato", "grey66"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("tomato3", "grey33"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)
graphics.off()


#
# For future reference- VENNs of hpk-1 DE genes and:
# * Hobert "neuronal genome"
# * Kaletsky et al 2018 neuron-specific genes
# * Kaletsky et al 2018 neuron-enriched genes
#


pdf("venn- hpk-1 DE up vs hobert neuronal genome vs kaletsky unique_enriched.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 upregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "Neuronal genome" = neuron_gene_classes$gene,
         "Kaletsky neuron unique" = kaletsky_unique_neurons,
         "kaletsky neuron enriched" =kaletsky_enriched_neurons),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold", "grey66"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4", "grey33"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()

pdf("venn- hpk-1 DE down vs hobert neuronal genome vs kaletsky unique_enriched.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 downregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0],
         "Neuronal genome" = neuron_gene_classes$gene,
         "Kaletsky neuron unique" = kaletsky_unique_neurons,
         "kaletsky neuron enriched" =kaletsky_enriched_neurons),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold", "grey66"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4", "grey33"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()

pdf("venn- hpk-1 DE up vs hobert neuronal genome vs kaletsky unique OR enriched.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 upregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "Neuronal genome" = neuron_gene_classes$gene,
         "Kaletsky neuron unique/enriched union" = union(kaletsky_unique_neurons, kaletsky_enriched_neurons)),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()


# neuron expressed/enriched and upregulated with age
summary(aging_atlas_DE_up_down_count_all$num_clusters_aging_up[aging_atlas_DE_up_down_count_all$count_sum_cat == "all_up"])
summary(aging_atlas_DE_up_down_count_all$num_clusters_aging_down[aging_atlas_DE_up_down_count_all$count_sum_cat == "all_down"])


pdf("venn- hpk-1 DE up vs aging upregulated (only upreg any cluster count) vs kaletsky unique OR enriched.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 upregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "aging atlas- genes only upregulated with age in any number of clusters" = rownames(aging_atlas_DE_up_down_count_all[
           aging_atlas_DE_up_down_count_all$count_sum_cat == "all_up" , ]),
         "Kaletsky neuron unique/enriched union" = union(kaletsky_unique_neurons, kaletsky_enriched_neurons)),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()

pdf("venn- hpk-1 DE up vs aging downregulated (only downreg any cluster count) vs kaletsky unique OR enriched.pdf", width = 10, height = 8)
plot.new()
grid.draw(venn.diagram(
  x=
    list("hpk-1 upregulated" = DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0],
         "aging atlas- genes only downregulated with age in any number of clusters" = rownames(aging_atlas_DE_up_down_count_all[
           aging_atlas_DE_up_down_count_all$count_sum_cat == "all_down" , ]),
         "Kaletsky neuron unique/enriched union" = union(kaletsky_unique_neurons, kaletsky_enriched_neurons)),
  filename = NULL,
  fill = c("cornflowerblue", "tomato", "gold"),
  alpha = 0.50,
  fontface = "bold",
  cex = 2,
  cat.cex=1,
  cat.col = c("darkblue",  "tomato3", "gold4"),
  margin = 0.1, print.mode = c("raw"), cat.dist = 0.1)
)

graphics.off()


# also want to know hypergeometric p-value for overlap between hobert list and up-regulated genes
#  "universe": # of genes expressed in RNA-Seq dataset (union of RNA-Seq and the hobert list)
#  Set A: DE genes in hpk-1 DE
#  Set B: Genes in hobert list
#  specified as formulated for phyper with lower.tail as false: hitInSample-1 (q), hitInPop (m), failInPop (n), sampleSize (k)
#  m = the # of genes of Set A 
#  n = length(universe) - length(set A)
#  k = length(set B)
#  q = length(intersect(Set A, Set B)) - 1

# up genes
m = length(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0])
n = length(union(neuron_gene_classes$gene, DE_hpk1.all$GeneSymbol)) - m
k = length(unique(neuron_gene_classes$gene))
q = length(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0], unique(neuron_gene_classes$gene)))

phyper(q-1, m, n, k, lower.tail = FALSE)

# down genes
m = length(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0])
n = length(union(neuron_gene_classes$gene, DE_hpk1.all$GeneSymbol)) - m
k = length(unique(neuron_gene_classes$gene))
q = length(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC < 0], unique(neuron_gene_classes$gene)))

phyper(q-1, m, n, k, lower.tail = FALSE)

rm(m, n, k, q)

# what about with the neuron unique/enriched list?
# up genes
m = length(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0])
n = length(union(union(kaletsky_enriched_neurons, kaletsky_unique_neurons),
                 DE_hpk1.all$GeneSymbol)) - m
k = length(union(kaletsky_enriched_neurons, kaletsky_unique_neurons))
q = length(intersect(DE_hpk1$GeneSymbol[DE_hpk1$DESeq2_log2FC > 0], union(kaletsky_enriched_neurons, kaletsky_unique_neurons)))

phyper(q-1, m, n, k, lower.tail = FALSE)

rm(m, n, k, q)

#
# Sidequest- GATK annotated variants- genes with some kind of predicted impact
#


path_varFile <- file.path("E:/project_data_ssd/myc_project/GATK/snpeff_annotated_variants/summary")

varFiles <- list("HPK1_US_64" =  read.table(file.path(path_varFile, "clt_HPK-1_US_64_summary.genes.txt"), sep = "\t", comment.char = "", skip = 1, header = TRUE),
                 "HPK1_US_68" =  read.table(file.path(path_varFile, "clt_HPK-1_US_68_summary.genes.txt"), sep = "\t", comment.char = "", skip = 1, header = TRUE),
                 "HPK1_US_72" =  read.table(file.path(path_varFile, "clt_HPK-1_US_72_summary.genes.txt"), sep = "\t", comment.char = "", skip = 1, header = TRUE),
                 "N2_US_63" = read.table(file.path(path_varFile, "clt_N2_US_63_summary.genes.txt"), sep = "\t", comment.char = "", skip = 1, header = TRUE),
                 "N2_US_67" = read.table(file.path(path_varFile, "clt_N2_US_67_summary.genes.txt"), sep = "\t", comment.char = "", skip = 1, header = TRUE),
                 "N2_US_71" = read.table(file.path(path_varFile, "clt_N2_US_71_summary.genes.txt"), sep = "\t", comment.char = "", skip = 1, header = TRUE))

varFiles <- lapply(varFiles, function(x){
  x$gene_transcript = paste(x[,"X.GeneName"], x[,"TranscriptId"], sep = "_")
  return(x)
})

# gene/transcript IDs are DIFFERENT for each file!
union_gene_transcript <- unique(unlist(sapply(names(varFiles), function(currName){varFiles[[currName]][,"gene_transcript"]})))
varSum_impact_high <- data.frame(row.names = union_gene_transcript, mat.or.vec(length(union_gene_transcript), length(varFiles)))
colnames(varSum_impact_high) <- names(varFiles)

varSum_impact_moderate <- varSum_impact_high # since these start off from the same state.
varSum_impact_low <- varSum_impact_high
# relevant columns in original matrix- variants_impact_HIGH, variants_impact_MODERATE, variants_impact_LOW

for(currName in names(varFiles)){
  
  currMat = varFiles[[currName]][  ,
                                   c("gene_transcript", "variants_impact_HIGH", "variants_impact_MODERATE", "variants_impact_LOW")]
  currMat = merge(union_gene_transcript, currMat, by.x = 1, by.y = "gene_transcript", all.x = TRUE)
  rownames(currMat) = currMat[,1]
  currMat = currMat[,-1]
  # reorder
  currMat = currMat[ union_gene_transcript , ]
  
  varSum_impact_high[,currName] <- currMat[,"variants_impact_HIGH"]
  varSum_impact_moderate[,currName] <- currMat[,"variants_impact_MODERATE"]
  varSum_impact_low[,currName] <- currMat[,"variants_impact_LOW"]
}
rm(currMat, currName)

# replace NAs with 0
varSum_impact_low[is.na(varSum_impact_low)] <- 0
varSum_impact_moderate[is.na(varSum_impact_moderate)] <- 0
varSum_impact_high[is.na(varSum_impact_high)] <- 0

# remove rows that are all 0s

varSum_impact_low <- varSum_impact_low[rowSums(varSum_impact_low == 0) < ncol(varSum_impact_low) , ]
varSum_impact_moderate <- varSum_impact_moderate[rowSums(varSum_impact_moderate == 0) < ncol(varSum_impact_moderate) , ]
varSum_impact_high <- varSum_impact_high[rowSums(varSum_impact_high == 0) < ncol(varSum_impact_high) , ]

nrow(varSum_impact_low)
nrow(varSum_impact_moderate)
nrow(varSum_impact_high)

# write out results
wb <- createWorkbook()
addWorksheet(wb, "high_impact")
addWorksheet(wb, "moderate_impact")
addWorksheet(wb, "low_impact")
writeDataTable(wb, sheet = "high_impact", x =  varSum_impact_high, rowNames = TRUE)
writeDataTable(wb, sheet = "moderate_impact", x =  varSum_impact_moderate, rowNames = TRUE)
writeDataTable(wb, sheet = "low_impact", x =  varSum_impact_low, rowNames = TRUE)
saveWorkbook(wb, file = "Hpk-1 project unstressed- variant analysis affected genes by predicted effect.xlsx")
rm(wb)

#
# add WBGeneID to hpk-1 DE result table for paper
#
wb_gene_info_all <- read.csv(file.path(paths$wormbase, "c_elegans.PRJNA13758.WS284.geneIDs.txt"), header = FALSE)[ , 2:6 ]
colnames(wb_gene_info_all) <- c("WBGeneID", "gene_name", "molecular_name", "status", "gene_biotype")
# if gene name is NA, and molecular name is not NA, assign the molecular name as the gene name
wb_gene_info_all[ wb_gene_info_all == "" ] <- NA
wb_gene_info_all$gene_name[is.na(wb_gene_info_all$gene_name)] <- wb_gene_info_all$molecular_name[is.na(wb_gene_info_all$gene_name)]

paper_SF06_1 <- merge(DE_hpk1.all[,c(1:3, 8,9, 11,12)], wb_gene_info_all, by.x = "GeneSymbol", by.y = "gene_name", all.x = TRUE)
# for cases where the WBGeneID is still NA, try matching to the molecular name
# pass 1
paper_SF06_1_unmapped <- paper_SF06_1[ is.na(paper_SF06_1$WBGeneID) , 1:7 ]
paper_SF06_1 <- paper_SF06_1[ !is.na(paper_SF06_1$WBGeneID) , ]
paper_SF06_1_unmapped <- merge(paper_SF06_1_unmapped, wb_gene_info_all, by.x = "GeneSymbol", by.y = "molecular_name", all.x = TRUE)
paper_SF06_1_unmapped_2 <- paper_SF06_1_unmapped[ is.na(paper_SF06_1_unmapped$WBGeneID) , 1:7 ]
paper_SF06_1_unmapped <- paper_SF06_1_unmapped[ !is.na(paper_SF06_1_unmapped$WBGeneID) , ]
colnames(paper_SF06_1_unmapped)[colnames(paper_SF06_1_unmapped) == "GeneSymbol"] <- "molecular_name"
colnames(paper_SF06_1_unmapped)[colnames(paper_SF06_1_unmapped) == "gene_name"] <- "GeneSymbol"

paper_SF06_1 <- rbind(paper_SF06_1, paper_SF06_1_unmapped[ , colnames(paper_SF06_1) ])
# pass 2- just use old annotation

wb_gene_info_OLD <- read.csv(file.path("E:/public_data_annotation/wormbase/WS254",
                                       "c_elegans.PRJNA13758.WS254.geneIDs.txt"), header = FALSE)[ , 2:4 ]
colnames(wb_gene_info_OLD) <- c("WBGeneID", "gene_name", "molecular_name")
# if gene name is NA, and molecular name is not NA, assign the molecular name as the gene name
wb_gene_info_OLD[ wb_gene_info_OLD == "" ] <- NA
wb_gene_info_OLD$gene_name[is.na(wb_gene_info_OLD$gene_name)] <- wb_gene_info_OLD$molecular_name[is.na(wb_gene_info_OLD$gene_name)]

paper_SF06_1_unmapped_2 <- merge(paper_SF06_1_unmapped_2, wb_gene_info_OLD, by.x = "GeneSymbol", by.y = "gene_name", all.x = TRUE)
colnames(paper_SF06_1_unmapped_2)

paper_SF06_1 <- paper_SF06_1[,1:9]
colnames(paper_SF06_1)
sum(is.na(paper_SF06_1_unmapped_2$WBGeneID)) # still a few that seem to not be dead, but have symbols that changed in 2014ish- 
#  which was before WS254- this must have been from the ensemble annotation in the GTF I was using being a little bit behind wormbase then

paper_SF06_1 <- rbind(paper_SF06_1, paper_SF06_1_unmapped_2)

rm(paper_SF06_1_unmapped, paper_SF06_1_unmapped_2)

# merge in homology information
paper_SF06_1 <- merge(paper_SF06_1, ortholist_collapsed[,1:3], by.x = "WBGeneID", by.y = 1, all.x = TRUE)

paper_SF06_2 <- paper_SF06_1[ (abs(paper_SF06_1$DESeq2_log2FC) >= 1) & (paper_SF06_1$DESeq2_fdrpval < 0.05)  , ]

write.xlsx(paper_SF06_1,file = file.path("E:/project_data_ssd/hpk-1_dataset", "SF06_S1_all_genes_FC_hpk1.xlsx"), asTable = TRUE)
write.xlsx(paper_SF06_2,file = file.path("E:/project_data_ssd/hpk-1_dataset", "SF06_S2_DE_genes_hpk1.xlsx"), asTable = TRUE)


#
# Read in original count data-
#  this was assembled with Featurecounts on ALL of the RNA-Seq samples, so we need to re-subset it and probably also
# keep only the genes we called expressed
#

origCounts <- read.table(file.path("F:/project_data/myc project/myc and hpk-1 project/hpk-1 GEO submission 2022",
                                   "STARaligned_Featurecounts_gene.txt"), sep = "\t", quote = "", header = TRUE, row.names = 1)[,c(-1,-2,-3,-4,-5)]
origCounts <- origCounts[ , grepl("*US*", colnames(origCounts)) ]

# rename samples appropriately, using original IDs
colnames(origCounts)[grepl("64", colnames(origCounts))] <- "hpk-1_basal_d2_1"
colnames(origCounts)[grepl("68", colnames(origCounts))] <- "hpk-1_basal_d2_2"
colnames(origCounts)[grepl("72", colnames(origCounts))] <- "hpk-1_basal_d2_3"

colnames(origCounts)[grepl("63", colnames(origCounts))] <- "N2_basal_d2_1"
colnames(origCounts)[grepl("67", colnames(origCounts))] <- "N2_basal_d2_2"
colnames(origCounts)[grepl("71", colnames(origCounts))] <- "N2_basal_d2_3"

# Need original GTF I used from ensembl 82-
# ...which I already have loaded.
colnames(gtf)

origCounts$WBGeneID <- rownames(origCounts)
# This must be true
identical(origCounts$WBGeneID, gtf$gene_id)
origCounts$GeneSymbol <- gtf$gene_name

# subset count rows

origCounts_filt <- origCounts[ origCounts$GeneSymbol %in% DE_hpk1.all$GeneSymbol , ]
nrow(origCounts_filt)
nrow(DE_hpk1.all)

# and that's basically all I need I think, so I could just write this out then.
origCounts_filt <- origCounts_filt[,c(7,8,1:6)]

write.table(origCounts_filt, file.path("F:/project_data/myc project/myc and hpk-1 project/hpk-1 GEO submission 2022",
                                       "hpk-1 and N2 D2 RNA-Seq- filtered raw counts.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# Also write out the TPM matrix, again, after filtering in the same manner
tpm_filt <- tpm[ origCounts_filt$WBGeneID , ]
identical(rownames(tpm_filt), rownames(origCounts_filt)) # should be true!

# rename samples appropriately, using original IDs
colnames(tpm_filt)[grepl("64", colnames(tpm_filt))] <- "hpk-1_basal_d2_1"
colnames(tpm_filt)[grepl("68", colnames(tpm_filt))] <- "hpk-1_basal_d2_2"
colnames(tpm_filt)[grepl("72", colnames(tpm_filt))] <- "hpk-1_basal_d2_3"

colnames(tpm_filt)[grepl("63", colnames(tpm_filt))] <- "N2_basal_d2_1"
colnames(tpm_filt)[grepl("67", colnames(tpm_filt))] <- "N2_basal_d2_2"
colnames(tpm_filt)[grepl("71", colnames(tpm_filt))] <- "N2_basal_d2_3"

tpm_filt <- as.data.frame(tpm_filt)
tpm_filt$WBGeneID <- origCounts_filt$WBGeneID
tpm_filt$GeneSymbol <- origCounts_filt$GeneSymbol
tpm_filt <- tpm_filt[,c(7,8,1:6)]

write.table(tpm_filt, file.path("F:/project_data/myc project/myc and hpk-1 project/hpk-1 GEO submission 2022",
                                       "hpk-1 and N2 D2 RNA-Seq- filtered TPM.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
