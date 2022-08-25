#Day 28 -------------
#Load Packages #Enhanced Volcano

#library(readr)
#library(dplyr) has conflict with tidyverse
#library(magrittr)
library(BiocManager)
library(tximport)
library(DESeq2)
library(plotly)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

##NOT data.frame BUT tibble
#sample_table = read_csv("SraRunTable2.txt") %>% #View() #pipe symbol= %>%
#  select(`Sample Name`, source_name, Treatment,
#         Cell_Line, Cell_type, Time_point) %>%
#  slice(seq(1, 48, by=4)) #%>% View()
#, df[!duplicated(df$column),]
#install('tximport')

sample_table = read_csv("SraRunTable2.txt") %>% #select is dplyr function
  dplyr::select(`Sample Name`,
                source_name,
                Treatment,
                Cell_Line,
                Cell_type,
                Time_point) %>%
  unique()

#sample_table[sample_table$Cell_Line=='A549',]
st_a549 = filter(sample_table, Cell_Line=='A549')
'GSM4432378/quant.sf'
loc = paste0(pull(st_a549, `Sample Name`), '/quant.sf')
names(loc) = pull(st_a549, `Sample Name`)
gene_map = read_csv("gene_map.csv", col_names = c('esntid', 'ensgid'))

count_data_a549 = tximport(files = loc, type = "salmon", 
                           tx2gene = gene_map, ignoreTxVersion = TRUE)

st_a549$conditions = factor(rep(c('mock', 'infected'), each =3),
                                      levels = c('mock', 'infected'))
dds_a549 = DESeqDataSetFromTximport(txi = count_data_a549,
                                    colData = st_a549,
                                    design =  ~ conditions)
#DESeq does: SizeFactors, Dispersiosns, nbinomWaldTest
dds_a549 = DESeq(dds_a549)
vst_a549 = varianceStabilizingTransformation(dds_a549)
plotPCA(vst_a549, intgroup = 'conditions') + theme_bw()
pca_a549  = prcomp(t(assay(vst_a549)))
df = as.data.frame(pca_a549$x)
df$condition = st_a549$conditions
pve = round(pca_a549$sdev^2/sum(pca_a549$sdev^2) * 100, 2)
ggplot(df, aes(x=PC1, y = PC2)) + geom_point(aes(colour = condition)) + 
  theme_bw() + xlab(label = paste0("PC1 (", pve[1], "%)")) + 
  ylab(label = paste0("PC1 (", pve[2], "%)"))

#Retrieve and filter results #contrast used to set control and test condition
results_a549 = results(dds_a549, contrast = c('conditions', 'infected', 'mock'))
results_a549df = as.data.frame(results_a549)
#filter_a549df = filter(as.data.frame(results_a549, complete.cases))
filter_a549df = as.data.frame(results_a549[complete.cases(results_a549),])
filter_a549df1 = filter_a549df[filter_a549df$padj < 0.05, ]
filter_a549df1 = filter_a549df1[abs(filter_a549df1$log2FoldChange) > 1, ]

#Plots
plotMA(results_a549)

filter_a549df$test = filter_a549df$padj < 0.05 & abs(filter_a549df$log2FoldChange) > 1
filter_a549df = rownames_to_column(filter_a549df, var = 'ensgene')
g = ggplot(filter_a549df, aes(x = log2FoldChange, y = -log10(padj), name = ensgene )) +
  geom_point(aes(color=test), size = 1, alpha = 0.3) + 
  scale_color_manual(values = c('black', 'red')) +
  geom_vline(xintercept = 1, color = 'green', linetype=3) +
  geom_vline(xintercept = -1, color = 'green', linetype=3) +
  geom_hline(yintercept = -log10(0.05), color = 'blue', linetype=3)+
  theme_bw() + theme(legend.position = 'none')
ggplotly(g)

#Annotation #If data not available on Biomart use readGFF()
ensembl = useEnsembl(biomart = 'ensembl')
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
a549_anno = getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                  'start_position', 'end_position', 'strand',
                                  'gene_biotype', 'description',
                                  'external_gene_name'), 
                   filters = c('ensembl_gene_id'),
                   values = filter_a549df$ensgene, mart = ensembl)

a549_anno_df = left_join(filter_a549df, a549_anno, 
                         by=c('ensgene' = 'ensembl_gene_id'))

g = ggplot(a549_anno_df, aes(x = log2FoldChange, y = -log10(padj), name = external_gene_name )) +
  geom_point(aes(color=test), size = 1, alpha = 0.3) + 
  geom_vline(xintercept = 1, color = 'green', linetype=3) +
  geom_vline(xintercept = -1, color = 'green', linetype=3) +
  geom_hline(yintercept = -log10(0.05), color = 'blue', linetype=3)+
  scale_color_manual(values = c('black', 'red')) + theme_bw() +
  theme(legend.position = 'none')
ggplotly(g)

a549_anno1 = filter(a549_anno_df, padj < 0.05)
a549_anno2 = filter(a549_anno1, abs(log2FoldChange) > 1)
a549_degs = a549_anno2$ensgene
vst_a549_mat = assay(vst_a549)
hm_a549_data = vst_a549_mat[a549_degs,]
rownames(hm_a549_data) = a549_anno2$external_gene_name
heatmap(hm_a549_data)
pheatmap(hm_a549_data, fontsize_row = 4, scale = 'row', cutree_cols = 2, 
         cutree_rows = 2) #, color = blue)

#GO enrichment, Kegg ------------------
ent_gene_a549 = getBM(attributes = c('entrezgene_id'), 
                 filters = c('ensembl_gene_id'),
                 values = a549_anno2$ensgene, mart = ensembl)
ent_gene_a549 =ent_gene_a549$entrezgene_id
ent_gene_a549 = as.character(ent_gene_a549)
ent_uni_a549 = getBM(attributes = c('entrezgene_id'), 
                filters = c('ensembl_gene_id'),
                values = a549_anno_df$ensgene, mart = ensembl)
ent_uni_a549 =ent_uni_a549$entrezgene_id
ent_uni_a549 = as.character(ent_uni_a549)
ego = enrichGO(gene = ent_gene_a549, OrgDb = org.Hs.eg.db, 
               ont = "BP", universe = ent_uni_a549, readable = TRUE)
View(summary(ego))
barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)
fold_change = a549_anno2$log2FoldChange
names(fold_change) = a549_anno2$external_gene_name
cnetplot (ego, showCategory = 20, foldChange = fold_change)
goplot(ego, showCategory = 20)
keg = enrichKEGG(gene = ent_gene, universe = ent_uni)
View(summary(keg))