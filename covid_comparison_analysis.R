#Day40 Comparing two analyses ------------------
#Load Packages
library(tidyverse)
library(VennDiagram)
library(pheatmap)

source('covid_rnaseq_analysis.R', echo = TRUE)
source('covid_a549_analysis.R')

#nhbe 1st filter results
#annotated_df
#a549 1st filter results
#a549_anno_df

nhbe_ensgene= annotated_df$ensgene
a549_ensgene= a549_anno_df$ensgene

common_ensgene  = intersect(nhbe_ensgene, a549_ensgene)
length(common_ensgene)
nhbe_commomn = annotated_df[annotated_df$ensgene %in% common_ensgene, 
             c('ensgene', 'log2FoldChange')]
a549_commomn = a549_anno_df[a549_anno_df$ensgene %in% common_ensgene, 
                            c('ensgene', 'log2FoldChange')]
common_fc = left_join(nhbe_commomn,a549_commomn, by = c ('ensgene' = 'ensgene'))

ggplot(common_fc, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = 0.3) + geom_abline(slope = 1, colour ='red') #+
  #xlim(min=-5, max =5) + ylim(min=-5, max =5)

cor.test(common_fc$log2FoldChange.x, common_fc$log2FoldChange.y,
         method = 'spearman')

nhbe_deg= annotated_df1.2$ensgene
a549_deg= a549_anno2$ensgene
common_deg  = intersect(nhbe_deg, a549_deg)
length(common_deg)
nhbe_commomn_deg = annotated_df1.2[annotated_df1.2$ensgene %in% common_deg, 
                            c('ensgene', 'log2FoldChange')]
a549_commomn_deg = a549_anno2[a549_anno2$ensgene %in% common_deg, 
                            c('ensgene', 'log2FoldChange')]
common_fc_deg = left_join(nhbe_commomn_deg,a549_commomn_deg, 
                          by = c ('ensgene' = 'ensgene'))

ggplot(common_fc_deg, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = 0.3) + geom_abline(slope = 1, colour ='red') #+
#xlim(min=-5, max =5) + ylim(min=-5, max =5)

cor.test(common_fc_deg$log2FoldChange.x, common_fc_deg$log2FoldChange.y,
         method = 'spearman')
nhbe_only = setdiff(nhbe_deg, a549_deg)
a549_only = setdiff(a549_deg, nhbe_deg)

#venn.diagram(x=list(nhbe_deg,a549_deg)) #Does not work
plot.new()
draw.pairwise.venn(area1 = length(nhbe_deg), area2 = length(a549_deg), 
                   cross.area = length(common_deg), scaled = TRUE, 
                   fill =c('red', 'blue'), alpha = 0.5)
union_deg = union(nhbe_deg, a549_deg)
union_fc_hm = filter(common_fc, ensgene %in% union_deg)
union_fc_hm_mat = as.matrix(union_fc_hm[,2:3])

my_colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
#my_breaks = c(seq(min(union_fc_hm_mat), -0.01, length.out=50), 0,
              #seq(0.01, max(union_fc_hm_mat), length.out=50))
my_breaks = c(seq(-4, -0.01, length.out=50), 0, 
              seq(0.01, 4, length.out=50))
pheatmap(union_fc_hm_mat, breaks = my_breaks, color = my_colours)
