#Install packages#

#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#BiocManager::install("tidyverse") # includes ggplot2 and dplyr
#BiocManager::install("EnhancedVolcano")
#install.packages("plotly")
#BiocManager::install("ComplexHeatmap") #alternative = pheatmap
#BiocManager::install( "AnnotationDbi" )
#install.packages("devtools")
#install.packages("Rtools")
#BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install( "clusterProfiler" )
# BiocManager::install( "enrichplot" )
# BiocManager::install("DOSE")
#BiocManager::install("ReactomePA", force = TRUE)
#BiocManager::install("msigdbr", force = TRUE)
#BiocManager::install("KEGGREST", force = TRUE)
# BiocManager::install("BiocUpgrade") ## you may need this
# BiocManager::install("GOSemSim", force = TRUE)
# install.packages("remotes")
# library(devtools)
# devtools::install_github("GuangchuangYu/GOSemSim")
#install.packages(c("ggraph", "ggnetwork")) #altering ggplots 
#devtools::install_github("datapplab/pathview")
#update.packages(c("lattice", "spatial"))
#install.packages("igraph")

#Load libraries#
#library(DESeq2)
library(AnnotationDbi)
library(devtools)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(clusterProfiler) # for PEA analysis
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(ReactomePA)
library(msigdbr)
library(KEGGREST)
library(httr2)
library(DBI)
library(GOSemSim)
library(lattice)
library(spatial)
library(pathview)
library(cowplot)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
library(igraph)

#Set directory, input and output path
getwd()
in_path <- paste0(getwd(),"/Input/") # input path, where your input data is located
out_path <- paste0(getwd(),"/ORA/") # output path, where you want your results exported to
out_heat <- paste0(getwd(),"/ORA-heatmap/") # output path for heatmap results

#Read in the raw gene expression data
df <- read.csv(paste0(in_path, 'DEgenes-sig_ens.csv'), row.names = 1)
#Annotate genes based on differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2FoldChange > 1 & padj < 0.05 ~ 'UP',
  log2FoldChange < -1 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))


unique(df$diffexpressed)
na.omit(df)

##Creating a list of up or down DEGs for analysis
# Remove non-significant genes
df_O <- df[df$diffexpressed != 'NO', ]

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
ORA_deg_list <- split(df_O, df_O$diffexpressed)

#saveRDS(ORA_deg_list, file = paste0(in_path,"ORA_deg_list.RDS"))

#Split the genelist into up and down list vectors after converting the IDs to entrez and selecting entrez column

edb <- EnsDb.Hsapiens.v86 #alternative would be org.Hs.eg.db
ORA_up <- bitr(rownames(ORA_deg_list$UP), fromType = "GENEID", toType = "ENTREZID", OrgDb = edb)$ENTREZID

ORA_down <- bitr(rownames(ORA_deg_list$DOWN), fromType = "GENEID", toType = "ENTREZID", OrgDb = edb)$ENTREZID

head(ORA_up)

head(ORA_down)
#save datasets
# write.csv(ORA_up, file = paste0(in_path,"ORA_UP_input.csv"))
# write.csv(ORA_down, file = paste0(in_path,"ORA_DN_input.csv"))

##ORA GO analysis#
##up genes##
ego_BP_UP <- enrichGO(ORA_up,
         org.Hs.eg.db,
         keyType = "ENTREZID",
         ont = "BP",
         pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

#Save enrichGO output
# saveRDS(ego_BP_UP, file = paste0(out_path,"GO/ego_BP_UP.RDS"))
# ego_BP_UP <- readRDS(file = "ego_BP_UP.RDS")

#Extract the result table(top 10) and create a new column saying up genes
BP_UP <- as.data.frame(ego_BP_UP@result) %>% 
  slice_head(n = 10) %>%
  mutate(FC = "UP")

slice_head(BP_UP, n =10)
#Save All or top 10 up enrichment results
# write.csv(ego_BP_UP, file = paste0(out_path, "GO/EGO_BP_UP_all.csv"))
# write.csv(BP_UP, file = paste0(out_path, "GO/EGO_BP_UP_top10.csv"))

#goplot Up genes
#pdf(file = paste0(out_path, "GO/EGO_BP.pdf"), width = 10, height = 8, bg = "white")
gpu <- goplot(ego_BP_UP, showCategory = 10)
print(gpu)

#barplot Up genes
#pdf(file = paste0(out_path, "GO/EGO_bar_BP.pdf"), width = 6, height = 6, bg = "white")
bpu <- barplot(ego_BP_UP, showCategory=10, title = "Biological Processes(up)")

print(bpu)


##ORA Down genes##
ego_BP_DN <- enrichGO(ORA_down,
         org.Hs.eg.db,
         keyType = "ENTREZID",
         ont = "BP",
         pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

# Save the RDS file
#saveRDS(ego_BP_DN, file = paste0(out_path,"GO/ego_BP_DN.RDS"))

#Load data from RDS file
#ego_BP_DN <- readRDS("ego_BP_DN.RDS")

#Extract the result table(top 10) and create a new column saying down genes
BP_DN <- ego_BP_DN@result %>% 
  slice_head(n=10) %>%
  mutate(FC = "DOWN")

slice_head(BP_DN, n=10)
#Save All or top 10 down enrichment results
# write.csv(ego_BP_DN, file = paste0(out_path, "GO/EGO_BP_DN_all.csv"))
# write.csv(BP_DN, file = paste0(out_path, "GO/EGO_BP_DN_top10.csv"))

# GO plot down genes
#pdf(file = paste0(out_path, "GO/EGO_BP_DN.pdf"), width = 10, height = 8, bg = "white")
gpd <- goplot(ego_BP_DN, color = "p.adjust", title = "Biological Processes(down") 
print(gpd)


# For barplot:
#pdf(file = paste0(out_path, "GO/EGO_bar_DN.pdf"), width = 6, height = 6, bg = "white")
bpd <- barplot(ego_BP_DN, showCategory=10, fill = "p.adjust", title = "Biological Processes(down)") 

print(bpd)



#Combine UP and down results then make GO Plots with ggplot##
#Combine up and down top 10 plots
BP_comb_top10 <- rbind(BP_UP, BP_DN) 

slice_head(BP_comb_top10, n =20)
#save combined data frame
# write.csv(BP_comb_top10, file = paste0(out_path, "GO/GO_BP_top10_combined.csv"))

#import the file
# BP_comb_top10 <- read.csv(file = "GO_BP_top10_combined.csv")

# Create combined ordering column (most robust) for indexing purposes:
BP_comb_top10 <- BP_comb_top10 %>%
  mutate(order_combined = ifelse(FC == "UP", paste0("0_", Description), paste0("1_", Description))) %>% # Add prefix for correct sorting
  arrange(order_combined, -abs(Count)) # Sort by the combined column

# Make counts negative for DOWN-regulated pathways
BP_comb_top10$Count[BP_comb_top10$FC == "DOWN"] <- -BP_comb_top10$Count[BP_comb_top10$FC == "DOWN"]
BP_comb_top10 <- BP_comb_top10 %>% 
  group_by(FC) %>%
      arrange(desc(Count)) %>%
  ungroup()

# 4. Convert Description to factor with correct levels (CRUCIAL):
BP_comb_top10$Description <- factor(BP_comb_top10$Description, levels = rev(unique(BP_comb_top10$Description)))

#Create a ggplot of combined bargraph and save as pdf
#pdf(file = paste0(out_path, "GO/EGO_BP_combo_bar.pdf"), width = 8, height = 6, bg = "white")

GBP <- ggplot(BP_comb_top10, aes(x = Description, y = Count, fill = factor(FC, levels = c("UP", "DOWN")))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  scale_y_continuous(breaks = seq(min(BP_comb_top10$Count), max(BP_comb_top10$Count), by = 25), # Set breaks as needed
                     labels = abs(seq(min(BP_comb_top10$Count), max(BP_comb_top10$Count), by = 25))) + # Show absolute values on axis
  labs(title = "GO BP Pathways", x = NULL, y = "Gene Count", fill = "Fold Change") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.text.y = element_text(size = 11), plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = 0, color = "brown")  # Dashed line

print(GBP)


#Save the Figure object
#saveRDS(GBP, file = paste0(out_path, "GO/EGO_BP_combo_gg.RDS"))

#Repeat the same thing with other ORA analysis(Reactome, Wikipathways...)


#Create dotplots with the top gene lists
#import the corrected geneRatio top10 file
BP_comb_top10 <- read.csv(file = paste0(out_path,"GO/GO_BP_top10_combined.csv"), row.names = 1)

#Function to Split GR and return a ratio
Ratio_GR <- function(y){# Function to separate the gene names in the list, input geneList and Description
  z <- strsplit(y, "/")[[1]] #split and return the geneIds as a vector
  A <- as.numeric(z[1])
  B <- as.numeric(z[2])    
  return(A/B)
  }
BP_comb_top10$GeneRatio <- sapply(BP_comb_top10$GeneRatio, Ratio_GR)  #Ratio changed to numeric and properly

#save
#write.csv(BP_comb_top10, file = paste0(out_path, "GO/BP_top10_correctfinal.csv", row.names = F))

#read in the corrected file
# BP_comb_top10 <- read.csv(file = paste0(out_path,"GO/BP_top10_correctfinal.csv"), row.names = 1)

BP_comb_top10$FC <- factor(BP_comb_top10$FC) 
levels(BP_comb_top10$FC) <- c("UP", "DOWN")

BP_comb_top10 <- BP_comb_top10 %>%
  arrange(-GeneRatio)

#Create ggplot dotplot
#pdf(file = paste0(out_path, "GO/EGO_BP_combo_dot.pdf"), width = 8, height = 6, bg = "white")
dot_BP <- ggplot(BP_comb_top10, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust, shape = FC)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "Gene Ratio", y = NULL, size = "Gene Count", shape = "Fold Change", color = "p.adjust", title = "Top GO pathways") +
  theme_bw(base_size = 14, base_line_size = 0.2) +
  theme( plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 13)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))

print(dot_BP)


#Save the Figure object
#saveRDS(dot_BP, file = paste0(out_path, "GO/EGO_BP_combo_dot.RDS"))


##Heatplot with GO UP and Down results

#Read the enrichGO results
ego_BP_UP <- readRDS(paste0(out_path,"GO/ego_BP_UP.RDS"))

## convert gene ID to Symbol
edou <- setReadable(ego_BP_UP, edb, 'ENTREZID')

#Read the input geneList with fold change(GSEA input list)
df_RL_E <- read.csv(file = paste0(in_path,"Ranked_ENTREZ_GSEA.csv")) 

#Create a vector with logFC with names as ENTREZ from the gene List
#Function to create geneList with FC vector from unfiltered genes with FC list
geneList_vec_E <- function(x){
  logFC <- x$log2FC
  names(logFC) <- x$ENTREZID
  return(logFC)
}

#Convert genelist into symbol list from entrez and extract into a vector(all unfiltered genes)
#Check keytypes in org.Hs.eg.db- keytypes(org.Hs.eg.db)
geneList_vec_S <- function(x){
  EID <-as.character(x$ENTREZID) #Character vector with ENTREZID for mapping
  x$SYMBOL <-mapIds(edb, keys = EID, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first") 
  logFC <- x$log2FC
  names(logFC) <- x$SYMBOL
  logFC <- logFC[!is.na(names(logFC))]
  return(logFC)
}

geneList_E <- geneList_vec_E(df_RL_E) #Vector with log2FC and ENTREZID as names
geneList_S <- geneList_vec_S(df_RL_E) # Vector with Log2FC and Symbol as names

#write.csv(geneList_S, file = paste0(in_path, "Ranked_ENTREZ_GSEA_Symbol.csv"), row.names = F) #Save the ranked genelist with symbols
slice_head(edou@result, n=10)

#heatplot of up genes
p7 <- heatplot(edou, foldChange=geneList_S, symbol = "dot", showCategory=5) + viridis::scale_fill_viridis(name = "log2FC") + ggtitle("Upregulated") + theme(plot.title = element_text(hjust = 0.5, vjust = -3))


##Heatplot with GO Down results###
#Read the enrichGO results and convert it to gene symbols
ego_BP_DN <- readRDS(paste0(out_path,"GO/ego_BP_DN.RDS"))

## convert gene ID to Symbol
edox <- setReadable(ego_BP_DN, edb, 'ENTREZID')

slice_head(edox@result, n=10)
#Heatplot of down genes
p8 <- heatplot(edox, foldChange=geneList_S, symbol = "dot", showCategory=5) + viridis::scale_fill_viridis(name = "log2FC") + ggtitle("Downregulated") + theme(plot.title = element_text(hjust = 0.5, vjust = -3))


#Visualize and save both up and down heatplot figures
#pdf(file = paste0(out_path, "GO/EGO_BP_heatplot_combo.pdf"), width = 36, height = 14, bg = "white")
combo <- cowplot::plot_grid(p7,p8, ncol=1, labels=letters[1:2])
print(combo)



#cnet analysis - interaction network
#Import unfiltered genelist vector with logFC 
#Read the enrichGO results and convert it to gene symbols
#ego_BP_DN <- readRDS("ego_BP_DN.RDS")

## convert gene ID to Symbol
#edox <- setReadable(ego_BP_DN, edb, 'ENTREZID')

#Up genes
## categorySize can be scaled by 'pvalue' or 'geneNum'
p1 <- cnetplot(edou, categorySize="pvalue", foldChange=geneList_S) + ggtitle("Upregulated") + theme(plot.title = element_text(hjust = 0.5))
p2 <- cnetplot(edou, foldChange=geneList_S, igraph::layout_in_circle, showCategory = 2, colorEdge = TRUE) + ggtitle("Upregulated") + theme(plot.title = element_text(hjust = 0.5))

#Down genes
## categorySize can be scaled by 'pvalue' or 'geneNum'
p3 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList_S) + ggtitle("Downregulated") + theme(plot.title = element_text(hjust = 0.5))
p4 <- cnetplot(edox, foldChange=geneList_S, igraph::layout_in_circle, showCategory = 2, colorEdge = TRUE) + ggtitle("Downregulated") + theme(plot.title = element_text(hjust = 0.5))

#Save the combined up and down cnet plots
#pdf(file = paste0(out_path, "GO/EGO_BP_cnet_up or down.pdf"), width = 25, height = 30, bg = "white")
cn <- cowplot::plot_grid(p1, p3, p2, p4, ncol=2, labels=LETTERS[1:4])
print(cn)


#enrichment map
#Up genes
edou3 <- pairwise_termsim(edou)
pe1 <- emapplot(edou3) + theme(plot.margin = margin(2,2,2,2)) + ggtitle("Upregulated") + theme(plot.title = element_text(hjust = 0.5))

#Down genes
edox3 <- pairwise_termsim(edox)
pe2 <- emapplot(edox3) + theme(plot.margin = margin(2,2,2,2)) + ggtitle("Downregulated") + theme(plot.title = element_text(hjust = 0.5))

#pdf(file = paste0(out_path, "GO/EGO_BP_emap_up or down.pdf"), width = 16, height = 7, bg = "white")
ep <- cowplot::plot_grid(pe1, pe2, ncol=2, labels=letters[1:2], rel_widths = c(0.5, 0.5))

print(ep)

#Treeplot

#Up and down genes
edox2 <- pairwise_termsim(edox)
edou2 <- pairwise_termsim(edou)

pt1 <- treeplot(edox2, foldChange = geneList_S, hclust_method = "average") + ggtitle("Downregulated") + theme(plot.title = element_text(hjust = 0.3, vjust = -1))
pt2 <- treeplot(edou2, foldChange = geneList_S, hclust_method = "average") + ggtitle("Upregulated") + theme(plot.title = element_text(hjust = 0.3, vjust = -1))

# pdf(file = paste0(out_path, "GO/EGO_BP_tree_up or down_A.pdf"), width = 22, height = 18, bg = "white")
# tpa <- aplot::plot_list(pt2, pt1, tag_levels='a')
# print(tpa)
# 

#pdf(file = paste0(out_path, "GO/EGO_BP_tree_up or down.pdf"), width = 16, height = 18, bg = "white")
tp <- cowplot::plot_grid(pt2, pt1, ncol=1, labels=letters[1:2], rel_widths = c(0.5, 0.5))
print(tp)



sessionInfo() 

