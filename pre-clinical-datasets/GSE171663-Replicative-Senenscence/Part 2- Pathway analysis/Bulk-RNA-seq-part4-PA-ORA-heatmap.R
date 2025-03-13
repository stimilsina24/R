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

#Load libraries and set directory paths#
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
GO_BP <- paste0(getwd(),"/ORA/GO/") # Path for counts table
out_heat_GO <- paste0(getwd(),"/ORA-heatmap/GO/") # output path for heatmap results


#Heatmap of pathways using gene symbols#
#load counts table, Ensgene names as rownames
norm_counts <- read.csv(file= paste0(in_path,"norm_counts.csv")) %>% column_to_rownames(var = "X")

#import the top10 combined pathways 
BP_comb_top10 <- read.csv(file = paste0(GO_BP,"BP_top10_correctfinal.csv"))

#create functions to extract data from the enrichGO object, create a list of ensembl genes, and 

slice_head(BP_comb_top10, n=20)

#Extract ORA list from results data frame
extract_ORA <- function(x) { #Function to extract geneID and description from enrichGO as a list
  y <- x$geneID #extract geneIds from enrichGO results
  z <- x$Description #extract Description from enrichGO results
  return(list(geneIDs = y, Description = z))
}

#Split gene names from 
Split_ORA_l <- function(y,z){# Function to separate the gene names in the list, input geneList and Description
  x <- lapply(seq_along(y), function(i){ # y is the geneList from 1st step
    unlist(strsplit(y[i], "/")) #split and return the geneIds as a vector
  })
    names(x) <- z #name of each pathway/list
    return(x) #return the function as a list
}
#Convert ORA list from ENTREZ to gene symbols
Sym_ORA <- function(x) { #list of split genes
  y <- lapply(x, function(i){ # create a list
    z <- mapIds(org.Hs.eg.db, keys = i, keytype = "ENTREZID", column = "SYMBOL") #Use annotationDbi to return the ENSEMBL Ids from Split list
    return(z[!is.na(z)])
  })
  return(y)
}

#Annotate Counts table with gene symbol, remove rownames then convert symbol col to rownames
Counts_S <- norm_counts
Counts_S$symbol <- mapIds(org.Hs.eg.db, keys = rownames(norm_counts), keytype = "ENSEMBL", column = "SYMBOL")

rownames(Counts_S) <- NULL 

Counts_S <- Counts_S[!duplicated(Counts_S$symbol),] %>% na.omit(Counts_S) 

rownames(Counts_S) <- NULL 

Counts_S <- Counts_S %>% column_to_rownames(var = "symbol")

slice_head(Counts_S, n=20)

#Create gene list for heatmap for each pathway by matching ORA list with counts table
ORA_hm_l <- function(x,y) {#x is the ens enrichGO list(from previous step) and y is the counts list
  z <- lapply(x, function(i){
    k <- y[rownames(y) %in% i,] #k is the list of common genes between the counts table and ens list
    return(k)
  })
  return(z) # z is the whole list
}

#Heatmap function using gene symbols
ORA_hm_S <- function(x, top_n = 25){
  k <- lapply(seq_along(x), function(i){ #whole heatmap function
   # print(nrow(x[[i]]))
    if (nrow(x[[i]]) > 0) { #Make sure there are genes in each list
      # print(sapply(x[[i]], is.numeric))
       # Subset to top 20 genes
      if (nrow(x[[i]]) > top_n) {
        w <- x[[i]][1:top_n, ] # Select top 20 rows
      } else {
        w <- x[[i]] # If less than 20 rows, use the whole dataframe
      }

      y <- w[, sapply(w, is.numeric)]
      z <- scale(y)
      # print("Heatmap is being called")
      hm <- Heatmap(z, name = "z-score", #heatmap settings and input
                    row_labels = rownames(w),
                    column_labels = colnames(y),
                    column_title = names(x)[i],
                    column_title_gp = gpar(fontsize = 16), #Bigger title
                    cluster_rows = T, cluster_columns = T,
                    col = colorRamp2(c(min(z), 0, max(z)), c("blue", "white", "red")))      # Save heatmap as PNG
        pdf(file = paste0(out_heat_GO, "Symbol_", names(x)[i], ".pdf")) # Adjust width and height as needed
        draw(hm)
        dev.off()
    return(hm)
    } else { #make sure there are genes in the list or no problems, otherwise it will have no hm and thus return NULL
      return(NULL)
      }
      })
 print(k)
     }
  
 
#Perform all data extraction and heatmap functions
ORA1 <- extract_ORA(BP_comb_top10) #Extract the data from the combined top 10 pathways
ORA2 <- Split_ORA_l(ORA1$geneIDs, ORA1$Description) # Split each gene into a separate character
ORA3a <- Sym_ORA(ORA2) # Convert gene names to symbol
ORA4a <- ORA_hm_l(ORA3a, Counts_S) # Create the gene list with counts using the genes present in normalized counts table
ORA5a <- ORA_hm_S(ORA4a) # Visualize the heatmap
print(ORA5a) #Visualize the heatmap


#Heatmap of pathways with ensgene#
# #load counts table, Ensgene names as rownames
# norm_counts <- read.csv(file= paste0(in_path,"norm_counts.csv")) %>% column_to_rownames(var = "X")
# 
# #import the top10 combined pathways 
# BP_comb_top10 <- read.csv(file = paste0(GO_BP,"BP_top10_correctfinal.csv"))
# 
# #create functions to extract data from the enrichGO object, create a list of ensembl genes, and 
# 
# slice_head(BP_comb_top10, n=10)
# 
# #Extract ORA list from results data frame
# extract_ORA <- function(x) { #Function to extract geneID and description from enrichGO as a list
#   y <- x$geneID #extract geneIds from enrichGO results
#   z <- x$Description #extract Description from enrichGO results
#   return(list(geneIDs = y, Description = z))
# }
# 
# #Split gene names from 
# Split_ORA_l <- function(y,z){# Function to separate the gene names in the list, input geneList and Description
#   x <- lapply(seq_along(y), function(i){ # y is the geneList from 1st step
#     unlist(strsplit(y[i], "/")) #split and return the geneIds as a vector
#   })
#     names(x) <- z #name of each pathway/list
#     return(x) #return the function as a list
# }
# 
# #Convert ORA list from ENTREZ to ENSEMBL ID
# Ens_ORA <- function(x) { #list of split genes
#   y <- lapply(x, function(i){ # create a list
#     z <- mapIds(org.Hs.eg.db, keys = i, keytype = "ENTREZID", column = "ENSEMBL") #Use annotationDbi to return the ENSEMBL Ids from Split list
#     return(z[!is.na(z)])
#   })
#   return(y)
# }
# 
# #Create gene list for heatmap for each pathway by matching ORA list with counts table
# ORA_hm_l <- function(x,y) {#x is the ens enrichGO list(from previous step) and y is the counts list
#   z <- lapply(x, function(i){
#     k <- y[rownames(y) %in% i,] #k is the list of common genes between the counts table and ens list
#     return(k)
#   })
#   return(z) # z is the whole list
# }
# 
# #Function to create heatmap for all pathways.. Loop through each pathway and save as a new file
# ORA_hm <- function(x, top_n = 20){
#    k <- lapply(seq_along(x), function(i){ #whole heatmap function
#     if (nrow(x[[i]]) > 0) { #Make sure there are genes in each list
#       print(sapply(x[[i]], is.numeric))
#       print(nrow(x[[i]]))
#        # Subset to top 20 genes
#       if (nrow(x[[i]]) > top_n) {
#         w <- x[[i]][1:top_n, ] # Select top 20 rows
#       } else {
#         w <- x[[i]] # If less than 20 rows, use the whole dataframe
#       }
# 
#       y <- w[, sapply(w, is.numeric)]
#       z <- scale(y)
#       print("Heatmap is being called")
#       hm <- Heatmap(z, name = "z-score", #heatmap settings and input
#                     row_labels = rownames(w),
#                     column_labels = colnames(y),
#                     column_title = names(x)[i],
#                     column_title_gp = gpar(fontsize = 14), #Bigger title
#                     cluster_rows = T, cluster_columns = T,
#                     col = colorRamp2(c(min(z), 0, max(z)), c("blue", "white", "red")))      # Save heatmap as PNG
#         pdf(file = paste0(out_heat_GO, paste0(names(x)[i], ".pdf"))) # Save each file as a pathway name
#         draw(hm)
#         dev.off()
#     return(hm)
#     } else { #make sure there are genes in the list or no problems, otherwise it will have no hm and thus return NULL
#       return(NULL)
#       }
#       })
#  print(k)
#      }
#   
#   
# #Perform data extraction and heatmap
# ORA1 <- extract_ORA(BP_comb_top10) #Extract the data from the combined top 10 pathways
# ORA2 <- Split_ORA_l(ORA1$geneIDs, ORA1$Description) # Split each gene into a separate character
# ORA3 <- Ens_ORA(ORA2) # Convert gene names to ensembl ID
# ORA4 <- ORA_hm_l(ORA3, norm_counts) # Create the gene list with counts using the genes present in normalized counts table
# ORA5 <- ORA_hm(ORA4) # Visualize the heatmap


sessionInfo() 

