# This script is used for plotting GLM results

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))      
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(fmsb))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(qdap))
suppressPackageStartupMessages(library(Logolas))

# ==============================================================================
# 0. Set up utilities
# ==============================================================================
# Set up DNA shape featrues
all_shape <- c("MGW", "HelT","ProT","Roll","EP",  
               "Stretch","Tilt","Buckle","Shear","Opening","Rise","Shift","Stagger","Slide")

long_shape <- c("HelT","Roll","Tilt","Rise","Shift","Slide")
short_shape <- setdiff(all_shape, long_shape)

# ==============================================================================
# 1. plot gradient predictors
# ==============================================================================
plot_gradient_predictors <- function(files){
  
  plot_mat <- as.data.frame(matrix(0, ncol = length(files), nrow = 5))
  colnames(plot_mat) <- files
  
  for (file in files){
    res <- read.table(file)
    colnames(res) <- c("AUROC","Sensitivity","Specificity","AUPRC","Precision",
                       "Recall","F-score","Accurancy","plan")
    plot_mat[,file] <- res$Accurancy
  }
  
  plot_mat <- round(plot_mat,4)
  rownames(plot_mat) <- c("motif",gsub(".*shuffled_|1-", "", res$plan)[2:5])
  rownames(plot_mat) <- sub("shape_","",sub("_shape_"," + ",rownames(plot_mat)))
  rownames(plot_mat) <- sub("MGW,ProT,Roll,HelT","4shapes",rownames(plot_mat))
  plot_mat$shape <- rownames(plot_mat)
  
  colnames(plot_mat) <- gsub("_trainning_ext25bp|_matrice.*", "", colnames(plot_mat))
  #colnames(plot_mat) <- mgsub(l_pattern, l_replacement, colnames(plot_mat))
  
  write.table(plot_mat,"Metrice_data_matrix_nakedDNA.txt",sep = "\t", quote = F, row.names = T, col.names = T)
  
  #Select predictors for plot
  plot_mat_final <- plot_mat
  plot_mat_melt <- reshape2::melt(plot_mat_final,c("shape"))
  colnames(plot_mat_melt) <- c("shape", "sample","Accurancy")
  
  plot_mat_melt <- plot_mat_melt[plot_mat_melt$shape %in% 
                                   c("motif","4shapes","14shapes","motif + 14shapes","motif + 4shapes"),]
  #seperate original and shuffled samples
  plot_mat_melt_true <- plot_mat_melt[-grep("_shuffled", plot_mat_melt$sample),]
  plot_mat_melt_true$species <- gsub("_.*","",plot_mat_melt_true$sample)
  
  #replace sample names for better visualization
  meta_patern <- as.vector(unique(plot_mat_melt_true$sample))
  meta_replacement <- c("HIV-1 HEK293T","HIV-1 PBMC","HIV-1 T cell","HTLV-1 hemato.","HTLV-1 T cell",
                        "MLV HepG2","MLV K562","MMTV","MoMLV","SIV")
  plot_mat_melt_true$sample <- mgsub(meta_patern, meta_replacement,plot_mat_melt_true$sample)
  plot_mat_melt_true$sample <- factor(plot_mat_melt_true$sample, 
                                          levels = unique(plot_mat_melt_true$sample)[c(1:7,10,9,8)])
  
  plot_mat_melt_shuffled <- plot_mat_melt[grep("_shuffled", plot_mat_melt$sample),]
  plot_mat_melt_shuffled$sample <- sub("_shuffled", "", plot_mat_melt_shuffled$sample)
  plot_mat_melt_shuffled$species <- gsub("_.*","",plot_mat_melt_shuffled$sample)
  plot_mat_melt_shuffled$sample <- mgsub(meta_patern, meta_replacement,plot_mat_melt_shuffled$sample)
  plot_mat_melt_shuffled$sample <- factor(plot_mat_melt_shuffled$sample, 
                                          levels = unique(plot_mat_melt_shuffled$sample)[c(1:7,10,9,8)])
  
  pdf(paste0("GLM performance on different predictors.pdf"), height = 5, width = 7)
  ggplot(plot_mat_melt_true[-grep("^14shapes",plot_mat_melt_true$shape),], aes(sample, Accurancy, color=shape)) + 
    geom_point(size = 2.5, alpha = 0.9) + 
    geom_point(data = plot_mat_melt_shuffled[grep(" 14shapes",plot_mat_melt_shuffled$shape),], 
               aes(sample, Accurancy), color="grey", size = 2.5, alpha = 0.5) +
    theme_bw() + ylim(0.5, 0.75) +
    scale_colour_manual(values = brewer.pal(n = 8, name = "Dark2")[c(2,3,4,1)]) + 
    labs(x="", y = "Accuracy (ACC)") +
    #facet_grid(~ group, scales = "free_x", space='free') +
    theme(
      plot.title = element_text(color="black", size=20, face="bold"),
      axis.title.x = element_text(color="black", size=18,face="bold"),
      axis.text.x = element_text(color="black", size=11,face="bold", angle=70, hjust=1),
      #axis.text.x = element_blank(),
      axis.title.y = element_text(color="black", size=18, face="bold"),
      axis.text.y = element_text(color="black", size=20,face="bold"),
      legend.title = element_blank(), 
      legend.text = element_text(color="black", size=12, face="bold")
    )
  #plot(p)
  dev.off()
}

files <- list.files("./",pattern = "matrice")
plot_gradient_predictors(files)

# ==============================================================================
# 2. Plot model performance (radar plot)
# ==============================================================================
plot_radar <- function(files){
  files <- files[-grep("shuffled",files)][-c(1,4)]
  pdf("predictor_gradient_matrices_radar.pdf", width = 16, height = 8)
  layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE), widths=rep(3,4), heights=rep(3,4))
  
  for (file in files){
    #file <- files[6]
    res <- read.table(file)
    file <- gsub("_trainning_ext25bp|_matrice.*","",file)
    # file <- mgsub(meta_patern, meta_replacement,file)
    # file <- sub(" .*","",file)
    
    colnames(res) <- c("AUROC","Sensitivity","Specificity","AUPRC",
                       "Precision","Recall","F-score","accuracy","plan")
    rownames(res) <- c("motif",gsub("bp_|1-","",str_extract(res$plan,"bp_.*$"))[2:5])
    rownames(res) <- sub("_shape","",sub("_shape_"," + ",rownames(res)))
    rownames(res) <- sub("MGW,ProT,Roll,Tilt","4shapes",rownames(res))
    write.table(res,paste0(file,".txt"),sep = "\t",row.names = T,col.names = T,quote = F)
    
    res <- res[c(1,2,3,5),c(1,4,5,6,7)]
    res <- rbind(rep(0.8,7) , rep(0.5,7) , res)
    
    colors_border <- brewer.pal(n = 8, name = "Dark2")[c(3,2,1,4)]
    #colors_in <- alpha(colors_border, 0.3)
    radarchart(res, axistype=1, pcol = colors_border, plwd=2, plty=1, title = file, cglcol="grey", 
               cglty=2, axislabcol="grey", caxislabels=seq(0.6, 1, length = 5), cglwd=1,vlcex=1.3)
  }
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x=0.6, y=1.2, legend = rownames(res[-c(1,2),]), bty = "n", pch=20, 
         col=colors_border, text.col = "black", cex=2, pt.cex=3)
  
  dev.off()
}

files <- list.files("./",pattern = "matrice.*txt")
plot_radar(files)
