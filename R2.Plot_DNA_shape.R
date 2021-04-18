setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

Plot_DNA_shape <- function(flank_window = 25){
  files <- list.files("./", "ext100bp_DNAshape.rds$")
  
  for (file in files){

    cat("Processing ",file,"...\n")
    tmp <- readRDS(file)
    shapes <- c("MGW", "HelT","ProT","Roll","EP", "Stretch","Tilt","Buckle","Shear","Opening","Rise","Shift","Stagger","Slide")
    long_shape <- c("HelT","Roll","Tilt","Rise","Shift","Slide")
    short_shape <- setdiff(shapes, long_shape)
    
    NUMC = ncol(tmp[[shapes[1]]])
    NUMR = nrow(tmp[[shapes[1]]])
    
    mean_shape <- matrix(0, nrow = length(shapes), ncol = NUMC ) %>% as.data.frame()
    rownames(mean_shape) <- shapes
    
    for (shape in shapes){
      if (shape %in% short_shape){
        mean_shape[shape,] <- colMeans(tmp[[shape]], na.rm = T)
      } else {
        mean_shape[shape,1:(NUMC-1)] <- colMeans(tmp[[shape]], na.rm = T)
      }
    }
    rownames(mean_shape) <- paste0(rownames(mean_shape), "_" ,file)
    assign(paste0(file,"_mean_shape"), mean_shape)
  }
  
  v <- paste0(files,"_mean_shape")
  final_matrix <- do.call(rbind,mget(v))
  write.table(final_matrix,"shape_value_for_plot2.txt", sep = "\t",quote = F)

  final_matrix <- read.table("shape_value_for_plot2.txt")
  final_matrix <- final_matrix[grep("MGW|HelT|Roll|ProT",rownames(final_matrix)),]
  
  n <- (ncol(final_matrix)-1)/2 + 1
  final_matrix <- as.data.frame(t(final_matrix[,(n-flank_window):(n+flank_window)]))
  final_matrix$position <- -flank_window:flank_window
  
  final_matrix_melt <- reshape2::melt(final_matrix,c("position"))
  
  final_matrix_melt$variable <- gsub(".*_mean_shape.|_ext.*","",final_matrix_melt$variable)
  final_matrix_melt$shape <- gsub("_.*","",final_matrix_melt$variable)
  final_matrix_melt$shape <- factor(final_matrix_melt$shape, levels = c("MGW","HelT","ProT","Roll"))
  final_matrix_melt$sample <- sub(".*?_","",final_matrix_melt$variable)
  final_matrix_melt$value <- round(as.numeric(final_matrix_melt$value),2)
  final_matrix_melt$vt <- ""
  final_matrix_melt[grep("MoMLV_",final_matrix_melt$sample),]$vt <- "MoMLV"
  final_matrix_melt[grep("_MLV_",final_matrix_melt$sample),]$vt <- "MLV"
  final_matrix_melt[grep("HIV",final_matrix_melt$sample),]$vt <- "HIV-1"
  final_matrix_melt[grep("MMTV",final_matrix_melt$sample),]$vt <- "MMTV"
  final_matrix_melt[grep("HTLV",final_matrix_melt$sample),]$vt <- "HTLV-1"
  final_matrix_melt[grep("SIV",final_matrix_melt$sample),]$vt <- "SIV"
  final_matrix_melt[grep("danRer11",final_matrix_melt$sample),]$vt <- "danRer11"
  final_matrix_melt[grep("hg38",final_matrix_melt$sample),]$vt <- "hg38"
  final_matrix_melt[grep("mm10",final_matrix_melt$sample),]$vt <- "mm10"
  final_matrix_melt$vt <- factor(final_matrix_melt$vt, levels = unique(final_matrix_melt$vt)[c(4,5,6,7,2,9,1,3,8)])
  
  for (i in c("MGW","HelT","Roll","ProT")){
    pdf(paste0("DNA_shapes_around_DNA_invaders_",i,".pdf"), height = 3, width = 18)
    p <- ggplot(final_matrix_melt[grep(i,final_matrix_melt$shape),], aes(position, value, color = sample)) + 
      geom_smooth(method="loess", span=0.3, se=F, fullrange=F, level=0.95) + 
      facet_wrap(~vt, ncol = 9) +
      theme_bw() +
      scale_x_continuous(breaks=seq(-25,25,25)) + 
      labs(x="Relative position to insertion (bp)", y = "Mean DNA shape values") +
      scale_color_manual(values = c("gray","#1F78B4","gray","#FB9A99","#E7298A","#E31A1C",
                                    "#B2DF8A","#33A02C","#CAB2D6","#6A3D9A","#FF7F00","gray","#B15928")) +
      theme(
        plot.title = element_text(color="black", size=20, face="bold"),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.text.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        axis.text.y = element_text(color="black", size=20, face="bold"),
        legend.title = element_text(color="black", size=14, face="bold"), 
        legend.text = element_text(color="black", size=12, face="bold"),
        strip.text.x = element_text(size = 20, face="bold")
      ) 
    plot(p)
    dev.off()
  }
}

for (i in c(25,50,100,200,300,498)[6]){
  Plot_DNA_shape(path, flank_window = i)
}

#Do statistics
stat <- final_matrix
colnames(stat) <- gsub(".*_mean_shape.|_ext.*","",colnames(stat))
stat <- stat[,grep("MGW|HelT|ProT|Roll",colnames(stat))]

fcv <- cbind(combn(colnames(stat)[grep("MGW",colnames(stat))],2)[,c(1,24,25,26,27,28,29,30,31,78)],
      combn(colnames(stat)[grep("HelT",colnames(stat))],2)[,c(1,24,25,26,27,28,29,30,31,78)],
      combn(colnames(stat)[grep("ProT",colnames(stat))],2)[,c(1,24,25,26,27,28,29,30,31,78)],
      combn(colnames(stat)[grep("Roll",colnames(stat))],2)[,c(1,24,25,26,27,28,29,30,31,78)])
sink("test.txt")
for (i in 1:ncol(fcv)){
  x = fcv[,i][1]
  y = fcv[,i][2]
  cc = wilcox.test(as.numeric(stat[,x]), as.numeric(stat[,y]))
  if(cc$p.value > 0.05){cat(x,y,"-->",cc$p.value,"\n")}
  
}
sink()
