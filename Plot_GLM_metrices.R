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

# ==============================================================================
# 3. Batch plot position importance
# ==============================================================================
plot_comprehensive_position <- function(position_matrix_scaled, plot_range = 23, plot_abs = F){
  #select range for plot
  plot_range <- plot_range
  position_matrix_scaled <- as.data.frame(position_matrix_scaled[,-c(1:(25-plot_range),(25+plot_range+2):51)])
  
  #The predictor_weight is calculated by central +-9 bp
  center <- (ncol(position_matrix_scaled) - 1 )/2
  position_matrix_scaled$predictor_weight <- rowSums(abs(position_matrix_scaled[,c((center-9):(center+9))]))
  
  position_matrix_scaled <- rbind(position_matrix_scaled, colSums(abs(position_matrix_scaled)))
  rownames(position_matrix_scaled)[16] <- "position_weight"
  
  #Sort the matrix according to predictor_weight
  ord <- order(position_matrix_scaled[-nrow(position_matrix_scaled),]$predictor_weight, decreasing = T)
  position_matrix_scaled <- position_matrix_scaled[c(ord,16),]
  
  a <- position_matrix_scaled[nrow(position_matrix_scaled), -ncol(position_matrix_scaled)]
  ha1 = HeatmapAnnotation(Position_weight = as.numeric(a),
                          annotation_name_side = "left",
                          show_annotation_name = T,
                          annotation_legend_param = list(Position_weight = list(title = "Position\nweight")),
                          annotation_name_gp = gpar(fontsize = 14, fontface = "bold", col = "green4"),
                          show_legend = T,col = list(Position_weight = colorRamp2(c(0,2), c("white","green4"))),
                          #show_legend = T,col = list(Position_weight = colorRamp2(c(0,max(a)), c("white","green4"))),
                          gp = gpar(col = "grey",lwd = 1.5)
  )
  if (plot_abs){
    to_plot_mat <- abs(position_matrix_scaled[-nrow(position_matrix_scaled),-ncol(position_matrix_scaled)])
    color <- colorRamp2(c(0,1), c("white","red"))} 
  else {
    to_plot_mat <- position_matrix_scaled[-nrow(position_matrix_scaled),-ncol(position_matrix_scaled)]
    color <- colorRamp2(c(-1,0,1), c("blue","white","red"))
  }
  
  p1 <- Heatmap(to_plot_mat,
                col = color,
                name = "Relative\nimportance",
                row_names_side = "left",row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                column_names_gp = gpar(fontsize = 11,fontface = "bold"),
                column_names_rot = 0, column_names_centered = T,
                rect_gp = gpar(col = "grey", lwd = 1.5),
                cluster_columns = F,cluster_rows = F,
                column_title = "Relative position to integration center (bp)",
                column_title_side = "bottom", column_title_gp = gpar(fontsize = 14,fontface = "bold"),
                bottom_annotation = ha1
  )
  #p1
  b <- as.data.frame(position_matrix_scaled[-nrow(position_matrix_scaled), ncol(position_matrix_scaled)])
  p2 <- Heatmap(b,
                column_title = "Predictor\nweight",
                column_title_gp = gpar(fontsize = 12, fontface = "bold", col = "orange"),
                cluster_columns = F, cluster_rows = F,
                show_heatmap_legend = F, rect_gp = gpar(type = "none"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                  grid.circle(x = x, y = y, r = b[i, j] * min(unit.c(width, height)), 
                              gp = gpar(fill = "orange",col = NA))}
  )
  draw(p1 + p2)
}

# start running all files
files <- list.files("./","bp_varImp_list.rds$")
for (file in files){
  # file <- files[5]
  cat("Process",file,"...\n")
  GLM_varImp_list <- readRDS(file)
  
  #only plot models with DNA shape
  shape_lists <- names(GLM_varImp_list)[grep("ext25bp_motif$|ext25bp_shape_14shapes", names(GLM_varImp_list))]

  model_Imp1 <- GLM_varImp_list[[shape_lists[1]]]
  model_Imp2 <- GLM_varImp_list[[shape_lists[2]]]
  position_matrix <- matrix(0, nrow = length(all_shape) + 1, ncol = 51)
  rownames(position_matrix) <- c("Motif",all_shape)
  colnames(position_matrix) <- c(-25:25)
  
  position_matrix["Motif",] <- model_Imp1$Overall
  st <- 1
  for (i in c(all_shape)){
    ext = 47; fill_range <- 3:49
    if (i %in% long_shape){ext = 48; fill_range <- 3:50}
    position_matrix[i,fill_range] <- model_Imp2[st:(st+ext-1),]
    st <- st + ext
  }
  write.table(position_matrix, paste0(file,"_positionWeight.txt"), col.names = T, row.names = T,quote = F, sep = "\t")
  #Normalize values to the biggest one
  position_matrix_scaled <- position_matrix/max(abs(position_matrix))
  #plot_comprehensive_position(position_matrix_scaled, plot_range = 9)

  plot_range = 23; plot_abs = F
  pdf(paste(file,"position_matrix_scale.pdf", sep = "_"), width = 10, height = 6)
  plot_comprehensive_position(position_matrix_scaled, plot_range = plot_range, plot_abs = plot_abs)
  dev.off()
}

# ==============================================================================
# 4. Combine all values for a consensus matrix
# ==============================================================================
combined_Imp_list <- list()
files <- list.files("./","bp_varImp_list.rds")
for (file in files){
  # file <- "ce_input_shuf20000_ext25bp_GLM_varImp_list.rds" #for test only
  GLM_varImp_list <- readRDS(file)
  shape_list <- names(GLM_varImp_list)[grep("14shapes_shape_motif",names(GLM_varImp_list))]
  combined_Imp_list[[shape_list]] <- GLM_varImp_list[[shape_list]]
}

combined_Imp <- matrix(0, ncol = length(names(combined_Imp_list)), 
                       nrow = nrow(combined_Imp_list[[1]])) %>% as.data.frame()
colnames(combined_Imp) <- names(combined_Imp_list)

for (Imp in names(combined_Imp_list)){combined_Imp[,Imp] <- combined_Imp_list[[Imp]]}

combined_Imp <- scale(combined_Imp,center = F,scale = T)
concensus <- apply(combined_Imp, 1, mean) %>% as.matrix()

position_matrix <- matrix(0, nrow = length(all_shape) + 1, ncol = 51)
rownames(position_matrix) <- c("Motif",all_shape)
colnames(position_matrix) <- c(-25:25)

st <- 1
for (i in c("Motif",all_shape)){
  ext = 47; fill_range <- 3:49
  if (i %in% long_shape){ext = 48; fill_range <- 3:50}
  if (i == "Motif"){ ext = 19; fill_range <- 17:35}
  position_matrix[i,fill_range] <- concensus[st:(st+ext-1),]
  st <- st + ext
}

#Center values
position_matrix_scaled <- position_matrix/max(abs(position_matrix))
write.table(position_matrix_scaled,"Consensus_position_matrix_scaled.txt",sep = "\t",quote = F)

plot_range = 23; plot_abs = F
pdf(paste(scale_bg,plot_range,plot_abs,"position_matrix_consensus.pdf", sep = "_"), 
    width = unit(14, "cm"), height = unit(5, "cm"))
plot_comprehensive_position(position_matrix_scaled, plot_range = plot_range, plot_abs = plot_abs)
dev.off()

# ==============================================================================
# 5. Shuffled DNA shapes
# ==============================================================================
files <- list.files("./",pattern = "matrice_2020")
plot_mat <- as.data.frame(matrix(0, ncol = length(files), nrow = 2))
colnames(plot_mat) <- files

for (file in files){
  res <- read.table(paste0("./",file))
  colnames(res) <- c("AUROC","Sensitivity","Specificity","AUPRC","Precision",
                     "Recall","F-score","Accurancy","plan")
  plot_mat[,file] <- res$Accurancy
}

plot_mat <- round(plot_mat,4)
rownames(plot_mat) <- gsub(".*ext25bp_", "", res$plan)
rownames(plot_mat) <- sub("shape_","",sub("_shape_"," + ",rownames(plot_mat)))
rownames(plot_mat) <- sub("MGW,ProT,Roll,HelT","4shapes",rownames(plot_mat))
plot_mat$shape <- rownames(plot_mat)

colnames(plot_mat) <- gsub("_trainning_ext25bp|_matrice.*|_", "", colnames(plot_mat))
#colnames(plot_mat) <- mgsub(l_pattern, l_replacement, colnames(plot_mat))


ss <- read.table("../Metrice_data_matrix.txt",sep = "\t")[c(4,5),]
iii <- intersect(colnames(plot_mat),colnames(ss))[-32]
final <- t(rbind(plot_mat[,iii],ss[,iii])) %>% as.data.frame()
final$species <- gsub("_.*","",rownames(final))
write.table(plot_mat,"Metrice_data_matrix_shuffledshape.txt",sep = "\t", quote = F, row.names = T, col.names = T)


pdf(paste0("GLM performance on shuffled shapes.pdf"), height = 4, width =5)
final
ggplot(final[-nrow(final),], aes(as.numeric(`14shapes1`), as.numeric(`14shapes`), color=species)) + 
  geom_point(size = 2) + 
  theme_bw() + ylim(0.5, 0.8) + xlim(0.5, 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype="solid", size=0.8)+
  scale_colour_manual(values = brewer.pal(n = 8, name = "Dark2")) + 
  labs(x="Accuracy (14 true shapes)", y = "Accuracy (14 shuffled shapes)") +
  theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=14,face="bold"),
    axis.text.x = element_text(color="black", size=11,face="bold"),
    #axis.text.x = element_blank(),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.y = element_text(color="black", size=11,face="bold"),
    legend.title = element_text(color="black", size=14, face="bold"), 
    legend.text = element_text(color="black", size=12, face="bold")
  )
dev.off()

# ==============================================================================
# 6. Cross-species validation
# ==============================================================================
plot_ML_matrix <- function(mat=AUC_matrix, scheme="all", metric = "AUC"){
  plot_mat <- mat
  if(scheme == "naked"){
    plot_ids <- grep("input", colnames(mat))
    plot_mat <- mat[plot_ids, plot_ids]
  } else if(scheme == "chromatin"){
    plot_ids <- grep("input", colnames(mat))
    plot_mat <- mat[-plot_ids, -plot_ids]
  }
  
  colnames(plot_mat) <- gsub("_|_trainning_ext25bp_splits", "", colnames(plot_mat))
  rownames(plot_mat) <- gsub("_|_trainning_ext25bp_splits", "", rownames(plot_mat))
  plot_mat$species <- sub("_.*","",colnames(plot_mat))
  
  remove <- grep("MA9|ZHBTC4",colnames(plot_mat))
  plot_mat <- plot_mat[-remove,-remove]
  
  ha = rowAnnotation(species = plot_mat$species,
                     gp = gpar(col = "#FFFFFF"), 
                     col = list(
                       species = structure(brewer.pal(n = 8, name = "Dark2"), names = unique(plot_mat$species))
                     )
  )
  
  p <- Heatmap(plot_mat[,-ncol(plot_mat)],
               show_row_names = F, show_column_names = F,
               row_names_gp = gpar(fontsize = 14),
               column_names_gp = gpar(fontsize = 14),
               #column_title = paste0(metric, " of cross-test"),
               cluster_columns = F, cluster_rows = F,
               col = colorRamp2(c(0.5,0.85), c("white","black")),
               rect_gp = gpar(col = "black", lwd = 1),
               right_annotation = ha,
               name = metric
  )
  pdf(paste0(metric, " of cross-test", "_", scheme,".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

AUC_matrix <- read.csv("CorssValidationForSpeciesTrainned_AUC_matrix.csv", row.names = 1)
ROC_matrix <- read.csv("CorssValidationForSpeciesTrainned_ROC_matrix.csv", row.names = 1)
ACC_matrix <- read.csv("CorssValidationForSpeciesTrainned_ACC_matrix.csv", row.names = 1)
for (i in c("naked","chromatin","all")[1]){
  plot_ML_matrix(mat=AUC_matrix, scheme=i, metric = "AUPRC")
  plot_ML_matrix(mat=ROC_matrix, scheme=i, metric = "AUROC")
  plot_ML_matrix(mat=ACC_matrix, scheme=i, metric = "ACC")
}
# ==============================================================================
# 5. Structure motif
# ==============================================================================
position_matrix_scaled <- read.table("Consensus_position_matrix_scaled.txt", sep = "\t")
colnames(position_matrix_scaled) <- -25:25

position_matrix_scaled_logo <- position_matrix_scaled[-which(rownames(position_matrix_scaled) == "Motif"),]
rownames(position_matrix_scaled_logo) <- c("M","O","B","R","E","S","H","I","T","F","A","L","P","V")

pdf("position_matrix_scaled_logo.pdf",width = unit(14, "cm"), height = unit(3, "cm"))
logomaker(abs(position_matrix_scaled_logo), type = "Logo")
dev.off()

