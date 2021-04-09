# This script is the Step1 for QRIS framework that for 
# extracting model coefficients and plotting

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))      
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(fmsb))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(qdap))
suppressPackageStartupMessages(library(Logolas))

# ==============================================================================
# Batch plot position importance
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

  cat("Process",file,"...\n")
  GLM_varImp_list <- readRDS(file)
  
  #only plot models with DNA shape
  shape_lists <- names(GLM_varImp_list)[grep("ext25bp_shape_14shapes",names(GLM_varImp_list))]
  
  for (shape_list in shape_lists){
    #shape_list <- shape_lists[1]
    model_Imp <- GLM_varImp_list[[shape_list]]
    
    #Build matrix for each model and assign values
    position_matrix <- matrix(0, nrow = length(all_shape) + 1, ncol = 51)
    rownames(position_matrix) <- c("Motif",all_shape)
    colnames(position_matrix) <- c(-25:25)
    
    st <- 1
    for (i in c("Motif",all_shape)){
      
      ext = 47; fill_range <- 3:49
      
      if (i %in% long_shape){ext = 48; fill_range <- 3:50}
      if (i == "Motif"){ ext = 19; fill_range <- 17:35}
      position_matrix[i,fill_range] <- model_Imp[st:(st+ext-1),]
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

  GLM_varImp_list <- readRDS(file)
  shape_list <- names(GLM_varImp_list)[grep("ext25bp_shape_14shapes",names(GLM_varImp_list))]
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