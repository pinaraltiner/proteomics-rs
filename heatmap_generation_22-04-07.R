library("pheatmap")
library("stringr")
library("RColorBrewer")

# List of all spectrum matched result
## setwd does not recognize once! It keeps previous directory that's why you need to rerun the same line!
## Find a reasonable solution!

# set_dir <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/trial/07_04_22_count")
# set_dir <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/trial/07_04_22_count")


c <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/results_eyers_correct_list/")
c <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/results_eyers_correct_list/")
list_all_PEMs <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/results_eyers_correct_list/")
main_dir <- list.files(list_all_PEMs)

for(i in 1:length(main_dir)) {
  
  set_dir <- assign(paste0("set_dir", i),setwd((paste0(c,"/",main_dir[i],"/"))))
  set_dir <- assign(paste0("set_dir", i),setwd((paste0(c,"/",main_dir[i],"/"))))
  
  temp <- str_split_fixed(set_dir,"_",7)
  cont_count <- grep("count", temp)
  int_count <- grep("sum", temp)
  # List of all files in this directory
  
  d <- list.files(set_dir)
  
  # To save all figures into one pdf
  
  pdf(paste0(set_dir,"_correct_list_pheatmap_Minkowski_based_distance_clustering.pdf"),width=15, height = 15)
  
  # Adjustment of color palette of heatmaps
  
  paletteFunc <- colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929","#CC4C02"))
  palette     <- paletteFunc(100)
  
  
  ### Minkowski based distance clustering (both considered row and column for clustering)
  adjusted_hm <- function(df,break_values) {
    hm <- pheatmap(df,color=palette,border_color = "grey90",main=paste("HEATMAP of", d[i]),
                   fontsize = 10, fontsize_row = 12, fontsize_col = 12,
                   #breaks = seq(0, max(df), length.out = 5),
                   cellwidth = 18, cellheight = 10, clustering_distance_cols = dcols, angle_col = "45", width=1096, height = 609,cluster_cols  = TRUE)
    return(hm)
    
  }
  
  
  # Iteration of heatmap generation
  
  for(i in 1:length(d)) {
    
    # Adjustment of data to fit a suitable format for heatmap generation
    
    data <- assign(paste0("data", i),read.delim(paste0(d[i]),header=TRUE))
    data_sep_titles <- str_split_fixed(data$SPECTRUM, ";",7)
    data_refined_titles <- as.data.frame(paste(data_sep_titles[,1],data_sep_titles[,7], sep = "_"))
    tmp <- cbind(data_refined_titles, data[,-1])
    rownames(data) <- tmp[,1]
    df <- data[,-1]
    dcols = dist(t(df), method = "minkowski")
    
    # This condition eliminates undesirable format of data (in case of lack info of both spec_id and position)
    df_sum <- sum(df)
    
    if (sum(df) > 0 && max(df) > 1 && length(rownames(df))>1 && length(colnames(df)) > 2) { ### sum(df) > 0 && max(df) > 1 && 
      
      df <- as.matrix(df)
      
      # Condition discrimination btw count and intensity
      
      if (length(cont_count) == 0 && length(int_count)==1){
        
        df[df == 0] <- 1
        df <- log10(df) ##use only intensity summation
        adjusted_hm(df)
        
      }else{
        #print("correct file_name does not find! Please make sure file name either contain int_sum or count")
      }
      
      if (length(cont_count) == 1 && length(int_count)==0){
        
        adjusted_hm(as.matrix(df))    
      }
      
      
    }else {}
    
    if (sum(df) > 0 && max(df) <= 1 && length(rownames(df))>1 && length(colnames(df)) > 2){
      break_values = seq(0, max(df), length.out = 10)
      adjusted_hm(df,break_values)
      # pheatmap(df, color=palette,border_color = "grey90",main=paste("HEATMAP of", d[i]),
      #          fontsize = 10, fontsize_row = 12, fontsize_col = 12,
      #          breaks = seq(0, max(df), length.out = 10),
      #          cellwidth = 18, cellheight = 10, clustering_distance_cols = dcols, angle_col = "45", width=1096, height = 609,cluster_rows = TRUE)
    }
    
    # To save separate file
    
    #tiff(paste0("../figures/pheatmap_minkowski_dist_clust_method_row_col/",d[i],"rplot.tiff"),res = 100)
    #dev.off()
    
    rm(df)
  }
  dev.off()
  
  setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/results_eyers_correct_list/")
  setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/results_eyers_correct_list/")
  
  }



# List of all phospho evidence matrices 

#setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-03-30_monophospho_count")
#setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-03-30_monophospho_intensity_sum")

#set_dir <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-04-05_twophospho_intensity_sum")
#set_dir <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-04-05_twophospho_count")

#setwd("D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-02-24_monophospho_intensity_sum")
#setwd("D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-03-30_mono_phospho_count")

#setwd("D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-03-30_two_phospho_intensity_sum")
#setwd("D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-03-17_twophospho_count")

# For run all PEMs at the same time 
# Implement it to the rest of the code

# list_all_PEMs <- c("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-03-30_monophospho_count",
#                    "D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-03-30_monophospho_intensity_sum",
#                    "D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-02-24_monophospho_intensity_sum",
#                    "D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-03-30_mono_phospho_count",
#                    "D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-03-30_two_phospho_intensity_sum",
#                    "D:/dev/rust/proteomics-rs/data/rims_dataset_phospho_evidenice_matrices_list/2022-03-17_twophospho_count",
# "D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-04-05_twophospho_intensity_sum",
# "D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/2022-04-05_twophospho_count")
#                   
#                   
# for(i in 1:length(list_all_PEMs)) {
#   set_dir <- assign(paste0("set_dir", i),setwd(paste0(d[i])))
# }