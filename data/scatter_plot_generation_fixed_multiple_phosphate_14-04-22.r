library("stringr")
library("pheatmap")
library("dplyr")
library("ggplot2")
library("ggalluvial")
a<- list_all_PEMs <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/final_results/08_04_22_twophospho_int_count/")
#08_04_22_twophospho_int_count
aa <- list.files(a)

list_all_PEMs <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/final_results/")
list_all_PEMs <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/final_results/")
main_dir <- list.files(list_all_PEMs)

intentisty_mod <- c("int_count","int_sum")

tp_final  <- NULL
#tp_final_sum  <- NULL
fp_final <- NULL
#fp_final_res_sum <- NULL

# Change index of main_dir for other phospho_types 
  # for monophospo main_dir[1,2]
  # for twophospho main_dir[5,6]
  # for threephospho main_dir[3,4], etc.

for(i in 1:length(aa)) {
  
  set_dir <- assign(paste0("set_dir", i), setwd((paste0(list_all_PEMs,"/",main_dir[5],"/"))))
  set_dir <- assign(paste0("set_dir", i), setwd((paste0(list_all_PEMs,"/",main_dir[5],"/"))))
  list_all_pep <- list.files(set_dir)
  # Read the tsv
  assign(paste0(intentisty_mod[1]),read.delim(paste0(list_all_pep[i]),header=TRUE))
  
  # Remove spectrum ID
  int_count_wo_spec <- int_count[,-1]
  
  # Do a summation column-wise
  data_int_count_agg <- t(as.data.frame(colSums(int_count_wo_spec)))
  
  # Indicate its peptide info and type of phospho evidence calculation (either count or sum)
  #rownames(data_int_count_agg) <- paste(list_all_pep[1], sep = "_")
  
  rownames(data_int_count_agg) <- paste(intentisty_mod[1], list_all_pep[i], sep = "-")
  
  # # # Change directory to intensity summation
  set_dir <- assign(paste0("set_dir", i), setwd((paste0(list_all_PEMs,"/",main_dir[6],"/"))))
  set_dir <- assign(paste0("set_dir", i), setwd((paste0(list_all_PEMs,"/",main_dir[6],"/"))))
  
  # Read the tsv
  assign(paste0(intentisty_mod[2]),read.delim(paste0(list_all_pep[i]),header=TRUE))
  
  # Remove spectrum ID
  int_sum_wo_spec <- int_sum[,-1]
  
  # Do a summation column-wise
  data_int_sum_agg <- t(as.data.frame(colSums(int_sum_wo_spec)))
  
  # Indicate its peptide info and type of phospho evidence calculation (either count or sum)
  #rownames(data_int_sum_agg) <- paste(list_all_pep[i], sep = "_")
  
  rownames(data_int_sum_agg) <- paste(intentisty_mod[2], list_all_pep[i], sep = "-")
  
  int_both <- rbind(data_int_count_agg,data_int_sum_agg)
  
  # It is still true for monophosphate calculation 
  #true_pos <- data_sep_titles %>% str_match_all("[0-9]+") %>% unlist %>% as.character
  
  # This is for two or more different phosphorlation states
  # Extracting position information from the file name
  list_all_pep_pos <- str_extract_all(list_all_pep[i], ' \\(([A-Z]\\d+)\\)')
  pep_pos <- str_extract_all(list_all_pep_pos, "\\d+")
  
  # Do unlist
  pep_pos1 <- unlist(pep_pos)
  
  if (length(pep_pos1) > 1){
    
    # Add same separator like column name of the tsv file to match between them
    true_pos <- paste(pep_pos1[1], pep_pos1[2], sep = "..")
    
  }else{
    
    true_pos <- paste(pep_pos1[1])
  }
  
  # Search the correct positions in the all PEMs
  tmp <- as.data.frame(int_both) %>% select(contains(true_pos))
  ftmp <- as.data.frame(int_both) %>% select(!contains(true_pos))
  
  #tmp1 <- as.data.frame(int_both) %>% select(contains(true_pos))
  #ftmp1 <- as.data.frame(int_both) %>% select(!contains(true_pos))
  
  # Collect data for each peptide
  
  # The huge difference btw intensity sum and count, they stored separate object 
  tp_final <- bind_rows(as.data.frame(tp_final), as.data.frame(tmp))
  
  #fp_final_res_count <- bind_rows(as.data.frame(fp_final_res_count), as.data.frame(ftmp))
  
  fp_final <- bind_rows(as.data.frame(fp_final), as.data.frame(ftmp))
  
  #fp_final_res_sum <- bind_rows(as.data.frame(fp_final_res_sum), as.data.frame(ftmp1))
  
  
}

setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/")
setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/final_results/")

# Write output as tsv
setwd("D:/dev/rust/proteomics-rs/data/benchmarks/phospho_evidence_matrices_list/final_results")

write.table(tp_final, file = "TP_pos_values_Eyers_twophosphate.tsv", sep = "\t")
write.table(fp_final, file = "FP_pos_values_Eyers_twophosphate.tsv", sep = "\t")

########################## DATA TRANSFORMATION ###############################
  

# Convert all NA to 1
tp_final[is.na(tp_final)] <- 1
tp_final[tp_final == 0] <- 1

fp_final[is.na(fp_final)] <- 1
fp_final[fp_final == 0] <- 1


# Log transformation ( do this before int_xx col addition )
tp_final_log <- log10(tp_final)
fp_final_log <- log10(fp_final)

# To reduce dimension take the row sum of individual Ps and NPs
tp_final_log_Ps <-as.data.frame(tp_final_log) %>% select(!contains("NP"))
tp_final_log_P_rowsum1 <- as.data.frame(rowSums(tp_final_log_Ps))

tp_final_log_NPs <-as.data.frame(tp_final_log) %>% select(matches("NP"))
tp_final_log_NP_rowsum1 <- as.data.frame(rowSums(tp_final_log_NPs))

tp_final_log_rowsum1 <- cbind(tp_final_log_P_rowsum1,tp_final_log_NP_rowsum1)
colnames(tp_final_log_rowsum1) <- c("P","NP")

fp_final_log_Ps <-as.data.frame(fp_final_log) %>% select(!matches("NP"))
fp_final_log_P_rowsum1 <- as.data.frame(rowSums(fp_final_log_Ps))

fp_final_log_NPs <-as.data.frame(fp_final_log) %>% select(matches("NP"))
fp_final_log_NP_rowsum1 <- as.data.frame(rowSums(fp_final_log_NPs))

fp_final_log_rowsum1 <- cbind(fp_final_log_P_rowsum1,fp_final_log_NP_rowsum1)
colnames(fp_final_log_rowsum1) <- c("P","NP")

# Rowname extraction
rname_tp_final <- as.data.frame(rownames(tp_final))

# Extract int_xx info and add it to the final result as a new column
sep_int_mod <- as.data.frame(sapply(strsplit(rname_tp_final$`rownames(tp_final)`, "-"), "[", 1))

tp_final_int_mod <- cbind(tp_final_log_rowsum1, sep_int_mod,"true_positive")

fp_final_int_mod <- cbind(fp_final_log_rowsum1, sep_int_mod,"false_positive")

# Change the colname of newly added col
last_col <- length(tp_final_int_mod)
last_col1 <- length(fp_final_int_mod)

names(tp_final_int_mod)[last_col-1] <- paste("intensity_mod")
names(fp_final_int_mod)[last_col1-1] <- paste("intensity_mod")

names(tp_final_int_mod)[last_col] <- paste("type")
names(fp_final_int_mod)[last_col1] <- paste("type")

#######################################
############ Mono Phospho ############
# Re-run the code by changing index of main_dir 
monophosphate_all <- rbind(tp_final_int_mod,fp_final_int_mod)
monophosphate_all_p_state <- cbind(monophosphate_all,"monophosphate")
names(monophosphate_all_p_state)[last_col+1] <- paste("phospho_state")

#monophosphate_all_p_state contains 61 int_sum + 61 int_count for each of them has also TP and FP =244

# Before merging btw mono and two phosphate, we collect FP and TP in the separate column 
# It would make easier to generate regression model btw these features.
dimm <- dim(monophosphate_all_p_state)
only_tp_mono <- filter(monophosphate_all_p_state, between(row_number(),1,dimm[1]/2))
only_fp_mono <- filter(monophosphate_all_p_state, between(row_number(),dimm[1]/2 +1 ,dimm[1]))

comb_tp_fp_mono <- cbind(only_tp_mono,only_fp_mono)
comb_tp_fp_mono <- comb_tp_fp_mono[,c(1:7,9)]
colnames(comb_tp_fp_mono) <- c("TP","TNP","intensity_mod","type_TP","phospho_state","FP","FNP","type_FP")

#comb_tp_fp_mono has the same data as monophosphate_all_p_state but fp (int_sum and int_count) was bound to col-wise. 

# You have to clear the environment to eliminate any overwrite issue
rm(list= ls()[!(ls() %in% c('monophosphate_all_p_state','comb_tp_fp_mono'))])


#######################################
############ Two Phospho ############
# Re-run the code by changing index of main_dir 
twophosphate_all <- rbind(tp_final_int_mod,fp_final_int_mod)
twophosphate_all_p_state <- cbind(twophosphate_all,"twophosphate")
names(twophosphate_all_p_state)[last_col+1] <- paste("phospho_state")

# Before merging btw mono and two phosphate, we collect FP and TP in the separate column 
# It would make easier to generate regression model btw these features.
dimm <- dim(twophosphate_all_p_state)
only_tp_two <- filter(twophosphate_all_p_state, between(row_number(),1,dimm[1]/2))
only_fp_two <- filter(twophosphate_all_p_state, between(row_number(),dimm[1]/2 +1 ,dimm[1]))

comb_tp_fp_two <- cbind(only_tp_two,only_fp_two)
comb_tp_fp_two <- comb_tp_fp_two[,c(1:7,9)]
colnames(comb_tp_fp_two) <- c("TP","TNP","intensity_mod","type_TP","phospho_state","FP","FNP","type_FP")
# Merge mono and two
mono_two_all_states <- rbind(monophosphate_all_p_state,twophosphate_all_p_state)

# Remove the rest of the object except for main working ones
rm(list= ls()[!(ls() %in% c('monophosphate_all_p_state','comb_tp_fp_mono','mono_two_all_states', 'twophosphate_all_p_state','comb_tp_fp_two'))])

comb_tp_fp_mono_TP <- comb_tp_fp_mono$TP
comb_tp_fp_two_TP<-comb_tp_fp_two$TP
comb_tp_mono_two <- cbind2(comb_tp_fp_mono_TP,comb_tp_fp_two_TP)
## Basic Scatter Plot ##
 ### For TP and TNP of all MONO and TWO
par(mar = c(5, 4, 4, 8),                                
    xpd = TRUE)
plot(comb_tp_fp_mono$TP, xlim=c(0,122), ylim=c(0,10),pch = 19, xlab = "indeces of each score", 
     ylab="scores", main = "Distribution of TP and TNP values", sub = "mono and two phosphate")
par(new=TRUE)
plot(comb_tp_fp_two$TP, col="red",xlim=c(0,122), ylim=c(0,10),pch = 19,xlab = "", ylab="")
par(new=TRUE)
plot(comb_tp_fp_two$TNP, col="blue",xlim=c(0,122), ylim=c(0,10),pch = 19,xlab = "", ylab="")
par(new=TRUE)
plot(comb_tp_fp_mono$TNP, col="#69b3a2", xlim=c(0,122), ylim=c(0,10),pch = 19,xlab = "", ylab="")
coord <- par("usr")
legend(x = coord[2] * 1.05, y = coord[4],
       pch = 19, c("mono_TP","two_TP","two_TNP","mono_TNP"),col = c("black","red","blue","#69b3a2"))


### For FP and FNP of all MONO and TWO
## Basic Scatter Plot ##
par(mar = c(5, 4, 4, 8),                                
    xpd = TRUE)
plot(comb_tp_fp_mono$FP, xlim=c(0,122), ylim=c(0,110),pch = 19, xlab = "indeces of each score", 
     ylab="scores", main = "Distribution of FP and FNP values", sub = "mono and two phosphate")
#axis(side = 2, at=seq(0,110,10))
par(new=TRUE)
plot(comb_tp_fp_two$FP, col="red",xlim=c(0,122), ylim=c(0,110),pch = 19,xlab = "", ylab="")
#axis(side = 2, at=seq(1,110,10))
par(new=TRUE)
plot(comb_tp_fp_two$FNP, col="blue",xlim=c(0,122), ylim=c(0,110),pch = 19,xlab = "", ylab="")
#axis(side = 2, at=seq(1,110,10))
par(new=TRUE)
plot(comb_tp_fp_mono$FNP, col="#69b3a2", xlim=c(0,122), ylim=c(0,110),pch = 19,xlab = "", ylab="")
#axis(side = 2, at=seq(1,110,10))
coord <- par("usr")
legend(x = coord[2] * 1.05, y = coord[4],
       pch = 19, c("mono_FP","two_FP","two_FNP","mono_FNP"),col = c("black","red","blue","#69b3a2"))

# ggplot2 scatter plot
library(hrbrthemes)
library(reshape2)

melted <- melt(mono_two_all_states)
dimmelt <- dim(melted)
melted <- cbind(melted,1:dimmelt[1])

ggplot(melted, aes(x=melted[,"1:dimmelt[1]"], y=value, alpha=melted[,"intensity_mod"] )) + geom_point(size=2,shape=19, color=("black")) + 
  theme_ipsum() +
  labs(x="indeces of each scores for all peptides",y="scores", 
       title="Distribution of TP and FP scores on Eyers dataset",subtitle="mono and two phosphate")

  
## FIXME: change dimension of the data to become suitable for ggplot
# ggplot(melted, aes(x=melted[,"1:664"], y=value )) + geom_point(size=2,shape=19) +
#   
#   geom_point(y=melted[1:122,"value"], color="#69b3a2",size=2,shape=19) + #TP mono
#   geom_point(y=melted[123:244,"value"],color="#FFFFB3",size=2,shape=19) + #FP mono
#   
#   geom_point(y=melted[245:288,"value"],color="#ABB065",size=2,shape=7) + #TP two
#   geom_point(y=melted[245:332,"value"],color="#ABB065",size=2,shape=7) + #FP two
#   
#   geom_point(y=melted[333:454,"value"],color="#ABB065",size=2,shape=7) + #TNP mono
#   geom_point(y=melted[455:576,"value"],color="#ABB065",size=2,shape=7) + #FNP mono
#   
#   geom_point(y=melted[577:620,"value"],color="#ABB065",size=2,shape=7) + #TNP two
#   geom_point(y=melted[621:664,"value"],color="#ABB065",size=2,shape=7) + #FNP two
#   
#   theme_ipsum() 
#    

comb_tnps_tps <- cbind(comb_tp_fp_mono$TP,comb_tp_fp_two$TP, comb_tp_fp_mono$TNP, comb_tp_fp_two$TNP)

dim_comb_tnps_tps <- dim(comb_tnps_tps)[1]

comb_tnps_tps[45:dim_comb_tp[1],2] <- 0
comb_tnps_tps[45:dim_comb_tp[1],4] <- 0
colnames(comb_tnps_tps) <- c("TP_mono","TP_two","TNP_mono","TNP_two")

# <- seq(comb_tp_fp_two$TP, 166

 # All variation of comb_ objects are for ggplot generation ##
## Linear Regression ##
LM_comb_tp_fp_mono <- lm(TP~TNP,data=comb_tp_fp_mono)
plot(LM_comb_tp_fp_mono)


ggplot(as.data.frame(comb_tnps_tps), aes(x=TP_mono, y=TP_two)) + 
  geom_point( color="#69b3a2") +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum() +labs(x="Mono TP of Phospho Evidence scores",y=" Two TP of Phospho Evidence scores", 
                      title="Distribution of TP scores on Eyers dataset",subtitle="mono and two phosphate")


########## Figure genenation ##########
ggplot(as.data.frame(melted),
       aes(y = value, axis1 = phospho_state, axis2 = melted[,"variable"], axis3=melted[,"intensity_mod"])) +
       geom_alluvium(aes(fill = type), width = 1/12) +
       theme_ipsum()  +
  
       geom_stratum(width = 1/12, fill = "black", color = "grey") +
       geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
       scale_x_discrete(limits = c("phospho_state", "variable","intensity_mod"), expand = c(.05, .05)) +
       #scale_y_continuous(breaks = seq(0,200,10)) +
       scale_fill_brewer(type = "qual", palette = "Set1") +
       ggtitle("Distribution of N/PE Values of Mono-Two Phosphate from Eyer's Data")



#HEATMAP generation

paletteFunc <- colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929","#CC4C02"))
palette     <- paletteFunc(100)


pdf("tp_sum_two_phospho.pdf",width=50, height = 50)

pheatmap(as.matrix(comb_tp_fp_mono[,c(1,2,6,7)]),color=palette,  length.out = 50, show_rownames = F ,border_color="white", main = "Distribution of N/PE Values of Mono-Phosphate from Eyer's Data")
dev.off()

pheatmap(as.matrix(comb_tp_fp_two[,c(1,2,6,7)]),color=palette,  length.out = 50, show_rownames = F ,border_color="white", main = "Distribution of N/PE Values of Two-Phosphate from Eyer's Data")


#,cellwidth = 20, cellheight = 20,width = 7, height=10)









