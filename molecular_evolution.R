#------------------------------
# Molecular Evolution Analysis
# Augusto Franco
#------------------------------

library("Biostrings")
library(seqinr)
library(stringr)
library(dplyr)
library(phytools)
library(DECIPHER)
library(ape)
library(adegenet)
library(pegas)


path = "Genetic_diversity_metrics/"
setwd(path)

genes = list("GENE_LIST")
d = NULL

for (i in genes){
  
  gene_seq <- readDNAStringSet(paste('Genes_sequences/Seqs_',i,'.fasta',sep = ""), format = "fasta")
  gene_size <- width(gene_seq)
  
  #Read hit table + description table
  
  alig_hit = read.csv(paste("Hit_tables/",i,"-HitTable.csv",sep = ""),header = FALSE)
  names(alig_hit) <- c("query","subject","identity","align_length","mismatches","gaps_opens","q_start","q_end","s_start","s_end","evalue","bit_score")
  
  alig_hit$Query_cover <- (alig_hit$align_length/gene_size)*100
  #alig_hit$align_subject <- paste0(alig_hit$s_start,"-",alig_hit$s_end, sep = "")
  alig_hit$seq_num <- 1:nrow(alig_hit) 
  
  # filter step Query_Cover, Length, identity
  filter_data <- alig_hit %>% filter(evalue < 0.001 & Query_cover>90 & identity>90)
  
  #Count sequences
  count_seqs<- filter_data %>% group_by(subject) %>% summarise(count=n())
  
  #Filter fasta sequences
  fastaFile <- readDNAStringSet(paste("blast_alignments/",i,".txt",sep = ""), format = "fasta")
  names(fastaFile) <- alig_hit$seq_num
  filtered_seqs <- fastaFile[c(which(names(fastaFile) %in% filter_data$seq_num))]
  names(filtered_seqs) <- filter_data$subject
  
  
  # Multiple sequence alignment
  writeXStringSet(c(gene_seq, filtered_seqs), paste("Seqs_for_align/Finalseqs_",i,".fasta",sep = ""),format = "fasta")
  seq_x_align <- readDNAStringSet(paste("Seqs_for_align/Finalseqs_",i,".fasta",sep = ""),format = "fasta")
  alignment <- AlignSeqs(seq_x_align,anchor=NA,processors = 6)
  writeXStringSet(alignment,filepath = paste("alignments/align_",i,".fasta", sep = ""),format = "fasta")
  align <- as.DNAbin(alignment)
  
  #Genetic diversity metrica
  nuc_diversity <- nuc.div(x = align)
  tajD <- tajima.test(align)
  d = rbind(d, data.frame(i,nuc_diversity,tajD))
}

names(d)[names(d)== 'i'] <- 'Gene_id'
write.csv(d,file = "Nuc_div_Results.csv")