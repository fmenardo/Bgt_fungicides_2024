
setwd("~/projects/project_fungicides/analysis/coverage/")

path_bam = "../../../alignments/july2024/"

files <- list.files(path=path_bam,pattern =".cov$")

genome_av_cov <- data.frame(matrix(ncol = 2, nrow = 0))


for (i in 1:length(files)){
#for (i in 1:5){
  table<-read.table(paste0(path_bam,files[i]))
  bases=0
  length=0
  for (r in 1:11){
    bases=bases+(table[r,3]*table[r,7])
    length=length+table[r,3]
  }
  av_cov=bases/length
  genome_av_cov[i,1] <- files[i]
  genome_av_cov[i,2] <- av_cov
}

genome_av_cov[,1] <- gsub("_compiled_marked_dup.bam.cov","",genome_av_cov[,1])
genome_av_cov[,1] <- gsub("aln_","",genome_av_cov[,1])
genome_av_cov[,1] <- gsub("renamed_","",genome_av_cov[,1])


colnames(genome_av_cov) <- c("Sample.Name","Average_coverage")

files <- list.files(pattern ="raw$")

files

names=c("Btub","cyp51","cytB","erg2","erg24","mit","sdhB","sdhC","sdhD" )

list_df <- list()

for (f in 1:9){

  lines <- readLines(files[f])
  head(lines)
  if (lines[1] == "../../alignments/july2024/*.bam"){lines <- lines[-1]}

  col1 <- c()
  col2 <- c()

  # Iterate over the lines in steps of three
  for (i in seq(1, length(lines), by = 3)) {
    # Ensure there are at least three lines remaining
    if ((i + 2) <= length(lines)) {
      # Add the first and third lines to the respective columns
      col1 <- c(col1, lines[i])
      col2 <- c(col2, lines[i + 2])
    }
  }

  # Create a data frame with the two columns
  list_df[[f]] <- data.frame(Column1 = col1, Column2 = col2)
  
  colnames(list_df[[f]]) <- c("Sample.Name",paste0(names[[f]],"_cov"))

  list_df[[f]]$Sample.Name <- gsub("../../../alignments/july2024/aln_","",list_df[[f]]$Sample.Name)
  list_df[[f]]$Sample.Name <- gsub("_compiled_marked_dup.bam","",list_df[[f]]$Sample.Name)
  list_df[[f]]$Sample.Name <- gsub("renamed_","",list_df[[f]]$Sample.Name)
  

}  

all_targets <- Reduce(function(x, y) merge(x, y, by= "Sample.Name"), list_df)


#gw_stats <- read.csv("~/projects/project_data_prep/analysis/pl_stats/2022+before2022+2023+ncsu_737_gatkpl_stats_new.csv")
#gw_covg <- gw_stats[,c(1,5)]
#gw_covg$Isolate <- gsub("renamed_","",gw_covg$Isolate)
#covg <- merge(all_targets,gw_covg,by.x = "Sample.Name", by.y = "Isolate")

covg <- merge(all_targets,genome_av_cov,by = "Sample.Name")


for (i in 2:11) {
  covg[, i] <- as.numeric(as.character(covg[, i]))
}

ratios <- data.frame(matrix(nrow = nrow(covg), ncol = 0))

#caalculate ratios for nuclear targets
for (i in c(1:2,4:5,7:9)){
  ratios[[paste0("ratio_", names[i])]] <- covg[, i + 1] / covg[, 11]
#  ratio <- covg[,i+1]/covg[,11]
}

#caalculate ratios for cytb target
  ratios[[paste0("ratio_", names[3])]] <- covg[, 3 + 1] / covg[, 7]

table_with_covg<-cbind(covg,ratios)

samples <- read.csv("~/projects/vcf_project_tritici/tritici_extended_europe_2022+before2022+2023+ncsu.args",header=FALSE)

setdiff(table_with_covg$Sample.Name,samples[,1])
setdiff(samples[,1],table_with_covg$Sample.Name)

final_ds <- merge(table_with_covg,samples, by.x="Sample.Name",by.y= "V1")
nrow(final_ds)


write.csv(final_ds,file="fungicides_coverage_table.csv")

labels=c("Btub","cyp51","erg2","erg24","sdhB","sdhC","sdhD","cytB")
pdf("coverage_ratio_plot.pdf")
boxplot(final_ds[,12:19],names=labels, ylab = "Coverage ratio" )
abline(h=1,lty=2)
dev.off()
