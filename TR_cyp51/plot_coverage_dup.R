library(vcfR)
vcf <- read.vcfR( "tritici_ext_eur_cyp51_combined_genotyped_filtered.vcf.gz" )
library(isoRelate,lib = "~/data/R_lib")

gt <- extract.gt(vcf, element = "AD")

total_coverage_matrix <- t(apply(gt, 1, function(x) {
  sapply(strsplit(x, ","), function(y) sum(as.numeric(y)))  # Sum of reference and alternate alleles
}))


samples_names_matrix<-colnames(total_coverage_matrix)
cov_genome <- read.csv("../phasing/cyp51_cnv.csv")

cov_genome_filtered <- subset(cov_genome, Isolate %in% samples_names_matrix)
cov_genome_filtered <- cov_genome_filtered[match(samples_names_matrix, cov_genome_filtered$Isolate), ]

ratio_matrix<-total_coverage_matrix

for (i in 1:ncol(total_coverage_matrix)) {
  ratio_matrix[,i] <- total_coverage_matrix[,i] / cov_genome_filtered$genome_wide_coverage[i]
}


##### plot function
plot_dup <- function(samples_to_plot,ratio_matrix,pos_label){
  positions_range <- 4200:4300

    # Extract the coverage data for the specified range and samples
  coverage_data <- ratio_matrix[positions_range, samples_to_plot]
  
  # Create a sequence for the x-axis (positions), shifted so that 5000 = 0
  positions <- positions_range - 5000
  
  par(mar = c(1.5, 2.5, 0, 0.5))
  if (pos_label == "Position"){
    par(mar = c(2, 2.5, 0, 0.5))  
  }

  plot(positions, coverage_data[,1], type = "l",col=alpha("gray", 0.4), xlab=pos_label,ylab = "Coverage ratio",cex.lab=0.7,cex.axis=0.6,#axes=FALSE,
       main = "", ylim = c(0,8))
  # Add lines for the other samples
  for (i in 2:ncol(coverage_data)) {
    lines(positions, coverage_data[,i],col=alpha("gray", 0.4))
  }
  
  positions_range <- 6700:6800
  
  # Extract the coverage data for the specified range and samples
  coverage_data <- ratio_matrix[positions_range, samples_to_plot]

  # Create a sequence for the x-axis (positions), shifted so that 5000 = 0
  positions <- positions_range - 5000

    par(mar = c(1.5, 1, 0, 0.5))
    if (pos_label == "Position"){
      par(mar = c(2, 1, 0, 0.5))  
    }
  # Plot the coverage for each sample in a loop
  plot(positions, coverage_data[,1], type = "l",col=alpha("gray", 0.4), xlab=pos_label, ylab = "",cex.axis=0.6,cex.lab=0.7,
       main = "", ylim = c(0,8))
  # Add lines for the other samples
  for (i in 2:ncol(coverage_data)) {
    lines(positions, coverage_data[,i],col=alpha("gray", 0.4))
  }
 
  
}  


## read isorelate clusters

load("../isorelate_cyp51/BgtE+r_cyp51_geno.RData")

my_i_clusters <- getIBDiclusters(ped.genotypes = my_genotypes, 
                                 ibd.segments = my_ibd, 
                                 interval = c("8", 5728014, 5729706), 
                                 prop=1, 
                                 hi.clust = FALSE)

get_first_occurrence <- function(x) {
  split_parts <- strsplit(x, "/")[[1]]
  return(split_parts[1])
}




# plot 8 main clusters

pdf("duplication_8_cl.pdf")
par(mfrow=c(8,2))
par(oma = c(1, 1, 1, 1))
par(mgp = c(1, 0.2, 0),tck=-0.04)
pos_label=c("","","","","","","","Position")
for (i in 1:8){
  samples_cl_1 <- sapply(my_i_clusters$clusters[[i]], get_first_occurrence)
  
  samples_to_plot <- which(cov_genome_filtered$Isolate %in% as.vector(samples_cl_1))
  
  plot_dup(samples_to_plot = samples_to_plot, ratio_matrix = ratio_matrix,pos_label = pos_label[i])
}

dev.off()

