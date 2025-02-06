#if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
#BiocManager::install(version = "3.16")
BiocManager::install("ggtree")

library(ggtree)
library(ape)
library(ggplot2)
#library(TreeTools)
#library(TDbook)
library(dplyr)

#head(df_tip_data)

#tree_boots



tree <- read.tree("Europe_recent_Bgs.raxml.support")



meta<-read.csv("~/projects/project_fungicides/analysis/cyp51/isorelate/Minadakis_et_al_2024_Supplementary_Data_S1_rev.csv")


table <- meta[c("Sample_ID","cytb_G143A")]
table
table <- table[1:415,]

p<-ggtree(tree)+ geom_treescale()
p
  
tree$tip.label
is.rooted(tree)

secalis=c("S-1459","S-1203","S-1391","S-1201","S-1400")

is.monophyletic(tree, secalis)
getMRCA(tree,secalis)


tree_rooted<- root(tree,node=getMRCA(tree,secalis))

# Modify node labels based on their values
tree_rooted$node.label <- as.numeric(tree_rooted$node.label) # Ensure numeric


num_tips <- length(tree_rooted$tip.label)
num_nodes <- length(tree_rooted$node.label)

nodes_of_interest <- which(tree_rooted$node.label <= 50)
nodes_of_interest <- nodes_of_interest + num_tips

tree_data <- fortify(tree_rooted)

tree_tips <- data.frame(Sample_ID = tree_rooted$tip.label)
merged_data <- left_join(tree_tips, table, by = "Sample_ID")

p <- ggtree(tree_rooted, layout = "rectangular", size = 0.5) +  # Reduced line size
  geom_treescale(
    fontsize = 3,    # Increase label font size
    linesize = 0.8,  # Increase scale bar thickness
    offset = 2,      # Adjust vertical offset
    label = "Scale Bar"  # Add a label to the scale bar
  ) +
#  geom_nodelab(size=1)+
  geom_nodepoint(
    data = tree_data[tree_data$node %in% nodes_of_interest, ],  # Filter nodes of interest
    color = "purple3",  # Set color for the points
    size = 2,  # Adjust the size of the points
    shape = 15  # Use filled circle as the symbol
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(), # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    legend.position = "none"  # Remove the legend
  )

p

# Display the plot
print(p)


meta1 <- merged_data
rownames(meta1) <- meta1$Sample_ID
head(meta1)


meta1[1] <- list(NULL)

gheatmap(p,meta1,width=0.2,colnames=FALSE) +  
  scale_fill_manual(breaks=c("G", "A"), values=c("#ff8577","#1a82d2"), name=NULL)+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title =element_blank(), #change legend title font size
        legend.text = element_text(size=8)) #change legend text font size



#p <- p %<+% table + geom_tippoint(aes(color = as.factor(cytb_G143A)),shape=20,size=1)
#p<-p+ scale_color_manual(
#  values = c("G" = "#ff8577", "A" = "#1a82d2") )

p
ggsave("tree_europe_support.pdf", width=10, height= 6)
