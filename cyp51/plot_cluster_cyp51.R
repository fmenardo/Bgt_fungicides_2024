setwd("~/projects/project_fungicides/analysis/cyp51/isorelate/")
library(isoRelate,lib = "~/data/R_lib")
load("BgtE+r_cyp51_geno.RData")
library(ggplot2, lib.loc="/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2")
library(maps, lib.loc="/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2")
library(dplyr, lib.loc="/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2")
library(igraph, lib.loc="/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2")
library(viridis, lib.loc="/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2")


## make clusters with isoRelate
my_i_clusters <- getIBDiclusters(ped.genotypes = my_genotypes, 
                                 ibd.segments = my_ibd, 
                                 interval = c("8", 5728014, 5729706), 
                                 prop=1, 
                                 hi.clust = FALSE)


g<-my_i_clusters$i.network


Samples <- my_genotypes$pedigree$fid
samples <- paste(Samples, Samples, sep = "/")


## plot all clusters
# first attach metainformation to graph

my_groups <- my_genotypes[[1]][,1:3]
my_groups


## identify samples not in the graph (not clustered) and add them
vertices<-V(g)$name

unique_elements <- setdiff(samples, vertices)
length(unique_elements)
g <- add_vertices(g, length(unique_elements), attr = list(name = unique_elements))
vertices<-V(g)$name


# read meta, extract samples and reorder them based on graph vertex order
meta<-read.csv("~/projects/nikos/fungicide/resistance/2_output/final_metadata_20_8_24.csv")
meta<-read.csv("~/projects/project_fungicides/analysis/cyp51/isorelate/Minadakis_et_al_2024_Supplementary_Data_S1_rev.csv")

meta1<-data.frame(meta$Sample_ID,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$cyp51_Y136F,meta$cyp51_cn_interpretation,meta$cyp51_S509T,meta$cyp51_S79T,meta$cyp51_K175N,meta$cyp51_Y136F_het,meta$cyp51_S509T_het,meta$cyp51_S79T_het,meta$cyp51_K175N_het)

colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Y136F","cnv","S509T","S79T","K175N","Y136F_het","S509T_het","S79T_het","K175N_het")

groups<- merge(my_groups,meta1,by.x="fid",by.y="Sample.Name")
groups
groups$merged <- paste(groups$fid, groups$iid, sep = "/")
colnames(groups) <- c("fid","iid","pid","Country","Population","Longitude","Latitude","Y136F","cnv","S509T","S79T","K175N","Y136F_het","S509T_het","S79T_het","K175N_het","ID")

vertices<-V(g)$name
vertices

groups_sorted <- groups[match(vertices, groups$ID), ]



## attach attributes to graph
vertex_attr(g, "Population", index = V(g)) <- groups_sorted$Population
vertex_attr(g, "Y136F_het", index = V(g)) <- groups_sorted$Y136F_het
vertex_attr(g, "cnv", index = V(g)) <- groups_sorted$cnv
vertex_attr(g, "S509T_het", index = V(g)) <- groups_sorted$S509T_het
vertex_attr(g, "S79T_het", index = V(g)) <- groups_sorted$S79T_het
vertex_attr(g, "K175N_het", index = V(g)) <- groups_sorted$K175N_het
#vertex_attr(g, "cyp51_hap", index = V(g)) <- groups_sorted$cyp51_ps_haplotype

#Define a color mapping
color_mapping <- c("N_EUR" = "#377EB8", "S_EUR1" = "#EA9999","S_EUR2" = "#E41A1C", "TUR" = "#E5B110", "ME" = "#984EA3")
shape_mapping <- c("Y" = "circle", "F" = "square")
color_mapping_cnv <- c("1" = "#FDE725FF", "2" = "#8FD744FF","3" = "#35B779FF", "4" = "#21908CFF", "5" = "#31688EFF","6" = "#443A83FF", "7" = "#440154FF")
color_mapping_Y136F <- c("Y" = "#ff8577", "F" = "#1a82d2", "X" = "#9F79EE")
color_mapping_S509T <- c("S" = "#ff8577", "T" = "#1a82d2", "X" = "#9F79EE")
color_mapping_S79T <- c("S" = "#ff8577", "T" = "#1a82d2", "X" = "#9F79EE")
color_mapping_K175N <- c("K" = "#ff8577", "N" = "#1a82d2", "X" = "#9F79EE")
#color_mapping_cyp51_hap <- c("a" = "#FCFDBFFF", "b" = "#FEBA80FF", "d" = "#F8765CFF", "e" = "#D3436EFF", "n" = "#982D80FF", "o" = "#5F187FFF",
#                             "p" = "#231151FF")

#color_mapping_cyp51_hap <- c("a" = "black", "b" = "#443A83FF", "d" = "#31688EFF", "e" = "#21908CFF", "n" = "#35B779FF", "o" = "#8FD744FF",
#                             "p" = "#FDE725FF")
#color_mapping_cyp51_hap <- c("a" = "#76EEC6", "b" = "#556B2F", "d" = "#8B4500", "e" = "#EE9A49", "n" = "#EE3A8C", "o" = "#008B8B",
#                             "p" = "#F5DEB3")





#############   CUSTOMIZE LAYOUT 
#layout= layout_with_kk(g)

#layout <- layout_nicely(g)


#cl1 <- which(vertices %in% my_i_clusters$clusters[[1]])
#cl2 <- which(vertices %in% my_i_clusters$clusters[[2]])


#layout1 <- layout


##### jitter cl 1
#for (i in 1:length(cl1)){
  
#  layout1[cl1[i],1] <- jitter(layout[cl1[i],1],factor=1,amount=4)
#  layout1[cl1[i],2] <- jitter(layout[cl1[i],2],factor=1,amount=4)
#}

#layout2 <- layout1
###### jitter cl 2
#for (i in 1:length(cl2)){
  
#  layout2[cl2[i],1] <- jitter(layout1[cl2[i],1],factor=1,amount=2)
#  layout2[cl2[i],2] <- jitter(layout1[cl2[i],2],factor=1,amount=2)
#}

#save(layout2,file="layout.Rdata")
##############################################

load("layout.Rdata")
layout <- layout2




###################################################
#################### plot combined

vertex_colors <- color_mapping[V(g)$Population]
pdf("Cyp51_ibd_clusters_part1.pdf")
par(mfrow = c(2, 2))
par(oma=c(0,0,0,0))

par(mar = c(0,0,0,0))
plot.igraph(g,
            layout=layout,
            #     vertex.shape=vertex_shape,
            vertex.size=3.2,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray85"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)
#legend(x=-0.8,y=-1, legend = names(color_mapping), fill = color_mapping, cex = 0.8, bty="n",horiz= "True")
legend(x=-1.11,y=1.1, legend = names(color_mapping), fill = color_mapping, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("Population"))),title.adj = 0.2,title.cex=0.7)#,horiz= "True")

vertex_colors <- color_mapping_cnv[V(g)$cnv]

par(mar = c(0,0,0,0))
plot.igraph(g,
            layout=layout2,
            #            vertex.shape=vertex_shape,
            vertex.size=3.2,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)
#legend(x=-0.7,y=-1, legend = names(color_mapping_cnv), fill = color_mapping_cnv, cex = 0.8, bty="n",horiz= "True")
legend(x=-1.11,y=1.1, legend = names(color_mapping_cnv), fill = color_mapping_cnv, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("Copy number"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")

vertex_colors <- color_mapping_Y136F[V(g)$Y136F_het]

par(mar = c(0,0,0,0))
plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=3.2,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1.1,y=1.1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("Y136F"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")

vertex_colors <- color_mapping_S509T[V(g)$S509T_het]

plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=3.2,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1.1,y=1.1, legend = names(color_mapping_S509T), fill = color_mapping_Y136F, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("S509T"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")

dev.off()
pdf("Cyp51_ibd_clusters_part2.pdf")
par(mfrow = c(2, 2))
par(oma=c(0,0,0,0))

par(mar = c(0,0,0,0))
vertex_colors <- color_mapping_S79T[V(g)$S79T_het]


plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=3.2,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1,y=1.1, legend = names(color_mapping_S79T), fill = color_mapping_S79T, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("S79T"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")

vertex_colors <- color_mapping_K175N[V(g)$K175N_het]


plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=3.2,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1,y=1.1, legend = names(color_mapping_K175N), fill = color_mapping_K175N, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("K175N"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")


dev.off()








###########################
#######age of clusters

my_ibd_cyp <- subset(my_ibd,chr == "8" & start_position_bp < 5728014 & end_position_bp > 5729706 ) 


nrow(my_ibd_cyp)

my_ibd_cl <- my_ibd[0,]

for (i in 1:10){
  is <-my_i_clusters$clusters[[i]]
  for (t in 1:length(is)) {
    is_name<-strsplit(is[t],"/")
    print(is_name[[1]][1])
    temp<-subset(my_ibd_cyp,fid1==is_name[[1]][1] | fid2 ==is_name[[1]][1])
    #print(nrow(temp))
    temp$cluster <- i
    my_ibd_cl <- rbind(my_ibd_cl,temp)
  }
}

my_ibd_cl_u <- unique(my_ibd_cl)

par(mar = c(4, 5, 2, 2)) # Adjust as needed (bottom, left, top, right)

pdf("age_of_cluster_cyp51.pdf")
boxplot(length_M*100 ~ cluster, data = my_ibd_cl_u,
        xlab = "Cluster",
        ylab = "Length of IBD segment (cM)")
dev.off()

medians <- my_ibd_cl_u %>%
  group_by(cluster) %>%
  summarize(median_length = median(length_M * 100))

medians$gen <- (100/medians$median_length)/2
