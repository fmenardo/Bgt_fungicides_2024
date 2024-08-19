setwd("~/projects/project_fungicides/analysis/isoRelate_erg24/")
library(isoRelate,lib = "~/data/R_lib")
load("BgtE+r_erg24_geno.RData")
library(ggplot2)
library(maps)
library(dplyr)
library(igraph)
library(viridis)


## make clusters with isoRelate
my_i_clusters <- getIBDiclusters(ped.genotypes = my_genotypes, 
                                 ibd.segments = my_ibd, 
                                 interval = c("7", 6622512, 6624033), 
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
meta<-read.csv("~/projects/nikos/fungicide/resistance/2_output/final_metadata_9_7_24.csv")
meta1<-data.frame(meta$Sample.Name,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$erg24_Y165F,
                  meta$erg24_F289H,meta$erg24_V295L,meta$erg24_D137E, meta$erg24_D291N)

colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Y165F","F289H","V295L","D137E","D291N")

groups<- merge(my_groups,meta1,by.x="fid",by.y="Sample.Name")
groups
groups$merged <- paste(groups$fid, groups$iid, sep = "/")
colnames(groups) <- c("fid","iid","pid","Country","Population","Longitude","Latitude","Y165F","F289H","V295L","D137E","D291N","ID")

vertices<-V(g)$name
vertices

groups_sorted <- groups[match(vertices, groups$ID), ]

groups_sorted <- groups_sorted %>%
  mutate(haplotype = paste(Y165F, V295L, F289H, D137E, D291N, sep = ""))

## attach attributes to graph
vertex_attr(g, "Population", index = V(g)) <- groups_sorted$Population
vertex_attr(g, "Y165F", index = V(g)) <- groups_sorted$Y165F
vertex_attr(g, "F289H", index = V(g)) <- groups_sorted$F289H
vertex_attr(g, "V295L", index = V(g)) <- groups_sorted$V295L
vertex_attr(g, "D137E", index = V(g)) <- groups_sorted$D137E
vertex_attr(g, "D291N", index = V(g)) <- groups_sorted$D291N
vertex_attr(g, "hap", index = V(g)) <- groups_sorted$haplotype



#Define a color mapping
color_mapping <- c("N_EUR" = "#377EB8", "S_EUR1" = "#EA9999","S_EUR2" = "#E41A1C", "TUR" = "#E5B110", "ME" = "#984EA3")
color_mapping_Y165F <- c("Y" = "#ff8577", "F" = "#1a82d2")
color_mapping_F289H <- c("F" = "#ff8577", "H" = "#1a82d2")
color_mapping_V295L <- c("V" = "#ff8577", "L" = "#1a82d2")
color_mapping_D137E <- c("D" = "#ff8577", "E" = "#1a82d2")
color_mapping_D291N <- c("D" = "#ff8577", "N" = "#1a82d2")
color_mapping_hap <- c("FLFDD" = "cadetblue1", "FVFDD" = "indianred4", "YLFDD" = "#1a82d2",
                       "YLHDD" = "darkblue", "YVFDD" = "#ff8577",
                       "YVFDN" = "cadetblue4", "YVFED" = "firebrick1", "YVHDD" = "blue")



#layout=layout_nicely(g)

#save(layout, file="layout.RData")

load("layout.RData")

###################################################
#################### plot combined

vertex_colors <- color_mapping[V(g)$Population]
pdf("erg24_ibd_clusters.pdf")
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


vertex_colors <- color_mapping_V295L[V(g)$V295L]

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

legend(x=-1.1,y=1.1, legend = names(color_mapping_V295L), fill = color_mapping_V295L, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("V295L"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")


vertex_colors <- color_mapping_Y165F[V(g)$Y165F]

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

legend(x=-1.1,y=1.1, legend = names(color_mapping_Y165F), fill = color_mapping_Y165F, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("Y165F"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")



vertex_colors <- color_mapping_F289H[V(g)$F289H]

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
legend(x=-1.1,y=1.1, legend = names(color_mapping_F289H), fill = color_mapping_F289H, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("F289H"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")

dev.off()

pdf("erg24_ibd_clusters_part2.pdf")


par(mfrow = c(1, 2))
par(oma=c(0,0,0,0))
vertex_colors <- color_mapping_D137E[V(g)$D137E]

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
legend(x=-1.11,y=1.1, legend = names(color_mapping_D137E), fill = color_mapping_D137E, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("D137E"))),title.adj = 0.2,title.cex=0.7)#,horiz= "True")


vertex_colors <- color_mapping_D291N[V(g)$D291N]

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

legend(x=-1.1,y=1.1, legend = names(color_mapping_D291N), fill = color_mapping_D291N, cex = 0.65, bty="n",
       title = as.expression(bquote(bold("D291N"))),title.adj = 0.5,title.cex=0.7)#,horiz= "True")

dev.off()

#################################################
####### plot hap

pdf("erg24_ibd_clusters_hap.pdf")
vertex_colors <- color_mapping[V(g)$Population]

par(mfrow = c(1, 2))
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

vertex_colors <- color_mapping_hap[V(g)$hap]

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

dev.off()


############### plot only one cluster at the time
world <- map_data("world")

for (i in 1:10){

# identify vertices in cluster
subcluster_vertices <- which(V(g)$name %in% my_i_clusters$clusters[[i]])

# Create a subgraph containing only the vertices of the subcluster
subgraph <- induced_subgraph(my_i_clusters$i.network, subcluster_vertices)

vertices<-V(subgraph)$name


# fetch metainformations and make attributes for subgraph
sub<-groups_sorted[groups_sorted$ID %in% vertices,]

sub_sorted <- sub[match(vertices, sub$ID), ]

vertex_attr(subgraph, "Population", index = V(subgraph)) <- sub_sorted$Population
vertex_attr(subgraph, "Longitude", index = V(subgraph)) <- sub_sorted$Longitude
vertex_attr(subgraph, "Latitude", index = V(subgraph)) <- sub_sorted$Latitude

# Step 2: Assign colors to the vertices based on their population
vertex_colors <- color_mapping[V(subgraph)$Population]

pdf(paste0("cyp51_ibd_cluster_",i,".pdf"))

plot(subgraph,
     edge.width=0.2,
     vertex.size=2.5,
#     vertex.label.cex=0.5,
#     vertex.label.dist=1,
#     vertex.label.color="black",
     vertex.label=NA,
     vertex.color=vertex_colors
)
legend("bottomleft", legend = names(color_mapping), fill = color_mapping, cex = 0.5)
dev.off()

## plot on map with ggplot
# layout for map plotting
lo <- layout.norm(as.matrix(sub_sorted[,6:7]))
#convert nodes and edges to df
nodes <- data.frame(
  id = V(subgraph)$name,
  lat = V(subgraph)$Latitude,
  lon = V(subgraph)$Longitude
)

edges <- as.data.frame(get.edgelist(subgraph))
colnames(edges) <- c("from", "to")
edges <- edges %>%
  left_join(nodes, by = c("from" = "id")) %>%
  rename(lat_from = lat, lon_from = lon) %>%
  left_join(nodes, by = c("to" = "id")) %>%
  rename(lat_to = lat, lon_to = lon)


### make plot with nodes
p <- ggplot() +
#  geom_point(data = nodes, aes(x = lon, y = lat), color = vertex_colors, size = 2.5)
  geom_jitter(data = nodes, aes(x = lon, y = lat), color = vertex_colors, size = 1.5,width=0.6, height=0.6)+
  geom_point(data = nodes, aes(x = lon, y = lat), color = "black", size = 1)

##add edges

p <- p +
  geom_segment(data = edges, aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to), color = "gray50",linewidth=0.15,alpha=0.4)
## add map
p +  geom_path(data = world, aes(x = long, y = lat, group = group), color = "gray50",linewidth=0.2)+
  xlim(c(-12, 45)) +  # Set specific x-axis limits
  ylim(c(20, 65))+
  coord_fixed(ratio = 1)+
  theme_void()

ggsave(paste0("erg24_ibd_cluster_",i,"_map.pdf"),width=10,height=10)

}
  


###########################
#######age of clusters

my_ibd_cyp <- subset(my_ibd,chr == "7" & start_position_bp < 6622512 & end_position_bp > 6624033 ) 


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

pdf("age_of_cluster_erg24.pdf")
boxplot(length_M*100 ~ cluster, data = my_ibd_cl_u,
        xlab = "Cluster",
        ylab = "Length of IBD segment (cM)")
dev.off()

medians <- my_ibd_cl_u %>%
  group_by(cluster) %>%
  summarize(median_length = median(length_M * 100))

medians$gen <- (100/medians$median_length)/2


##############################
##### test association muttaion / cluster, cluster with 6 or more (13 clusters)

clustered_samples <- unlist(my_i_clusters$clusters[1:13])

n_ncl <- length(unique_elements) ### these are singleton
n_cl <- length(clustered_samples)

V295L_tot <- subset(groups_sorted, groups_sorted$V295L =="L")

nV295L_cl<- length(intersect(V295L_tot$ID,clustered_samples))
nV295L_ncl <- length(intersect(V295L_tot$ID,unique_elements))

V295L_freq<-matrix(c(nV295L_cl, n_cl-nV295L_cl, nV295L_ncl, n_ncl-nV295L_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("V295L", "not_V295L"),
                                   period = c("clustered", "not clustered")))

fisher.test(V295L_freq, alternative = "t")


Y165F_tot <- subset(groups_sorted, groups_sorted$Y165F =="F")

nY165F_cl<- length(intersect(Y165F_tot$ID,clustered_samples))
nY165F_ncl <- length(intersect(Y165F_tot$ID,unique_elements))

Y165F_freq<-matrix(c(nY165F_cl, n_cl-nY165F_cl, nY165F_ncl, n_ncl-nY165F_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("Y165F", "not_Y165F"),
                                   period = c("clustered", "not clustered")))

fisher.test(Y165F_freq, alternative = "t")



Y165F_noV295L_tot <- subset(groups_sorted, groups_sorted$Y165F =="F" & groups_sorted$V295L =="V")

nY165F_noV295L_cl<- length(intersect(Y165F_noV295L_tot$ID,clustered_samples))
nY165F_noV295L_ncl <- length(intersect(Y165F_noV295L_tot$ID,unique_elements))

Y165F_noV295L_freq<-matrix(c(nY165F_noV295L_cl, n_cl-nY165F_noV295L_cl, nY165F_noV295L_ncl, n_ncl-nY165F_noV295L_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("Y165F_noV295L", "not_Y165F_noV295L"),
                                   period = c("clustered", "not clustered")))

fisher.test(Y165F_noV295L_freq, alternative = "t")


F289H_tot <- subset(groups_sorted, groups_sorted$F289H =="H")

nF289H_cl<- length(intersect(F289H_tot$ID,clustered_samples))
nF289H_ncl <- length(intersect(F289H_tot$ID,unique_elements))

F289H_freq<-matrix(c(nF289H_cl, n_cl-nF289H_cl, nF289H_ncl, n_ncl-nF289H_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("F289H", "not_F289H"),
                                   period = c("clustered", "not clustered")))

fisher.test(F289H_freq, alternative = "t")

D137E_tot <- subset(groups_sorted, groups_sorted$D137E =="E")

nD137E_cl<- length(intersect(D137E_tot$ID,clustered_samples))
nD137E_ncl <- length(intersect(D137E_tot$ID,unique_elements))

D137E_freq<-matrix(c(nD137E_cl, n_cl-nD137E_cl, nD137E_ncl, n_ncl-nD137E_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("D137E", "not_D137E"),
                                   period = c("clustered", "not clustered")))

fisher.test(D137E_freq, alternative = "t")


D291N_tot <- subset(groups_sorted, groups_sorted$D291N =="N")

nD291N_cl<- length(intersect(D291N_tot$ID,clustered_samples))
nD291N_ncl <- length(intersect(D291N_tot$ID,unique_elements))

D291N_freq<-matrix(c(nD291N_cl, n_cl-nD291N_cl, nD291N_ncl, n_ncl-nD291N_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("D291N", "not_D291N"),
                                   period = c("clustered", "not clustered")))

fisher.test(D291N_freq, alternative = "t")


wild_tot <- subset(groups_sorted, groups_sorted$hap =="YVFDD")

nwild_cl<- length(intersect(wild_tot$ID,clustered_samples))
nwild_ncl <- length(intersect(wild_tot$ID,unique_elements))

wild_freq<-matrix(c(nwild_cl, n_cl-nwild_cl, nwild_ncl, n_ncl-nwild_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("wild", "not_wild"),
                                   period = c("clustered", "not clustered")))

fisher.test(wild_freq, alternative = "t")


D291N_tot <- subset(groups_sorted, groups_sorted$D291N =="N")

nD291N_cl<- length(intersect(D291N_tot$ID,clustered_samples))
nD291N_ncl <- length(intersect(D291N_tot$ID,unique_elements))

D291N_freq<-matrix(c(nD291N_cl, n_cl-nD291N_cl, nD291N_ncl, n_ncl-nD291N_ncl),
                   nrow = 2,
                   dimnames = list(Var = c("D291N", "not_D291N"),
                                   period = c("clustered", "not clustered")))

fisher.test(D291N_freq, alternative = "t")
