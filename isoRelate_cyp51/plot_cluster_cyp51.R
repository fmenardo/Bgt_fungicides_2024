setwd("~/projects/project_fungicides/analysis/cyp51/isorelate/")
library(isoRelate,lib = "~/data/R_lib")
load("BgtE+r_cyp51_geno.RData")
library(ggplot2)
library(maps)
library(dplyr)
library(igraph)

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

meta<-read.csv("~/projects/nikos/fungicide/resistance/2_output/cyp51_metadata.csv")
meta1<-data.frame(meta$Sample.Name,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$cyp51_Y136F,meta$cnv,meta$cyp51_S509T,meta$cyp51_S79T, meta$cyp51_K175N)
colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Y136F","cnv","S509T","S79T","K175N")

groups<- merge(my_groups,meta1,by.x="fid",by.y="Sample.Name")
groups
groups$merged <- paste(groups$fid, groups$iid, sep = "/")
colnames(groups) <- c("fid","iid","pid","Country","Population","Longitude","Latitude","Y136F","cnv","S509T","S79T","K175N","ID")

vertices<-V(g)$name
vertices

groups_sorted <- groups[match(vertices, groups$ID), ]



## attach attributes to graph
vertex_attr(g, "Population", index = V(g)) <- groups_sorted$Population
vertex_attr(g, "Y136F", index = V(g)) <- groups_sorted$Y136F
vertex_attr(g, "cnv", index = V(g)) <- groups_sorted$cnv
vertex_attr(g, "S509T", index = V(g)) <- groups_sorted$S509T
vertex_attr(g, "S79T", index = V(g)) <- groups_sorted$S79T
vertex_attr(g, "K175N", index = V(g)) <- groups_sorted$K175N


#Define a color mapping
color_mapping <- c("N_EUR" = "#377EB8", "S_EUR1" = "#EA9999","S_EUR2" = "#E41A1C", "TUR" = "#E5B110", "ME" = "#984EA3")
shape_mapping <- c("Y" = "circle", "F" = "square")
color_mapping_cnv <- c("1" = "#FDE725FF", "2" = "#8FD744FF","3" = "#35B779FF", "4" = "#21908CFF", "5" = "#31688EFF","6" = "#443A83FF", "7" = "#440154FF")
color_mapping_Y136F <- c("Y" = "#ff8577", "F" = "#1a82d2")
color_mapping_S509T <- c("S" = "#ff8577", "T" = "#1a82d2")
color_mapping_S79T <- c("S" = "#ff8577", "T" = "#1a82d2")
color_mapping_K175N <- c("K" = "#ff8577", "N" = "#1a82d2")


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




vertex_colors <- color_mapping[V(g)$Population]
pdf("Cyp51_ibd_population.pdf")
par(mfrow = c(1, 1))
par(oma=c(0,0,0,0))

par(mar = c(0,0,0,0))
plot.igraph(g,
     layout=layout,
#     vertex.shape=vertex_shape,
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.color=vertex_colors,
     edge.color="gray85"
#     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)
#legend(x=-0.8,y=-1, legend = names(color_mapping), fill = color_mapping, cex = 0.8, bty="n",horiz= "True")
legend(x=0.75,y=1.1, legend = names(color_mapping), fill = color_mapping, cex = 0.8, bty="n",
       title = as.expression(bquote(bold("Population"))),title.adj = 0.2,title.cex=1)#,horiz= "True")
dev.off()

vertex_colors <- color_mapping_cnv[V(g)$cnv]



pdf("Cyp51_cnv.pdf")
par(mfrow = c(1, 1))
par(oma=c(0,0,0,0))

par(mar = c(0,0,0,0))
plot.igraph(g,
            layout=layout2,
#            vertex.shape=vertex_shape,
            vertex.size=2.5,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)
#legend(x=-0.7,y=-1, legend = names(color_mapping_cnv), fill = color_mapping_cnv, cex = 0.8, bty="n",horiz= "True")
legend(x=-1.1,y=1.1, legend = names(color_mapping_cnv), fill = color_mapping_cnv, cex = 0.8, bty="n",
       title = as.expression(bquote(bold("Copy number"))),title.adj = 0.5,title.cex=1)#,horiz= "True")

dev.off()

vertex_colors <- color_mapping_Y136F[V(g)$Y136F]

pdf("Cyp51_Y136F.pdf")
par(mfrow = c(1, 1))
par(oma=c(0,0,0,0))

par(mar = c(0,0,0,0))
plot.igraph(g,
            layout=layout,
#            vertex.shape=vertex_shape,
            vertex.size=2.5,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1,y=1.1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n",
       title = as.expression(bquote(bold("Y136F"))),title.adj = 0.5,title.cex=1)#,horiz= "True")
dev.off()

vertex_colors <- color_mapping_S509T[V(g)$S509T]

pdf("Cyp51_S509T.pdf")
par(mfrow = c(1, 1))
par(oma=c(0,0,0,0))
par(mar = c(0,0,0,0))

plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=2.5,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1,y=1.1, legend = names(color_mapping_S509T), fill = color_mapping_Y136F, cex = 0.8, bty="n",
       title = as.expression(bquote(bold("S509T"))),title.adj = 0.5,title.cex=1)#,horiz= "True")

dev.off()

vertex_colors <- color_mapping_S79T[V(g)$S79T]

pdf("Cyp51_S79T.pdf")
par(mfrow = c(1, 1))
par(oma=c(0,0,0,0))
par(mar = c(0,0,0,0))

plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=2.5,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1,y=1.1, legend = names(color_mapping_S79T), fill = color_mapping_S79T, cex = 0.8, bty="n",
       title = as.expression(bquote(bold("S79T"))),title.adj = 0.5,title.cex=1)#,horiz= "True")


dev.off()


vertex_colors <- color_mapping_K175N[V(g)$K175N]plot.new()

pdf("Cyp51_K175N.pdf")
par(mfrow = c(1, 1))
par(oma=c(0,0,0,0))
par(mar = c(0,0,0,0))

plot.igraph(g,
            layout=layout,
            #            vertex.shape=vertex_shape,
            vertex.size=2.5,
            vertex.label.cex=0.5,
            vertex.label.dist=0.4,
            vertex.label.color="black",
            vertex.label=NA,
            vertex.color=vertex_colors,
            edge.color="gray80"
            #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Population"))
)

#legend(x=-0.6,y=-1, legend = names(color_mapping_Y136F), fill = color_mapping_Y136F, cex = 0.8, bty="n,horiz= "True")
legend(x=-1,y=1.1, legend = names(color_mapping_K175N), fill = color_mapping_K175N, cex = 0.8, bty="n",
       title = as.expression(bquote(bold("K175N"))),title.adj = 0.5,title.cex=1)#,horiz= "True")


dev.off()

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

vertex_colors <- color_mapping_Y136F[V(g)$Y136F]

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

vertex_colors <- color_mapping_S509T[V(g)$S509T]

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
vertex_colors <- color_mapping_S79T[V(g)$S79T]


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

vertex_colors <- color_mapping_K175N[V(g)$K175N]


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





############### plot only one cluster at the time
world <- map_data("world")

for (i in 1:1){

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

ggsave(paste0("cyp51_ibd_cluster_",i,"_map.pdf"),width=10,height=10)

}
  




