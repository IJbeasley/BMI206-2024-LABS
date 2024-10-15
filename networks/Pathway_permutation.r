rm(list=ls())
library("plyr")
library("RColorBrewer")
library("ggplot2")

##### A function to read SIF files
## PPI reading and cleaning function

SIF.reader <- function(filename,verbose=F){
  PPI <- read.table(as.is=T,sep='\t',fill=T, file = filename,
                    colClasses = c("character","character","character"))
  ## optionally get rid of unconnected genes
  # PPI <- PPI[PPI[,3] != "",]
  ## make a two column matrix
  PPI <- PPI[,c(1,3)]
  ## change the order of columns so edges are alphabetical
  ## (this removes any directedness of edges if any exists)
  PPI <- cbind(GeneA=pmin(PPI[,1],PPI[,2]),
               GeneB=pmax(PPI[,1],PPI[,2]))
  ## sort the edges
  PPI <- PPI[order(PPI[,1],PPI[,2]),]
  ## make sure there are no duplicates
  tmp  <- apply(PPI[,1:2],1,paste,collapse='--')
  if(verbose) print(summary(!duplicated(tmp))) ## all true
  ## end
  return(PPI)
}

## neighbor function and giant Net function
## arguments:
##   x   = any given subnetwork
##   PPI = a set of edges (must be a two column matrix)
##   min.fraction = the minimum fraction of the PPI that the
##                  giant net must cover
##   total.nodes  = total.nodes in PPI (do not change)
## value:
##   neighbor.func returns all neighbors of x in the PPI
##   giantNet.func finds the giant net in a given PPI that is
##     at least as large as min.fraction * total.nodes

get.nodes <- function(PPI) setdiff(unique(unlist(c(PPI))),"")
get.edges <- function(PPI) PPI[PPI[,1]!="",,drop=F]
n.nodes <- function(PPI) length(get.nodes(PPI))
n.edges <- function(PPI) nrow(get.edges(PPI))

## subset function
subset.func <- function(x,PPI) PPI[PPI[,1] %in% c(x,"") & PPI[,2] %in% c(x,""),,drop=F]

## neighbor function
neighbor.func <- function(x,PPI) get.nodes(PPI[PPI[,1] %in% x | PPI[,2] %in% x,])

## giant Net function
giantNet.func <- function(PPI,min.fraction=0.50,total.nodes = n.nodes(PPI)){  
  cat("Finding Giant net")
  PPI <- get.edges(PPI)
  running.giantNet <- c()
  for(loop.count in 1:nrow(PPI)){    
    giantNet <- PPI[loop.count,1]
    size <- length(giantNet)
    while(T){      
      giantNet <- neighbor.func(giantNet,PPI)      
      if(length(giantNet) == size) ## if no growth        
        break else size <- length(giantNet)      
    }    
    size <- length(giantNet)    
    if(size >= min.fraction * total.nodes){      
      ## if you have achieved the goal      
      running.giantNet <- giantNet      
      break      
    } else {      
      ## if you have not acheived the goal      
      if(size > length(running.giantNet))        
        ## but you have acheived a relative goal        
        running.giantNet <- giantNet      
      next      
    }  
  }
  cat(" Done\n")
  return(running.giantNet)
}

## read in parent network
backgroundPPI <- SIF.reader("parent_PPI.sif")
set.seed(1980)

## read in GWAS profiles
VEGAS1 <- read.table("HT.pvals.out", header=T, as.is=T)
VEGAS2 <- read.table("MS.pvals.out", header=T, as.is=T)
p.BlockDefine.threshold <- 0.05
PPInodes <- get.nodes(backgroundPPI)
VEGAS1 <- subset(VEGAS1, GenePvalue<p.BlockDefine.threshold)
## get the number of genes passint threshold
nodes1 <- VEGAS1$Gene
print(length(nodes1))
nodes1 <- unique(nodes1)
print(length(nodes1))
nodes1 <- intersect(nodes1, PPInodes)
print(length(nodes1))
subNet1 <- subset.func(nodes1,backgroundPPI)

##Obaserved network size
actual.size <- data.frame("Region" = "HT", 
               "Extracted_nodes" = length(nodes1), 
               "Edges" = length(subNet1), 
               "largest_nodes" = length(giantNet.func(subNet1)), stringsAsFactors = F)

VEGAS2 <- subset(VEGAS2, GenePvalue<p.BlockDefine.threshold)
## get the number of genes passint threshold
nodes2 <- VEGAS2$Gene
print(length(nodes2))
nodes2 <- unique(nodes2)
print(length(nodes2))
nodes2 <- intersect(nodes2, PPInodes)
print(length(nodes2))
subNet2 <- subset.func(nodes2,backgroundPPI)

actual.size <-rbind(actual.size, c("Region"="MS", 
                                 "Extracted_nodes" = length(nodes2), 
                                 "Edges" =length(subNet2), 
                                 "largest_nodes" = length(giantNet.func(subNet2)), stringsAsFactors = F))

my.cut <- function(x,bin.size=10){round(x/bin.size)*bin.size}


rand.full.net.sizes <- rep(seq(300,2000,by=100),each = 100)  ## as they exist in the folders
rand.giant.net.sizes <- c()
rand.full.net.connectivity <- c()
for(rand.count in 1:length(rand.full.net.sizes)){
  cat(rand.count,'\n')
  rand.genelist <- sample(get.nodes(backgroundPPI),
                          rand.full.net.sizes[rand.count])
  PPI.subset <- subset.func(rand.genelist,backgroundPPI)
  rand.full.net.connectivity[rand.count] <- n.edges(PPI.subset)
  if(n.edges(PPI.subset) == 0)
    rand.giant.net.sizes[rand.count] <- 0 else
      rand.giant.net.sizes[rand.count] <- n.nodes(giantNet.func(PPI.subset))
  rm(rand.count,PPI.subset,rand.genelist)
}

full.edge <- c()
big.size <- c()

for (i in c(0.50, 0.75, 0.90, 0.95, 0.99)){
  full.edge <- rbind(full.edge, cbind(percentile=i, edge=apply(matrix(rand.full.net.connectivity,100),2,function(x){quantile(x,p=i)})))
  big.size <- rbind(big.size, cbind(percentile=i, Max_size=apply(matrix(rand.giant.net.sizes,100),2,function(x){quantile(x,p=i)})))
}

full.edge <- as.data.frame(full.edge)
big.size <- as.data.frame(big.size)

full.edge$Nodes <- rep(unique(rand.full.net.sizes), 5)
big.size$Nodes <- rep(unique(rand.full.net.sizes), 5)
full.edge$percentile <- as.factor(full.edge$percentile)
big.size$percentile <- as.factor(big.size$percentile)
big.full <- data.frame(edge=rand.full.net.connectivity,giantSize=rand.giant.net.sizes)

big.full <- tapply(rand.giant.net.sizes, my.cut(rand.full.net.connectivity,bin.size=10), quantile,p=c(0.50, 0.75, 0.90, 0.95, 0.99))
big.full <- ldply(big.full)
big.full <- as.data.frame(rbind(cbind(percentile=0.50, Edge=as.numeric(big.full[,1]), giantSize=as.numeric(big.full[,2])),
                                cbind(percentile=0.75, Edge=as.numeric(big.full[,1]), giantSize=as.numeric(big.full[,3])),
                                #cbind(percentile=0.80, Edge=as.numeric(big.full[,1]),giantSize=as.numeric(big.full[,4])),
                                cbind(percentile=0.90, Edge=as.numeric(big.full[,1]), giantSize=as.numeric(big.full[,4])), 
                                cbind(percentile=0.95, Edge=as.numeric(big.full[,1]), giantSize=as.numeric(big.full[,5])),
                                cbind(percentile=0.99, Edge=as.numeric(big.full[,1]), giantSize=as.numeric(big.full[,6]))))

qcolours <- rep(c(brewer.pal(n = 7, name = "Set1")),6)

pdf("pathway_permutation.pdf",width=9, height=9, onefile=T)
fill.data <- full.edge[full.edge$percentile==0.5,]
fill.loess <- loess(edge~Nodes,fill.data,span=1)
actual.size <- as.data.frame(actual.size)
print(actual.size)
print(
    ggplot(data= full.edge, aes(x=Nodes, y=edge))+geom_smooth(aes(colour=percentile),method = "loess", size = 0.5,se=F, span=1)+
       geom_area(aes(x,y),data.frame(x=fill.data$Nodes,y=predict(fill.loess)),stat="identity", fill="gray",alpha=0.4)+
       scale_colour_manual(aes(percentile), values = qcolours)+
       guides(colour=guide_legend(reverse=TRUE))+
       geom_point(data=actual.size,aes(x=as.numeric(Extracted_nodes), y=as.numeric(Edges), label=Region), colour="red", size=4)+       
#       geom_text(data=actual.size,aes(x=as.numeric(Extracted_nodes), y=as.numeric(Edges), label=Region), size=5, hjust=0, vjust=1.5)+
    theme_bw()+
    scale_y_continuous("Total number of edges")+
     xlab("Number of significant genes")+
    coord_cartesian(xlim=c(300,1000), ylim=c(1,950))+
    theme(plot.title = element_text(face="bold", size=16), 
          axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16, angle=90),
          axis.text = element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size=14),
          legend.text=element_text(size=14),
          legend.position=c(0,1), legend.justification=c(0,1),
          legend.key=element_blank())
)
fill.data <- big.size[big.size$percentile==0.5,]
fill.loess <- loess(Max_size~Nodes,fill.data,span=1)
print(
  ggplot(data=big.size, aes(x=Nodes, y=Max_size))+geom_smooth(aes(colour=percentile),method = "loess", size = 0.5,se=F,span=1)+
    geom_area(aes(x,y),data.frame(x=fill.data$Nodes,y=predict(fill.loess)),stat="identity", fill="gray",alpha=0.4)+
    scale_colour_manual(aes(percentile), values = qcolours)+
    guides(colour=guide_legend(reverse=TRUE))+
    geom_point(data=actual.size,aes(x=as.numeric(Extracted_nodes), y=as.numeric(largest_nodes), label=Region))+

  theme_bw()+
  coord_cartesian(xlim=c(300,1000), ylim=c(1,300))+
   
  scale_y_continuous("Size of largest\nconnected component")+
    xlab("Number of significant gene")+
  theme(plot.title = element_text(face="bold", size=16), 
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16, angle=90),
        axis.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=14),
        legend.text=element_text(size=14),
        legend.position=c(0,1), legend.justification=c(0,1),
        legend.key=element_blank())
)
fill.data <- big.full[big.full$percentile==0.5,]
fill.loess <- loess(giantSize~Edge,fill.data,span=0.5)
print(
    ggplot(data=big.full, aes(x=Edge, y=giantSize))+geom_smooth(aes(colour=as.factor(percentile)), method = "loess", size = 0.5,se=F, span=1)+
      geom_area(aes(x,y),data.frame(x=fill.data$Edge,y=predict(fill.loess)),stat="identity", fill="gray",alpha=0.4)+
      scale_colour_manual(aes(percentile), values = qcolours)+
      guides(colour=guide_legend(reverse=TRUE))+
      geom_point(data=actual.size,aes(x=as.numeric(Edges), y=as.numeric(largest_nodes), label=Region))+
            coord_cartesian(xlim=c(50,1000), ylim=c(1,300))+
            theme_bw()+
      scale_y_continuous("Size of largest\nconnected component")+
      xlab("Total number of edges")+
      theme(plot.title = element_text(face="bold", size=16), 
            axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16, angle=90),
            axis.text = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size=14),
            legend.text=element_text(size=14),
            legend.position=c(0,1), legend.justification=c(0,1),
            legend.key=element_blank())
    )
dev.off()
