#!/usr/bin/env Rscript
library(vegan)
library(expss)
library(parallel)

nameDistToMatrix <- function(x){
  x.names <- sort(unique(c(x[[1]], x[[2]])))
  x.dist <- matrix(0, length(x.names), length(x.names))
  dimnames(x.dist) <- list(x.names, x.names)
  x.ind <- rbind(cbind(match(x[[1]], x.names), match(x[[2]], x.names)), cbind(match(x[[2]], x.names), match(x[[1]], x.names)))
  x.dist[x.ind] <- rep(x[[3]], 2)
  return(x.dist)
}

readInMatrix <- function(file_a){
 if(grepl("otu_table",file_a)){
   otu.table <- read.table(file_a, as.is=TRUE, header=TRUE, row.names = 1)
   otu.dist <- vegdist(t(otu.table), method="jaccard", diag=TRUE, upper=TRUE)
   return(otu.dist)
 } else {
   mash.table <- read.table(file_a, as.is=TRUE, header=FALSE)
   mash.dist <- as.dist(1-nameDistToMatrix(mash.table), diag =TRUE, upper =TRUE)
   return(mash.dist)
 }
}

readAndMantel <- function(file_a, file_b){
  file_a.dist <- readInMatrix(file_a)
  file_b.dist <- readInMatrix(file_b)
  a_b.mantel <- mantel(file_a.dist, file_b.dist, method="pearson")$statistic
  print(c(file_a,file_b,a_b.mantel))
  #return(c(file_a,file_b,a_b.mantel))
}

list_of_dist <- function(file_a){
  file_b <- sub("_otu_table","_mash_dists", file_a)
  otu.table <- read.table(file_a, as.is=TRUE, header=TRUE, row.names = 1)
  otu.dist <- vegdist(t(otu.table), method="jaccard", diag=TRUE, upper=TRUE)
  otu.dist.mean <- mean(otu.dist)
  mash.table <- read.table(file_b, as.is=TRUE, header=FALSE)
  mash.dist <- as.dist(1-nameDistToMatrix(mash.table), diag=TRUE, upper=TRUE)
  mash.dist.mean <- mean(mash.dist)
  return(c(otu.dist.mean,mash.dist.mean))
}

color_me_not_suprised <- function(a_suprise,meta.data,colours){
  return(colours[vlookup_df(a_suprise, meta.data, lookup_column = 'File')$Id])
}

all_distances_plot <- function(file_a,meta.data,colours,cores){
  file_b <- sub("_otu_table","_mash_dists", file_a)
  otu.table <- read.table(file_a, as.is=TRUE, header=TRUE, row.names = 1)
  otu.dist <- vegdist(t(otu.table), method="jaccard", diag=TRUE, upper=TRUE)
  otu_matrix <- as.matrix(otu.dist)
  mash.table <- read.table(file_b, as.is=TRUE, header=FALSE)
  mash.dist <- as.dist(1-nameDistToMatrix(mash.table), diag =TRUE, upper =TRUE)
  mash_matrix <- as.matrix(mash.dist)
  mclapply(simplify2array(rownames(otu_matrix)), 
           function(x) sapply(colnames(otu_matrix), 
                              function(y) points(otu_matrix[x,y],mash_matrix[x,y],col=color_me_not_suprised(file_a,meta.data,colours)))
           , mc.cores = cores)
}


meta.data <- read.table( "meta.data.file.txt", as.is=TRUE, header=TRUE)
meta.data <- transform(meta.data,Id=as.numeric(factor(Type)))
meta.file.list <- list.files(path='.', pattern="_run_")
meta.file.list.otu <- list.files(path='.', pattern="_otu_table")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
  args[2] = 2
}

if( args[1] == "mantel"){
  #meta.cor.list <- 
  mclapply(simplify2array(meta.file.list), function(x) sapply(meta.file.list, function(y) readAndMantel(x,y)), mc.cores = args[2])
  #print(meta.cor.list)
} else if ( args[1] == "distances"){
  distances <- mclapply(simplify2array(meta.file.list.otu), function(x){
    list_of_dist(x)
  }, mc.cores = args[2])
  meta.data.distances <- vlookup_df(colnames(distances), meta.data, lookup_column = 'File')
  plot(distances[1,],distances[2,], col=factor(meta.data.distances$Type), xlab="OTU Method Distance", ylab="K-mer Method Distance", main="Mash and Qiime Avg Sample to Sample Distance of Simulated Datasets")
  legend("bottomleft", 
         legend=levels(factor(meta.data.distances$Type)), 
         col=factor(levels(factor(meta.data.distances$Type))),
         lwd=1, 
         lty=NA, 
         pch=1, 
         cex=0.8,
         y.intersp=0.5)
} else if ( args[1] == "all_distances"){
  png(filename = "all_distances.png", width = 1024, height = 1024)
  plot(1, type="n", xlab="OTU Method Distance", ylab="K-mer Method Distance", xlim=c(0, 1), ylim=c(0, 1))
  colours <- terrain.colors(25, alpha = 1)
  distances <- sapply(meta.file.list.otu, function(x){
    all_distances_plot(x,meta.data,colours,args[2])
  })
  dev.off()
} else {
  stop("Command line argument not recognized.", call.=FALSE)
}
