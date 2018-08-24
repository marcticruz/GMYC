# Load packages   

library(ape)
library(splits)
library(spider)
library(scales)  
library(phangorn) 
library(grDevices)  
library(base)     

# Set directory 

setwd("C:\\Users\\SAMSUNG\\Dropbox\\Oligo_barcode_manuscript\\analise2\\gmyc")

# Import newick file 

trcoalescent <- read.tree(file = "oligo_gmyc_tree.nwk") 

# Remove outgroups  

trcoalescent$tip.label

tip <- c(148, 149, 150)
trcoalescent <- drop.tip(trcoalescent, tip, trim.internal = TRUE)

# Ladderize the tree 

trcoalescent <- ladderize(trcoalescent , right = TRUE)

   
# Function to correct rounding erros for branch lengths values

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}   

trcoalescent <- force.ultrametric(trcoalescent)  

# Run gmyc analysis on the ultrametric tree   

test6 <- gmyc(trcoalescent, method= "single", interval = c(0, 10))

# A more flexible function to plot and edit the ultrametric tree that passed through gmyc analysis

plot.result <- function(test6, cex = 0.5, edge.width = 1, no.margin = F, show.tip.label = T, label.offset = 0){  
  plot.tree <- function(trcoalescent, MRCA, cex = 0.5, edge.width = 1, no.margin = F, show.tip.label = T, label.offset = 0){ 
    traverse.subtree <-  function(trcoalescent, n = 1){    
      numtip <- length(trcoalescent$tip.label)  
      sub.node <- trcoalescent$edge[which(trcoalescent$edge[,1] == n + numtip), 2] 
      
      test6.node <- c()
      for(sn in sub.node){  
        test6.node <- c(test6.node, sn)  
        test6.node <- c(test6.node, traverse.subtree(trcoalescent, sn-numtip))  
      }   
      return(test6.node)  
    }  
    
    numtip <- length(trcoalescent$tip.label)  
    br.col <- rep(1, length(trcoalescent$edge.length))  
    
    for (i in MRCA){  
      for (j in traverse.subtree(trcoalescent, i-numtip)){  
        br.col[which(trcoalescent$edge[,2] == j)] <- "blue"
      }   
    }
    plot(trcoalescent, edge.color = br.col, show.tip.label = show.tip.label, cex = cex, edge.width = edge.width, no.margin = no.margin, label.offset = label.offset) 
  }   
  
  plot.tree(test6$tree, test6$MRCA[[which.max(test6$likelihood)]] + length(test6$tree$tip.label), cex = cex, edge.width = edge.width, no.margin = no.margin, show.tip.label = show.tip.label, label.offset = label.offset) 
} 




# Plot the ultrametric tree with the clustes identified by GMYC  

windows(width = 85, height = 120, rescale = "fixed")

plot.result(test6, cex = 4, edge.width = 8, no.margin = FALSE) 


# Represent posterior probability > 0.95 as black circles

trcoalescent$node.label[trcoalescent$node.label >= 0.95] <- "black"  
trcoalescent$node.label[trcoalescent$node.label < 0.95 & trcoalescent$node.label >=0.75] <- alpha("grey", alpha = 0)  
trcoalescent$node.label[trcoalescent$node.label < 0.75] <- alpha("blue", alpha = 0) 

nodelabels(pch = 21, bg = trcoalescent$node.label, cex = 8, frame = "n", lwd = 0, col =  alpha("white", alpha = 0))

# Add a scale

axisPhylo()

# Identify the maximum branch time value

maxbranchtime <- max(branching.times(test6$tree))

# Subtract the maximum branch time value by the threshold time identified by GMYC 

lineposition <- maxbranchtime + ((test6$threshold.time[which.max(test6$likelihood)]))

# Plot a vertical line that represents the threshold time infered by gmyc

lines(x = c(lineposition, lineposition), y = c(0,273), col = "red", lwd = 8, lty = "dashed")   

# Subtract the maximum branch time value by the minimum value for threshold time error     

lineposition1 <- maxbranchtime - 0.0069

# Plot a vertical line representing the minimum value for the threshold time error   

lines(x = c(lineposition1, lineposition1), y = c(0,273), col = "black", lwd = 8)

# Subtract the maximum brach time value by the maximum value for the threshold time error

lineposition2 <- maxbranchtime - 0.0149

# Plot a vertical line representing the maximum value for the threshold time error  

lines(x = c(lineposition2, lineposition2), y = c(0,273), col = "black", lwd = 8)
