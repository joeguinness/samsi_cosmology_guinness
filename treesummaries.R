
rm(list=ls())
# I use rhdf5 from bioconductor.
library(rhdf5)

# Look at the contents of an HDF5 file.
# These files are structured like a file system.
# At the lowest level (I think) are datasets.
# These can be put into groups, which can also be put into groups, up to /.
h5ls('planck_p100_1000trees.hdf5')

# Read in the whole thing as a list of lists/
planck.data <- h5read('planck_p100_1000trees.hdf5', '/')

# get (0 indexed) indices of first nodes of each tree
idxvec <- planck.data$forestIndex$firstNode

# get number of trees and number of nodes of each tree
nnodes <- planck.data$forestIndex$numberOfNodes
ntrees <- length(nnodes)

# tree masses
treemass <- planck.data$forestHalos$nodeMass[idxvec+1]

# position, velocity, redshift, expansion factor
treeposition <- matrix(NA,ntrees,3)
velocity     <- matrix(NA,ntrees,3)
redshift     <- rep(NA,ntrees)
avgexpfac    <- rep(NA,ntrees)

for( t in 1:ntrees ){
    inds <- idxvec[t]+(1:nnodes[t])
    treeposition[t,] <-  rowMeans(planck.data$forestHalos$position[,inds])
    redshift[t]      <-  mean(planck.data$forestHalos$redshift[inds])
    velocity[t,]     <-  rowMeans(planck.data$forestHalos$velocity[,inds])
    avgexpfac[t]     <-  mean(planck.data$forestHalos$expansionFactor[inds])
}

treesummaries <- data.frame(
    numnodes = nnodes,
    mass = treemass,
    pos1 = treeposition[,1],
    pos2 = treeposition[,2],
    pos3 = treeposition[,3],
    vel1 = velocity[,1],
    vel2 = velocity[,2],
    vel3 = velocity[,3]
) 


write.csv(treesummaries,file="treedataframe.csv",row.names=FALSE)

