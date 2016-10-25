
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

# maximum redshift
maxredshift     <- rep(NA,ntrees)
for( t in 1:ntrees ){
    inds <- idxvec[t]+(1:nnodes[t])
    maxredshift[t] <-  max(planck.data$forestHalos$redshift[inds])
}

# number of branches
# look for nodes with no descendents
nbranches <- rep(NA,ntrees)
for( t in 1:ntrees ){
    inds <- idxvec[t]+(1:nnodes[t])
    nbranches[t] <-  sum(planck.data$forestHalos$descendentIndex[inds] == 0)
}

masslabs = seq(11,14)
redlabs = seq(-1,1)
branchlabs = seq(0,3)
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(log10(treemass),log10(maxredshift),pch=16,cex=0.5,xlim=c(11,14),ylim=c(-1,1),
     axes=FALSE,xlab="Tree Mass",ylab="Max. Redshift")
box()
axis(1,at=masslabs,labels=c(expression(10^11),expression(10^12),
                            expression(10^13),expression(10^14)))
axis(2,at=redlabs,labels=c(expression(10^{-1}),expression(10^0),
                            expression(10^1)))

plot(log10(treemass),log10(nbranches),pch=16,cex=0.5,xlim=c(11,14),ylim=c(0,3),
     axes=FALSE,xlab="Tree Mass",ylab="# of Branches")
box()
axis(1,at=masslabs,labels=c(expression(10^11),expression(10^12),
                            expression(10^13),expression(10^14)))
axis(2,at=branchlabs,labels=c(expression(10^{0}),expression(10^1),
                           expression(10^2),expression(10^3)))



# Let n be the total number of trees. I am thinking that we want to run 
# galacticus with M different parameter settings, where M is large but less than n.
# One approach is to split the merger trees into M representative
# groups and run galacticus with the same setting on the approximately
# n/M trees within each group. Assume for now that (n/M) is an integer.

# there are many ways to split the trees into groups, but here are two
# based on tree mass and max redshift.

# 1. Simplest version: randomize within blocks defined by tree mass
# Divide the tree mass domain into (n/M) subintervals, where the i'th breakpoint 
# is the i/(n/M) quantile of the merger tree mass distribution. This gives approximately
# M trees in each subinterval. Take one tree randomly from each subinterval.

M = round(ntrees/32)
nbins = ceiling(ntrees/M)
ord = order(treemass)
ordplusNA = c(ord,rep(NA,M*nbins - length(ord)))
binmatrix = matrix(ordplusNA,nbins,M,byrow=TRUE)
for(j in 1:nbins){  # reorder columns within each row
    binmatrix[j,] = sample(binmatrix[j,])
}
# now each of M columns is a representative sample with respect to treemass

# plot masses as a function of column
plot(0,xlim=c(1,M),ylim=c(11,14),xlab="group",ylab="log 10 tree mass")
for(j in 1:M) points(rep(j,nbins),log10(treemass[binmatrix[,j]]))


# 2. Slightly more complicted. Latin Square in tree mass and redshift.
# After splitting the trees into mass bins, order trees within each subinterval by 
# redshift. Then force each group to take one tree from each subinterval,
# subject to the constraint that you cannot select the same redshift rank twice.
# Doing this exactly requires M = sqrt(n), but you can create heuristics to 
# mimic this design when M \neq sqrt(n).

# haven't written anything for this yet



























# other stuff, not used yet for anything

# position, velocity, redshift, expansion factor
treeposition <- matrix(NA,ntrees,3)
velocity     <- matrix(NA,ntrees,3)
avgexpfac    <- rep(NA,ntrees)

for( t in 1:ntrees ){
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

