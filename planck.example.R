rm(list=ls())
# I use rhdf5 from bioconductor.
# I just barely know how to do this, so the help is a good place to start.
library(rhdf5)
library(igraph)

# Look at the contents of an HDF5 file.
# These files are structured like a file system.
# At the lowest level (I think) are datasets.
# These can be put into groups, which can also be put into groups, up to /.
h5ls('planck_p100_1000trees.hdf5')

# Read in the whole thing as a list of lists/
planck.data <- h5read('planck_p100_1000trees.hdf5', '/')
# Like with any list, use $ to access a component, e.g.
hist(log(planck.data$forestHalos$nodeMass), xlab='Mass', main='')
# Can save as an R data file.
#save(planck.data, file='planck.RData')

# Some things have their own attributes, even if they don't have data.
csm <- h5readAttributes('planck_p100_1000trees.hdf5', '/cosmology')
units <- h5readAttributes('planck_p100_1000trees.hdf5', '/units')
sim <- h5readAttributes('planck_p100_1000trees.hdf5', '/simulation')
grp <- h5readAttributes('planck_p100_1000trees.hdf5', '/groupFinder')

# Don't have to read in the whole thing.
# The second argmuent in h5read tells R which piece of the structure to grab.
halos <- h5read('planck_p100_1000trees.hdf5', '/forestHalos')
index <- h5read('planck_p100_1000trees.hdf5', '/forestIndex')
# Again, usual list rules apply
hist(log(index$numberOfNodes), xlab='# of nodes', main='')

# Read in individual pieces as whatever they are.
# When you read in a dataset, it should show up as the correct type of object.
position <- h5read('planck_p100_1000trees.hdf5', '/forestHalos/position')
pairs(t(position[,1:1000]), labels=c('x','y','z'), gap=0)

# Show how the index info is used to access the node info
# Plot one of the smaller trees with igraph.
# Node 0 is the all of the others coalesce into by step 499, z=0.
tree.number <- 636
tree.index <- (1:index$numberOfNodes[tree.number]) + index$firstNode[tree.number]
gdf <- cbind(halos$nodeIndex[tree.index], halos$descendentIndex[tree.index])
gdf <- graph.data.frame(gdf[-1,])
plot(gdf)

# Stuff about time using the first tree
tree.number <- 1
tree.index <- (1:index$numberOfNodes[tree.number]) + index$firstNode[tree.number]
# Timestep, redshift, and expansion factor are related.
# redshift decreases as time moves forward
plot(halos$timestep[tree.index], halos$redshift[tree.index], 
     xlab='time', ylab='z', main='Redshift is like moving back in time.')

# redshift and expansion factor are a transformation of each other
plot(halos$redshift[tree.index], -1 + 1/halos$expansionFactor[tree.index], 
     xlab='z', ylab='1/a - 1', main='Redshift and expansionFactor are sort of inverses.')
abline(0,1)

# Look at the mass of a node changes over time as it merges with others
# Let's try a different tree
tree.number <- 10
tree.index <- (1:index$numberOfNodes[tree.number]) + index$firstNode[tree.number]
with(halos, {
    # Use the descendentIndex variable to get the bits we want.
    node <- nodeIndex[tree.index[index$numberOfNodes[tree.number]]]
    merge.index <- c()
    while(node > -1) {
        merge.index <- c(node+1, merge.index)
        node <- descendentIndex[tree.index[node+1]]
    }
    plot(redshift[tree.index[merge.index]], nodeMass[tree.index[merge.index]], 
         xlab='z', ylab='mass')
})
