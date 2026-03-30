library(data.table)
library(tidyverse)
library(igraph)
library(ergm)
library(ergm.tapered)

# Load in the 'test' set 

library(readxl)
heldout <- read_xlsx('test_data.xlsx')
glimpse(heldout)

# Mask the test edges in the adjacency matrix to create training set 

adj <- as.matrix.network(g_net) # adjacency matrix for the full network
adj_hold <- adj

for(k in seq_len(nrow(heldout))) {
  i <- heldout$i[k]
  j <- heldout$j[k]
  adj_hold[i,j] <- NA
  adj_hold[j,i] <- NA
} # mask

g_hold <- network(adj_hold, directed=FALSE, matrix.type="adjacency") # create the training graph

# Copy vertex attributes 
for(a in list.vertex.attributes(g_net)) {
  g_hold %v% a <- g_net %v% a
}

# Fit the ERGM

m <- ergm.tapered(
  g_hold ~ edges
   + gwdegree(1.25, fixed=TRUE)
   # gwesp(0.1, fixed = TRUE) 
  + degree(1:4) 
  + gwnsp(0.5, fixed=TRUE)
    + nodematch('nucleus') +
    nodematch('membrane') +
    nodematch('endomembrane_secretory') +
    nodematch('complexes') +
    nodematch('cytoplasm') +
    nodematch('mitochondrion') +
    nodematch('other')
  ,
  control=control.ergm.tapered(
    MCMLE.maxit=60,
    MCMC.samplesize=20000,
    MCMC.interval=1000,
    MCMC.burnin=10000,
    parallel=4
  ),
  verbose=TRUE,
  estimate = 'MPLE'
) # Best model. THe other specifications were simply added into the m object when computing their corresponding predictions

summary(m) # parameter estimates

# Batch simulation (500 per batch)

get_presence <- function(net, heldout) {
  el <- as.edgelist(net) # convert net to edgelist 
  key <- paste(el[,1], el[,2], sep="_") # convert each to a string with _ seperating the 2 nodes
  key <- c(key, paste(el[,2], el[,1], sep="_")) # handle the reverse (because net is undirected)
  query <- paste(heldout$i, heldout$j, sep="_") # converts heldout edges to same format 
  as.integer(query %in% key)  # check if the heldout edges are 0 or 1 in the network
}


batch_size <- 500
checkpoint_file <- "hope_cumulative_counts_m3.csv"

if(file.exists(checkpoint_file)) {    
  cat("Load previous cumulative counts\n") 
  checkpoint <- read.csv(checkpoint_file) # load in previous batch if not at iteration 0 
  
  presence_sum <- checkpoint$presence_sum 
  sims_done <- unique(checkpoint$sims_done) # compute predictions up until current iteration 
  
} else {
  
  cat("Start first batch\n")
  
  presence_sum <- rep(0, nrow(heldout))
  sims_done <- 0 # initialize in the case of iteration 0
}

cat("Run 500 new simulations\n")

set.seed(as.integer(Sys.time())) # important to make sure we have truly different seed (previous versions did not work because the seed was the same at each batch, leading to degenerate batches)

sims_batch <- simulate(
  m,
  nsim = batch_size,
  output = "network",
  constraints = ~ observed,
  control = control.simulate.ergm(parallel = 4)) # simulate 500 networks (to be more precise conditional simulation on observed 'training' edges)

sim_mat_batch <- sapply(sims_batch, get_presence, heldout = heldout) # compute predictions

presence_sum <- presence_sum + rowSums(sim_mat_batch) # add to existing predictions 
sims_done <- sims_done + batch_size # update the true number of simulations


write.csv(
  data.frame(
    i = heldout$i,
    j = heldout$j,
    label = heldout$label,
    presence_sum = presence_sum,
    sims_done = sims_done
  ),
  checkpoint_file,
  row.names = FALSE
) # save to XL file

cat("Batch complete. Total simulations so far:",sims_done, "\n")

rm(sims_batch, sim_mat_batch) # free up memory (I personally found that this still led to slower batches after the first batch so I recommend closing R and running everything from scratch)

