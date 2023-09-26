############The method of network interactions' stability based on motif analysis under different attacks in R
######1 motif under random attacks
library(parallel)
cl <- makeCluster(10) #defined the number of CPUs
clusterEvalQ(cl, {
  library(igraph)
  fun <- function(g) {    # 'g' is a network
    a1 <- numeric(0)
    for (i in 1:vcount(g)) {
      s1 <- as.data.frame(motifs(g,size=4)) #motifs: calculate the numbers of different subgraph; size=3 :3-nodes subgraph,generally size=3 or 4.
      s1 <- s1[5,1]              #!!! to calcualte  the numbers of different subgraphs; '[4,1]' should be altered when calculate different types of subgraphs.
      a1 <- c(a1, s1)
      node_to_remove <- sample(1:vcount(g), 1)
      g <- delete_vertices(g, node_to_remove)
    }
    return(a1)
  }
})
clusterExport(cl, "ina")    # ina: your networks
results <- parLapply(cl, 1:1000, function(i) {
  fun(ina)
})
#stopCluster(cl)
averages <- numeric(length = 547)     # length: the number of nodes present in the network
system.time(for (i in 1:1000) {
  values <- sapply(results, function(x) x[i])
  averages[i] <- mean(values) 
})
ss <- data.frame(averages)
write.table(ss,'rand_in_m41.txt',sep="\t") 
averages <- numeric(length = 547)
for (i in 1:1000) {
  values <- sapply(results, function(x) x[i])
  averages[i] <- sd(values) 
}
sssd <- data.frame(averages)
write.table(sssd,'rand_in_m41_sd.txt',sep="\t")

##2——motif under betweenness attacks
#3-nodes
motif_bet <- function(g) {                       # 'g' is a network
  natcon <- function(g) {
    s1 <- as.data.frame(motifs(g,size=3))   #motifs: calculate the numbers of different subgraph; size=3 :3-nodes subgraph,generally size=3 or 4.
    s1[4,1]                                #!!! to calcualte  the numbers of different subgraphs '[4,1]' should be altered
  }
  motif.attack <- function(g) {
    hubord <- order(rank(betweenness(g)), decreasing=TRUE)          #betweenness(g):under betweenness attack
    sapply(1:round(vcount(g)*.8), function(i) {
      ind <- hubord[1:i]
      tmp <- delete_vertices(g, V(g)$name[ind])
      natcon(tmp)
    })
  }
  motif_bet <- motif.attack(g)
  motif_bet
}
in_bet <- motif_bet(ina)             # ina: your networks
dat1 <- data.frame(
  network = c(rep('ina', length(in_bet)), rep('outa', length(out_bet))),
  'Proportion_of_removes_nodes' = c((1:length(in_bet))/length(in_bet), (1:length(out_bet))/length(out_bet)),
  'Natural_Connectivity' = c(in_bet, out_bet), check.names = FALSE
)
write.table(dat1,'m32_intentional_degree.txt',sep="\t")

##2——motif under degree attacks
#3-nods
nc_degree <- function(g) {
  natcon <- function(g) {
    s1 <- as.data.frame(motifs(g,size=3))
    s1[4,1]
  }
  nc.attack <- function(g) {
    hubord <- order(rank(degree(g)), decreasing=TRUE)
    sapply(1:round(vcount(g)*.8), function(i) {
      ind <- hubord[1:i]
      tmp <- delete_vertices(g, V(g)$name[ind])
      natcon(tmp)
    })
  }
  nc <- nc.attack(g)
  nc
}
in_bet <- nc_degree(ina)
out_bet <- nc_degree(outa)
dat1 <- data.frame(
  network = c(rep('ina', length(in_bet)), rep('outa', length(out_bet))),
  'Proportion_of_removes_nodes' = c((1:length(in_bet))/length(in_bet), (1:length(out_bet))/length(out_bet)),
  'Natural_Connectivity' = c(in_bet, out_bet), check.names = FALSE
)
write.table(dat1,'m32_degree.txt',sep="\t")







