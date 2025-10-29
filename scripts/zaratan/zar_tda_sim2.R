library(TDA)

rips_sim <- function (i) {
  for (mya in 10:0) {
    data <- read.csv(paste0("data/simulate/sim", i, "/trait_", mya, "mya.csv"))
    tic <- proc.time() # record time
    tda <- ripsDiag(X = data, maxdimension = 1, maxscale = 10, dist = "euclidean", library = "Dionysus", location = T, printProgress = T)
    saveRDS(tda, file = paste0("data/simulate/sim", i, "/tda_", mya, "mya.rds"))
    print(paste0("sim", i, ", ", mya, " Mya, Elapsed time ", (proc.time() - tic)[3])) # print runtime
  }
}
rips_sim(2)
