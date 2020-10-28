rm(list = ls())
library(reshape2)
library(ggplot2)
source("R/.Source.R")

# Default Setting
N <- c(40,70)
T0 <- c(40,70)
T1 <- 10
BS <- 1:2000
nthreads <- 2
seed <- 1
Process <- 1:9

for(i in Process){
  file <- paste0("DGP", i, "_data", "(BS=", max(BS), ").Rdata")
  if(!file %in% dir(path = "data")){
    object <- getQCM(getDGP(
      Process = i,
      nthreads = nthreads,
      seed = seed,
      N = N,
      T0 = T0,
      T1 = T1,
      BS = BS),nthreads = nthreads, seed = seed)
    save(object, file = paste0("data/", file))
  }else load(paste("data/", file, sep = ""))
  
  result <- getIndicator(
    data = object, 
    Process = i, 
    N = N, 
    T0 = T0, 
    T1 = T1,
    BS = BS)
  openxlsx::write.xlsx(result, file = paste0("DGP", i, ".xlsx"))
  
  temp <- melt(result, id.vars = c("Process", "N", "T0", "Time","Indicator"))
  temp <- dcast(temp, Process + N + T0 + Time ~ Indicator + variable)
  openxlsx::write.xlsx(temp, file = paste0("table/DGP", i, "_table.xlsx"))
  
  for(n in N) for(t0 in T0){
    getPlot(result, N = n, T0 = t0, T1 = T1)
    ggsave(paste0("photo/DGP", i, ",N=", n, ",T0=", t0, ".png"),
           width = 8, height = 4.5)}
}
