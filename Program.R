rm(list = ls())
source("R/.Source.R")

# Default Setting
N = seq(45, 90, 15)
T0 = seq(30, 90, 10)
T1 = 10
BS = 1:20

DGP1_data <- getDGP(
  Process = 1,
  nthreads = 4,
  seed = 1,
  N = N,
  T0 = T0,
  T1 = T1,
  BS = BS
)
DGP1_data <- getQCM(data = DGP1_data, nthreads = 4, seed = 1)
DGP1_data <- getIndicator(data = DGP1_data, N = N, T0 = T0, T1 = T1, BS = BS)
getPlot(DGP1_data, N = 90, T0 = 80, T1 = 10)
