parallelLapply <- function(nthreads, seed, ...){
  if(nthreads <= 1){
    if(!is.null(seed)) set.seed(seed)
    return(lapply(...))
  }else{
    cl <- parallel::makeCluster(nthreads)
    if(!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    ans <- parallel::parLapply(cl, ...)
    parallel::stopCluster(cl)
    return(ans)
  }
}

getQCM <- function(data, nthreads, seed){
  return(
    parallelLapply(nthreads, seed, data, function(x){
      lapply(x, function(x){
        lapply(x, function(x){
          object <- x
          N <- object$N; T0 <-object$T0;
          cat("N = ", N, ", T0 = ", T0, "\n")
          QCM.x <- t(object$y.ctfl)[,-1]
          QCMX.x <- QCM.x
          for(k in 1:object$K) QCMX.x <- cbind(QCMX.x, t(object$x[k,,]))
          y <- object$y.ctfl[1,]
          QCM.fit <- quantregForest::quantregForest(x = QCM.x[1:T0,], y = y[1:T0], nthreads = nthreads)
          QCM.y <- predict(QCM.fit, QCM.x, what = c(0.025, 0.5, 0.975))
          QCM.y <- cbind(QCM.y, mean = predict(QCM.fit, QCM.x, what = mean))              
          QCMX.fit <- quantregForest::quantregForest(x = QCMX.x[1:T0,], y = y[1:T0], nthreads = nthreads)
          QCMX.y <- predict(QCMX.fit, QCMX.x, what = c(0.025, 0.5, 0.975))
          QCMX.y <- cbind(QCMX.y, mean = predict(QCMX.fit, QCMX.x, what = mean)) 
          object$QCM.y <- QCM.y
          object$QCMX.y <- QCMX.y
          object
        })
      })
    }))
}
getDGP <- function(Process, nthreads, seed, N, T0, T1, BS){
  source(paste("R/DGP", Process, ".R", sep = ""))
  FUN <- paste("dataGeneratingProcess", Process, sep = "")
  return(parallelLapply(nthreads, seed, N, function(x, FUNC, T0, T1, BS){
    n <- x
    lapply(T0, function(x){
      t0 <- x
      lapply(BS, function(x){
        object <- do.call(
          FUNC,
          list(N = n, T0 = t0, T1 = T1)
        )
        object
      })
    })
  }, FUNC = eval(parse(text = paste("dataGeneratingProcess", Process, sep = ""))), 
  T0 = T0,
  T1 = T1,
  BS = BS))
}
# DGP1_data <- DGP(
#   Process = 1,
#   nthreads = 4,
#   seed = 1,
#   N = seq(60, 90, 10),
#   T0 = seq(10, 90, 10),
#   BS = 1:10
# )
# str(DGP1_data[[4]][[9]][[10]])
CP <- function(actl, pred, T0){
  return(as.numeric((actl[1, -1:-T0] - pred[-1:-T0, 3] <= 1) & 
               (actl[1, -1:-T0] - pred[-1:-T0, 1] >= 1)))
}
WCI <- function(actl, pred, T0){
  return((actl[1, -1:-T0] - pred[-1:-T0, 1]) - (actl[1, -1:-T0] - pred[-1:-T0, 3]))
}
# Squard Error
SE <- function(actl, pred, T0){
  return((actl[1, -1:-T0] - pred[-1:-T0, 4])^2)
}
# Absolute Deviation
AD <- function(actl, pred, T0){
  return(abs(actl[1, -1:-T0] - pred[-1:-T0, 2]))
}
getIndicator <- function(data, Process, N, T0, T1, BS){
  ans <- NULL
  for(n in 1:length(N)){
    for(t0 in 1:length(T0)){
      QCM.CP <- NULL; QCM.WCI <- NULL; QCM.MSE <- NULL; QCM.MAD <- NULL;
      QCMX.CP <- NULL; QCMX.WCI <- NULL; QCMX.MSE <- NULL; QCMX.MAD <- NULL;
      for(bs in BS){
        object <- data[[n]][[t0]][[bs]]
        QCM.CP <- rbind(QCM.CP, CP(object$y.actl, object$QCM.y, T0[t0]))
        QCM.WCI <- rbind(QCM.WCI, WCI(object$y.actl, object$QCM.y, T0[t0]))
        QCM.MSE <- rbind(QCM.MSE, SE(object$y.actl, object$QCM.y, T0[t0]))
        QCM.MAD <- rbind(QCM.MAD, AD(object$y.actl, object$QCM.y, T0[t0]))
        QCMX.CP <- rbind(QCMX.CP, CP(object$y.actl, object$QCMX.y, T0[t0]))
        QCMX.WCI <- rbind(QCMX.WCI, WCI(object$y.actl, object$QCMX.y, T0[t0]))
        QCMX.MSE <- rbind(QCMX.MSE, SE(object$y.actl, object$QCMX.y, T0[t0]))
        QCMX.MAD <- rbind(QCMX.MAD, AD(object$y.actl, object$QCMX.y, T0[t0]))
      }
      ans <- rbind(ans, data.frame(
        Process = rep(paste0("DGP", Process), 4*T1),
        N = rep(N[n], 4*T1), 
        T0 = rep(T0[t0], 4*T1), 
        Time = rep(1:T1, 4),
        Indicator = c(rep("CP", 10), rep("WCI", 10), rep("MSE", 10), rep("MAD", 10)),
        QCM = c(
          apply(QCM.CP, 2, mean),
          apply(QCM.WCI,2, mean),
          apply(QCM.MSE, 2, mean),
          apply(QCM.MAD, 2, median)),
        QCMX = c(
          apply(QCMX.CP, 2, mean),
          apply(QCMX.WCI, 2, mean),
          apply(QCMX.MSE, 2, mean),
          apply(QCMX.MAD, 2, median)
        )))
    }
  }
  return(ans)
}
getPlot <- function(data, N, T0, T1){
  ggplot(data = data[data$N == N & data$T0 == T0, ][,-1:-2]) +
    geom_line(aes(x = Time, y = QCM, group = Indicator, colour = "QCM")) +
    geom_point(aes(x = Time, y = QCM, group = Indicator, colour = "QCM")) +
    geom_line(aes(x = Time, y = QCMX, group = Indicator, colour = "QCMX")) + 
    geom_point(aes(x = Time, y = QCMX, group = Indicator, colour = "QCMX")) +
    facet_wrap(.~Indicator, scales = "free_y") +
    xlab("Post-treatment period") + 
    ylab("")+
    scale_colour_hue("")+
    scale_x_continuous(breaks = 1:T1)
}