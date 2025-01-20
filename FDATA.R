## Andi Tenri Ola

## EWMA
## ============================================================================================== ##
setClass("EWMA", slots = list(output = "ANY"))
setMethod("initialize", "EWMA", function(.Object, data, lambda, L, plot.ewma, out.control){
  EWMA <- c()
  UCL  <- c()
  LCL  <- c()
  miu  <- mean(data)
  CL   <- miu
  sig  <- sd(data)
  lmd  <- lambda
  data <- as.vector(data)
  n    <- length(data)
  
  EWMA[1] <- (lmd * data[1]) + (1 - lmd)*miu 
  UCL[1]  <- miu + (L * sig)*sqrt((lmd*(1 - ((1-lmd)^(2*1))))/(2 - lmd))
  LCL[1]  <- miu - (L * sig)*sqrt((lmd*(1 - ((1-lmd)^(2*1))))/(2 - lmd))
  
  for (i in 2:n) {
    EWMA[i] <- (lmd * data[i]) + (1 - lmd)*EWMA[i-1] 
    UCL[i]  <- miu + (L * sig)*sqrt((lmd*(1 - ((1-lmd)^(2*i))))/(2 - lmd))
    LCL[i]  <- miu - (L * sig)*sqrt((lmd*(1 - ((1-lmd)^(2*i))))/(2 - lmd))
  }
  
  ooc <- 0
  ooc.index <- c()
  tmp.index <- 0
  for (i in 1:n) {
    if(EWMA[i] > UCL[i] || EWMA[i] < LCL[i]){
      tmp.index <- tmp.index + 1
      ooc <- ooc + 1
      ooc.index[tmp.index] <- i
    }
  }
  
  if(out.control == T){
    message("==================================================================")
    message(" EWMA")
    message("------------------------------------------------------------------")
    message("Mean UCL   : ",mean(UCL))
    message("Mean LCL   : ",mean(LCL))
    message("Out-control: ",ooc)
    message("index      :\n", paste0(ooc.index, sep = " "))
    message("==================================================================")
  }
  
  hasil.ewma <- data.frame(
    index = c(1:n),
    EWMA,
    UCL,
    CL,
    LCL
  )
  
  hasil.ewma <- hasil.ewma %>%
    mutate(
      Cp_color = ifelse((EWMA > UCL | EWMA < LCL), "out control", "EWMA")
    )
  
  ggplot(hasil.ewma) +
    geom_line(aes(x = index, y = UCL, color = "UCL"), size = 0.8) +
    geom_line(aes(x = index, y = EWMA, color = "EWMA"), size = 0.8) +
    geom_point(aes(x = index, y = EWMA, color = "EWMA")) +
    geom_point(aes(x = index, y = EWMA, color = Cp_color)) +
    geom_line(aes(x = index, y = CL, color = "CL"), size = 0.8) +
    geom_line(aes(x = index, y = LCL, color = "LCL"), size = 0.8)+
    labs(
      title = "Plot EWMA",
      subtitle = paste0("lambda = ",lmd),
      y = "EWMA",
      x = "index"
    ) +
    scale_color_manual(
      name = "Keterangan:",
      values = c(
        "EWMA" = "#000000",
        "UCL" = "#f35e5a",
        "CL" = "#17b5b3",
        "LCL" = "#0000F0",
        "out control" = "red"
      ),
      breaks = c("EWMA", "UCL", "CL", "LCL", "out control")
    )->pewma
  
  if(plot.ewma == T){
    print(pewma)
  }
  
  output <- NULL
  output$EWMA <- hasil.ewma
  output$plot <- pewma
  output$total.out <- ooc
  output$index.out <- ooc.index
  
  .Object@output <- output
  
  return(.Object)
  
})
rewma <- function(data, lambda = 0.1, L = 3, plot.ewma = F, out.control = T){
  new("EWMA", data, lambda, L, plot.ewma, out.control)@output
}
## ============================================================================================== ##

## DEWMA
## ============================================================================================== ##
setClass("DEWMA", slots = list(output = "ANY"))
setMethod("initialize", "DEWMA", function(.Object, data, lambda, L, sig.value, konstant, plot.dewma, out.control){
  DEWMA<- c()
  miu  <- mean(data)
  CL   <- miu
  data <- as.vector(data)
  n    <- length(data)
  z    <- rewma(data, lambda, out.control = F)$EWMA$EWMA
  lmd  <- lambda
  sig  <- sd(data)
  
  if(sig.value == T){
    sig2  <- sig^2
    var.z <- c();for (i in 1:n) {var.z[i] <- ((lmd^4)*(1+((1-lmd)^2)-(((i+1)^2)*((1-lmd)^(2*i)))+(((2*(i^2))+(2*i)-1)*((1-lmd)^((2*i)+2)))-((i^2)*((1-lmd)^((2*i)+4))))*(sig2))/((1-(1-lmd)^2)^3)}
    std.z <- sqrt(var.z)
    UCL   <- miu + (L * std.z)
    LCL   <- miu - (L * std.z)
  } else {
    sig2   <- 1
    var.z <- c();for (i in 1:n) {var.z[i] <- ((lmd^4)*(1+((1-lmd)^2)-(((i+1)^2)*((1-lmd)^(2*i)))+(((2*(i^2))+(2*i)-1)*((1-lmd)^((2*i)+2)))-((i^2)*((1-lmd)^((2*i)+4))))*(sig2))/((1-(1-lmd)^2)^3)}
    std.z <- sqrt(var.z)
    UCL   <- miu + (L * sig * std.z)
    LCL   <- miu - (L * sig * std.z)
  }
  
  if(konstant == T){
    ## untuk varinsi konstant
    var.c <- (lmd * (lmd^2 - 2*lmd + 2) *sig2) / ((2-lmd)^3)
    std.c <- sqrt(var.c)
    UCL   <- miu + (L * std.c)
    LCL   <- miu - (L * std.c)
  }
  
  DEWMA[1] <- (lmd * z[1]) + (1 - lmd)*miu
  for (i in 2:n) {
    DEWMA[i] <- (lmd * z[i]) + (1 - lmd)*DEWMA[i-1]
  }
  
  ooc <- 0
  ooc.index <- c()
  tmp.index <- 0
  for (i in 1:n) {
    if(DEWMA[i] > UCL[i] || DEWMA[i] < LCL[i]){
      tmp.index <- tmp.index + 1
      ooc <- ooc + 1
      ooc.index[tmp.index] <- i
    }
  }
  
  if(out.control){
    message("==================================================================")
    message(" DEWMA")
    message("------------------------------------------------------------------")
    message("Mean UCL   : ",mean(UCL))
    message("Mean LCL   : ",mean(LCL))
    message("Out-control: ",ooc)
    message("index      :\n", paste0(ooc.index, sep = " "))
    message("==================================================================")
  }
  
  hasil.dewma <- data.frame(
    index = c(1:n),
    DEWMA,
    UCL,
    CL,
    LCL
  )
  
  hasil.dewma <- hasil.dewma %>%
    mutate(
      Cp_color = ifelse((DEWMA > UCL | DEWMA < LCL), "out control", "DEWMA")
    )
  
  ggplot(hasil.dewma) +
    geom_line(aes(x = index, y = UCL, color = "UCL"), size = 0.8) +
    geom_line(aes(x = index, y = DEWMA, color = "DEWMA"), size = 0.8) +
    geom_point(aes(x = index, y = DEWMA, color = Cp_color)) +
    geom_line(aes(x = index, y = CL, color = "CL"), size = 0.8) +
    geom_line(aes(x = index, y = LCL, color = "LCL"), size = 0.8)+
    labs(
      title = "Plot DEWMA",
      subtitle = paste0("lambda = ",lmd),
      y = "DEWMA",
      x = "index"
    ) +
    scale_color_manual(
      name = "Keterangan:",
      values = c(
        "DEWMA" = "#000000",
        "UCL" = "#f35e5a",
        "CL" = "#17b5b3",
        "LCL" = "#0000F0",
        "out control" = "red"
      ),
      breaks = c("DEWMA", "UCL", "CL", "LCL", "out control")
    )->pdewma
  
  if(plot.dewma == T){
    print(pdewma)
  }
  
  # output
  output <- NULL
  output$DEWMA <- hasil.dewma
  output$plot  <- pdewma
  
  output$total.out <- ooc
  output$index.out <- ooc.index
  
  .Object@output <- output
  
  return(.Object)
})
rdewma.sign <- function(data, lambda = 0.1, L = 3, sig.value = T, konstant = F, plot.dewma = F, out.control = T){
  new("DEWMA",data, lambda, L, sig.value, konstant, plot.dewma, out.control)@output
}
## ============================================================================================== ##


## CUSUM ##
## ============================================================================================== ##
setClass("CUSUM", slots = list(output = "ANY"))
setMethod("initialize", "CUSUM", function(.Object, data, lambda, h, theta, plot.cusum, out.control){
  data<- as.vector(data)
  miu <- mean(data)
  sig <- sd(data)
  lmd <- lambda
  n   <- length(data)
  UCL <- H <- sig*h
  
  k     <- theta/2
  K     <- k * sig
  miu0p <- (miu + K)
  miu0m <- (miu - K)
  
  Cp <- c()
  Cm <- c()
  
  Cp[1] <- max(0, data[1]-miu0p+0)
  Cm[1] <- max(0, miu0m-data[1]+0)
  
  for (i in 2:n) {
    Cp[i] <- max(0, data[i]-miu0p+Cp[i-1])
    Cm[i] <- max(0, miu0m-data[i]+Cm[i-1])
  }
  
  ooc.Cp <- 0
  ooc.Cm <- 0
  ooc.index.Cp <- c()
  tmp.index.Cp <- 0
  ooc.index.Cm <- c()
  tmp.index.Cm <- 0
  
  for (i in 1:n) {
    if(Cp[i] > UCL){
      tmp.index.Cp <- tmp.index.Cp + 1
      ooc.Cp <- ooc.Cp + 1
      ooc.index.Cp[tmp.index.Cp] <- i
    }
    
    if(Cm[i] > UCL){
      tmp.index.Cm <- tmp.index.Cm + 1
      ooc.Cm <- ooc.Cm + 1
      ooc.index.Cm[tmp.index.Cm] <- i
    }
  }
  
  if(out.control){
    message("==================================================================")
    message("CUSUM")
    message("------------------------------------------------------------------")
    message("Mean UCL   : ",mean(UCL))
    message("------------------------------")
    message("Sum of Out-control: ")
    message("------------------------------")
    message("C+   : ",ooc.Cp)
    message("index:\n", paste0(ooc.index.Cp," "))
    message("------------------------------------------------------------------")
    message("C- : ",ooc.Cm)
    message("index:\n", paste0(ooc.index.Cm," "))
    message("------------------------------------------------------------------")
    message("total: ",(ooc.Cp+ooc.Cm))
    message("==================================================================")
  }
  
  
  hasil.cusum <- data.frame(
    index = c(1:n),
    Cp,
    Cm,
    UCL
  )
  
  hasil.cusum <- hasil.cusum %>%
    mutate(
      Cp_color = ifelse(Cp > UCL, "out control", "C+"),
      Cm_color = ifelse(Cm > UCL, "out control", "C-")
    )
  
  ggplot(hasil.cusum) +
    geom_line(aes(x = index, y = UCL, color = "UCL"), size = 0.8, linetype = 5) +
    geom_line(aes(x = index, y = Cp, color = "C+"), size = 0.8) +
    geom_point(aes(x = index, y = Cp, color = Cp_color)) +
    geom_line(aes(x = index, y = Cm, color = "C-"), size = 0.8) +
    geom_point(aes(x = index, y = Cm, color = Cm_color)) +
    labs(
      title = "Plot CUSUM",
      subtitle = paste0("lambda = ", lmd, " h = ", h, " theta = ", theta),
      y = "CUSUM",
      x = "index"
    ) +
    scale_color_manual(
      name = "Keterangan:",
      values = c(
        "UCL" = "#f35e5a",
        "C+" = "#17b5b3",
        "C-" = "#808080",
        "out control" = "red"
      ),
      breaks = c("UCL", "C+", "C-", "out control")
    ) -> pcusum
  
  if(plot.cusum == T){
    print(pcusum)
  }
  
  # output
  output <- NULL
  output$CUSUM <- hasil.cusum
  output$plot  <- pcusum
  
  output$total.out.Cp <- ooc.Cp
  output$index.out.Cp <- ooc.index.Cp
  output$total.out.Cm <- ooc.Cm
  output$index.out.Cm <- ooc.index.Cm
  
  .Object@output <- output
  
  return(.Object)
})
rcusum.sign <- function(data, lambda = 0.2, h = 4, theta = 0.1, plot.cusum = T, out.control = T){
  new("CUSUM",data, lambda, h, theta, plot.cusum, out.control)@output
}
## ============================================================================================== ##


## MIX DEWMA CUSUM ##
## ============================================================================================== ##
setClass("MDC", slots = list(output = "ANY"))
setMethod("initialize", "MDC", function(.Object, data, lambda, h, L, theta, plot.cusum, out.control){
  data<- as.vector(data)
  miu <- mean(data)
  sig <- sd(data)
  lmd <- lambda
  n   <- length(data)
  UCL <- H <- sig*h
  
  k <- theta/2
  K <- k * sig
  
  DEWMA <- rdewma.sign(data, lambda = lmd, L = L, sig.value = T, konstant = F, plot.dewma = F, out.control = F)
  dewma <- DEWMA$DEWMA
  miu0p <- (miu + K)
  miu0m <- (miu - K)
  
  Cp <- c()
  Cm <- c()
  
  Cp[1] <- max(0, data[1]-miu0p+0)
  Cm[1] <- max(0, miu0m-data[1]+0)
  
  for (i in 2:n) {
    Cp[i] <- max(0, data[i]-miu0p+Cp[i-1])
    Cm[i] <- max(0, miu0m-data[i]+Cm[i-1])
  }
  
  ooc.Cp <- 0
  ooc.Cm <- 0
  ooc.index.Cp <- c()
  tmp.index.Cp <- 0
  ooc.index.Cm <- c()
  tmp.index.Cm <- 0
  
  for (i in 1:n) {
    if(Cp[i] > UCL){
      tmp.index.Cp <- tmp.index.Cp + 1
      ooc.Cp <- ooc.Cp + 1
      ooc.index.Cp[tmp.index.Cp] <- i
    }
    
    if(Cm[i] > UCL){
      tmp.index.Cm <- tmp.index.Cm + 1
      ooc.Cm <- ooc.Cm + 1
      ooc.index.Cm[tmp.index.Cm] <- i
    }
  }
  
  if(out.control){
    message("==================================================================")
    message("MIX DEWMA CUSUM")
    message("------------------------------------------------------------------")
    message("Mean UCL   : ",mean(UCL))
    message("------------------------------")
    message("Sum of Out-control: ")
    message("------------------------------")
    message("MDC+   : ",ooc.Cp)
    message("index:\n", paste0(ooc.index.Cp," "))
    message("------------------------------------------------------------------")
    message("MDC- : ",ooc.Cm)
    message("index:\n", paste0(ooc.index.Cm," "))
    message("------------------------------------------------------------------")
    message("total: ",(ooc.Cp+ooc.Cm))
    message("==================================================================")
  }
  
  
  hasil.cusum <- data.frame(
    index = c(1:n),
    Cp,
    Cm,
    UCL
  )
  
  hasil.cusum <- hasil.cusum %>%
    mutate(
      Cp_color = ifelse(Cp > UCL, "out control", "MDC+"),
      Cm_color = ifelse(Cm > UCL, "out control", "MDC-")
    )
  
  ggplot(hasil.cusum) +
    geom_line(aes(x = index, y = UCL, color = "UCL"), size = 0.8, linetype = 5) +
    geom_line(aes(x = index, y = Cp, color = "MDC+"), size = 0.8) +
    geom_point(aes(x = index, y = Cp, color = Cp_color)) +
    geom_line(aes(x = index, y = Cm, color = "MDC-"), size = 0.8) +
    geom_point(aes(x = index, y = Cm, color = Cm_color)) +
    labs(
      title = "Plot MIX DEWMA CUSUM",
      subtitle = paste0("lambda = ", lmd, " h = ", h, " theta = ", theta),
      y = "MDC",
      x = "index"
    ) +
    scale_color_manual(
      name = "Keterangan:",
      values = c(
        "UCL" = "#f35e5a",
        "MDC+" = "#17b5b3",
        "MDC-" = "#808080",
        "out control" = "red"
      ),
      breaks = c("UCL", "MDC+", "MDC-", "out control")
    ) -> pcusum
  
  if(plot.cusum == T){
    print(pcusum)
  }
  
  # output
  output <- NULL
  output$CUSUM <- hasil.cusum
  output$plot  <- pcusum
  
  output$total.out.Cp <- ooc.Cp
  output$index.out.Cp <- ooc.index.Cp
  output$total.out.Cm <- ooc.Cm
  output$index.out.Cm <- ooc.index.Cm
  
  .Object@output <- output
  
  return(.Object)
})
rmdc.sign <- function(data, lambda = 0.2, h = 4, L = 3,theta = 0.1, plot.cusum = T, out.control = T){
  new("MDC", data, lambda, h, L, theta, plot.cusum, out.control)@output
}
## ============================================================================================== ##



## ARL EWMA, DEWMA, TEWMA, CUSUM
## ============================================================================================== ##
setClass("ARL", slots = list(output="ANY"))
setMethod("initialize", "ARL", function(.Object, data, lambda, theta, L, h, info, max.iterasi, method, konstant.tewma){
  data <- as.vector(data)
  miu  <- mean(data)
  sig  <- sd(data)
  lmd  <- lambda
  
  n  <- length(data)
  mx <- max(data)
  pr <- mean(data)/n
  mi <- max.iterasi
  RL <- c()
  
  if(method[1] == "EWMA"){
    for (i in 1:mi) {
      set.seed(i)
      
      # tmpd  <- rnorm(n, miu, sig)
      tmpx <- rnorm(n, miu, sig)
      mode <- rewma(tmpx, lambda = lmd, L = L, out.control = F, plot.ewma = F)
      data <- mode$EWMA$EWMA
      outc <- mode$EWMA$Cp_color
      tmp0 <- 0
      RL[i]<- 0
      
      in.control <- TRUE
      while(in.control){
        tmp0 <- tmp0 + 1
        if(outc[tmp0] == "out control"){
          RL[i] <- tmp0 - 1
          in.control <- FALSE
        }
        if(tmp0 == n){in.control <- FALSE}
      }
    }
  }
  if(method[1] == "DEWMA"){
    for (i in 1:mi) {
      set.seed(i)
      
      # tmpd <- rnorm(n, miu, sig)
      tmpx <- rnorm(n, miu, sig)
      mode <- rdewma.sign(tmpx, lambda = lmd, L = L, out.control = F, plot.dewma = F)
      data <- mode$DEWMA$DEWMA
      outc <- mode$DEWMA$Cp_color
      tmp0 <- 0
      RL[i]<- 0
      
      in.control <- TRUE
      while(in.control){
        tmp0 <- tmp0 + 1
        if(outc[tmp0] == "out control"){
          RL[i] <- tmp0 - 1
          in.control <- FALSE
        }
        if(tmp0 == n){in.control <- FALSE}
      }
    }
  }
  if(method[1] == "CUSUM"){
    for (i in 1:mi) {
      set.seed(i)
      
      # tmpd  <- rnorm(n, miu, sig)
      tmpx  <- rnorm(n, miu, sig)
      mode  <- rcusum.sign(tmpx, lambda = lmd, h = h, theta = theta, plot.cusum = F, out.control = F)
      Cp    <- mode$CUSUM$Cp
      Cm    <- mode$CUSUM$Cm
      outCp <- mode$CUSUM$Cp_color
      outCm <- mode$CUSUM$Cm_color
      
      tmp0 <- 0
      RL[i]<- 0
      
      in.control <- TRUE
      while(in.control){
        tmp0 <- tmp0 + 1
        if(outCp[tmp0] == "out control" || outCm[tmp0] == "out control" ){
          RL[i] <- tmp0 - 1
          in.control <- FALSE
        }
        if(tmp0 == n){in.control <- FALSE}
      }
    }
  }
  
  ARL <- mean(RL);ARL
  if(info){
    message("Metode: ", method)
    message("ARL   : ",ARL)
  }
  
  .Object@output <- ARL
  
  return(.Object)
})
rarl <- function(data, lambda = 0.5, theta = 0.25, L = 3, h = 4, info = T, max.iterasi = 1000, method = c("EWMA", "DEWMA", "TEWMA", "CUSUM"), konstant.tewma = F){
  new("ARL", data, lambda, theta, L, h, info, max.iterasi, method, konstant.tewma)@output
}


## ARL MIX DEWMA CUSUM ##
## ============================================================================================== ##
setClass("ARLMDC", slots = list(output = "ANY"))
setMethod("initialize", "ARLMDC", function(.Object, data, lambda, theta, L, h, info, version, max.iterasi){
  data <- as.vector(data)
  miu  <- mean(data)
  sig  <- sd(data)
  lmd  <- lambda

  n  <- length(data)
  mx <- max(data)
  pr <- mean(data)/n
  mi <- max.iterasi
  RL <- c()
  
  for (i in 1:mi) {
    set.seed(i)
    
    tmpd  <- rnorm(n, miu, sig)
    tmpx  <- ifelse(tmpd > L,0,1)
    mode  <- rmdc.sign(tmpx, lambda = lmd, h = h, L = L, theta = theta, plot.cusum = F, out.control = F)
    MCEp  <- mode$CUSUM$Cp
    MCEm  <- mode$CUSUM$Cm
    outCp <- mode$CUSUM$Cp_color
    outCm <- mode$CUSUM$Cm_color
    tmp0  <- 0
    RL[i] <- 0
    
    in.control <- TRUE
    while(in.control){
      tmp0 <- tmp0 + 1
      if(outCp[tmp0] == "out control" || outCm[tmp0] == "out control" ){
        RL[i] <- tmp0 - 1
        in.control <- FALSE
      }
      if(tmp0 == n){in.control <- FALSE}
    }
  }
  
  ARL <- mean(RL);ARL
  if(info){
    message("Metode: MIX DEWMA CUSUM")
    message("ARL   : ",ARL)
  }
  
  .Object@output <- ARL
  
  return(.Object)
})
arl_mdc <- function(data, lambda = 0.5, theta = 0.25, L = 3, h = 4, info = T, version = 1, max.iterasi = 1000){
  new("ARLMDC", data, lambda, theta, L, h, info, version, max.iterasi)
}
## ============================================================================================== ##
