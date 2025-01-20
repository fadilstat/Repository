
## ============================================================================================== ##
## FISKA
## ============================================================================================== ##

## Normalitas
## ============================================================================================== ##
setClass("BinomialTest",slots = list(output = "ANY"))
setMethod("initialize", "BinomialTest", function(.Object, data, alpha){
  x <- as.vector(data)
  
  # motode
  km <- ks.test(x,'pbinom') # Kolmogorov-Smirnov
  
  ket     <- c()
  p_value <- c(km$p.value)
  stat    <- c(km$statistic, sp$statistic)
  ket[1]  <- ifelse(p_value[1] < alpha,'Tolak H0','Terima H0')
  
  hasil <- data.frame(
    uji = c('Kolmogorov-Smirnov', 'Shapiro-Wilk'),
    stat,
    p_value,
    ket
  )
  
  colnames(hasil) <- c('Motode', 'Statistik','p-value','Keptusan')
  
  message('Hipotesis:')
  cat('H0 : Data mengikuti distribusi normal\n')
  cat('H1 : Data tidak mengikuti distribusi normal\n')
  message('\nStatistik Uji:')
  print(hasil)
  message('\nKeputusan:')
  
  cat('Kolmogorov-Smirnov: ')
  if(p_value[1] < alpha){cat('Data tidak mengikuti distribusi normal\n')}else{cat('Data mengikuti distribusi normal\n')}
  cat('Shapiro-Wilk      : ')
  if(p_value[2] < alpha){cat('Data tidak mengikuti distribusi normal\n')}else{cat('Data mengikuti distribusi normal\n')}
  
  .Object@output <- hasil
  
  return(.Object)
})
binomial.test <- function(data, alpha = 0.05){
  new("BinomialTest", data = data, alpha = alpha)@output
}
## ============================================================================================== ##

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
    message(" EWMA Poisson")
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
      title = "Plot EWMA Poisson",
      subtitle = paste0("lambda = ",lmd),
      y = "EWMA Poisson",
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

## DEWMA Poisson
## ============================================================================================== ##
setClass("DEWMA", slots = list(output = "ANY"))
setMethod("initialize", "DEWMA", function(.Object, data, lambda, L, plot.dewma, out.control){
  data <- as.vector(data)
  miu  <- mean(data)
  std  <- miu
  CL   <- miu
  n    <- length(data)
  Y    <- rewma(data, lambda, out.control = F)$EWMA$EWMA
  lmd  <- lambda
  sig  <- sd(data)
  Z    <- c();Z[1]<-miu
  t    <- c()
  UCL  <- LCL <- c()
  
  for (i in 1:n) {
    if(i == 1){
      Z[i]<-miu
    } else {
      Z[i] <- (lmd*Y[i])+((1-lmd)*Z[i-1])
    }
    
    t_1  <- (1-lmd)^2
    t_2  <- ((i+1)^2)*((1-lmd)^(2*i))
    t_3  <- (((2*(i^2))+(2*i)-1)*((1-lmd)^((2*i)+2)))
    t_4  <- ((i^2)*((1-lmd)^((2*i)+4)))
    t[i] <- (lmd^4)*(1+(t_1-t_2+t_3-t_4))/(1-((1-lmd)^2)^3)
    
    UCL[i] <- miu + L * sqrt(std * (t[i]))
    LCL[i] <- miu - L * sqrt(std * (t[i]))
  }

  
  DEWMA <- Z
  
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
    message(" DEWMA Poisson")
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
      title = "Plot DEWMA Poisson",
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
rdewmap <- function(data, lambda = 0.1, L = 3, plot.dewma = F, out.control = T){
  new("DEWMA",data, lambda, L, plot.dewma, out.control)@output
}
## ============================================================================================== ##

## DEWMA Poisson
## ============================================================================================== ##
setClass("DEWMA.ISRT", slots = list(output = "ANY"))
setMethod("initialize", "DEWMA.ISRT", function(.Object, data, lambda, L, plot.dewma, out.control){
  da   <- as.vector(data)
  data <- sqrt(da)
  miu  <- mean(data)
  std  <- miu
  CL   <- miu
  n    <- length(data)
  Y    <- rewma(data, lambda, out.control = F)$EWMA$EWMA
  lmd  <- lambda
  sig  <- sd(data)
  Z    <- c();Z[1]<-miu
  t    <- c()
  UCL  <- LCL <- c()
  
  for (i in 1:n) {
    if(i == 1){
      Z[i]<-miu
    } else {
      Z[i] <- (lmd*Y[i])+((1-lmd)*Z[i-1])
    }
    
    t_1  <- (1-lmd)^2
    t_2  <- ((i+1)^2)*((1-lmd)^(2*i))
    t_3  <- (((2*(i^2))+(2*i)-1)*((1-lmd)^((2*i)+2)))
    t_4  <- ((i^2)*((1-lmd)^((2*i)+4)))
    t[i] <- (lmd^4)*(1+(t_1-t_2+t_3-t_4))/(1-((1-lmd)^2)^3)
    
    UCL[i] <- miu + L*(1/2 - (1/6 * miu)) * sqrt(std * (t[i]))
    LCL[i] <- miu - L*(1/2 - (1/6 * miu)) * sqrt(std * (t[i]))
  }
  
  
  DEWMA <- Z
  
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
    message("DEWMA ISRT Poisson")
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
      title = "Plot DEWMA ISRT Poisson",
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
rdewmap.isrt <- function(data, lambda = 0.1, L = 3, plot.dewma = F, out.control = T){
  new("DEWMA.ISRT",data, lambda, L, plot.dewma, out.control)@output
}
## ============================================================================================== ##



## ============================================================================================== ##
## ARL
## ============================================================================================== ##

## ARL EWMA, DEWMA, TEWMA, CUSUM
## ============================================================================================== ##
setClass("ARLP", slots = list(output="ANY"))
setMethod("initialize", "ARLP", function(.Object, data, lambda, L, info, max.iterasi, method){
  data <- as.vector(data)
  miu  <- mean(data)
  sig  <- mean(data)
  lmd  <- lambda
  
  n  <- length(data)
  mi <- max.iterasi
  RL <- c()
  
  if(method[1] == "EWMAP"){
    for (i in 1:mi) {
      set.seed(i)
      
      tmpx <- rpois(n, miu)
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
  if(method[1] == "DEWMAP"){
    for (i in 1:mi) {
      set.seed(i)
      
      tmpx <- rpois(n, miu)
      mode <- rdewmap(tmpx, lambda = lmd, L = L, out.control = F, plot.dewma = F)
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
  if(method[1] == "DEWMAP.ISRT"){
    for (i in 1:mi) {
      set.seed(i)
      
      tmpx <- rpois(n, miu)
      mode <- rdewmap.isrt(tmpx, lambda = lmd, L = L, out.control = F, plot.dewma = F)
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
  
  
  ARL <- mean(RL);ARL
  if(info){
    message("Metode: ", method)
    message("ARL   : ",ARL)
  }
  
  .Object@output <- ARL
  
  return(.Object)
})
rarl.p <- function(data, lambda = 0.5, L = 3, info = T, max.iterasi = 1000, method = c("EWMAP", "DEWMAP", "DEWMAP.ISRT")){
  new("ARLP", data, lambda, L, info, max.iterasi, method)@output
}
## ============================================================================================== ##


