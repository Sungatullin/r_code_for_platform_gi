Sys.setenv(LANG = "en")

suppressPackageStartupMessages({
library(stringr)
library(pracma)
library(dplyr)
library(magrittr)
})


fmax <- function(t0,lam,dat,zk,sigk,ro){
  sum(sqrt(dat[t0+ro])) - lam*sum(((zk-diff(t0+ro))^2)/(sigk+1))
}

fmovement <- function(t0,lam,dat,zk,sigk,ro,j){
  ro[j] <- t0 + ro[j]
  if(ro[j] < 0){
    ro[j] <- 1
  }
  return(sum(sqrt(dat[ro])) - lam*sum(((zk-diff(ro))^2)/(sigk+1)))
}

turnUp <- function(dat){
  a1 <- diff(dat)[-1]
  a2 <- rev(diff(rev(dat))[-1])
  a1 <- c(0,a1,0)
  a2 <- c(0,a2,0)
  a3 <- which(a1 > mean(a1) + 0.5*sd(a1) & a2 > mean(a2) + 0.5*sd(a2) & dat < 1.5e+6)
  dat[a3] <- (dat[a3-1]+dat[a3+1])/2
  return(dat)
}

turnUpNull <- function(dat){
  a3 <- which(dat == 0)
  n3 <- length(a3)
  a3_1 <- a3
  if (a3[1] == 1) a3_1[1] <- 2
  if (a3[n3] == length(dat)) a3_1[n3] <- length(dat)-1
  dat[a3] <- (dat[a3_1-1]+dat[a3_1+1])/2
  return(dat)
}

calculate_area <- function(idpeak,text){
  sCoef <- 1.016112e-07
  table_area_neeed <- idpeak
  table_area_neeed[,"area"] <- 0
  for(j in 1:nrow(table_area_neeed)){
    len_peak <- table_area_neeed[j,3]:table_area_neeed[j,4]
    
    table_area_neeed[j,"area"] <- 
      sCoef*(trapz(len_peak,text[len_peak])-
               (min(text[table_area_neeed[j,4]],text[table_area_neeed[j,3]])*
                  (table_area_neeed[j,4]-table_area_neeed[j,3])))
  }
  return(table_area_neeed)
}

`%p%` <- paste0

mu_index <- function(dat,LA=50){
  
  ro <- c(7280, 9018, 9120, 9405, 9489)
  
  mean_z <- c(1738, 101, 286, 84)
  
  var_sig <- c(18, 11, 23, 22)
  
  NNN <- 200
  datT0 <- data.frame("t0" = 1:NNN,
                      "sum" = rep(0,length(NNN)))
  datT0$sum<- sapply(datT0[,"t0"],
                     fmax, lam=LA, dat=dat,ro = ro,
                     zk = mean_z,sigk = var_sig)
  # print(datT0[which.max(datT0[,"sum"]),"t0"])
  pams_ro <- datT0[which.max(datT0[,"sum"]),"t0"]
  if(isempty(pams_ro)){
    pams <- ro
  }else{
    pams <- ro + pams_ro
    pams[1] <- ro[1]#pams_ro
  }
  
  
  lociloc <- pams
  prom_lociloc <- pams #+jghjgh
  
  datT0 <- data.frame("t0" = -50:50,
                      "sum" = rep(0,length(-50:50)))
  
  prom_lociloc <- prom_lociloc+105
  while(sum(abs(lociloc - prom_lociloc))>=10){
    j<-1
    for(i in 1:length(pams)){
      prom_lociloc[i] <- lociloc[i]
      
      datT0$sum<- sapply(datT0[,"t0"],
                         fmovement, lam=LA, dat=dat,ro = lociloc,
                         zk = mean_z,sigk = var_sig,j=i)
      prom_del <- datT0[which.max(datT0[,"sum"]),"t0"]
      if(isempty(prom_del)){
        next
      }else{
        if((prom_del + lociloc[i]) < 1){
          lociloc[i] <- pams[i]+50
        }else{
          lociloc[i] <- prom_del + lociloc[i]
        }
      }
      
    }
    j <- j+1
  }
  
  pams <- lociloc
  return(pams)
}

create_list<- function(len_vec=4){
  df <- list("Phen"=rep(0,len_vec),
             "X3MPh"=rep(0,len_vec),
             "X2MPh"=rep(0,len_vec),
             "X9MPh"=rep(0,len_vec),
             "X1MPh"=rep(0,len_vec),
             "MPI" = 0,
             "MPI1" = 0,
             "MPI2" = 0)
  
  return(df)
}

search_start_end_extr <- function(dat,smo){
  #for smooth data
  data <- smo$y
  eps <- 0.001
  N <- length(smo$y)
  log_vector <- smo$y[2:(N-1)]-smo$y[1:(N-2)] < eps & 
    smo$y[2:(N-1)]-smo$y[3:N] < eps
  
  min_of_peak  <- which(log_vector)
  
  ind_max_data <- mu_index(smo$y)
  
  ind_peak <- data.frame("Value_max" = data[ind_max_data],
                         "ind_Max" = ind_max_data, 
                         "ind_Start" = rep(0,length(ind_max_data)), 
                         "idn_End" = rep(0,length(ind_max_data)))
  
  for(i in 1:length(ind_max_data)){
    if(i == 1){
      max_peak <- sort(c(min_of_peak,ind_max_data[i]))
      index_start_end <- max_peak[which(max_peak == ind_max_data[i]) + c(-1,1)]
      ind_peak[i,3] <- index_start_end[1]
      ind_peak[i,4] <- index_start_end[2]
    }else{
      max_peak <- sort(c(min_of_peak,ind_max_data[i]))
      index_start_end <- max_peak[which(max_peak == ind_max_data[i]) + c(-1,1)]
      ind_peak[i,"Value_max"] <- data[ind_max_data[i]]
      ind_peak[i,3] <- index_start_end[1]
      ind_peak[i,4] <- index_start_end[2]
    }
  }
  
  vector_eps_for_peak <- matrix(0,ncol = nrow(ind_peak),nrow=2)
  vector_eps_for_peak[1,] <- ind_peak[,2]
  for(i in 1:nrow(ind_peak)){
    if(ind_peak[i,1]>1000000){
      vector_eps_for_peak[2,i]<- 500
    }
    if(ind_peak[i,1]<1000000){
      vector_eps_for_peak[2,i]<- 0.01
    }
    
  }
  ind_max_data <- vector_eps_for_peak[1,] #+"'khgjh"
  
  for(i in which(vector_eps_for_peak[2,i]!=0)){
    log_vector <- data[2:(N-1)]-data[1:(N-2)] < vector_eps_for_peak[2,i] & 
      data[2:(N-1)]-data[3:N] < vector_eps_for_peak[2,i]
    min_peak  <- which(log_vector)
    max_peak <- sort(c(min_peak,ind_max_data[i]))
    index_start_end <- max_peak[which(max_peak == ind_max_data[i]) + c(-1,1)]
    ind_peak[i,3] <- index_start_end[1]
    ind_peak[i,4] <- index_start_end[2]
  }
  
  # table_peak <- search_start_end_extr_old(smo$y)#+jhg
  
  return(ind_peak)
}

area_percent_peak <- function(text, name_well, file_error){
  
  vec_error_output <- c()
  vec_error_wells <- c()
  
  listFiltred <- list("m/z" = "178+192",
                      "intensity" = 0)
  listRes <- create_list()
  
  
      tryCatch(
        error = function(cnd) {
          vec_error_wells <<- c(vec_error_wells,i)
          
          cl1 <- cnd$call %>% deparse %>% paste0(collapse="") %>% gsub("  ", "", .) %>%
            gsub("\n", " ", .) %>% gsub("\r", " ", .)
          msg <- "Error in " %p% cl1 %p% ": " %p% cnd$message
          print(msg)
          vec_error_output <<- c(vec_error_output,paste("OIL 178,192mz",name_well, msg, sep = ": "))
        },
        {
          
          while(!isempty(which(text == 0))){
            text <- turnUpNull(text)
            text <- turnUp(text)
          }
          t<-1:length(text)
          
          listFiltred$intensity <- text
          
          smooth <- ksmooth(t,text,kernel = "normal",bandwidth = 50)
          idpeak <- search_start_end_extr(text,smooth)
          
          table_area <- calculate_area(idpeak,text)
          
          sumArea <- sum(table_area[,"area"])
          table_area[,"area_ratio"] <- (table_area[,"area"]/sumArea)*100
          nameL <- names(listRes)
          
          for(l in 1:5){
            listRes[[nameL[l]]][1] <- table_area[l,"area_ratio"]
            listRes[[nameL[l]]][2] <- table_area[l,"ind_Start"]
            listRes[[nameL[l]]][3] <- table_area[l,"idn_End"]
            listRes[[nameL[l]]][4] <- table_area[l,"ind_Max"]
          }
          
        }
      )

  
  listRes$MPI[1] <- 
    (listRes$X3MPh[1]+listRes$X2MPh[1])/(listRes$X9MPh[1]+listRes$X1MPh[1])
  listRes$MPI1[1] <- 
    1.5*(listRes$X3MPh[1]+listRes$X2MPh[1])/(listRes$X9MPh[1]+listRes$X1MPh[1]+listRes$Phen[1])
    
  listRes$MPI2[1] <- 
    3*listRes$X2MPh[1]/(listRes$X9MPh[1]+listRes$X1MPh[1]+listRes$Phen[1])
  

  
  write_error_file <- file(file_error)  
  
  if(!isempty(vec_error_output)){
    if(file.size(file_error) == 0){
      writeLines(vec_error_output, write_error_file)
    }else{
      write(vec_error_output, file = write_error_file, append = T)
    }
  }
  
  close(write_error_file)
  
  listR <<- listRes
  listF <<- listFiltred
  
}

area_percent_peak(text, name_well, file_error)

# start.time <- Sys.time()
# area_percent_peak(text, name_well, file_error)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# print(time.taken)


