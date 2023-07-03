Sys.setenv(LANG = "en")

suppressPackageStartupMessages({
library(stringr)
library(pracma)
library(dplyr)
library(magrittr)
})

fmax <- function(inds,valueTarget,distMatrix,diffPenalty,lam1,lam2){
  return(
    sum(valueTarget[inds]) - 
      lam2*sum(distMatrix[inds,inds]) - 
      lam1*sum(diffPenalty[inds])
  )
}

OptimizeOnF <- function(i0,N,M,valueTarget,distMatrix,diffPenalty,lam1,lam2){
  eps <- 1
  while(eps > 0)
  {
    i1 <- i0
    for (m in 1:M)
    {
      i2 <- i1
      f1 <- sapply(1:N, function(i) {
        i2[m] <- i
        fmax(i2,valueTarget,distMatrix,diffPenalty,lam1,lam2)
      })
      i1[m] <- which.max(f1)
    }
    eps <- sum(abs(i1 - i0))
    i0 <- i1
  }
  i0
}

mu_index <- function(peak, datP, M = 27, deltaD = 1700, a1 = 600, 
                     lam1 = 1e-0, lam2 = 1e+5){
  #параметры штрафных функций
  #deltaD - Расстояние поиска максимального значения кривой слева и справа
  #a1 - Расстояние квадратичного штрафа на сближение
  #lam1 - Штраф на слишком большую депрессию по сравнению с соседними пиками
  #lam2 - Штраф на близость пиков
  tt <- 1:length(peak)
  delta <- 0
  indsLocMax <- diff(tail(peak, -1)) < delta & diff(head(peak, -1)) > delta
  indsLocMax %<>% which
  indsTarget <- tt[1+indsLocMax]
  valueTarget <- datP[1+indsLocMax] %>% sqrt
  #Количество точек, среди которых ищем
  N <- length(indsTarget)

  #Вычисляем штрафы.
  #Заранее вычисляем для каждого пика его депрессию.
  lrMaxes <- sapply(1:N, function(i) {
    c(
      max(valueTarget[indsTarget <= indsTarget[i] & indsTarget >= indsTarget[i] - deltaD]),
      max(valueTarget[indsTarget >= indsTarget[i] & indsTarget <= indsTarget[i] + deltaD])
    )
  }) %>% t
  diffPenalty <- apply(lrMaxes - valueTarget > 0, 1, all) * rowSums(lrMaxes - valueTarget)

  #Штраф за расстояние между пиками.
  distMatrix <- as.matrix(dist(indsTarget)^2)
  distMatrix <- 1 - distMatrix/(a1^2)
  distMatrix[distMatrix < 0] <- 0
  distMatrix <- lower.tri(distMatrix, TRUE)*distMatrix
  
  #Вычисляем оптимум
  i0 <- round((1:M)*(N/M))
  i0 <- OptimizeOnF(i0,N,M,valueTarget,distMatrix,diffPenalty,lam1,lam2)
  
  #Подтянуть до истинных локальных максимумов
  result <- c()
  L <- 30
  for (i in 1:M){
    result[i] <- indsTarget[i0[i]]-L-1 + which.max(datP[indsTarget[i0[i]] + seq(-L,L)])
  }
  result %<>% sort
  
  # browser() 
  #Найдём PrPh
  #`%p%` <- paste0
  L1 <- 150
  L2 <- 2*L1+1
  Ls <- seq(-L1,L1)
  allVals <- c()
  allPos <- c()
  for (i in 1:M){
    idResult <- result[i]+Ls
    if(!isempty(which(idResult<=0))){
      idResult[which(idResult<=0)] <- 1
      dat0 <- datP[idResult]
    }else{
      dat0 <- datP[idResult]
    } 
    i1 <- L1+1
    targetVal <- sapply(1:L2, function(x) (dat0[x] - min(dat0[i1:x]))*(dat0[x]/dat0[i1] > 0.1)*((dat0[x] - min(dat0[i1:x]))/dat0[x] > 0.1))
    i2 <- targetVal %>% which.max
    allVals[i] <- targetVal[i2]/dat0[i1]
    allPos[i] <- i2
  }
  
  is2 <- (order(allVals[4:10], decreasing = TRUE)[1:2] + 3) %>% sort
  PrPh <- result[is2] + Ls[allPos[is2]]
  

  resu <- rep(0,M)
  if((is2[1]==6)){
    resu[2:27] <- result[-length(result)]
  }else{
    if(is2[1]==5){
      resu[3:27] <- result[-c(length(result)-1,length(result))]
    }else{
      resu <- result
    }
  }

  idPoint <- c(resu,PrPh)
  
  
  return(idPoint)
}

strange_spline <- function(data,index){
  l <- length(index)
  difs <- index[2:l]-index[1:(l-1)]
  if(all(difs == 1)){
    ends <- which(difs==1)
    q <- c(1,ends)
    plus <- 10
    SplineFun <- splinefun(y = c(data[(index[q[1]]-plus):index[q[1]]],
                                 data[index[q[length(q)]]:(index[q[length(q)]]+plus)]), 
                           x =c((index[q[1]]-plus):index[q[1]],
                                index[q[length(q)]]:(index[q[length(q)]]+plus)),method = "natural")
    SplineFit_2 <- SplineFun(index[(1:q[length(q)])])
    data[index[(1:q[length(q)])]] <- SplineFit_2
    return(data)
  }else{
    ends <- which(difs>1)
    q <- c(1,ends)
    if(length(q) == 2){
      #first
      SplineFun <- splinefun(y = c(data[(index[q[1]]-del):index[q[1]]],
                                   data[index[q[2]]:(index[q[2]]+del)]), 
                             x =c((index[q[1]]-del):index[q[1]],
                                  index[q[2]]:(index[q[2]]+del)),method = "natural")
      SplineFit_2 <- SplineFun(index[(1:q[2])])
      data[index[(1:q[2])]] <- SplineFit_2
      #last
      SplineFun <- splinefun(y = c(data[(index[q[length(q)]+1]-del):index[q[length(q)]+1]],
                                   data[index[length(index)]:(index[length(index)]+del)]), 
                             x =c((index[q[length(q)]+1]-del):index[q[length(q)]+1],
                                  index[length(index)]:(index[length(index)]+del)),method = "natural")
      SplineFit_5 <- SplineFun(index[(q[length(q)]+1):(length(index))])
      
      data[index[(q[length(q)]+1):(length(index))]] <- SplineFit_5
      return(data)
    }else{
      #first
      SplineFun <- splinefun(y = c(data[(index[q[1]]-del):index[q[1]]],
                                   data[index[q[2]]:(index[q[2]]+del)]), 
                             x =c((index[q[1]]-del):index[q[1]],
                                  index[q[2]]:(index[q[2]]+del)),method = "natural")
      SplineFit_2 <- SplineFun(index[(1:q[2])])
      data[index[(1:q[2])]] <- SplineFit_2
      
      for(i in 2:(length(q)-1)){
        SplineFun <- splinefun(y = c(data[(index[q[i]+1]-del):index[q[i]+1]],
                                     data[index[q[i+1]]:(index[q[i+1]]+del)]), 
                               x = c((index[q[i]+1]-del):index[q[i]+1],
                                     index[q[i+1]]:(index[q[i+1]]+del)),method = "natural")
        SplineFit_3 <- SplineFun(index[(q[i]+1):q[i+1]])
        data[index[(q[i]+1):q[i+1]]] <- SplineFit_3
      }
      
      #last
      SplineFun <- splinefun(y = c(data[(index[q[length(q)]+1]-del):index[q[length(q)]+1]],
                                   data[index[length(index)]:(index[length(index)]+del)]), 
                             x =c((index[q[length(q)]+1]-del):index[q[length(q)]+1],
                                  index[length(index)]:(index[length(index)]+del)),method = "natural")
      SplineFit_5 <- SplineFun(index[(q[length(q)]+1):(length(index))])
      
      data[index[(q[length(q)]+1):(length(index))]] <- SplineFit_5
      return(data)
    }
  }
}

search_start_end_extr <- function(dat){
   #for smooth data
  smooth <- ksmooth(1:length(dat),dat,kernel = "normal",bandwidth = 50)
  peak <- smooth$y
  eps <- 0.001
  N <- length(peak)
  smooth2 <- ksmooth(1:length(dat),dat,kernel = "normal",bandwidth = 15)
  peak2 <- smooth2$y
  NN <- length(peak2)
  log_vector <- peak2[2:(N-1)]-peak2[1:(N-2)] < eps & 
    peak2[2:(N-1)]-peak2[3:N] < eps
  
  min_of_peak  <- which(log_vector)
  # browser()
  ind_max_data <- mu_index(peak,dat)
  
  ind_peak <- data.frame("Value_max" = rep(0,length(ind_max_data)),
                         "ind_Max" = ind_max_data, 
                         "ind_Start" = rep(0,length(ind_max_data)), 
                         "idn_End" = rep(0,length(ind_max_data)))
  
  for(i in 1:length(ind_max_data)){
    if(ind_max_data[i] == 0){
      next
    }
      ind_peak[i,1] <- dat[ind_max_data[i]]
      max_peak <- sort(c(min_of_peak,ind_max_data[i]))
      index_start_end <- max_peak[which(max_peak == ind_max_data[i]) + c(-1,1)]
      ind_peak[i,3] <- index_start_end[1]
      ind_peak[i,4] <- index_start_end[2]
  }
  
  vector_eps_for_peak <- matrix(0,ncol = nrow(ind_peak),nrow=2)
  vector_eps_for_peak[1,] <- ind_peak[,2]
  for(i in 1:nrow(ind_peak)){
    if(ind_peak[i,1]>30000000){
      vector_eps_for_peak[2,i]<- 5000
    }
    if(ind_peak[i,1]<30000000 & ind_peak[i,1]>4500000){
      vector_eps_for_peak[2,i]<- 9000
    }
    if(ind_peak[i,1]<9000000 & ind_peak[i,1]>4500000){
      vector_eps_for_peak[2,i]<- 6000
    }
  }
  ind_max_data <- vector_eps_for_peak[1,] #+"'khgjh"
  
  for(i in which(vector_eps_for_peak[2,i]!=0)){
    if(ind_max_data[i] == 0){
      next
    }
    log_vector <- peak2[2:(N-1)]-peak2[1:(N-2)] < vector_eps_for_peak[2,i] & 
      peak2[2:(N-1)]-peak2[3:N] < vector_eps_for_peak[2,i]
    min_peak  <- which(log_vector)
    max_peak <- sort(c(min_peak,ind_max_data[i]))
    index_start_end <- max_peak[which(max_peak == ind_max_data[i]) + c(-1,1)]
    ind_peak[i,3] <- index_start_end[1]
    ind_peak[i,4] <- index_start_end[2]
  }
  return(ind_peak)
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

calculate_area <- function(idpeak,text){
  sCoef <- 1.016112e-07
  # browser()
  table_area_neeed <- idpeak
  table_area_neeed[,"area"] <- 0
  for(j in 1:nrow(table_area_neeed)){
    len_peak <- table_area_neeed[j,3]:table_area_neeed[j,4]
    
    if((j == 1 || j == 2) & (table_area_neeed[j,2] == 0)){
      next
    }
    table_area_neeed[j,"area"] <- 
      sCoef*(trapz(len_peak,text[len_peak])-
               (min(text[table_area_neeed[j,4]],text[table_area_neeed[j,3]])*
                  (table_area_neeed[j,4]-table_area_neeed[j,3])))
    if((j == 1 || j == 2) & (table_area_neeed[j,"area"] < 9)){
      table_area_neeed[j,] <- 0
    }
  }
  wLess <- which(table_area_neeed[20:27,"area"]<0.8)
  if(!isempty(wLess)){
    vvv <- c(20:27)
    table_area_neeed[vvv[wLess],] <- 0
  }
  # browser()
  return(table_area_neeed)
}

create_list<- function(len_vec = 4){
  ls <- list("C11"=rep(0,len_vec),
                   "C12"=rep(0,len_vec),
                   "C13"=rep(0,len_vec),
                   "C14"=rep(0,len_vec),
                   "C15"=rep(0,len_vec),
                   "C16"=rep(0,len_vec),
                   "C17"=rep(0,len_vec),
                   "C18"=rep(0,len_vec),
                   "C19"=rep(0,len_vec),
                   "C20"=rep(0,len_vec),
                   "C21"=rep(0,len_vec),
                   "C22"=rep(0,len_vec),
                   "C23"=rep(0,len_vec),
                   "C24"=rep(0,len_vec),
                   "C25"=rep(0,len_vec),
                   "C26"=rep(0,len_vec),
                   "C27"=rep(0,len_vec),
                   "C28"=rep(0,len_vec),
                   "C29"=rep(0,len_vec),
                   "C30"=rep(0,len_vec),
                   "C31"=rep(0,len_vec),
                   "C32"=rep(0,len_vec),
                   "C33"=rep(0,len_vec),
                   "C34"=rep(0,len_vec),
                   "C35"=rep(0,len_vec),
                   "C36"=rep(0,len_vec),
                   "C37"=rep(0,len_vec),
                   "Pr"=rep(0,len_vec),
                   "Ph"=rep(0,len_vec),
                   "Pr/Ph" = 0,
                   "TAR" = 0,
                   "CPI" = 0,
                   "Pr/C17" = 0,
                   "Ph/C18" = 0,
                   "Ki" = 0,
                   "KVN" = 0,
                   "C27/C17" = 0,
                   "2C29/(C28+C30)" = 0)
  return(ls)
}

`%p%` <- paste0

area_percent_peak <- function(text, name_well, file_error){
  
  vec_error_wells <- c()
  vec_error_output <- c()
  
  listFiltred <- list("m/z" = "71",
                      "intensity" = 0)
  listRes <- create_list()
  
    tryCatch(
        error = function(cnd) {
          vec_error_wells <<- c(vec_error_wells,i)

          cl1 <- cnd$call %>% deparse %>% paste0(collapse="") %>% gsub("  ", "", .) %>%
            gsub("\n", " ", .) %>% gsub("\r", " ", .)
          msg <- "Error in " %p% cl1 %p% ": " %p% cnd$message
          print(msg)
          vec_error_output <<- c(vec_error_output, paste("OIL 71mz", name_well, msg, sep = ": "))
        },
        {
          text <- turnUp(text)
          text[text == 0] <- 0.1
          t<-1:length(text) 
          
          listFiltred$intensity <- text

          idpeak <- search_start_end_extr(text)
		  
		  table_area <- calculate_area(idpeak,text)
        
      sumArea <- sum(table_area[,"area"])
      table_area[,"area_ratio"] <- (table_area[,"area"]/sumArea)*100
      nameL <- names(listRes)
      
      for(l in 1:29){
        listRes[[nameL[l]]][1] <- table_area[l,"area_ratio"]
        listRes[[nameL[l]]][2] <- table_area[l,"ind_Start"]
        listRes[[nameL[l]]][3] <- table_area[l,"idn_End"]
        listRes[[nameL[l]]][4] <- table_area[l,"ind_Max"]
      }
          
        }
      )

    listRes$`Pr/Ph` <- listRes$Pr[1]/listRes$Ph[1]
    listRes$TAR <- (listRes$C27[1] + listRes$C29[1] +
                      listRes$C31[1])/(listRes$C15[1] +
                                          listRes$C17[1] + listRes$C19[1])
    listRes$CPI <-
      2*(listRes$C23[1] +
           listRes$C25[1] + listRes$C27[1]+
           listRes$C29[1])/(listRes$C22[1]+
                                2*(listRes$C24[1] + listRes$C26[1] + listRes$C28[1]) + listRes$C30[1])
    listRes$`Pr/C17` <- listRes$Pr[1]/listRes$C17[1]
    listRes$`Ph/C18` <- listRes$Ph[1]/listRes$C18[1]
    listRes$Ki <- (listRes$Pr[1] + listRes$Ph[1])/(listRes$C17[1] + listRes$C18[1])
    listRes$KVN <-
      (listRes$C27[1] + listRes$C28[1] + listRes$C29[1] +
         listRes$C30[1] + listRes$C31[1])/(listRes$C15[1] +
                                              listRes$C16[1] + listRes$C17[1] +
                                              listRes$C18[1] + listRes$C19[1])
    listRes$`C27/C17`<- listRes$C27[1]/listRes$C17[1]
    listRes$`2C29/(C28+C30)` <- 2*listRes$C29[1]/(listRes$C28[1] + listRes$C30[1])

#####    
        
#   }else{
#     
#     dataArea <- read.xlsx(paste(path_to,"dataPercent_with_coefs.xlsx",sep="//"),1)
#     dataArea_Area <- read.xlsx(paste(path_to,"dataArea_with_sum.xlsx",sep="//"),1)
#     dataPeaks <- read.xlsx(paste(path_to,"id_peaks.xlsx",sep="//"),"id_peak")
#     dataStart <- read.xlsx(paste(path_to,"id_peaks.xlsx",sep="//"),"id_start")
#     dataEnd <- read.xlsx(paste(path_to,"id_peaks.xlsx",sep="//"),"id_end")
#     
#     rem <- nrow(dataArea)
# 
#     for(i in 1:length_file_vector){
#       tryCatch(
#         error = function(cnd) {
#           #здесь обработка ошибок
#           vec_error_wells <<- c(vec_error_wells,i)
# 
#           cl1 <- cnd$call %>% deparse %>% paste0(collapse="") %>% gsub("  ", "", .) %>%
#             gsub("\n", " ", .) %>% gsub("\r", " ", .)
#           msg <- "Error in " %p% cl1 %p% ": " %p% cnd$message
#           print(msg)
#           vec_error_output <<- c(vec_error_output,paste("OIL 71mz",file_vector[i], msg, sep = ": "))
#         },
#         {
#           del_well <- c()
#           if((str_remove(strsplit(file_vector[i],"\\\\"), ".txt")) %in% dataArea$name){
#             del_well <- c(del_well,i)
#             next
#           }
#           #здесь код, на котором ловится ошибка
#           text <- read_data(paste(path_from,file_vector[i],sep= "/"))
#           text <- turnUp(text)
#           text[text == 0] <- 0.1
#           t<-1:length(text)
# 
#           wb_id1 <- createWorkbook()
#           addWorksheet(wb_id1,"data")
#           writeData(
#             wb_id1,
#             "data",
#             text,
#             startCol = 1,
#             startRow = 1,
#             colNames = TRUE
#           )
#           saveWorkbook(wb_id1,
#                        file = paste(paste(path_to,"clear data",sep="\\"),"\\",str_remove(strsplit(file_vector[i],"\\\\"), ".txt"),".xlsx",sep=""),
#                        overwrite = TRUE)
# 
#           #поиск вылетов в правой части данных
#           idpeak <- search_start_end_extr(text)
# 
# 		  # if(is.na(idpeak[28,1])){
#             # idpeak[28,] <- 0
#             # vec_error_output <- c(vec_error_output,paste("OIL 71mz",file_vector[i], "Pr is empty. Please, mark them on the picture", sep = ": "))
#             # vec_pr_ph_wells <- c(vec_pr_ph_wells,i)
#           # }
#           # if(is.na(idpeak[29,1])){
#             # idpeak[29,] <- 0
#             # vec_error_output <- c(vec_error_output,paste("OIL 71mz",file_vector[i], "Ph is empty. Please, mark them on the picture", sep = ": "))
#             # vec_pr_ph_wells <- c(vec_pr_ph_wells,i)
#           # }
# 		  
#           # if(is.na(idpeak[28,1])){
#             # idpeak[28,] <- 0
#             # vec_error_output <- c(vec_error_output,paste("OIL 71mz",file_vector[i],"Pr is empty. Please, mark them on the picture",sep = ": "))
#             # vec_pr_ph_wells <- c(vec_pr_ph_wells,i)
#           # }
#           # if(is.na(idpeak[29,1])){
#             # idpeak[29,] <- 0
#             # vec_error_output <- c(vec_error_output,paste("OIL 71mz",file_vector[i],"Ph is empty. Please, mark them on the picture", sep = ": "))
#             # vec_pr_ph_wells <- c(vec_pr_ph_wells,i)
#           # }
#           
#           # bOOOUT <- c()
#           # for(j in 15:27){
#             # if((idpeak[j,"Value_max"]>idpeak[14,"Value_max"])){bOOOUT <- c(bOOOUT,j)}
#           # }
#           # if(!isempty(bOOOUT)){
#             # text <- strange_spline(peak,idpeak[bOOOUT,"ind_Start"]:idpeak[bOOOUT,"idn_End"])
#             # t<-1:length(text)
#             # smooth <- ksmooth(t,text,kernel = "normal",bandwidth = 20)
#             # smooth1 <- ksmooth(t,text,kernel = "normal",bandwidth = 15)#25
#             # peak <- smooth$y
#             # tt <- smooth$x
#             # idpeak <- search_start_end_extr(peak,smooth1)
#             
#             # if(is.na(idpeak[28,1])){
#               # idpeak[28,] <- 0
#             # }
#             # if(is.na(idpeak[29,1])){
#               # idpeak[29,] <- 0
#             # }
#             
#           # }
# 		# if(table_area[1,"ind_Max"] == table_area[2,"ind_Max"]){
# 			# table_area[1,] <- 0
# 			# vec_error_output <- c(vec_error_output,paste("OIL 71mz",file_vector[i],"Please note: C11 is not defined correctly.",sep = ": "))		
# 		  # }
#           
#           dataArea_Area[nrow(dataArea_Area)+1,"name"] <- str_remove(strsplit(file_vector[i],"\\\\"), ".txt")
#           table_area <- calculate_area(idpeak,text)
# 
#           dataArea_Area[nrow(dataArea_Area),2:(ncol(dataArea_Area)-1)] <- unlist(table_area[,"area"])
#           
#    
#           dataArea[nrow(dataArea)+1,"name"] <- str_remove(strsplit(file_vector[i],"\\\\"), ".txt")
#           dataPeaks[nrow(dataPeaks)+1,"name"] <- str_remove(strsplit(file_vector[i],"\\\\"), ".txt")
#           dataEnd[nrow(dataEnd)+1,"name"] <- str_remove(strsplit(file_vector[i],"\\\\"), ".txt")
#           dataStart[nrow(dataStart)+1,"name"] <- str_remove(strsplit(file_vector[i],"\\\\"), ".txt")
#           
#           sumArea <- sum(table_area[,"area"])
#           dataArea[nrow(dataArea),2:30] <- (table_area[,"area"]/sumArea)*100
#           
#           # заполняем таблицы с id
#           dataPeaks[nrow(dataPeaks),2:ncol(dataPeaks)] <- unlist(table_area[,2])
#           dataEnd[nrow(dataEnd),2:ncol(dataEnd)] <- unlist(table_area[,4])
#           dataStart[nrow(dataStart),2:ncol(dataStart)] <- unlist(table_area[,3])
#           
#       }
#       )
#     }
#     
#     if(isempty(del_well)){
#       param <- rem:(rem+length_file_vector) 
#       
#       dataArea[param,"Pr/Ph"] <- dataArea[param,"Pr"]/dataArea[param,"Ph"]
#       dataArea[param,"TAR"] <- (dataArea[param,"C27"]+dataArea[param,"C29"]+
#                                   dataArea[param,"C31"])/(dataArea[param,"C15"]+
#                                                             dataArea[param,"C17"]+dataArea[param,"C19"])
#       dataArea[param,"CPI"] <-
#         2*(dataArea[param,"C23"]+
#              dataArea[param,"C25"]+dataArea[param,"C27"]+
#              dataArea[param,"C29"])/(dataArea[param,"C22"]+
#                                        2*(dataArea[param,"C24"]+dataArea[param,"C26"]+dataArea[param,"C28"])+dataArea[param,"C30"])
#       dataArea[param,"Pr/C17"] <- dataArea[param,"Pr"]/dataArea[param,"C17"]
#       dataArea[param,"Ph/C18"] <- dataArea[param,"Ph"]/dataArea[param,"C18"]
#       dataArea[param,"Ki"] <- (dataArea[param,"Pr"]+dataArea[param,"Ph"])/(dataArea[param,"C17"]+dataArea[param,"C18"])
#       dataArea[param,"KVN"] <-
#         (dataArea[param,"C27"]+dataArea[param,"C28"]+dataArea[param,"C29"]+
#            dataArea[param,"C30"]+dataArea[param,"C31"])/(dataArea[param,"C15"]+
#                                                            dataArea[param,"C16"]+dataArea[param,"C17"]+
#                                                            dataArea[param,"C18"]+dataArea[param,"C19"])
#       dataArea[param,"C27/C17"] <- dataArea[param,"C27"]/dataArea[param,"C17"]
#       dataArea[param,"2C29/(C28+C30)"] <- 2*dataArea[param,"C29"]/(dataArea[param,"C28"]+dataArea[param,"C30"])
#     }else{
#       param <- rem:(rem+length_file_vector-length(del_well)) 
#       
#       dataArea[param,"Pr/Ph"] <- dataArea[param,"Pr"]/dataArea[param,"Ph"]
#       dataArea[param,"TAR"] <- (dataArea[param,"C27"]+dataArea[param,"C29"]+
#                                   dataArea[param,"C31"])/(dataArea[param,"C15"]+
#                                                             dataArea[param,"C17"]+dataArea[param,"C19"])
#       dataArea[param,"CPI"] <-
#         2*(dataArea[param,"C23"]+
#              dataArea[param,"C25"]+dataArea[param,"C27"]+
#              dataArea[param,"C29"])/(dataArea[param,"C22"]+
#                                        2*(dataArea[param,"C24"]+dataArea[param,"C26"]+dataArea[param,"C28"])+dataArea[param,"C30"])
#       dataArea[param,"Pr/C17"] <- dataArea[param,"Pr"]/dataArea[param,"C17"]
#       dataArea[param,"Ph/C18"] <- dataArea[param,"Ph"]/dataArea[param,"C18"]
#       dataArea[param,"Ki"] <- (dataArea[param,"Pr"]+dataArea[param,"Ph"])/(dataArea[param,"C17"]+dataArea[param,"C18"])
#       dataArea[param,"KVN"] <-
#         (dataArea[param,"C27"]+dataArea[param,"C28"]+dataArea[param,"C29"]+
#            dataArea[param,"C30"]+dataArea[param,"C31"])/(dataArea[param,"C15"]+
#                                                            dataArea[param,"C16"]+dataArea[param,"C17"]+
#                                                            dataArea[param,"C18"]+dataArea[param,"C19"])
#       dataArea[param,"C27/C17"] <- dataArea[param,"C27"]/dataArea[param,"C17"]
#       dataArea[param,"2C29/(C28+C30)"] <- 2*dataArea[param,"C29"]/(dataArea[param,"C28"]+dataArea[param,"C30"])
#     }
#   }
#####
    
  
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


