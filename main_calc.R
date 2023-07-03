library(jsonlite)


go_sсript_calc <- function(val){
  
  switch(val,
         "71" = source("71.R"),
         "178+192" = source("178_192.R"),
         "184+198" = source("184_198.R"),
         "191" = source("191.R"),
         "gas" = source("gas.R"))
}


file_input <- "C:\\Users\\dlasp\\Desktop\\detection oil peak\\test\\89_input.json"

pathT <- strsplit(file_input,"\\\\")[[1]]
fileJ <- pathT[length(pathT)]

pathT <- paste(pathT[-length(pathT)], collapse = '\\')


id <- strsplit(fileJ, '_')[[1]][1]
# 
# fileOut <- paste(file_input, id, "_output.json", sep = '')

file_error <- paste(pathT, paste(id, "_error.txt", sep = ''), sep = '\\')

listWells <- fromJSON(file_input) #"id1_input.json")
wells <- names(listWells) #length(listWells)
listAns <- list()
listFil <- list()

for(i in 1:length(wells)){
  lp <- listWells[[i]]
  
  name_well <- wells[i]
  print(name_well)

  # 
  # toJSON(listFil, auto_unbox = T)
  
  val_array <- unlist(lp$`m/z`)
  
  for(j in 1:length(val_array)){
    val <- val_array[j]
    print(val)
    # listAns[[name_well]][[length(listAns)+1]]
    
    text <- unlist(lp$intensity)
    
    go_sсript_calc(val)
    
    promRes  <- list("m/z" = val,"elements" = listR)
    listAns[[name_well]] <- append(listAns[[name_well]],list(promRes))
    
    # promFil  <- list("m/z" = val,"elements" = listF)
    listFil[[name_well]] <- append(listFil[[name_well]],list(listF))
  }

}


# exportAnsJSON <- toJSON(listAns, pretty = T, auto_unbox=TRUE)
# exportFilJSON <- toJSON(listFil, pretty = T, auto_unbox=TRUE)

write_json(listAns, paste(pathT, paste(id, "_output.json", sep = ''), sep = '\\'), pretty = T, auto_unbox=TRUE)
write_json(listFil, paste(pathT, paste(id, "_filtered.json", sep = ''), sep = '\\'), pretty = T, auto_unbox=TRUE)







