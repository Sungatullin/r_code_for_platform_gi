options(encoding = "UTF-8")

args <- commandArgs(trailingOnly = T)
file_input <- args[2]

if(args[1] == 'calc'){
  source("main_calc.R")
}else{
  source("main_recalc.R")
}