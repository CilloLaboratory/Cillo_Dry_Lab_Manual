## Write an R script to fix code fencing on the website 
## Oct 30 2023 

## Read in arguments - input .md file to modify 
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("One input file must be specified", call.=FALSE)
}

library(tidyverse)

# Read in data 
dat <- trimws(readLines(args[1]))
dat_store <- dat

# Remove empty lines 
pos_lines <- !dat==""
dat <- dat[!dat==""]

# While loop to add code fencing
new_output <- vector()
i <- 1

while(i<length(dat)+1) {
  
  line_marker <- grepl("^##",dat[i])
  
  if(line_marker==T) {
    ind <- i+1
    while(i<length(dat)) {
      if(grepl("^##",dat[ind])) {
        ind <- ind+1
      } else {
        break
      } 
    }
    if (ind==(i+1)) {
      new_output[i] <- trimws(paste("```tpl\n",dat[[i]],"\n```",sep=""))
      i <- i + 1
      } else {
      new_output[i] <- trimws(paste("```tpl\n",dat[[i]],sep=""))
      new_output[ind-1] <- trimws(paste(dat[[ind-1]],"\n```",sep=""))
      i <- ind
      }
    } else {
    new_output[i] <- dat[i]
    i <- i + 1
  }
  
}

new_output[is.na(new_output)] <- dat[is.na(new_output)]

dat_final <- dat_store[pos_lines] <- new_output

writeLines(dat_final)
