list.of.packages <- c("arrow",
                      "compiler", 
                      "EML",
                      "mvtnorm",
                      "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
