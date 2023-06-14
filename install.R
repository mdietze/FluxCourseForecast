list.of.packages <- c("arrow",
                      "aws.s3",
                      "compiler", 
                      "EML",
                      "mvtnorm",
                      "remotes",
                      "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!("neon4cast" %in% installed.packages()[,"Package"])){
  remotes::install_github("eco4cast/neon4cast",upgrade="always")
}
