## Configuration
library(compiler)
library(tidyverse)
library(mvtnorm)
library(EML)
source("R/functions.R")
ne = 500 ## production run should be 200 - 5000, depending on what your computer can handle

timestep = 1800 #seconds

#today = Sys.time()

today = "2022-07-01 00:00:00"
jumpBack = 5  ## how many days do we go back to account for data latency?
horiz       = 35 #days, forecast horizon during forecast


## Get latest increment of data (flux & met)
target <- arrow::s3_bucket("targets/terrestrial_30min", endpoint_override="data.ecoforecast.org") %>% 
  dplyr::filter(site_id == "NIWO",time >= as.Date(today)-jumpBack) %>%
  dplyr::collect()

## build dat for Analysis  ******
#nep = -flux$NEE_VUT_REF[today]
#nep.qc = flux$NEE_VUT_REF_QC[today]
nep[nep.qc>0] = NA
dat = data.frame(nep = nep,
                 sigma.nep = flux$NEE_VUT_REF_JOINTUNC[today]
)

##  Get weather forecast
Sys.unsetenv("AWS_DEFAULT_REGION")
Sys.unsetenv("AWS_S3_ENDPOINT")
Sys.setenv(AWS_EC2_METADATA_DISABLED="TRUE")
s3 <- arrow::s3_bucket("drivers/noaa/neon/gefs", 
                       endpoint_override =  "js2.jetstream-cloud.org:8001",
                       anonymous=TRUE)
df <- arrow::open_dataset(s3)
dc <- df %>% 
  filter(site_id == "NIWO",start_time >= as.Date(today)-jumpBack, variable %in% c("TMP","DSWRF")) %>%
  collect()

## for reforecast, stitch together last observed met and first day of each forecast

## for forecast, grab matching start_time
met <- dc %>% filter(start_time == as.Date(today)) %>% pivot_wider(names_from = ensemble,values_from = predicted)
PAR = met %>% filter(variable == "DSWRF") 
inputs <- list(
  date = PAR$time,
  temp = temp = met %>% filter(variable == "TMP") %>% head(-1) %>% select(-(1:9)) %>% t(), ## air temperature (Celsius)
  PAR  = as.matrix(t(PAR[,-(1:9)])) / 0.486 ## gap filled PAR, conversion From Campbell and Norman p151
)

## Get previous Analysis
Analysis = readRDS("Analysis.RDS")
last = length(Analysis)

## Reforecast with updated drivers
out = ensemble_forecast(X = Analysis[[last]]$X,           ## initial conditions = yesterday's Analysis
                        params = Analysis[[last]]$params, ## today's parameters = yesterday's Analysis
                        inputs = inputs.reforecast)       ## today's subset of meteorology 

## Assimilate most recent observations
Analysis[[last+1]] = ParticleFilter(out,Analysis[[t]]$params, dat,wt = Analysis[[t]]$wt)

## Save analysis
saveRDS(Analysis,file="Analysis.RDS")

## Ensemble forecast into the future
out.f = ensemble_forecast(X = Analysis[[last+1]]$X,           ## initial conditions = yesterday's Analysis
                        params = Analysis[[last+1]]$params, ## today's parameters = yesterday's Analysis
                        inputs = inputs.forecast)           ## today's subset of meteorology 


## Convert Forecast to standard & save
dimnames(out.f) <- list(as.character(date[today]),as.character(1:ne),varnames) ## label array dimensions
fx = as.data.frame.table(out.f)                                                ## reorganize into long foremat
colnames(fx) = c("time","ensemble","variable","prediction")                    ## label columns
fx_file = file.path(outdir,
                    paste0("terrestrial_30min_",as.character(as.Date(today)),"_SSEM.csv")) ## output filename
write_csv(fx,fx_file)

## Submit
Sys.setenv("AWS_DEFAULT_REGION" = "data",
           "AWS_S3_ENDPOINT" = "ecoforecast.org")

#aws.s3::put_object(object = fx_file, bucket = "submissions")
#aws.s3::put_object(object = meta_file, bucket = "submissions")



