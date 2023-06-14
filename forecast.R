## Configuration
library(compiler)
library(tidyverse)
#library(mvtnorm)
#library(EML)
source("R/functions.R")
ne = 500 ## production run should be 200 - 5000, depending on what your computer can handle
timestep = 1800 #seconds

## Get previous Analysis
Analysis.files = dir("analysis")
last = length(Analysis.files)
last.date = substr(Analysis.files[last],1,10)
Analysis = readRDS(file.path("analysis",Analysis.files[last]))

#today = Sys.time()
today = "2022-07-01 00:00:00"
today_timestamp = strptime(today, "%Y-%m-%d %H:%M:%S",tz="UTC")
today_ch        = as.character(as.Date(today))
jumpBack = min(60,max(5,as.Date(today) - as.Date(last.date)))  ## how many days do we go back to account for data latency?
start_date = lubridate::as_datetime(today)-lubridate::days(jumpBack)
horiz       = 35 #days, forecast horizon during forecast

## Get latest increment of data (flux)
target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_30min/terrestrial_30min-targets.csv.gz", guess_max = 1e6) %>%
  dplyr::filter(site_id == "NIWO",datetime >= start_date)

## build dat for Analysis
nee = target |> filter(variable == "nee")
fu.params = flux.uncertainty(nee$observation)
unc = predict.flux.uncertainty(nee$observation,fu.params)
dat = data.frame(nep = -nee$observation,
                 sigma.nep = unc
)

##  Get weather forecast
dc = neon4cast::noaa_stage2() |>
  dplyr::filter(site_id == "NIWO",
                variable %in% c("air_temperature","surface_downwelling_shortwave_flux_in_air"),
                start_date == today_ch) |>
  dplyr::collect()
met <- dc  %>% pivot_wider(names_from = parameter,values_from = prediction)
PAR = met %>% filter(variable == "surface_downwelling_shortwave_flux_in_air") 
inputs <- list(
  date = PAR$datetime,
  temp = met %>% filter(variable == "air_temperature") %>% select(-(1:11)) %>% t()-273.15, ## air temperature (Celsius) ## prev included %>% head(-1)
  PAR  = as.matrix(t(PAR[,-(1:11)])) / 0.486 ## gap filled PAR, conversion From Campbell and Norman p151
)

## reforecast met: stitch together first day of each weather forecast
dr = neon4cast::noaa_stage3() |>
  dplyr::filter(site_id == "NIWO",
                  variable %in% c("air_temperature","surface_downwelling_shortwave_flux_in_air")) |>
  dplyr::collect()
dr = dr |> 
  dplyr::filter(between(datetime,lubridate::as_datetime(as.Date(today)-jumpBack),lubridate::as_datetime(today)))  ## couldn't get between to work in initial query
PAR = dr |>
  filter(variable == "surface_downwelling_shortwave_flux_in_air") |>
  na.omit() |>
  pivot_wider(names_from = parameter,values_from = prediction) 
inputs.reforecast <- list(
  date = PAR$datetime,
  temp = dr |> filter(variable == "air_temperature") |>
    na.omit() |>
    pivot_wider(names_from = parameter,values_from = prediction) |>
    select(-(1:11)) %>% t()-273.15, ## air temperature (Celsius) ## prev included %>% head(-1)
  PAR  = as.matrix(t(PAR[,-(1:11)])) / 0.486 ## gap filled PAR, conversion From Campbell and Norman p151
)

## Reforecast with updated drivers
out = ensemble_forecast(X = Analysis$X,           ## initial conditions = yesterday's Analysis
                        params = Analysis$params, ## today's parameters = yesterday's Analysis
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
fx = as.data.frame.table(out.f)                                                ## reorganize into long format
colnames(fx) = c("time","ensemble","variable","prediction")                    ## label columns
fx_file = file.path(outdir,
                    paste0("terrestrial_30min_",as.character(as.Date(today)),"_SSEM.csv")) ## output filename
write_csv(fx,fx_file)

## Submit
Sys.setenv("AWS_DEFAULT_REGION" = "data",
           "AWS_S3_ENDPOINT" = "ecoforecast.org")

aws.s3::put_object(object = fx_file, bucket = "submissions")
aws.s3::put_object(object = meta_file, bucket = "submissions")



