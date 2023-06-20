## Configuration
library(tidyverse)
source("R/functions.R")
n_submit = 50         ## sample size of ensemble members submitted to EFI
timestep = 3600       ## seconds, driven by timestep of met data
site_id  = "NIWO"     ## NEON site code
outdir   = "forecast" ## where should forecasts be saved locally before submitting
dir.create(outdir,showWarnings = FALSE)
dir.create("analysis",showWarnings = FALSE) ## should already exist
## ONE TIME: to initiate workflow, copy calibration analysis to day before start date
## e.g. cp analysis/2015-07-01.RDS analysis/2022-07-01.RDS

######## Get previous Analysis output ########
Analysis.files = dir("analysis")
last = length(Analysis.files)
last.date = substr(Analysis.files[last],1,10)
Analysis = readRDS(file.path("analysis",Analysis.files[last]))
ne = nrow(Analysis$X)      ## keep existing ensemble size
if(is.null(Analysis$met)){ ## met ensemble not resampled
  Analysis$met = sample(1:30,ne,replace=TRUE)
}

######### SET UP DATES ##########
today = Sys.time()
today_timestamp = strptime(today, "%Y-%m-%d %H:%M:%S",tz="UTC")
today_ch        = as.character(as.Date(today))
jumpBack = min(100,max(10,as.Date(today) - as.Date(last.date)-1))  ## how many days do we want go back to account for data latency?
## NOTE: max of 100 is arbitrary, system should happily jumpBack any amount of time
## TODO: improve min jumpBack (currently 10 days) by explicitly keeping track of last data assimilated
start_date = lubridate::as_datetime(today)-lubridate::days(jumpBack)
horiz       = 35 #days, forecast horizon during forecast

######## Get latest increment of data (flux) ########
target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_30min/terrestrial_30min-targets.csv.gz", guess_max = 1e6) |>
  dplyr::filter(site_id == site_id)

## build dat for Analysis
nee = target |> filter(variable == "nee")
fu.params = flux.uncertainty(nee$observation) ## use full dataset to calibrate uncertainty rel'n
nee.reforecast = nee |> 
  dplyr::filter(between(datetime,lubridate::as_datetime(as.Date(today)-jumpBack),lubridate::as_datetime(today))) |>
  dplyr::mutate(hour = lubridate::floor_date(datetime,"hour")) |>  ## aggregate fluxes to hourly
  dplyr::group_by(hour) |>
  dplyr::summarise(obs = mean(observation,na.rm = TRUE))
nee.reforecast$obs[is.nan(nee.reforecast$obs)] = NA
unc = predict.flux.uncertainty(nee.reforecast$obs,fu.params)
dat = data.frame(nep = -nee.reforecast$obs,
                 sigma.nep = unc
)

########  Get weather forecast ######
dc = neon4cast::noaa_stage2() |>
  dplyr::filter(site_id == site_id,
                variable %in% c("air_temperature","surface_downwelling_shortwave_flux_in_air"),
                start_date == today_ch) |>
  dplyr::collect()
met <- dc  |> pivot_wider(names_from = parameter,values_from = prediction)
PAR = met |> 
  filter(variable == "surface_downwelling_shortwave_flux_in_air") |> 
  select(-(1:11)) |> as.matrix()
temp = met |> filter(variable == "air_temperature")
fx.time = temp$datetime
temp = temp |> select(-(1:11)) |> as.matrix()
inputs = array(dim=c(nrow(PAR),ne,2))
dimnames(inputs)[[3]] = c("temp","PAR")
inputs[,,"temp"] = temp[,Analysis$met] - 273.15  ## air temperature (Celsius)
inputs[,,"PAR"]  = PAR[,Analysis$met] / 0.486 ## gap filled PAR, conversion From Campbell and Norman p151
## gap-fill
for(t in 1:nrow(inputs)){
  inputs[t,is.na(inputs[t, ,"temp"]),"temp"] = mean(inputs[t, ,"temp"],na.rm = TRUE)
  inputs[t,is.na(inputs[t, ,"PAR" ]),"PAR"]  = mean(inputs[t, ,"PAR"],na.rm = TRUE)
}

######## reforecast met: stitch together first day of each weather forecast ######
dr = neon4cast::noaa_stage3() |>
  dplyr::filter(site_id == site_id,
                  variable %in% c("air_temperature","surface_downwelling_shortwave_flux_in_air")) |>
  dplyr::collect()
dr = dr |> 
  dplyr::filter(between(datetime,lubridate::as_datetime(as.Date(today)-jumpBack),lubridate::as_datetime(today))) |> ## couldn't get between to work in initial query
  na.omit() |>
  pivot_wider(names_from = parameter,values_from = prediction)
PAR = dr |> 
  filter(variable == "surface_downwelling_shortwave_flux_in_air") |> 
  select(-(1:7)) |> as.matrix()
PAR[PAR < 0] = 0
temp = dr |> 
  filter(variable == "air_temperature") |> 
  select(-(1:7)) |> as.matrix()
inputs.reforecast = array(dim=c(nrow(PAR),ne,2))
dimnames(inputs.reforecast)[[3]] = c("temp","PAR")
inputs.reforecast[,,"temp"] = temp[, Analysis$met] - 273.15  ## air temperature (Celsius)
inputs.reforecast[,,"PAR"]  = PAR[, Analysis$met] / 0.486 ## gap filled PAR, conversion From Campbell and Norman p151

#########  REFORECAST  ############
date = nee.reforecast$hour
forecast <- array(NA,c(86400/timestep*(jumpBack),ne,12)) ## output storage [time, ensemble, state]
for(t in 1:jumpBack){
  # counter to help us know things are still running
  print(start_date + lubridate::days(t-1))
  
  # select input rows based on date
  now = which(date >= (start_date + lubridate::days(t-1)) & date < (start_date + lubridate::days(t)))
  
  # forecast
  out = ensemble_forecast(X = Analysis$X,           ## initial conditions = yesterday's Analysis
                          params = Analysis$params, ## today's parameters = yesterday's Analysis
                          inputs = inputs.reforecast[now, , ])       ## today's subset of meteorology     
  forecast[now,,] = out                           ## store today's forecast in overall array
  
  # Today's analysis
  newAnalysis = ParticleFilter(out,
                               params=Analysis$params, 
                               dat = dat[now,],
                               wt = Analysis$wt)
  newAnalysis$met = Analysis$met
  
  # Save the Analysis
  saveRDS(newAnalysis,file=file.path("analysis",paste0(start_date+lubridate::days(t),".RDS")))
  Analysis = newAnalysis
  
}

######## Ensemble forecast into the future ########
out.f = ensemble_forecast(X = Analysis$X,           ## initial conditions = yesterday's Analysis
                        params = Analysis$params,   ## today's parameters = yesterday's Analysis
                        inputs = inputs)            ## today's subset of meteorology 


######## Convert Forecast to standard & save ########
dimnames(out.f) <- list(as.character(fx.time),as.character(1:ne),varnames)      ## label array dimensions
fx = as.data.frame.table(out.f)                                                 ## reorganize into long format
colnames(fx) = c("datetime","parameter","variable","prediction")                ## label columns
fx2 = fx |> 
  dplyr::filter(variable == "NEP") |>                                           ## scored variables
  dplyr::mutate(dplyr::across(variable, ~ str_replace_all(.x,'NEP', 'nee'))) |> ## rename NEP -> NEE
  dplyr::mutate(prediction = -prediction) |>                                    ## change sign on NEE 
  dplyr::mutate(reference_datetime = today_timestamp) |> relocate(reference_datetime) |> ## add reference_datetime
  dplyr::mutate(site_id = site_id) |> relocate(site_id,.after=datetime) |>               ## add site_id
  dplyr::mutate(family = "ensemble") |> relocate(family,.before=parameter) |>            ## add family
  dplyr::mutate(datetime = lubridate::as_datetime(as.character(datetime))) |>            ## make sure datetime is formatted correctly
  dplyr::mutate(parameter = as.numeric(parameter))                                       ## make sure parameter is formatted correctly

## resample to make processing easier and eliminate the need to track weights
param = sample(seq_along(Analysis$wt),n_submit,replace = TRUE,prob = Analysis$wt)
fx3 = fx2 |> filter(parameter %in% param)

## save to file
setwd(outdir)
fx_file = paste0("terrestrial_30min-",today_ch,"-SSEM.csv") ## output filename
write_csv(fx2,fx_file)

######## Submit ########
neon4cast::submit(forecast_file = fx_file, metadata = NULL, ask = FALSE)
