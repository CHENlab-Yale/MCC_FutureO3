#########################################################################################
### Example R code for analysis in the MCC_FutureO3 paper
### Future ozone-related short-term excess mortality under climate and population 
### change scenarios:  a multi-location study in 407 cities in 20 countries
### Nina Domingo, october 2022
### Climate, Health, and Environment Nexus (CHEN) Lab, Yale School of Public Health
#########################################################################################

### load libraries

library(dplyr)
library(tidyverse)
library(tsibble)
library(stringr)

set.seed(123)

# load baseline mortality rates adjusted to present and future time periods
baseline_mortality <- read.csv('file path/file name.csv') # baseline mortality data

# get list of MCC cities included
mcc <- read.csv('file path/file name.csv') # load MCC data

# initialize model list
ensemble <- c() # ensemble list
realization <- c() # realization list
models <- data.frame(ensemble,realization)

# initialize number of Monte Carlo runs
n_runs = 1000

# create cbind fill function
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

for (e in unique(models$ensemble)){
  realization <- models$realization[models$ensemble==e]
  
  for (n in 1:n_runs){
    # read C-R function for sample n
    cr_n <- read.csv('file name/file path.csv')
    
    # get list of mcc cities available in present and future
    mcc.list<- list.files("file path") # file path with list of city files
    mcc.list <- str_extract(mcc.future, ".+?(?=_)") # clean file names to isolate city names
    
    # create empty dataframe where we can store attributable deaths
    add.df <- data.frame(city=character(),cityname=character(),countryname=character(),time_period=character(),add=numeric(),
                         baseline_mortality=numeric(),o3=numeric(),long=numeric(), lat=numeric(),days_greater_than_70=numeric())
    
    print(paste0("Computing attributable deaths with ",e))
    
    for (city in mcc.list){
      countryname <- unique(mcc$countryname_shp[mcc$city==city])
      
      print(paste0('Computing attributable deaths for ',city,'...'))
      
      # FUTURE PERIOD
      future_data <- data.frame()
      # load o3 concentration data across all realizations
      for (r in realization){
        # load bias-corrected future ozone concentrations
        load('file path/file name.RData')
        ts_futuresfo3 <- as.data.frame(future.downscale.fine$Data)
        # compute o3 concentrations above 70 microgram/m^3
        ts_futuresfo3 <- ts_futuresfo3-70
        future_data <- cbind.fill(future_data, ts_futuresfo3)
      }
      
      # compute model's o3 concentration data as the mean of the o3 concentrations across all realizations
      if (length(realization)>1){
        future_data <- as.data.frame(rowMeans(future_data))
      } else {
        future_data <- as.data.frame(future_data)
      }
      
      # select concentration data
      names(future_data) <- 'concentration'
      future_data <- future_data %>% dplyr::select(concentration)
      # select concentration data from 2050 to 2054
      future_data <- head(future_data, 1825)
      future_data$date <- seq(from = as.Date("2050-01-01"), to = as.Date("2054-12-30"), by = 'day')
      # remove variable for looping purposes
      remove(future.downscale.fine)
      
      # compute attributable fraction
      rr.future <- future_data %>%
        mutate(rr.base = cr_n$rr[cr_n$country==countryname],
               beta = cr_n$beta[cr_n$country==countryname],
               rr = exp(beta*concentration),
               af=ifelse(concentration < 0,0,((rr-1)/rr)), #only compute attributable fraction for days where o3 concentrations are at least 70 micrograms/m3
               month = as.numeric(format(as.Date(date,format="%Y-%m-%d"), format = "%m")),
               day = as.numeric(format(as.Date(date,format="%Y-%m-%d"), format = "%d")),
               year = as.numeric(format(as.Date(date,format="%Y-%m-%d"), format = "%y")))
      
      # compute annual attributable deaths
      baseline_mortality_c <- baseline_mortality[baseline_mortality$city==city & baseline_mortality$SSP=="SSP3",]
      add.future <- rr.future %>%
        left_join(., baseline_mortality_c %>%
                    dplyr::select(city, cityname, country, countryname, month, day, baseline_mortality_future)) %>%
        mutate(add=baseline_mortality_future*af,
               add=ifelse(add<0,0,add)) %>%
        dplyr:: group_by(city, cityname, country, countryname, month, day) %>%
        summarise(baseline_mortality=mean(baseline_mortality_future, na.rm=TRUE),
                  add=mean(add, na.rm=TRUE),
                  af=mean(af, na.rm=TRUE),
                  o3=mean(concentration,na.rm=TRUE)) %>%
        dplyr:: group_by(city, cityname, country, countryname) %>%
        summarise(baseline_mortality=sum(baseline_mortality, na.rm=TRUE),
                  add=sum(add, na.rm=TRUE),
                  days_greater_than_70 = length(o3[o3>0]),
                  o3=mean(o3,na.rm=TRUE)) %>%
        mutate(time_period = "2050 to 2054") %>%
        left_join(., mcc %>%
                    dplyr::select(city, long, lat) %>%
                    distinct())

      # add attributable deaths to dataframe
      add.df <- add.df %>%
        rbind(., as.data.frame(add.future))
      
      # PRESENT PERIOD
      # load bias-corrected present ozone concentrations
      present_data <- data.frame()
      for (r in realization){
        # load bias-corrected present ozone concentrations
        load('file path/file name.RData')
        ts_presentsfo3 <- as.data.frame(future.downscale.fine$Data)
        # compute o3 concentrations above 70 microgram/m^3
        ts_presentsfo3 <- ts_presentsfo3-70
        present_data <- cbind.fill(present_data, ts_presentfo3)
      }
      
      # select concentration data
      present_data$concentration <- rowMeans(present_data,na.rm=TRUE)
      present_data <- present_data %>% dplyr::select(concentration)
      # select concentration data from 2010 to 2014
      present_data <- head(present_data, 1825)
      present_data$date <- seq(from = as.Date("2010-01-01"), to = as.Date("2014-12-30"), by = 'day')
      # remove variable for looping purposes
      remove(future.downscale.fine)
        
      # compute attributable fraction
      rr.present <- present_data %>%
        mutate(rr.base = cr_n$rr.sample[cr_n$country==countryname],
               beta = cr_n$beta[cr_n$country==countryname],
               rr = exp(beta*concentration),
               af=ifelse(concentration < 0,0,((rr-1)/rr)), #only compute attributable fraction for days where o3 concentrations are at least 70 micrograms/m3
               month = as.numeric(format(as.Date(date,format="%Y-%m-%d"), format = "%m")),
               day = as.numeric(format(as.Date(date,format="%Y-%m-%d"), format = "%d")),
               year = as.numeric(format(as.Date(date,format="%Y-%m-%d"), format = "%y")))

      # compute annual attributable deaths
      baseline_mortality_c <- baseline_mortality[baseline_mortality$city==city & baseline_mortality$SSP=="SSP3",]
      add.present <- rr.present %>%
        left_join(., baseline_mortality_c %>%
                    dplyr::select(city, cityname, country, countryname, month, day, baseline_mortality_present)) %>%
        mutate(add=baseline_mortality_present*af,
               add=ifelse(add<0,0,add)) %>%
        dplyr:: group_by(city, cityname, country, countryname, month, day) %>%
        summarise(baseline_mortality=mean(baseline_mortality_present, na.rm=TRUE),
                  af=mean(af,na.rm=TRUE),
                  add=mean(add, na.rm=TRUE),
                  o3=mean(concentration,na.rm=TRUE)) %>%
        dplyr:: group_by(city, cityname, country, countryname) %>%
        summarise(baseline_mortality=sum(baseline_mortality, na.rm=TRUE),
                  add=sum(add, na.rm=TRUE),
                  days_greater_than_70 = length(o3[o3>0]),
                  o3=mean(o3,na.rm=TRUE)) %>%
        mutate(time_period = "2010 to 2014") %>%
        left_join(., mcc %>%
                    dplyr::select(city, long, lat) %>%
                    distinct())

      # add attributable deaths to main dataframe
      add.df <- add.df %>%
        rbind(., as.data.frame(add.present))
    }
    # write results to file
    write.csv(add.df,'file path/file name.csv')
  }
}