setwd("~/R_project/Extreme heat and cold")

library(tidyverse)
library(readr)
library(dlnm)
library(mgcv)
library(gnm)
library(splines)

first_stage_LAD_age_sepcific = readRDS("data/processed/model_data_lsoa11.rds")

#########################################################
#Set up data (Day of week + time trend)
df  = first_stage_LAD_age_sepcific %>%
  mutate(
    #systematic weekly patterns in death registrations, "weekend effect"
    dow = factor(as.character(lubridate::wday(date, label=TRUE, week_start=1)),
                 levels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),
    ym = format(date, "%Y-%m"),
    year = year(date),
    qtr = quarter(date),
    #this is the part to be conditioned out of the regression
    strata = factor(paste0(lsoa11,"-",ym)),
    #made age group as factor 
    # numeric calendar time for long-term trend
    doy = yday(date)
    
  ) %>% 
  group_by(date, dow,ym,year,doy) %>% 
  summarise(tasmean = mean(tasmean),
            deaths = sum(deaths))
  
########################################################
#define the spline for day of the year as long term seasonal trend

#calucalte the degree of freedom
df_per_year = 7
n_years = as.numeric(difftime(max(df$date), min(df$date), units = "days")) / 365.25
df_time = ceiling(df_per_year * n_years)

#the spline line
# spldoy = onebasis(df$doy, "ns", df = df_per_year)

########################################################
########################################################
#extremely omprtant 
#remove strate with 0 outcomes 
#if not removed, bias dispersion estimation, cause instability

argvar = list(fun = "bs", knots = quantile(df$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE))

# knots = quantile(df$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
# lag knots 
#Using logknots(lag_max, 3) is the gold standard for temperature.
#The log-scale puts more "statistical attention" (knots) on the early days (0, 1, 2) and spreads them out as you get toward day 21.
#puts more knots at the beginning (where the mortality effect changes quickly, like during a heatwave) and 
#fewer knots at the end (where the effect of a cold spell stays more constant)
arglag = list(fun = "ns", knots = logknots(21,3))

########################################################


cb = crossbasis(
  df$tasmean,
  lag= c(0,21), 
  argvar = argvar,
  arglag = arglag
)

fit = glm(deaths ~ cb +
            factor(dow) +
            ns(doy, df = 10):factor(year),
          family = quasipoisson,
          data = df)



cp = crosspred(cb,fit,at = -43:282/10,cen=10)

data.frame(
  temp = cp$predvar,
  rr   = cp$allRRfit
) %>% 
  slice_min(rr,n=1)


plot(cp, "overall", ylim = c(0.2, 10))




ref_temp <- 19.8

cp_cum <- crosspred(
  cb,
  fit,
  at    = seq(-4.3, 28.2, by = 0.1),
  cen   = ref_temp,
  cumul = TRUE
)

sim = plot(
  cp_cum,
  "overall",
  xlim = quantile(df$tasmean, c(0.001, 0.999)),
  ylim = c(0.5, 2)
)

sim = findmin(cb,fit,nsim = 5000)


mean(findmin(cb,fit,nsim = 5000, sim=TRUE))



median(findmin(cb,fit,nsim = 5000, sim=TRUE))


quantile(sim, c(0.1,0.5,0.975))



mean(quan)














