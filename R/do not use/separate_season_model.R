setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Extreme heat and cold")


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
    #this is the part to be conditioned out of the regression
    strata = factor(paste0(lsoa11,"-", ym)),
    #made age group as factor 
    age_group = factor(age_group, levels = c("0-64",
                                             "65-74",
                                             "75-84",
                                             "85+")),
    year = year(date),
    # numeric calendar time for long-term trend
    doy = yday(date)
  )


########################################################
#Add a season indicator
df = df %>%
  mutate(
    season = case_when(
      month(date) %in% c(12, 1, 2, 3) ~ "winter",
      month(date) %in% c(6, 7, 8, 9)  ~ "summer",
      TRUE                            ~ NA_character_
    )
  )

########################################################
#add covid effect
df <- df %>%
  mutate(
    date = as.Date(date),
    covid_period = date >= as.Date("2020-03-01") &
      date <= as.Date("2021-07-31"),
    covid_time = if_else(
      covid_period,
      as.numeric(date - as.Date("2020-03-01")),
      NA_real_
    )
  )

########################################################
#separate our dataset to winter and summer
df_winter = df %>% filter(season == "winter") %>% select(-tasmax)
df_summer = df %>% filter(season == "summer") %>% select(-tasmin)

########################################################
########################################################
########################################################
#winter modelling 

age_levels = levels(df$age_group)

fit_mod_winter = vector("list", length(age_levels))

names(fit_mod_winter) = age_levels


########################################################
#to create the time shifted copies for lags
cbs_winter  = vector("list", length(age_levels))  # store crossbasis used in each fit
names(cbs_winter) =  age_levels

########################################################
#extremely omprtant 
#remove strate with 0 outcomes 
#if not removed, bias dispersion estimation, cause instability

df_winter = df_winter %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()


########################################################
# we used a quadratic B-spline
# was used function with 3 df to model the non-linear
# effects of temperature on mortality with three equidistant knots from the temperature range. We specified a
# natural cubic spline with three equally spaced internal
# knots at logarithm scale for the lag periods, 0–21 days
# because it is flexible to handle endpoints of data where
# some degree of non-linearity is expected

#for winter range
tmin_range = range(df_winter$tasmin, na.rm = TRUE)

#three equidistant knots
tmin_knots = seq(
  from = tmin_range[1],
  to   = tmin_range[2],
  length.out = 5
)[2:4]

########################################################


########################################################


for (age in age_levels) {
  
  data <- df_winter %>%
    filter(age_group == age)
  
  # IMPORTANT: put crossbasis into the data frame
    cb_winter <- crossbasis(
    data$tasmin,
    lag = 21,
    argvar = list(fun = "bs", knots = tmin_knots, degree = 2),
    arglag = list(fun = "ns", knots = logknots(21, 3))
  )
    
    cbs_winter[[age]] <- cb_winter
  
  fit <- gnm(
    deaths ~ cb_winter +
      factor(dow) +
      splines::ns(doy, df = 2 * 11) +
      splines::ns(covid_time, k = 6),
    eliminate = factor(strata),
    family = quasipoisson,
    data = data
  )
  
  fit_mod_winter[[age]] <- fit
}

#############################################################

winter_df_0_64 = df_winter %>% 
  filter(age_group == "0-64") %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()

cb_winter_0_64 <- crossbasis(
  winter_df_0_64$tasmin,
  lag = 21,
  argvar = list(fun = "bs", knots = quantile(winter_df_0_64$tasmin, c(0.1,0.75,0.95), na.rm = TRUE)),
  arglag = list(fun = "ns", knots = logknots(21, 3))
)

fit_0_64  = gnm(
  deaths ~ cb_winter_0_64 +
    factor(dow) +
    splines::ns(doy, df = 2 * 11) +
    factor(covid_period),
  eliminate = factor(strata),
  family = quasipoisson,
  data = winter_df_0_64
)


temp_grid <- seq(
  from = quantile(winter_df_0_64$tasmin, 0.01, na.rm = TRUE),
  to   = quantile(winter_df_0_64$tasmin, 0.99, na.rm = TRUE),
  by   = 0.1
)



class(cb_winter_0_64)

cp_0_64 = crosspred(cb_winter_0_64,fit_0_64,cen=10, at=temp_grid,cumul = TRUE)



# Identify cumulative lag column
j= ncol(cp_0_64$cumfit)

keep = cp_0_64$predvar >= quantile(winter_df_0_64$tasmin, 0.01) &
  cp_0_64$predvar <= quantile(winter_df_0_64$tasmin, 0.99)


mmt = cp_0_64$predvar[keep][
  which.min(cp_0_64$cumRRfit[keep, j])
]

rr_df <- data.frame(
  temp    = cp_0_64$predvar,
  rr      = cp_0_64$cumRRfit[, j],
  rr_low  = cp_0_64$cumRRlow[, j],
  rr_high = cp_0_64$cumRRhigh[, j]
)

ggplot(rr_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high),
              fill = "grey80") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = mmt, linetype = "dashed") +
  labs(
    x = "Winter mean temperature (°C)",
    y = "Cumulative RR (lag 0–21)",
    title = "Winter temperature–mortality association (Age 60-64)"
  ) +
  theme_minimal()




plot(cp_0_64, "overall")


plot(cp_0_64, lag=7)

##########################################################################

winter_df_85 = df_winter %>% 
  filter(age_group == "85+") %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()

cb_winter_85 <- crossbasis(
  winter_df_85$tasmin,
  lag = 21,
  argvar = list(fun = "bs", knots = quantile(winter_df_85$tasmin, c(0.1,0.75,0.95), na.rm = TRUE)),
  arglag = list(fun = "ns", knots = logknots(21, 3))
)

fit_85  = gnm(
  deaths ~ cb_winter_85 +
    factor(dow) +
    splines::ns(doy, df = 2 * 11)+
    factor(covid_period),
  eliminate = factor(strata),
  family = quasipoisson,
  data = winter_df_85
)




class(cb_winter_85)


temp_grid <- seq(
  from = quantile(winter_df_85$tasmin, 0.01, na.rm = TRUE),
  to   = quantile(winter_df_85$tasmin, 0.99, na.rm = TRUE),
  by   = 0.1
)




cp_85 = crosspred(cb_winter_85,fit_85,cen=10, cumul = TRUE, at= temp_grid)

# Identify cumulative lag column
j= ncol(cp_85$cumfit)

keep = cp_85$predvar >= quantile(winter_df_85$tasmin, 0.01) &
  cp_85$predvar <= quantile(winter_df_85$tasmin, 0.99)


mmt = cp_85$predvar[keep][
  which.min(cp_85$cumRRfit[keep, j])
]

rr_df <- data.frame(
  temp    = cp_85$predvar,
  rr      = cp_85$cumRRfit[, j],
  rr_low  = cp_85$cumRRlow[, j],
  rr_high = cp_85$cumRRhigh[, j]
)

ggplot(rr_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high),
              fill = "grey80") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = mmt, linetype = "dashed") +
  labs(
    x = "Winter mean temperature (°C)",
    y = "Cumulative RR (lag 0–21)",
    title = "Winter temperature–mortality association (Age 85+)"
  ) +
  theme_minimal()


ggplot(winter_df_85, aes(y=tasmin))+
  geom_histogram()+
  coord_flip()


##################################################################


winter_df_65_74 = df_winter %>% 
  filter(age_group == "65-74") %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()

cb_winter_65_74 <- crossbasis(
  winter_df_65_74$tasmin,
  lag = 21,
  argvar = list(fun = "ns", knots = quantile(winter_df_65_74$tasmin, c(0.1,0.75,0.95), na.rm = TRUE)),
  arglag = list(fun = "ns", knots = logknots(21, 3))
)

fit_65_74  = gnm(
  deaths ~ cb_winter_65_74 +
    factor(dow) +
    splines::ns(doy, df = 2 * 11)+
    factor(covid_period),
  eliminate = factor(strata),
  family = quasipoisson,
  data = winter_df_65_74
)




class(cb_winter_65_74)


temp_grid <- seq(
  from = quantile(winter_df_65_74$tasmin, 0.01, na.rm = TRUE),
  to   = quantile(winter_df_65_74$tasmin, 0.99, na.rm = TRUE),
  by   = 0.1
)




cp_65_74 = crosspred(cb_winter_65_74,fit_65_74,cen=10, cumul = TRUE, at= temp_grid)

# Identify cumulative lag column
j= ncol(cp_65_74$cumfit)

keep = cp_65_74$predvar >= quantile(winter_df_65_74$tasmin, 0.01) &
  cp_65_74$predvar <= quantile(winter_df_65_74$tasmin, 0.99)


mmt = cp_65_74$predvar[keep][
  which.min(cp_65_74$cumRRfit[keep, j])
]

rr_df <- data.frame(
  temp    = cp_65_74$predvar,
  rr      = cp_65_74$cumRRfit[, j],
  rr_low  = cp_65_74$cumRRlow[, j],
  rr_high = cp_65_74$cumRRhigh[, j]
)

ggplot(rr_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high),
              fill = "grey80") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = mmt, linetype = "dashed") +
  labs(
    x = "Winter mean temperature (°C)",
    y = "Cumulative RR (lag 0–21)",
    title = "Winter temperature–mortality association (Age 65-74)"
  ) +
  theme_minimal()











##################################################################




winter_df_75_84 = df_winter %>% 
  filter(age_group == "75-84") %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()

cb_winter_75_84 <- crossbasis(
  winter_df_75_84$tasmin,
  lag = 21,
  argvar = list(fun = "ns", knots = quantile(winter_df_75_84$tasmin, c(0.1,0.75,0.95), na.rm = TRUE)),
  arglag = list(fun = "ns", knots = logknots(21, 3))
)

fit_75_84  = gnm(
  deaths ~ cb_winter_75_84 +
    factor(dow) +
    splines::ns(doy, df = 2 * 11)+
    factor(covid_period),
  eliminate = factor(strata),
  family = quasipoisson,
  data = winter_df_75_84
)




class(cb_winter_75_84)


temp_grid <- seq(
  from = quantile(winter_df_75_84$tasmin, 0.01, na.rm = TRUE),
  to   = quantile(winter_df_75_84$tasmin, 0.99, na.rm = TRUE),
  by   = 0.1
)




cp_75_84 = crosspred(cb_winter_75_84,fit_75_84,cen=10, cumul = TRUE, at= temp_grid)

# Identify cumulative lag column
j= ncol(cp_75_84$cumfit)

keep = cp_75_84$predvar >= quantile(winter_df_75_84$tasmin, 0.01) &
  cp_75_84$predvar <= quantile(winter_df_75_84$tasmin, 0.99)


mmt = cp_75_84$predvar[keep][
  which.min(cp_75_84$cumRRfit[keep, j])
]

rr_df <- data.frame(
  temp    = cp_75_84$predvar,
  rr      = cp_75_84$cumRRfit[, j],
  rr_low  = cp_75_84$cumRRlow[, j],
  rr_high = cp_75_84$cumRRhigh[, j]
)

ggplot(rr_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high),
              fill = "grey80") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = mmt, linetype = "dashed") +
  labs(
    x = "Winter mean temperature (°C)",
    y = "Cumulative RR (lag 0–21)",
    title = "Winter temperature–mortality association (Age 75-84)"
  ) +
  theme_minimal()









##################################################################
##################################################################
##################################################################


summer_df_65_74 = df_summer %>% 
  filter(age_group == "65-74") %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()

cb_summer_65_74 <- crossbasis(
  summer_df_65_74$tasmax,
  lag = 3,
  argvar = list(fun = "ns", knots = quantile(summer_df_65_74$tasmax, c(0.5,0.9), na.rm = TRUE)),
  arglag = list(fun = "ns", knots = 1)
)

fit_65_74  = gnm(
  deaths ~ cb_summer_65_74 +
    factor(dow) +
    splines::ns(doy, df = 2 * 11) +
    factor(covid_period),
  eliminate = factor(strata),
  family = quasipoisson,
  data = summer_df_65_74
)


temp_grid <- seq(
  from = quantile(summer_df_65_74$tasmax, 0.01, na.rm = TRUE),
  to   = quantile(summer_df_65_74$tasmax, 0.99, na.rm = TRUE),
  by   = 0.1
)



class(cb_summer_65_74)

cp_65_74 = crosspred(cb_summer_65_74,fit_65_74,cen=16, at=temp_grid,cumul = TRUE)



# Identify cumulative lag column
j= ncol(cp_65_74$cumfit)

keep = cp_65_74$predvar >= quantile(summer_df_65_74$tasmax, 0.1) &
  cp_65_74$predvar <= quantile(summer_df_65_74$tasmax, 0.9)


mmt = cp_65_74$predvar[keep][
  which.min(cp_0_64$cumRRfit[keep, j])
]

rr_df <- data.frame(
  temp    = cp_65_74$predvar,
  rr      = cp_65_74$cumRRfit[, j],
  rr_low  = cp_65_74$cumRRlow[, j],
  rr_high = cp_65_74$cumRRhigh[, j]
)

ggplot(rr_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high),
              fill = "grey80") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = mmt, linetype = "dashed") +
  labs(
    x = "summer mean temperature (°C)",
    y = "Cumulative RR (lag 0–21)",
    title = "summer temperature–mortality association (Age 65-74)"
  ) +
  theme_minimal()


##################################################################

summer_df_0_64 = df_summer %>% 
  filter(age_group == "0-64") %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()

cb_summer_0_64 <- crossbasis(
  summer_df_0_64$tasmax,
  lag = 3,
  argvar = list(fun = "ns", knots = quantile(summer_df_0_64$tasmax, c(0.5,0.9), na.rm = TRUE)),
  arglag = list(fun = "ns", knots = 1)
)

fit_0_64  = gnm(
  deaths ~ cb_summer_0_64 +
    factor(dow) +
    splines::ns(doy, df = 2 * 11) +
    factor(covid_period),
  eliminate = factor(strata),
  family = quasipoisson,
  data = summer_df_0_64
)


temp_grid <- seq(
  from = quantile(summer_df_0_64$tasmax, 0.01, na.rm = TRUE),
  to   = quantile(summer_df_0_64$tasmax, 0.99, na.rm = TRUE),
  by   = 0.1
)



class(cb_summer_0_64)

cp_0_64 = crosspred(cb_summer_0_64,fit_0_64,cen=16, at=temp_grid,cumul = TRUE)



# Identify cumulative lag column
j= ncol(cp_0_64$cumfit)

keep = cp_0_64$predvar >= quantile(summer_df_0_64$tasmax, 0.1) &
  cp_0_64$predvar <= quantile(summer_df_0_64$tasmax, 0.9)


mmt = cp_0_64$predvar[keep][
  which.min(cp_0_64$cumRRfit[keep, j])
]

rr_df <- data.frame(
  temp    = cp_0_64$predvar,
  rr      = cp_0_64$cumRRfit[, j],
  rr_low  = cp_0_64$cumRRlow[, j],
  rr_high = cp_0_64$cumRRhigh[, j]
)

ggplot(rr_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high),
              fill = "grey80") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = mmt, linetype = "dashed") +
  labs(
    x = "summer mean temperature (°C)",
    y = "Cumulative RR (lag 0–21)",
    title = "summer temperature–mortality association (Age 60-64)"
  ) +
  theme_minimal()





