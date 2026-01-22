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
    year = year(date),
    qtr = quarter(date),
    #this is the part to be conditioned out of the regression
    strata = factor(paste0(lsoa11,"-",ym)),
    #made age group as factor 
    age_group = factor(age_group, levels = c("0-64",
                                             "65-74",
                                             "75-84",
                                             "85+")),
    # numeric calendar time for long-term trend
    doy = yday(date)
    
  )



########################################################
#define the spline for day of the year as long term seasonal trend

#calucalte the degree of freedom
df_per_year = 7
n_years = as.numeric(difftime(max(df$date), min(df$date), units = "days")) / 365.25
df_time = ceiling(df_per_year * n_years)

#the spline line
# spldoy = onebasis(df$doy, "ns", df = df_per_year)

########################################################
#add a variable to take into account of covid period 
#check the covid period for birmingham

covid_check = df %>% 
  filter(date>= "2018-01-01" & date <="2024-01-01") %>% 
  group_by(date) %>% 
  summarise(deaths = sum(deaths))


#create a variable to df about the period of covid 
#must not use calendar time directly, otherwise the spline leaks outside COVID
df = df %>% 
  mutate(date = as.Date(date),  #  fix at the source
         covid_period = date >= as.Date("2020-03-01") &
           date <= as.Date("2021-07-31"),
         covid_time = if_else(
           covid_period,
           as.numeric(date - as.Date("2020-03-01")),
           NA_real_
         ))


#########################################################
#define the cross basis components for temp 
#data driven knots
#This gives the spline flexibility across Birmingham’s typical temperature range, including the cooler and warmer tails,
# which helps capture the non-linear risk shape without overfitting the very extremes
argvar = list(fun = "bs", knots = quantile(df$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE))

# knots = quantile(df$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
# lag knots 
#Using logknots(lag_max, 3) is the gold standard for temperature.
#The log-scale puts more "statistical attention" (knots) on the early days (0, 1, 2) and spreads them out as you get toward day 21.
#puts more knots at the beginning (where the mortality effect changes quickly, like during a heatwave) and 
#fewer knots at the end (where the effect of a cold spell stays more constant)
arglag = list(fun = "ns", knots = logknots(21,3))


##########################################################
# #made age group as factor 
# df = df %>% 
#   mutate(age_group = factor(age_group, levels = c("0-64",
#                                                   "65-74",
#                                                   "75-84",
#                                                   "85+"))


#########################################################
#initiates a list to store the model results
age_levels = levels(df$age_group)
fit_mod = vector("list", length(age_levels))

names(fit_mod) = age_levels

########################################################
#to create the time shifted copies for lags
cbs  = vector("list", length(age_levels))  # store crossbasis used in each fit
names(cbs) =  age_levels

########################################################
#calucalte the degree of freedom
# df_per_year = 7
# n_years = as.numeric(difftime(max(df$date), min(df$date), units = "days")) / 365.25
# df_time = ceiling(df_per_year * n_years)


########################################################

#data driven knots
#This gives the spline flexibility across Birmingham’s typical temperature range, including the cooler and warmer tails,
# which helps capture the non-linear risk shape without overfitting the very extremes

# temp_knots = quantile(df$tasmean, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)

# lag knots 
#Using logknots(lag_max, 3) is the gold standard for temperature.
#The log-scale puts more "statistical attention" (knots) on the early days (0, 1, 2) and spreads them out as you get toward day 21.
#puts more knots at the beginning (where the mortality effect changes quickly, like during a heatwave) and 
#fewer knots at the end (where the effect of a cold spell stays more constant)
# lag_knots = logknots(21, 3)

########################################################
#extremely omprtant 
#remove strate with 0 outcomes 
#if not removed, bias dispersion estimation, cause instability

df = df %>%
  group_by(strata) %>%
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()


########################################################
for (age in age_levels){
  data = df %>% 
    filter(age_group == age)
  
  #cross basis allows us to estimate non-linear effect of temp and delayed effect over time 
  cb = crossbasis(
    data$tasmean,
    lag= c(0,28), 
    argvar = argvar,
    arglag = arglag
  )
  
  fit = gnm(deaths ~ cb +
            factor(dow) +
            ns(doy, df = df_per_year):factor(year),
            eliminate = factor(strata),
            family = quasipoisson,
            data = data)
  
  fit_mod[[age]] = fit
  cbs[[age]] = cb
  
}

##############################################################################
#load the find min mmt fuinctiion 
source("R/find_min_MMT.R")
#######################################################################



res_0_64= find_MMT(
  ag      = "0-64",
  df      = df,
  cbs     = cbs,
  fit_mod = fit_mod,
  qrange  = c(0.02,0.98)
  
)

res_65_74= find_MMT(
  ag      = "65-74",
  df      = df,
  cbs     = cbs,
  fit_mod = fit_mod,
  qrange  = c(0.02,0.98)
  
)


res_75_84= find_MMT(
  ag      = "75-84",
  df      = df,
  cbs     = cbs,
  fit_mod = fit_mod,
  qrange  = c(0.02,0.98)
)


res_85= find_MMT(
  ag      = "85+",
  df      = df,
  cbs     = cbs,
  fit_mod = fit_mod,
  qrange  = c(0.02,0.98)
)

res_85$rr_curve


res_0_64$mmt
res_65_74$mmt
res_75_84$mmt
res_85$mmt


cp = crosspred(cb, fit_mod[["85+"]], cen=18)

plot(cp, "overall")


summary(fit_mod[["0-64"]])


ggarrange(res_0_64$plot, res_65_74$plot, res_75_84$plot, res_85$plot)

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
ag = "65-74"
d  = df[df$age_group == ag, ]
cb  = cbs[[ag]]          # use the crossbasis you fitted with
fit = fit_mod[[ag]]      # fitted glm

# Predict the effect relative to a reference temperature (e.g., 15°C)
pred = crosspred(cb, fit, by = 0.1, cen = median(d$tasmean, na.rm = TRUE), cumul = TRUE)

# lags to show
lags <- c(0, 1, 2, 3, 4, 5, 7, 14, 21)
lag_idx <- lags + 1   # R is 1-indexed

# temperature grid
temp <- pred$predvar

# restrict to supported temperature range
t_rng <- quantile(d$tasmean, c(0.01, 1), na.rm = TRUE)
keep  <- temp >= t_rng[1] & temp <= t_rng[2]

# build long dataframe
df_lag <- expand.grid(
  temp = temp[keep],
  lag  = lags
) %>%
  mutate(
    rr = as.vector(pred$matRRfit[keep, lag_idx])
  )


ggplot(df_lag, aes(x = temp, y = rr, colour = factor(lag))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_viridis_d(
    name = "Lag (days)",
    option = "plasma"
  ) +
  coord_cartesian(ylim = c(0.95, 100)) +
  labs(
    title = "Lag-specific temperature effects",
    x = "Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

#########################################################################
age_groups <- unique(df$age_group)
lags <- 0:21
lag_idx <- lags + 1

lag_response_all <- map_dfr(age_groups, function(ag) {
  
  d   <- df %>% filter(age_group == ag)
  cb  <- cbs[[ag]]
  fit <- fit_mod[[ag]]
  
  if (is.null(cb) || is.null(fit) || nrow(d) == 0)
    return(NULL)
  
  pred <- crosspred(
    cb,
    fit,
    by    = 0.1,
    cen   = median(d$tasmean, na.rm = TRUE),
    cumul = FALSE   # 🔴 REQUIRED
  )
  
  # temperatures to evaluate (percentiles)
  temps <- quantile(d$tasmean, c(0.0001), na.rm = TRUE)
  
  tempstable = data.frame(
    temp_lab = names(temps),
    temp = as.numeric(temps)
  )
  
  
  temp_idx <- sapply(temps, function(x)
    which.min(abs(pred$predvar - x))
  )
  
  expand.grid(
    age_group = ag,
    lag       = lags,
    temp_lab  = names(temps)
  ) %>%
    mutate(
      rr = as.vector(
        sapply(seq_along(temp_idx), function(i)
          pred$matRRfit[temp_idx[i], lag_idx]
        )
      ),
      temp_lab = factor(
        temp_lab,
        labels = paste0(tempstable$temp_lab, " percentile", "\n", round(tempstable$temp,2), "°C")
      )
    )
})





ggplot(lag_response_all,
       aes(x = lag, y = rr, colour = temp_lab)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ age_group, scales = "free_y") +
  scale_colour_viridis_d(name = "Temperature") +
  labs(
    title = "Lag–response relationship at high temperatures",
    x = "Lag (days)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )









cut(df$tasmean, breaks = quantile(df$tasmean, probs = c(0, 0.01, 0.05, 0.1, 0.9, 0.95, 0.99, 1), na.rm = TRUE), include.lowest = TRUE) %>% 
  table()

###########################################################################

cb  <- cbs[["75-84"]]
fit <- fit_mod[["75-84"]]


dev.off()


percentiles <- round(quantile(df$tasmean,c(0.0001,0.01, 0.99,0.9999)),1)
cp = crosspred(cb,fit,at = -43:282/10,cen=22)
par(mar = c(5, 4, 2, 1))
plot(cp,var=percentiles, lag=c(0,5,15,28))



plot(cp, "overall")













##########################################################################

plot(pred2, lag=c(21))






rr_curve = data.frame(
  temp   = pred_cum$predvar,
  rr     = pred_cum$cumRRfit[, j]  / pred_cum$cumRRfit[mmt_idx, j],
  rr_low = pred_cum$cumRRlow[, j]  / pred_cum$cumRRfit[mmt_idx, j],
  rr_high= pred_cum$cumRRhigh[, j] / pred_cum$cumRRfit[mmt_idx, j]
)


pred$cumRRfit[, j]


#########################################################################

ag = "65-74"
d  = df[df$age_group == ag, ]

summary(d$tasmean)
quantile(d$tasmean, c(.01,.05,.50,.95,.99), na.rm=TRUE)
table(cut(d$tasmean, breaks=c(-10,0,5,10,15,20,25,40)))
# 
# 


cb  = cbs[[ag]]          # use the crossbasis you fitted with
fit = fit_mod[[ag]]      # fitted glm

# Predict the effect relative to a reference temperature (e.g., 15°C)
pred = crosspred(cb, fit, by = 0.1, cen = median(d$tasmean, na.rm = TRUE), cumul = TRUE)

#Restrict to 1st–99th percentile temperature range (paper-style)
t_rng = quantile(d$tasmean, probs = c(0.01, 0.99), na.rm = TRUE)

keep = pred$predvar >= t_rng[1] & pred$predvar <= t_rng[2]

mmt = pred$predvar[keep][ which.min(pred$allRRfit[keep]) ]

j = ncol(pred2$cumfit)



rr = data.frame(
  temp = pred2$predvar,
  rr   = exp(pred2$cumfit[, j])
)

plot(
  rr$temp, rr$rr,
  type = "l",
  xlab = "Temperature",
  ylab = "RR",
  main = "Cumulative RR (0–21 days), 65-74"
)
abline(h = 1, lty = 2)
abline(v = mmt, col = "red", lty = 2)

pred2 = crosspred(cb, fit, by = 0.1, cen = mmt, cumul = TRUE)
# Plot the "Overall" effect (sum of all lags)
plot(pred2, "overall", xlab="Temperature", ylab="RR", main="Overall Effect for 65-74")

#the graph is not what we expected, this could due to uncontrained dlnm 
#many cross-basis and hard to follow every small wiggle
#highly correlated across several days, an unconstrained model might assign a negative coefficient to one lag and a positive to another 
#just to "fit" the noise, resulting in this non-sensical dip.
