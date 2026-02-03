#this is to find the MMT

find_MMT <- function(
    ag,                 # age group
    df,                 # model data
    cbs,                # crossbasis list
    fit_mod,            # fitted model list
    by      = 0.1,      
    qrange  = c(0.01, 0.99),
    cen     = 10,
    plot    = TRUE,
    ci_col  = "grey80"
) {
  
  # --------------------------------------------------
  # 0. Subset data
  # --------------------------------------------------
  d <- df[df$age_group == ag, ]
  if (nrow(d) == 0)
    stop("No data found for age group: ", ag)
  
  cb  <- cbs[[ag]]
  fit <- fit_mod[[ag]]
  if (is.null(cb) || is.null(fit))
    stop("Missing crossbasis or model for age group: ", ag)
  
  # Temperature range for MMT search
  t_rng <- quantile(d$tasmean, qrange, na.rm = TRUE)
  
  # --------------------------------------------------
  # 1. Cumulative prediction (lag 0–21)
  # --------------------------------------------------
  pred_cum <- crosspred(
    cb,
    fit,
    by    = by,
    cen   = cen,
    cumul = TRUE
  )
  
  # Identify cumulative lag column
  j <- ncol(pred_cum$cumfit)
  
  # Restrict MMT search to central range
  keep <- pred_cum$predvar >= t_rng[1] &
    pred_cum$predvar <= t_rng[2]
  
  mmt <- pred_cum$predvar[keep][
    which.min(pred_cum$cumRRfit[keep, j])
  ]
  mmt_idx <- which.min(abs(pred_cum$predvar - mmt))
  
  # --------------------------------------------------
  # 2. Correct log-scale re-centred RR + CI
  # --------------------------------------------------
  eta     <- pred_cum$cumfit[, j]      # log RR
  se      <- pred_cum$cumse[, j]       # SE on log scale
  eta_ref <- eta[mmt_idx]
  
  rr_curve <- data.frame(
    temp    = pred_cum$predvar,
    rr      = exp(eta - eta_ref),
    rr_low  = exp((eta - eta_ref) - 1.96 * se),
    rr_high = exp((eta - eta_ref) + 1.96 * se)
  )
  
  # --------------------------------------------------
  # 3. Plot
  # --------------------------------------------------
  gg <- NULL
  
  if (plot) {
    
    gg <- ggplot(rr_curve, aes(x = temp)) +
      
      geom_ribbon(
        aes(ymin = rr_low, ymax = rr_high),
        fill = ci_col,
        alpha = 0.6
      ) +
      
      geom_line(
        aes(y = rr),
        colour = "red",
        linewidth = 1
      ) +
      
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = mmt, linetype = "dashed") +
      
      labs(
        title = paste0("Age ", ag),
        x = "Temperature (°C)",
        y = "Cumulative relative risk (lag 0–21)"
      ) +
      
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    print(gg)
  }
  
  # --------------------------------------------------
  # 4. Return
  # --------------------------------------------------
  return(list(
    age_group = ag,
    mmt       = mmt,
    rr_curve  = rr_curve,
    plot      = gg
  ))
}
