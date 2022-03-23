globalVariables(c("Jahr", "prediction"))

##' Calculate the model
##' @param dat data.frame with data (Colnames Time, Jahr, ID must be available)
##' @param varname String with variablenames
##' @param alpha Alphavalue for area of rejection
##' @importFrom stats as.formula lm
##' @importFrom car dwt
##' @importFrom nlme corAR1 gls corARMA
##' @return lm or gls - model (linear, ML)
fit_trend <- function(dat, varname, alpha = 0.05){

  ## quadratic model formula
  mod_form_quadratic <- as.formula(sprintf("%s ~ Time + I(Time^2)", varname))
  ## linear model forumla
  mod_form_linear <- as.formula(sprintf("%s ~ Time", varname))

  ## Decide by number of observations whether a quadratic or linear term is used
  target_var <- dat[, varname]
  if (dim(dat)[1] > 12) {  # n > 12 => quadratic else linear
    mod_form <- mod_form_quadratic
    easy <- FALSE
  } else {
    mod_form <- mod_form_linear
    easy <- TRUE
  }

  ## Durbin-Watson-Test
  DW <- dwt(lm(mod_form, data = dat), max.lag = 2) # check if AR1 equal to AR2

  cor_dat <- NULL
  # In some extreme cases it is possible for all values to be equal (or. 0)
  # In this case the DWT and lm cant calculate a p-value
  # Therefore no autocorrelation is accepted and a linear model is used
  if (!is.na(DW$p[1])) {
    # Quadratic  term for n < 12 not allowed
    if (easy) {
      # significant autocorrelation 1. Order
      if (DW$p[1] < alpha ) {
        ## define correlation structure
        cor_dat <- corAR1(form = ~ Time|ID)
        mod <- gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
      } else {
        # significant autocorrelation 2. order
        if (DW$p[2] < alpha ) {
          # define correlation structure
          cor_dat <- corARMA(form = ~ Time|ID, p = 2)
          mod <- gls(mod_form_linear, data = dat, method="ML", correlation = cor_dat)
        } else {
          # no significant autocorrelation (up to 2nd order) can be shown
          mod <- lm(mod_form_linear, dat)
          cor_dat <- NULL
        }
      }
      used_formula <- mod_form_linear
    } else {
      # Quadratic term is possible beta2 is checked for significance
      # significant autocorrelation 2. order
      if (DW$p[1] < alpha ) {

        ## define correlation
        cor_dat <- corAR1(form = ~ Time|ID)

        ## AR1 fitted with ML
        mod <- gls(mod_form, data = dat, method = "ML", correlation = cor_dat)
        used_formula <- mod_form

        ## Check if beta2 is significantly different from 0
        smod <- summary(mod)
        beta2 <- smod$tTable["I(Time^2)", "Value"]
        p_betas <- smod$tTable[, "p-value"]
        p_beta2 <- p_betas["I(Time^2)"]
        beta2_zero <- p_beta2 > alpha

        ## if beta2 is not significant or not available fit without quadratic term
        if( beta2_zero) {
          mod <- gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
          used_formula <- mod_form_linear
        }

      } else {
        # significant autocorrelation 2. order
        if (DW$p[2] < alpha ) {
          ## define correlation structure
          cor_dat <- corARMA(form = ~ Time|ID, p = 2)

          ## AR1 fitted with ML
          mod <- gls(mod_form_quadratic, data = dat, method = "ML", correlation = cor_dat)
          used_formula <- mod_form_quadratic

          ## Check if beta 2 is significantly different form 0
          smod <- summary(mod)
          beta2 <- smod$tTable["I(Time^2)", "Value"]
          p_betas <- smod$tTable[, "p-value"]
          p_beta2 <- p_betas["I(Time^2)"]
          beta2_zero <- p_beta2 > alpha

          ## if beta2 is not significant fit without quadratic term
          if (beta2_zero) {
            mod <- gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
            used_formula <- mod_form_linear
          }

        } else {
          # no significant autocorrelation (up to 2. order) can be shown
          ### LiMo(quad)
          mod <- lm(mod_form_quadratic, dat)
          cor_dat <- NULL
          used_formula <- mod_form_quadratic
          mod_summary <- summary(mod)
          betaP <- mod_summary$coefficients["I(Time^2)","Pr(>|t|)"]
          ### LiMo(simple)
          if (betaP > alpha) {
            mod <- lm(mod_form_linear, dat)
            used_formula <- mod_form_linear
          }

        }
      }
    }
  } else {
    # If DW-Test cant be calculated fit a linear model
    mod <- lm(mod_form_linear, dat)
    used_formula <- mod_form_linear
    cor_dat <- NULL
  }
  return(list(mod = mod, used_formula = used_formula,
              correlation_structure = cor_dat))
}
