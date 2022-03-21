globalVariables(c("Year", "prediction"))

##' Calculate the model
##' @param dat data.frame with data (Colnames Time, Year, ID must be available)
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
  # es kann Extremfälle geben, bei denen alle Werte gleich sind (z.B. alle 0)
  # dann kann der Durbin-Watson-Test keinen p-Wert berechnen und das lm auch
  # nicht. In diesem Fall wird keine Autokorrelation und nur ein lineares
  # Modell angenommen
  if (!is.na(DW$p[1])) {
    if (easy) { # Quadratischer Term wegen n < 12 gar nicht erlaubt
      if (DW$p[1] < alpha ) { # d.h signifikante Autokorrelation 1. Ordnung
        ## Korrelationsstruktur festlegen
        cor_dat <- corAR1(form = ~ Time|ID)
        mod <- gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
      } else {
        if (DW$p[2] < alpha ) { # d.h signifikante Autokorrelation 2. Ordnung
          ## Korrelationsstruktur festlegen
          cor_dat <- corARMA(form = ~ Time|ID, p = 2)
          mod <- gls(mod_form_linear, data = dat, method="ML", correlation=cor_dat)
        } else { # d.h keine signifikante Autokorrelation(bis 2. Ordnung) nachweisbar
          mod <- lm(mod_form_linear, dat)
          cor_dat <- NULL
        }
      }
      used_formula <- mod_form_linear
    } else { # Quadratischer Term grundsaetzlich schon moeglich, beta2 wird aber auf signifikanz geprueft
      if (DW$p[1] < alpha ) { # d.h signifikante Autokorrelation 1. Ordnung

        ## Korrelationsstruktur festlegen
        cor_dat <- corAR1(form = ~ Time|ID)

        ## AR1 mit ML fitten
        mod <- gls(mod_form, data = dat, method = "ML", correlation = cor_dat)
        used_formula <- mod_form

        ## Pruefe ob beta2 signifikant von null verschieden
        smod <- summary(mod)
        beta2 <- smod$tTable["I(Time^2)", "Value"]
        p_betas <- smod$tTable[, "p-value"]
        p_beta2 <- p_betas["I(Time^2)"]
        beta2_zero <- p_beta2 > alpha

        ## falls beta2 nicht signifikant oder nicht vorhanden, fitte ohne quadratischen Term
        if( beta2_zero) {
          mod <- gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
          used_formula <- mod_form_linear
        }

      } else {

        if (DW$p[2] < alpha ) { # d.h signifikante Autokorrelation 2. Ordnung
          ## Korrelationsstruktur festlegen
          cor_dat <- corARMA(form = ~ Time|ID, p = 2)

          ## AR1 mit ML fitten
          mod <- gls(mod_form_quadratic, data = dat, method = "ML", correlation = cor_dat)
          used_formula <- mod_form_quadratic

          ## Pruefe ob beta2 signifikant von null verschieden
          smod <- summary(mod)
          beta2 <- smod$tTable["I(Time^2)", "Value"]
          p_betas <- smod$tTable[, "p-value"]
          p_beta2 <- p_betas["I(Time^2)"]
          beta2_zero <- p_beta2 > alpha

          ## falls beta2 nicht signifikant, fitte ohne quadratischen Term
          if (beta2_zero) {
            mod <- gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
            used_formula <- mod_form_linear
          }

        } else { # d.h keine signifikante Autokorrelation(bis 2. Ordnung) nachweisbar

          ### LiMo(quad)
          mod <- lm(mod_form_quadratic, dat)
          cor_dat <- NULL
          used_formula <- mod_form_quadratic
          mod_summary <- summary(mod)
          betaP <- mod_summary$coefficients["I(Time^2)","Pr(>|t|)"]

          if (betaP > alpha) { ### LiMo(einfach)
            mod <- lm(mod_form_linear, dat)
            used_formula <- mod_form_linear
          }

        }
      }
    }
  } else {
    # Der Fall wenn kein DW-Test berechnet werden kann
    mod <- lm(mod_form_linear, dat)
    used_formula <- mod_form_linear
    cor_dat <- NULL
  }
  return(list(mod = mod, used_formula = used_formula,
              correlation_structure = cor_dat))
}
