utils::globalVariables(c("Jahr", "prediction"))

##' Calculates the trend for data and variables
##' Die Funktion checkt ob genug Daten vorliegen und laesst dann das Modell rechnen.
##' @param dat dataset
##' @param varname colnames of dataset
##' @param log_trans defaults to FALSE and indicates whether variables of data
##' are log transformed or not
##' @param calc_infl_obs logical
##' @return Eine Liste mit Trendgrafike, Variablenname, Modell und Text zum Ergebnis
compute_trend <- function(dat, varname, log_trans = FALSE, calc_infl_obs = TRUE){

  ## encode variables numericly
  dat[, varname] <- as.numeric(dat[ ,varname])
  dat <- stats::na.omit(dat[, c("Zeit", "Jahr", "ID", varname)])

  ## check whether enough data points are avaliable
  enough_data <- check_n(dat[, varname, drop = TRUE])

  ## result, if not enough data points are avaliable
  if (!enough_data) {
    return(list(plot = ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes_string(x = "Jahr", y = varname)),
                varname = varname,
                mod = NULL, trend = "Keine Trendanalyse moeglich",
                phi = c(NA, NA),
                beta = list(beta0 = NA, beta1 = NA, beta2 = NA),
                rsq = NA,
                einfl_beob_index = NA,
                einfl_beob_value = NA,
                einfl_beob_cookd = NA,
                einfl_beob_jahr = NA))
  }
  ## fit the model
  model_data <- fitte_trend(dat, varname = varname)
  mod <- model_data[["mod"]]
  rsq <- round(1 - (sum(mod$residuals ^ 2) / sum(((mod$residuals + mod$fitted) - mean(mod$residuals + mod$fitted)) ^ 2)), 4)
  if (calc_infl_obs) {
  einfl_beobachtungen <- infl_obs(mod = mod,
                                                      used_formula =
                                                        model_data[["used_formula"]],
                                                      cor_struct =
                                                        model_data[["correlation_structure"]],
                                                      data = dat,
                                                      varname = varname)
  } else {
  einfl_beobachtungen <- list(einfl_beob_index = NA, einfl_beob_value = NA,
                              einfl_beob_cookd = NA,
                              einfl_beob_jahr = NA)
  }

  list(plot = plot_trend(mod, df = dat, log_trans = log_trans), varname = varname,
       mod = mod, trend = extrahiere_trend(mod)[4], phi = extrahiere_phi(mod),
       beta = extrahiere_trend(mod)[1:3], rsq = max(0,rsq),
       einfl_beob_index = paste0(einfl_beobachtungen[[1]], collapse = ", "),
       einfl_beob_value = paste(einfl_beobachtungen[[2]], collapse = ", "),
       einfl_beob_cookd = paste(einfl_beobachtungen[[3]], collapse = ", "),
       einfl_beob_jahr = paste(einfl_beobachtungen[[4]], collapse = ", "))
}



##' Checks if enough data points are avaliable
##'@param trend_variable vector with value
##'@param min_n min value for trendanalysis
##'@result TRUE / FALSE
check_n <- function(trend_variable, min_n = 7) {
  length(stats::na.omit(trend_variable)) >= min_n
}

##' Fit the model
##' @param dat data.frame (Column Zeit, Jahr, ID must be present)
##' @param varname string with variablename
##' @param alpha aplhavalue for region of rejection, standard is set to 0.05
##' @return lm or gls - model (linear, ML)
fitte_trend <- function(dat, varname, alpha = 0.05){

  ## quadratic model formula
  mod_form_quadratisch <- stats::as.formula(sprintf("%s ~ Zeit + I(Zeit^2)", varname))
  ## linear model formula
  mod_form_linear <- stats::as.formula(sprintf("%s ~ Zeit", varname))

  ## decide base on number of observations whether quadratic terms or linear
  ## terms are used
  Zielvar <- dat[, varname]
  if (dim(dat)[1] > 12) {  # n > 12 => quadratic, else linear
    mod_form <- mod_form_quadratisch
    easy <- FALSE
  } else {
    mod_form <- mod_form_linear
    easy <- TRUE
  }

  ## Durbin-Watson-Test
  # Check if AR1 und AR2 are equal
  DW <- car::dwt(stats::lm(mod_form, data = dat), max.lag = 2)

  cor_dat <- NULL
  # there can be extreme cases where all values are the same (e.g. all 0)
  # then the Durbin-Watson test cannot calculate a p-value and neither can the lm
  # In this case there is no autocorrelation and only a linear is model adopted
  if (!is.na(DW$p[1])) {
  if (easy) { # Quadratischer Term wegen n < 12 gar nicht erlaubt
    if (DW$p[1] < alpha ) { # d.h signifikante Autokorrelation 1. Ordnung
      ## Korrelationsstruktur festlegen
      cor_dat <- nlme::corAR1(form = ~ Zeit|ID)
      mod <- nlme::gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
    } else {
      if (DW$p[2] < alpha ) { # d.h signifikante Autokorrelation 2. Ordnung
        ## Korrelationsstruktur festlegen
        cor_dat <- nlme::corARMA(form = ~ Zeit|ID, p = 2)
        mod <- nlme::gls(mod_form_linear, data = dat, method="ML", correlation=cor_dat)
      } else { # d.h keine signifikante Autokorrelation(bis 2. Ordnung) nachweisbar
        mod <- stats::lm(mod_form_linear, dat)
        cor_dat <- NULL
      }
    }
    used_formula <- mod_form_linear
  } else { # Quadratischer Term grundsaetzlich schon moeglich, beta2 wird aber auf signifikanz geprueft
    if (DW$p[1] < alpha ) { # d.h signifikante Autokorrelation 1. Ordnung

      ## Korrelationsstruktur festlegen
      cor_dat <- nlme::corAR1(form = ~ Zeit|ID)

      ## AR1 mit ML fitten
      mod <- nlme::gls(mod_form, data = dat, method = "ML", correlation = cor_dat)
      used_formula <- mod_form

      ## Pruefe ob beta2 signifikant von null verschieden
      smod <- summary(mod)
      beta2 <- smod$tTable["I(Zeit^2)", "Value"]
      p_betas <- smod$tTable[, "p-value"]
      p_beta2 <- p_betas["I(Zeit^2)"]
      beta2_zero <- p_beta2 > alpha

      ## falls beta2 nicht signifikant oder nicht vorhanden, fitte ohne quadratischen Term
      if( beta2_zero) {
        mod <- nlme::gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
        used_formula <- mod_form_linear
      }

    } else {

      if (DW$p[2] < alpha ) { # d.h signifikante Autokorrelation 2. Ordnung
        ## Korrelationsstruktur festlegen
        cor_dat <- nlme::corARMA(form = ~ Zeit|ID, p = 2)

        ## AR1 mit ML fitten
        mod <- nlme::gls(mod_form_quadratisch, data = dat, method = "ML", correlation = cor_dat)
        used_formula <- mod_form_quadratisch

        ## Pruefe ob beta2 signifikant von null verschieden
        smod <- summary(mod)
        beta2 <- smod$tTable["I(Zeit^2)", "Value"]
        p_betas <- smod$tTable[, "p-value"]
        p_beta2 <- p_betas["I(Zeit^2)"]
        beta2_zero <- p_beta2 > alpha

        ## falls beta2 nicht signifikant, fitte ohne quadratischen Term
        if (beta2_zero) {
          mod <- nlme::gls(mod_form_linear, data = dat, method = "ML", correlation = cor_dat)
          used_formula <- mod_form_linear
        }

      } else { # d.h keine signifikante Autokorrelation(bis 2. Ordnung) nachweisbar

        ### LiMo(quad)
        mod <- stats::lm(mod_form_quadratisch, dat)
        cor_dat <- NULL
        used_formula <- mod_form_quadratisch
        mod_summary <- summary(mod)
        betaP <- mod_summary$coefficients["I(Zeit^2)","Pr(>|t|)"]

        if (betaP > alpha) { ### LiMo(einfach)
          mod <- stats::lm(mod_form_linear, dat)
          used_formula <- mod_form_linear
        }

      }
    }
  }
  } else {
    # Der Fall wenn kein DW-Test berechnet werden kann
    mod <- stats::lm(mod_form_linear, dat)
    used_formula <- mod_form_linear
    cor_dat <- NULL
  }
  return(list(mod = mod, used_formula = used_formula,
              correlation_structure = cor_dat))
}



##' Extrahiere phi aus gls - Modell
##' @param mod gls-Modell
##' @return phi Parameter der Korrelation
extrahiere_phi <- function(mod) {
  if (inherits(mod, "gls")) { # Die Funktion ist fuer gls-Objekt spezifiziert
    x <- summary(mod$modelStruct)$corStruct
    class(x) <- attr(x, "oClass")
    phi <- round(stats::coef(x, unconstrained = FALSE), 2)
  } else { # fuer lineares Modell betrdgt phi eben 0
    phi <- 0
  }
  if (length(phi) == 1) {phi[2] <- 0}
  return(as.character(phi))
}



##' Extrahiere Trend aus Modell
##' @param model gls oder lm Modell
##' @param alpha Signifikanzniveau
##' @param ... zur Zeit nicht benutzt
##' @return String mit Trendergebnis
extrahiere_trend <- function(model, alpha, ...) {
  UseMethod("extrahiere_trend")
}

# Fuer lm-Objekte
extrahiere_trend.lm <- function(model, alpha = 0.05) {
  ## extrahiere p-wert fuer linearen term
  smod <- summary(model)
  beta0 <- smod$coefficients["(Intercept)", "Estimate"]
  beta1 <- smod$coefficients["Zeit", "Estimate"]
  beta2 <- as.character(0)
  p_beta1 <- smod$coefficients["Zeit", "Pr(>|t|)"]

  ## enthaelt das Model einen quadratischen Term, so liegt ein quadratischer Trend vor
  quadratisch <- "I(Zeit^2)" %in% rownames(smod$coefficients)

  if (quadratisch){
    ## da in Funktion "fitte_trend" schon auf sign. geprueft wurde,
    ## kann man von einem quadratischen Trend ausgehen
    beta2 <- smod$coefficients["I(Zeit^2)", "Estimate"]
    extrem <- (-beta1 / 2 / beta2)
    maxi <- max(as.numeric(model$model$Zeit))
    mini <- min(as.numeric(model$model$Zeit))
    if (extrem >= mini && extrem <= maxi) {
      trend <- ifelse(beta2 > 0, "Quadratischer Trend(u)", "Quadratischer Trend(n)")
    } else {
      trend <- ifelse(model$fitted.values[1] < model$fitted.values[length(model$fitted.values)],
                      "Quadratischer Trend(+)", "Quadratischer Trend(-)")
    }
  } else {
    # überprüfe, ob überhaupt ein p-Wert berechnet werden konnte
    if (!is.na(p_beta1)) {
      if (p_beta1 < alpha) {
        trend <- ifelse(beta1 < 0, "Fallender Trend", "Steigender Trend")
      } else {
        trend <- "Kein signifikanter Trend"
      }
    } else {
      # wenn kein p-Wert berechnet werden konnte
      trend <- "es konnte kein p-Wert berechnet werden"
    }
  }
  return(list(beta0, beta1, beta2, trend))
}

# Fuer gls-Objekte
extrahiere_trend.gls <- function(model, alpha = 0.05){

  ## extrahiere p-wert fuer linearen term
  smod <- summary(model)
  beta0 <- smod$tTable["(Intercept)", "Value"]
  beta1 <- smod$tTable["Zeit", "Value"]
  beta2 <- as.character(0)
  p_beta1 <- smod$tTable["Zeit", "p-value"]

  ## enthaelt das Model einen quadratischen Term, so liegt ein quadratischer Trend vor
  quadratisch <- "I(Zeit^2)" %in% rownames(smod$tTable)

  if (quadratisch) {
    ## da in Funktion "fitte_trend" schon auf sign. geprueft wurde,
    ## kann man von einem quadratischen Trend ausgehen
    beta2 <- smod$tTable["I(Zeit^2)", "Value"]
    extrem <- (-beta1 / 2 / beta2)
    maxi <- max(as.numeric(names(model$fitted)))
    mini <- min(as.numeric(names(model$fitted)))
    if (extrem >= mini && extrem <= maxi) {
      trend <- ifelse(beta2 > 0, "Quadratischer Trend(u)", "Quadratischer Trend(n)")
    } else {
      trend <- ifelse(model$fitted[1] < model$fitted[length(model$fitted)],
                      "Quadratischer Trend(+)", "Quadratischer Trend(-)")
    }
  } else {
    if (p_beta1 < alpha) {
      trend <- ifelse(beta1 < 0, "Fallender Trend", "Steigender Trend")
    } else {
      trend <- "Kein signifikanter Trend"
    }
  }
  return(list(beta0, beta1, beta2, trend))
}



##' Plots for the trend model
##'
##' Plots data as Points and trends as lines
##'
##' @param mod lm or gls(nlme)- Modell
##' @param df Original data
##' @param log_trans gibt an, ob die modellierten Daten log-transformiert sind
##' @return plot - object
plot_trend <- function(mod, df, log_trans = FALSE) {
  rsq <- round(1 - (sum(mod$residuals ^ 2) / sum(((mod$residuals + mod$fitted) - mean(mod$residuals + mod$fitted)) ^ 2)), 4)

  trend <- extrahiere_trend(mod)
  phi <- extrahiere_phi(mod)
  if (log_trans) {
    # transformiere die log-transformierten Daten wieder auf die ursprüngliche
    # Skala zurück
    col_names <- colnames(df)
    index <- !col_names %in% "Jahr"
    df[, index] <- 10 ^ df[, index]
    # transformiere die Vorhersagen wieder auf die ursprüngliche Skala zurück
    df$prediction <- 10 ^ as.numeric(stats::predict(mod))
  } else {
    df$prediction <- as.numeric(stats::predict(mod))
  }
    target <- names(df)[4]
  ggplot2::ggplot(df, ggplot2::aes(x = Jahr)) +
    ggplot2::geom_point(ggplot2::aes_string(y = target)) +
    ggplot2::geom_line(ggplot2::aes_string(y = target), alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = prediction), color = "darkgreen", size = 1.3) +
    ggplot2::ggtitle(substitute(paste(trend, ", Phi = (", p1,", ", p2, "), Modell = ", modell, ", ", R^2, " = ", rsquared),
                       list(trend = trend[[4]], p1 = phi[1], p2 = phi[2], modell = class(mod)[1], rsquared = formatC(max(0, rsq), format = "f", digits = 4)))) +
    ggplot2::scale_x_continuous(breaks = round(seq(min(df$Jahr), max(df$Jahr), by = 1), 0)) +
    ggplot2::ylab(gsub("(\\_)+", "\\ ", target)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 24, vjust = 1.5),
          axis.title.x = ggplot2::element_text(size = 22, vjust = -.3),
          axis.title.y = ggplot2::element_text(size = 22, vjust = .3),
          axis.text.x = ggplot2::element_text(size = 18, angle = 90, vjust = 0.5),
          axis.text.y = ggplot2::element_text(size = 18))
}

#' Bestimme einflussreiche Beobachtungen nach Cook's Distance
#'
#' Alle Beobachtungen, deren Cook's distance größer als der Schwellenwert ist,
#' werden angegeben.
#'
#' @param mod lm oder gls Modell
#' @param threshold threshold at which influential observations are observed
#' @return Liste mit dem Index der Beobachtungen und den Beobachtungen selbst
#' @param ... variable number of further arguments
infl_obs <- function(mod, threshold = 0.5, ...) {
  UseMethod("infl_obs")
}

#' @param mod lm
#' @param threshold threshold at which influential observations are observed
#' @param data data to be fitted
#' @param ... variable number of further arguments
infl_obs.lm <- function(mod, threshold = 0.5, data, ...) {
  # Cook's distance berechnen
  cooks_distance <- stats::cooks.distance(mod)
  # Index der auffälligen Beobachtungen bestimmen
  index <- which(cooks_distance > threshold)
  attributes(index) <- NULL
  cooks_distance_selected <- cooks_distance[index]
  # die Werte der auffälligen Beobachtungen bestimmen
  values <- mod$model[index, 1]
  # das Jahr der auffälligen Beobachtung bestimmen
  jahr <- data[index, "Jahr"]

  list(einfl_beob_index = index, einfl_beob_value = values,
       einfl_beob_cookd = cooks_distance_selected,
       einfl_beob_jahr = jahr)
}

#' @param mod gls
#' @param used_formula Model formula which was used
#' @param cor_structure correlation structure which was used
#' @param threshold threshold at which influential observations are considered influential
#' @param data data used
#' @param varname colnames of observed variables
infl_obs.gls <- function(mod, used_formula, cor_structure,
                                             threshold = 0.5, data, varname) {
  # calculate cook's distance
  cooks_distance <- CookD_gls(mod, used_formula = used_formula,
                              cor_structure = cor_structure, data = data,
                              plot = FALSE)
  # Index der auffälligen Beobachtungen bestimmen
  index <- which(cooks_distance > threshold)
  attributes(index) <- NULL
  cooks_distance_selected <- cooks_distance[index]
  # erhalte Daten von genutztem Modell
  model_matrix <- data
  # die Werte der auffälligen Beobachtungen bestimmen
  values <- model_matrix[index, varname]
  # das Jahr der auffälligen Beobachtung bestimmen
  jahr <- data[index, "Jahr"]

  list(einfl_beob_index = index, einfl_beob_value = values,
       einfl_beob_cookd = cooks_distance_selected,
       einfl_beob_jahr = jahr)
}

#'
#' Copy of predictmeans::CookD, but with ability to include the used formula
#' and correlation structure so that the update function works. The model needs
#' to be fitted with ML.
#'
#' @inheritParams predictmeans::CookD
#' @param used_formula formula used when fitting the model
#' @param cor_structure correlation structure used when fitting the
#' @param data data used for the fit
#'
#' @details Rather quick and dirty method; not exported function from package
#' predictmeans is used. Maybe better implementation if scoping is adapted.
#'
#' @return vector with Cook's distance
CookD_gls <- function (model, group = NULL, plot = TRUE, idn = 3, newwd = TRUE,
                       used_formula, cor_structure, data)
{
  if (class(model)[1] != "gls") stop("The model has to be of class 'gls'")
  if (model$method != "ML") stop("The model has to be fitted with method 'ML'")

  mdf <- data

  mp <- predictmeans:::mymodelparm.default(model)
  beta0 <- mp$coef
  vcovb <- mp$vcov
  vb.inv <- solve(vcovb)
  if (is.null(group) || group %in% c("NULL", "")) {
    rn <- rownames(mdf)
    LOOmp <- lapply(rn, function(x) predictmeans:::mymodelparm.default(nlme::gls(model = used_formula,
                                                       data = mdf[rn != x, ],
                                                       correlation = cor_structure,
                                                       method = "ML")))
  }
  else {
    rn <- unique(mdf[, group])
    LOOmp <- lapply(rn, function(x) {
      rind <- mdf[, group] != x
      predictmeans:::mymodelparm.default(stats::update(model, data = mdf[rind, ]))
    })
  }
  LOObeta <- sapply(LOOmp, function(x) x$coef)
  rK <- t(LOObeta - beta0)
  CookD <- diag(rK %*% tcrossprod(vb.inv, rK)/length(beta0))
  names(CookD) <- rn
  if (plot) {
    if (newwd)
      grDevices::dev.new()
    outD <- CookD >= sort(CookD, decreasing = TRUE)[idn]
    labid <- names(CookD)
    plot(CookD, xlab = "Obs. number", col = "blue", ylim = c(0,
                                                             max(CookD) + 0.005), main = "Cook's Distance", ylab = "Cook's distance",
         type = "h")
    graphics::text((1:length(CookD))[outD], CookD[outD], labid[outD],
         pos = 3)
    graphics::points((1:length(CookD))[outD], CookD[outD], pch = 16,
           col = "blue")
  }
  return(invisible(CookD))
}
