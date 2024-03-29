##' Plot the trend model
##' @description Plots data as points and trends as lines
##' @param mod lm or gls model
##' @param df data
##' @param log_trans TRUE / FALSE for log-transformed data. If data was log transformed
##' in a preparation process beforehand
##' @import ggplot2
##' @importFrom stats predict
##' @return plot - object
plot_trend <- function(mod, df, log_trans = FALSE) {
  #return(df[2])
  f_r <- names(df[1])
  f_rn <- df[,1]
  divi <- round(nrow(df)*0.1)

  rsq <- round(1 - (sum(mod$residuals ^ 2) / sum(((mod$residuals + mod$fitted) - mean(mod$residuals + mod$fitted)) ^ 2)), 4)
  trend <- extract_trend(mod)
  phi <- extract_phi(mod)
  df[, "predi"] <- NULL
  if (log_trans) {
    # transform log-transformed data back to original scale
    col_names <- colnames(df)
    index <- !col_names %in% f_r
    df[, index] <- 10 ^ df[, index]
    # transforms the previous predictions back to the original scale

    df[, "predi"] <- 10 ^ as.numeric(stats::predict(mod))
  } else {
    df[, "predi"] <- as.numeric(stats::predict(mod))
  }
  target <- names(df)[4]
  ggplot(df, aes_string(x = f_rn)) +
    xlab(f_r) +
    geom_point(aes_string(y = target)) +
    geom_line(aes_string(y = target), alpha = 0.3) +
    geom_line(aes_string(y = "predi"), color = "darkgreen", size = 1.3) +
    ggtitle(substitute(paste(trend, ", Phi = (", p1,", ", p2, "), Model = ", model, ", ", R^2, " = ", rsquared),
                       list(trend = trend[[4]], p1 = phi[1], p2 = phi[2], model = class(mod)[1], rsquared = formatC(max(0, rsq), format = "f", digits = 4)))) +
    #scale_x_continuous(limits = c(min(df[f_r]), max(df[f_r])), n.breaks = 24) +
    #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    ylab(gsub("(\\_)+", "\\ ", target)) +
    theme(plot.title = element_text(size = 24, vjust = 1.5),
          axis.title.x = element_text(size = 22, vjust = -.3),
          axis.title.y = element_text(size = 22, vjust = .3),
          axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 18))
}

