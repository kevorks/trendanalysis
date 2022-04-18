##' Plot the trend model
##'
##' @description Plots data as points and trends as lines
##' @param mod lm or gls model
##' @param df data
##' @param log_trans TRUE / FALSE for log-transformed data
##' @import ggplot2
##' @importFrom stats predict
##' @return plot - object
plot_trend <- function(mod, df, log_trans = FALSE) {
  rsq <- round(1 - (sum(mod$residuals ^ 2) / sum(((mod$residuals + mod$fitted) - mean(mod$residuals + mod$fitted)) ^ 2)), 4)

  trend <- extract_trend(mod)
  phi <- extract_phi(mod)
  if (log_trans) {
    # transform log-transformed data back original scale
    col_names <- colnames(df)
    index <- !col_names %in% "Jahr"
    df[, index] <- 10 ^ df[, index]
    # transforms the previous predictions back to the original scale
    df$prediction <- 10 ^ as.numeric(predict(mod))
  } else {
    df$prediction <- as.numeric(predict(mod))
  }
  target <- names(df)[4]
  ggplot(df, aes(x = Jahr)) +
    geom_point(aes_string(y = target)) +
    geom_line(aes_string(y = target), alpha = 0.3) +
    geom_line(aes(y = prediction), color = "darkgreen", size = 1.3) +
    ggtitle(substitute(paste(trend, ", Phi = (", p1,", ", p2, "), Modell = ", modell, ", ", R^2, " = ", rsquared),
                       list(trend = trend[[4]], p1 = phi[1], p2 = phi[2], modell = class(mod)[1], rsquared = formatC(max(0, rsq), format = "f", digits = 4)))) +
    scale_x_continuous(breaks = round(seq(min(df$Jahr), max(df$Jahr), by = 1), 0)) +
    ylab(gsub("(\\_)+", "\\ ", target)) +
    theme(plot.title = element_text(size = 24, vjust = 1.5),
          axis.title.x = element_text(size = 22, vjust = -.3),
          axis.title.y = element_text(size = 22, vjust = .3),
          axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 18))
}
