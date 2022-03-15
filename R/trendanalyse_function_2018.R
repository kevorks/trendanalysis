#' Function to calculate trend analyses for Bosch & Partner
#'
#' @param data_dir directory where the data is stored
#' @param plot_dir directory where the generated plots are saved
#' @param result_name directory and name (without ending) where the result table
#' is stored
#' @param log_trans logical if log-transformed values are used (TRUE leads to
#' the adjustments for plot function)
#'
#' @return the function does not return a value but stores the plot and table in
#' the assigned directories
#' @export
trenda <- function(data_dir, plot_dir, result_name,
                                 log_trans = FALSE) {

  ### Benoetigt werden die Funktion in der Datei functions_quadreg_2018.R
  source("R/functions_quadreg_2018.R")

  ## set directory
  trend_files <- list.files(data_dir)

  ResTab <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(ResTab) <- c("File", "Index", "Beta0", "Beta1", "Beta2", "Phi1", "Phi2",
                        "Trend", "rQuadrat", "zu_entfernende_Beobachtungen_Index",
                        "zu_entfernende_Beobachtungen_Werte",
                        "Cooks_Distance_der_Beobachtungen",
                        "zu_entfernende_Beobachtungen_Jahr"
  )


  # set seed for durbin-watson-test
  set.seed(1)

  z <- 1

  output <- lapply(trend_files, function(trend_file) {

    print(sprintf("Data:----------%s", trend_file)) # Exceltitle

    ## Daten einlesen
    dat <- utils::read.csv2(sprintf("%s%s", data_dir, trend_file), header = TRUE, sep = ";")

    # zeige Datenstruktur zur Überprüfung an
    print(utils::str(dat))
    print(all(apply(dat, 2, is.numeric)))

    ## Wegen plot-Funktion muessen die Punkte in den Variablennamen ersetzt werden
    names(dat) <- gsub("\\.", "\\_", names(dat))

    ## Alle Variablen ausser ID, Zeit und Jahr sind Umweltindikatoren
    varnames <- setdiff(names(dat), c("ID", "Zeit", "Jahr"))

    ## ID notwendig fuer Korrelationsstruktur fuer gls-Funktion von
    # nlme package (alle Daten stammen jeweils von einer Beobachtungseinheit)
    dat$ID <- 1

    ## Zeit faengt bei 0 an
    dat$Zeit <- dat$Jahr - min(dat$Jahr)

    for(varname in varnames){
      print(varname)

      ## Trend berechnen
      res <- compute_trend(dat, varname, log_trans = log_trans)

      ## Ergebnisse speichern
      ResTab[z,"File"] <<- trend_file
      ResTab[z,"Index"] <<- varname
      ResTab[z,"Beta0"] <<- round(as.numeric(res$beta[1]),4)
      ResTab[z,"Beta1"] <<- round(as.numeric(res$beta[2]),4)
      ResTab[z,"Beta2"] <<- round(as.numeric(res$beta[3]),4)
      ResTab[z,"Phi1"] <<- res$phi[1]
      ResTab[z,"Phi2"] <<- res$phi[2]
      ResTab[z,"Trend"] <<- res$trend
      ResTab[z,"rQuadrat"] <<- res$rsq
      ResTab[z, "zu_entfernende_Beobachtungen_Index"] <<- res$einfl_beob_index
      ResTab[z, "zu_entfernende_Beobachtungen_Werte"] <<- res$einfl_beob_value
      ResTab[z, "Cooks_Distance_der_Beobachtungen"] <<- res$einfl_beob_cookd
      ResTab[z, "zu_entfernende_Beobachtungen_Jahr"] <<- res$einfl_beob_jahr

      ## Trendgrafik speichern
      if (log_trans) {
        name_file <- substr(trend_file, start = 6, stop = nchar(trend_file) - 4)
      } else {
        name_file <- abbreviate(trend_file, 8)
      }
      grDevices::jpeg(file = sprintf("%s%s_%s.jpeg", plot_dir, name_file, varname),
           width = 800, height = 800, quality = 640000)
      plot(res$plot)
      grDevices::dev.off()

      z <<- z + 1

    }

  })

  summary(ResTab)

  utils::write.table(ResTab,
              paste0(result_name, format(Sys.Date(), "%Y-%m-%d"), ".csv"),
              sep = ";", dec = ",", row.names = FALSE)

  # pruefen, ob alle Files benutzt wurden
  trend_files == unique(ResTab$File)

}
