#' Function to calculate trend analyses (without influential observations)
#'
#' @description Second main function of package trenda. Should only be executed
#' after the function trenda was performed on data and influential observations
#' were determined. trenda_obs removes those influential observations and runs
#' the analysis again
#'
#' @param data_dir directory where the data is stored
#' @param log_trans logical if log-transformed values were used (TRUE leads to
#' the adjustments for plot function)
#' @param res_tab_file result table of the same data generated by trenda;
#' needed for information about the influential observations which should are
#' left out in this analysis
#' @param calc_infl_obs logical; the influential observations should be
#' calculated again. Can be useful to turn off as this information is not
#' needed and sometimes the calculation of the influential observations for
#' GLS models does not work
#'
#' @return the function does not return a value but stores the plot and table in
#' the assigned directories
#'
#' @details Since this function needs information about influential
#' observations, the function trenda has to be run first on the same
#' dataset.
#' @importFrom grDevices dev.off jpeg
#' @export trenda_obs
trenda_obs <- function(data_dir, log_trans = FALSE, res_tab_file,
                       calc_infl_obs = TRUE) {
  if (!log_trans) {
    dir.create(paste0(data_dir, Sys.Date(), "results_standard_obs"))
    plot_dir <- paste0(data_dir, "results_standard_obs/")
    result_name <- paste0(plot_dir, "/result_rable_standard_obs_")
  } else {
    dir.create(paste0(data_dir, "results_log_obs"))
    plot_dir <- paste0(data_dir,"results_log_obs/")
    result_name <- paste0(plot_dir, "/Result_Table_log_obs")
  }
  trend_files <- list.files(data_dir)

  ResTab <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(ResTab) <- c("File", "Index", "Beta0", "Beta1", "Beta2", "Phi1", "Phi2",
                        "Trend", "rSquared", "observations_to_remove_index",
                        "observations_to_remove_value",
                        "cooks_distance_of_observation",
                        "observations_to_remove_year"
  )

  # Import the table generated by the trenda function to exclude influential observations
  result_infl <- read.csv2(res_tab_file,
                           stringsAsFactors = FALSE)

  # since generalized Durbin-Watson-Test is random -> set seed
  set.seed(1)

  z <- 1

  output <- lapply(trend_files, function(trend_file) {

    print(sprintf("Data:----------%s", trend_file)) # Exceltitle

    ## read dataset
    dat <- read.csv2(sprintf("%s%s", data_dir, trend_file), header = TRUE)
    print(str(dat))

    ## dots have to be replaced in the variablenames for the plot function
    names(dat) <- gsub("\\.", "\\_", names(dat))

    ## all variables except for ID, Time and Jahr are environmental indicators
    varnames <- setdiff(names(dat), c("ID", "Time", "Jahr"))

    ## ID is needed for the correlationstructure in the gls-function
    dat$ID <- 1

    ## Time begins at 0
    dat$Time <- dat$Jahr - min(dat$Jahr)

    for(varname in varnames){
      print(varname)

      # Remove observations which have been identified as influential in the previous
      # analysis
      # indices are saved as string, decimals are seperated by comma
      index <- result_infl[result_infl$File == trend_file &
                             result_infl$Index == varname,
                           "observations_to_remove_index"]
      if (!is.na(index)) {
        index <- as.character(index)
        index <- strsplit(index, ",")
        index <- as.vector(sapply(index, as.numeric))
      }

      if (!is.list(index) && !is.na(index)) {
        # remove rows with NA since the index referes to influential observations
        temp_dat <- na.omit(dat[, c("Time", "Jahr", "ID", varname)])
        temp_dat <- temp_dat[-index, ]

        print(str(temp_dat))

        ## calculate trend
        res <- compute_trend(temp_dat, varname, log_trans = log_trans,
                             calc_infl_obs = calc_infl_obs)

        ## save results
        ResTab[z,"File"] <<- trend_file
        ResTab[z,"Index"] <<- varname
        ResTab[z,"Beta0"] <<- round(as.numeric(res$beta[1]),4)
        ResTab[z,"Beta1"] <<- round(as.numeric(res$beta[2]),4)
        ResTab[z,"Beta2"] <<- round(as.numeric(res$beta[3]),4)
        ResTab[z,"Phi1"] <<- res$phi[1]
        ResTab[z,"Phi2"] <<- res$phi[2]
        ResTab[z,"Trend"] <<- res$trend
        ResTab[z,"rSquared"] <<- res$rsq

        ## save trend graphs
        if (log_trans) {
          name_file <- substr(trend_file, start = 6, stop = nchar(trend_file) - 4)
        } else {
          name_file <- abbreviate(trend_file, 8)
        }
        jpeg(file = sprintf("%s%s_%s.jpeg", plot_dir,
                            name_file,
                            varname),
             width = 800, height = 800, quality = 640000)
        plot(res$plot)
        dev.off()
        z <<- z + 1

      }
    }
  })

  summary(ResTab)

  write.table(ResTab,
              paste0(result_name, format(Sys.Date(), "%Y-%m-%d"), ".csv"),
              sep = ";", dec = ",", row.names = FALSE)

}