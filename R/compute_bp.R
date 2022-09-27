##' Wrapper for compute_bp
##'
##' @description Takes all observations from 2:n (since first column is taken as x-achis)
##'  of a data frame and computes the breakpoints for each observation
##' @param dat Data frame to be analysed
##' @param type Currently set to BIC only
##' @param control Argument used in segmented. Defaults to
##' seg.control(n.boot = 50, it.max = 1000, alpha = 0.1)
##' The alpha-value is set quit high to avoid finding break points near the edges
##' See also ?seg.control for more parameters
##' @param psi Psi can be a numeric vector of guesses on breakpoints passed on to the segmented function
##' Segmented takes the psi values as starting points to find segments in the data. Defaults to NULL
##' @param npsi A named vector or list meaning the number (and not locations) of breakpoints to be estimated.
##' The starting values will be internally computed via the quantiles or equally spaced values,
##' as specified in argument quant in seg.control. npsi can be missing and npsi=1
##' is assumed for all variables specified in seg.Z. If psi is provided, npsi is ignored.
##' @param msg logical. Indicates if results should be displayed
##' @param plot Plots the graph
##' @import dplyr
##' @import plotly
##' @return Returns a list containing the information
##' @export
trenda_bp <- function(dat, type = "BIC",
                      control = seg.control(n.boot = 50, it.max = 1000, alpha = 0.1),
                      psi = NULL,
                      npsi = 2,
                      msg = FALSE,
                      plot = TRUE,
                      return.fit = FALSE) {
  if (nrow(dat) < 30) {
    stop("Less than 30 observations")
  }
  if (length(psi) > 3 | npsi > 3) {
    warning("Psi and npsi should not be greater than 3 or psi longer than 3 elements")
  }
  varnames <- names(dat)
  varnames <- setdiff(names(dat), names(dat[1]))
  biglist <- list()
  for (varname in varnames) {
    #print(varname)
    if (return.fit) {
    t <- compute_bp(dat,
               varname,
               type = type,
               control = control,
               psi = psi,
               npsi = npsi,
               msg = msg,
               plot = plot,
               return.fit = return.fit)
    biglist[[varname]] <- t

    } else {
      compute_bp(dat,
                 varname,
                 type = type,
                 control = control,
                 psi = psi,
                 npsi = npsi,
                 msg = msg,
                 plot = plot,
                 return.fit = return.fit)
      }

  }
  if (length(biglist) > 0 | return.fit == TRUE) {
  return(biglist)
  }
}


##' Search for breakpoints using the segmented package. Picking the best model with BIC
##'
##' @description The package segmented is used to analyse data for break points.
##' Breakpoints are chosen by least BIC and are then, if present, are plotted in
##' a comprehensive graph using plotly
##' @param dat Data frame to be analysed
##' @param type Currently set to BIC only
##' @param control Argument used in segmented. Defaults to
##' seg.control(n.boot = 50, it.max = 1000, alpha = 0.1)
##' The alpha-value is set quit high to avoid finding break points near the edges
##' See also ?seg.control for more parameters
##' @param psi Psi can be a numeric vector of guesses on breakpoints passed on to the segmented function
##' Segmented takes the psi values as starting points to find segments in the data. Defaults to NULL
##' @param npsi A named vector or list meaning the number (and not locations) of breakpoints to be estimated.
##' The starting values will be internally computed via the quantiles or equally spaced values,
##' as specified in argument quant in seg.control. npsi can be missing and npsi=1
##' is assumed for all variables specified in seg.Z. If psi is provided, npsi is ignored.
##' @param msg logical. Indicates if results are displayed
##' @param plot Plots a graph using plotly
##' @details The package segmented provides a tool to determine segmented relationships in regression models
##' with breakpoints / changepoints (with possibly random effects) estimation
##' @import plotly
##' @import dplyr
##' @import segmented
##' @export
##'


compute_bp <- function(dat, varname, type = "BIC",
                        control = seg.control(n.boot = 50, it.max = 1000, alpha = 0.05),
                        psi = NULL,
                        npsi = 2,
                        msg = TRUE,
                        plot = TRUE,
                        return.fit = FALSE) {
  if (nrow(dat) < 30) {
    stop("Less than 30 observations")
  }
  if (length(psi) > 3 | npsi > 3) {
    warning("Psi and npsi should not be greater than 3")
  }

  mod_fit <- as.formula(sprintf("%s ~ %s", varname, names(dat[1])))
  mod_form <- lm(mod_fit, data = dat)

  nome <- all.vars(formula(mod_form))[[2]]
  seg.z <- as.formula(paste("~", nome))

  if (missing(type)) {
    type <- "BIC"
  }

  type <- match.arg(type)
  if(type %in% c("BIC")) {
    BIC.f <- if(type == "BIC")
      BIC
  }

  resi <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(resi) <- c("npsi", "Breakpoint", "BIC", "usedBP")

  if (is.null(psi) & npsi > 0) {
    npsi <- 1:npsi
    ris <- vector("list", length(npsi))
    BIC.values <- rep(NA, length(npsi))

    for (i in npsi) {
      ris[[i]] <- suppressWarnings(try(segmented(mod_form,
                                                 seg.z, npsi = i, control = control, silent = TRUE)))
      if (!inherits(ris[[i]], "segmented"))
        ris[[i]] <- suppressWarnings(try(segmented(mod_form,
                                                   seg.z, npsi = i, control = control, silent = TRUE)))
      if (inherits(ris[[i]], "segmented"))
        BIC.values[i] <- BIC.f(ris[[i]])
    }
    BIC.mod <- BIC.values
    BIC.values <- c(BIC.f(mod_form), BIC.mod)
    names(BIC.values) <- c("0", npsi)
    n.psi.ok <- as.numeric(names(BIC.values)[which.min(BIC.values)])

    if (n.psi.ok > 0) {
      model.ok <- ris[[n.psi.ok]]
      cuts <- data.frame(value = round(ris[[n.psi.ok]]$psi[,2]), BP = 1:n.psi.ok)
      row.names(cuts) <- 1:nrow(cuts)

      fit <- numeric(length(dat[varname])) * NA
      suppressWarnings(fit[complete.cases(rowSums(cbind(dat[varname], dat[1])))] <- broken.line(model.ok))

      dat1 <- data.frame(jahr = dat[1], varname = dat[varname], fit = fit)
      colnames(dat1) <- c("jahr", varname, "fit")

    } else if (n.psi.ok == 0) {
      model.ok <- mod_form
      cuts <- NULL
    }

    r <- list(BIC.values = BIC.values, n.psi.ok = n.psi.ok, breakpoints = cuts)


    res <- ris
    ris <- c(list(mod_form), ris)

    if (which.min(BIC.values) == 1) {
      n.psi.ok <- "0"
    } else {
      n.psi.ok <- which.min(BIC.values[2:length(BIC.values)])
    }
    #print(ris)
    # psi is predefined and passed on to segmented by a vector. this suggests a segmentation at the specific position
    # and is checked by segmented
    if (msg) {
      if(type == "BIC")
        cat("Auto-search BP with best BIC up to ", max(npsi), "breakpoints \n")
      cat("Variable: ", varname, "\n")
      cat("BIC values:\n")
      print(BIC.values)
      cat("Estimated Breakpoint(s)")
      print(cuts)
      #cat(paste("Lowest BIC among segmented: ", BIC.mod[which.min(BIC.mod)], "\n"))
      if (n.psi.ok > 0) {
        for (i in n.psi.ok) {
          #cat(paste("Estimated break-point(s): \n"))
          resi <- round(res[[i]]$psi[,2])
        }
      } else {
        cat(paste("No breakpoints found. Lowest BIC is lm \n"))
      }
    }
    if (plot) {
      ploti(dat = dat, varname = varname, cuts = cuts, n.psi.ok = n.psi.ok, model.ok = model.ok, dat1 = dat1)

    }
    # If npsi is empty but psi is given
  } else if (!is.null(psi)) {
    npsi <- NULL
    ris <- segmented(mod_form, control = control, psi = psi)
    BIC.values <- c(BIC(mod_form), BIC(ris))
    names(BIC.values) <- c("0", "1")
    n.psi.ok <- as.numeric(names(BIC.values)[which.min(BIC.values)])

    fit <- numeric(length(dat[varname])) * NA
    suppressWarnings(fit[complete.cases(rowSums(cbind(dat[varname], dat[1])))] <- broken.line(ris))

    dat1 <- data.frame(jahr = dat[1], varname = dat[varname], fit = fit)
    colnames(dat1) <- c("jahr", varname, "fit")

    #if (n.psi.ok > 0) {
    le <- as.numeric(length(psi))
    cuts <- data.frame(value = round(ris$psi[,2]), BP = 1:le)
    row.names(cuts) <- 1:nrow(cuts)
    #cuts <- data.frame(value = round(ris[[]]))
    #}
    if(msg) {
      cat(paste("Using selected psi", length(psi), "as starting point \n"))
      cat(paste("Variable: ", varname, "\n"))
      #print(cuts)
      cat("BIC Values: \n")
      print(BIC.values)
      if(BIC(mod_form) < BIC(ris)) {
        cat("No breakpoints found. Lowest BIC is lm \n")
      } else {
        cat("Estimated breakpoint(s): \n")
        print(cuts)
      }
      cat("_____________________________ \n")
    }
    if (plot) {
      #n.psi.ok = length(psi)
      ploti(dat, varname = varname, cuts = cuts, n.psi.ok = n.psi.ok, model.ok = model.ok, dat1 = dat1)
    }
  }
  if (return.fit) {
    return(r)
  }
}

##' Plots the function trenda_bp() and compute_bp()
##' @description Uses ggplot2 and plotly to plot the results from trenda_bp()
##' @param dat data frame with observations
##' @param varname variable name within the data frame dat
##' @param cuts breakpoints determined by segmented in trenda_bp()
##' @param n.psi.ok numeric value indicationg the best breakpoints according to BIC
##' @param model.ok segmented model
##' @param dat1 x and y values including the fitted segmented lines for plotting
##' @import plotly
##' @import ggplot2
##' @return interactive plot using plotly and ggplot2

# Function for plotting
ploti <- function(dat, varname, cuts, n.psi.ok, model.ok, dat1) {
  if (n.psi.ok > 0) {


  p <- ggplot(dat, aes_string(x = names(dat[1]), y = varname)) +
    geom_line() +
    theme_bw() + ylab(varname) + xlab(names(dat[1])) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    geom_vline(mapping = aes(xintercept = value),
               data = cuts,
               show.legend = FALSE, color = "red",
               linetype = "dotted") +
    geom_text(mapping = aes(x = value + 3,
                            y = max(dat[varname]) + 0.01*max(dat[varname]),
                            label = value),
              data = cuts) +
    geom_point(alpha = 0.5) +
    labs(title = paste(varname, "|", "Anz. BP: ", n.psi.ok)) +
    #geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
    #            data = dat) +
    geom_line(data = dat1, aes_string(x = names(dat1[1]), y = "fit"), color = "blue")
  print(ggplotly(p) %>% config(displayModeBar = FALSE))
  #plot.segmented(model.ok, add = TRUE)
  #p
  } else if (n.psi.ok == 0) {
    p <- ggplot(dat, aes_string(x = names(dat[1]), y = varname)) +
      geom_line() + theme_bw() + ylab(varname) + xlab(names(dat[1])) +
      #ggtitle(paste(varname, "Test", "stuff")) +
      suppressWarnings(geom_smooth(method = "loess", formula = y ~ x,  se = FALSE,
                 data = dat)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      geom_point(alpha = 0.5) +
      labs(title = paste("Varname: ", varname,"|", "Keine BP: ", n.psi.ok))
    #p <- p + labs(subtitle = "test sub")

    print(ggplotly(p) %>% config(displayModeBar = FALSE))
    #p + plot.segmented(model.ok, add = TRUE)
  }
}


