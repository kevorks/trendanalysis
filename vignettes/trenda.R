## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,  comment = "#>")
options(tibble.print_min = 6L, tibble.print_max = 6L, digits = 3)
library(trenda)

## ---- out.width="100%", echo=FALSE--------------------------------------------
knitr::include_graphics("diagram.JPG")

## -----------------------------------------------------------------------------

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  # load package
#  library(trenda)

## ---- echo = TRUE-------------------------------------------------------------
set.seed(1)
rnd_preci <- data.frame(Year = c(1991:2020),
                        precip_mm = rnorm(30, 770, 50),
                        height_m = c(rnorm(10, 10, 1),
                                     100,
                                     rnorm(9, 15, 1),
                                     rnorm(10, 20, 1)))

## -----------------------------------------------------------------------------
str(rnd_preci)

## ---- echo = TRUE-------------------------------------------------------------

trenda(rnd_preci, calc_infl_obs = TRUE)

## ---- out.width="70%", echo = FALSE-------------------------------------------
knitr::include_graphics("precip_mm.jpeg")
  
#Figure2 precip_mm
knitr::include_graphics("height_m.jpeg")
  
#Figure3 height_m
knitr::include_graphics("height_m_infl.jpeg")



## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  list.files(path = "./files_folder")

## ---- eval = FALSE, echo=TRUE-------------------------------------------------
#  trenda_dir("./files_folder/", log_trans = FALSE, create_dir = TRUE, calc_infl_obs = TRUE)

