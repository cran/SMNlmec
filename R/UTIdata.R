
#' @title Data set for Unstructured Treatment Interruption Study
#' @description Data set from a study of Unstructured Treatment Interruption in HIV-infected adolescents in four institutions in the US. The main outcome is the HIV-1 RNA viral load, which is subject to censoring below the lower limit of detection of the assay (50 copies/mL). The censored observations are indicated by the variable RNAcens.
#' @usage data(UTIdata)
#' @format
#' \describe{
#'  A data frame with 373 observations on the following 5 variables.
#'  \item{Patid}{patient ID}
#'  \item{Days.after.TI}{days after treatment interruption.}
#'  \item{Fup}{follow-up months}
#'  \item{RNA}{viral load RNA}
#'  \item{RNAcens}{censoring indicator for viral load}
#' }
#'
#' @references Saitoh, A., Foca, M, et al. (2008), Clinical outcome in perinatally acquired HIV-infected children and adolescents after unstructured treatment interruption, Pediatrics,121, e513-e521.
#' @docType data
#' @name UTIdata
NULL


