
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package EORF, developed by James Thorson for the Alaska Fisheries Science Center")
  packageStartupMessage("For details and citation guidance, please see http://github.com/james-thorson/EOFR/")
  packageStartupMessage("###########################################################################################")
  if( getOption("repos")["CRAN"] == "@CRAN@" ){
    options(repos = c("CRAN" = "http://cran.us.r-project.org"))
  }
  if( !"INLA" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: INLA...")
    #utils::install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
    utils::install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  }

  # Load `FishStatsUtils` via .onAttach because importFrom wasn't working
  # Also requries moving FishStatsUtils to SUGGESTS, so that it doesn't isntall main branch
  if( !"FishStatsUtils" %in% utils::installed.packages()[,1] || utils::packageVersion("FishStatsUtils") < numeric_version("1.1.0") ){
    packageStartupMessage("Updating package FishStatsUtils because previously using version < 1.1.0")
    devtools::install_github("james-thorson/FishStatsUtils", ref="development")
  }
  packageStartupMessage( "Loading package `FishStatsUtils`" )
  library(FishStatsUtils)
}

#' Copy of VAST::make_model
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?VAST::make_model} to see list of arguments and usage
#' @export
Build_TMB_Fn = function( ... ){
  .Deprecated( new="EOFR::make_model" )
  EOFR::make_model( ... )
}

#' Copy of VAST::make_data
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?VAST::make_data} to see list of arguments and usage
#' @export
Data_Fn = function( ... ){
  .Deprecated( new="EOFR::make_data" )
  EOFR::make_data( ... )
}

