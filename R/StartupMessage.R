################################################################################
# StartupMessage.R
################################################################################
# 2018-09-20
# Curtis Miller
################################################################################
# Package startup message functions, for fancy loading.
################################################################################

#' Create Package Startup Message
#'
#' Makes package startup message.
#'
#' @import utils
#' @examples
#' CPAT_startup_message()
CPAT_startup_message <- function() {
  c(paste0("      ________ _________ _________ __________\n     /       //  __",
           "_   //  ___   //         /\n    /   ____//  /  /  //  /  /  //___ ",
           "  ___/\n   /   /    /  /__/  //  /__/  /    /  /\n  /   /___ /  __",
           "____//  ___   /    /  /\n /       //  /      /  /  /  /    /  /\n/",
           "_______//__/      /__/  /__/    /__/        v. ",
           utils::packageVersion("CPAT")),
     "\nType citation(\"CPAT\") for citing this R package in publications")
}

#' Package Attach Hook Function
#'
#' Hook triggered when package attached
#'
#' @param lib a character string giving the library directory where the package
#'            defining the namespace was found
#' @param pkg a character string giving the name of the package 
#' @examples
#' .onAttach(.libPaths()[1], "CPAT")
.onAttach <- function(lib, pkg) {
  msg <- CPAT_startup_message()
  if (!interactive())
    msg[1] <- paste("Package 'CPAT' version", packageVersion("CPAT"))
  packageStartupMessage(msg)
  invisible()
}
