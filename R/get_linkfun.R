#' get link function, whether locations are lonlat and space time
#'
#' @param covfun_name string name of covariance function
get_linkfun <- function(covfun_name){
  
  link <- exp
  dlink <- exp
  invlink <- log
  lonlat <- FALSE
  space_time <- FALSE
  return(list(
    link = link, dlink = dlink, invlink = invlink,
    lonlat = lonlat, space_time = space_time
  ))
}

