.welcome <- function(text = NULL)
{
  if(is.null(text))
    text <- "Welcome to AgroR!"
  if(!inherits(text, "character") || length(text) != 1)
    stop("'text' must be a character vector of length 1!")
  vec <- strsplit(text, "")[[1]]
  lab <- c(vec, "\n")
  for(i in 1:length(lab)) {
    setTxtProgressBar(txtProgressBar(char = lab[i]), 0.01)
    Sys.sleep(0.01)
  }
}
.onAttach <- function(lib, pkg)
{
  vers <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage(.welcome(paste("------\nWelcome to AgroR! \n\nShimizu, G.D.; Marubayashi, R.Y.P.; Goncalves, L.S.A. (2021). Package AgroR version", vers)))
}
