## Source an R-code from GitHub
## This code was obtained from http://stackoverflow.com/questions/8229859/sourcing-an-r-script-from-github-for-global-session-use-from-within-a-wrapper
## Last update: 09/02/2016

source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
} 
