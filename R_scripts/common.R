########################################################################################################################
## common.R                                                                                                           ##
## ---------------------                                                                                              ##
## A script with commonly used functions to share them between scripts                                                ##
##                                                                                                                    ##
## A quite safe way to source the file (with respect to varying working directory; avoiding absolute paths to allow   ##
## moves of the files; and do not leaving any defined variables) is following (expecting it is in the same            ##
## directory):                                                                                                        ##
## (function(){                                                                                                       ##
##   paths = unlist(sapply(sys.frames(), function(f) f$ofile))                                                        ##
##   if (length(paths)) {                                                                                             ##
##     path = paths[length(paths)]                                                                                    ##
##   } else {                                                                                                         ##
##     args = commandArgs()                                                                                           ##
##     positions = which(startsWith(args, "--file="))                                                                 ##
##     if (length(positions)) {                                                                                       ##
##       path = substring(args[positions], 8)                                                                         ##
##     } else {                                                                                                       ##
##       path = "."                                                                                                   ##
##     }                                                                                                              ##
##   }                                                                                                                ##
##   path = file.path(dirname(path), "common.R")                                                                      ##
##   if (file.exists(path)) {                                                                                         ##
##     source(path)                                                                                                   ##
##   } else {                                                                                                         ##
##     stop("File 'common.R' is missing in the expected location, please repair the path or restore the file.")       ##
##   }                                                                                                                ##
## })()                                                                                                               ##
##                                                                                                                    ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                                                                 ##
## Last update: 2023-05-17                                                                                            ##
## Released under Apache License 2.0                                                                                  ##
########################################################################################################################


#### Prints help message and quit the script; error is in format c(error_message, args...)
help.common <- function(message, error = c()) {
  if (length(error)) {
    con     = stderr()
    status  = 1
    message = c(paste(error[1], paste0("'", error[-1], "'", collapse=", ")), "", message)
  } else {
    con     = stdout()
    status  = 0
  }
  writeLines(message, con = con)
  quit(status = status)
}


#### Checks whether manual page was requested and either shows it (and quit), or returns command-line arguments
get.arguments <- function(fun = help) {
  # Get arguments
  args = commandArgs(trailingOnly=TRUE)
  # Check whether manual page is requested
  if (length(args) == 0 || (length(args) == 1 && args[1] %in% c("-h","--help"))) {
    help()
  }
  # Return arguments
  return(args)
}


#### Checks whether number of arguments is even
check.parity <- function(args, even=T) {
  if (even && length(args) == 1) {
    help(error = c("Unrecognized argument:", args))
  }
  if (length(args) %% 2 == (if(even) 1 else 0)) {
    help(error = c(paste(if(even) "Even" else "Odd", "number of arguments expected:"), args))
  }
}


#### Parse argument and check whether it is a boolean
parse.boolean <- function(arg, name) {
  ret = as.logical(arg)
  if (is.na(ret)) stop(paste0("Unrecognized value of --",name," parameter, boolean expected: ",arg), call.=F)
  return(ret)
}


#### Check whether an extension is supported
check.extension <- function(extension) {
  if (!any(extension == c("pdf","svg","png"))) {
    stop(paste("Extension '",extension,"' is not supported, please use 'pdf' (default), 'svg', or 'png'"), call.=F)
  }
}


#### Generate an error that lengths of two list does not correspond
# An auxiliary function to format a list printing
pairwise.error.format <- function(name, list) {
  paste0(name, ": ", paste0(list, collapse="; "))
}
# The main function; <message> is an optional argument if the default error message is not suitable
pairwise.error <- function(name1, list1, name2, list2, message) {
  stop(paste(if(missing(message)) {
               paste0("The number of ", name1, " must be the same as the number of ", name2, ":")
             } else {
               message
             },
             pairwise.error.format(name1, list1),
             pairwise.error.format(name2, list2),
             sep="\n"),
       call.=F)
}


#### Checks whether libraries are installed
check.installed <- function(...) {
  # To save on c() when calling this function
  names = list(...)
  # Identify libraries that are not installed
  not.installed = names[!sapply(names, function(name) length(find.package(name, quiet=T)))]
  # If at least one library is not installed, then print the error message and exit
  # For security reasons, missing libraries are not installed within the script
  if (length(not.installed) == 1) {
    stop(paste("Library ", not.installed, " is not installed. Please, install it first - e.g. 'install.packages(",
               not.installed, ")'.", sep='"'),
         call.=F)
  } else if (length(not.installed) > 1) {
    stop(paste0("Libraries ", paste0("'", not.installed, "'", collapse=", "), " are not installed. ",
                "Please, install them first - e.g. 'install.packages(c(\"", paste0(not.installed, collapse='", "'),"\"))'."),
         call.=F)
  }
}
