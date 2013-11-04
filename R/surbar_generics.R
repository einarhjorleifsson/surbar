## Function that could go to generic packages, like fishvice

#' Shorthand for as.character
#' 
#' @export
#' @param x A vector
ch <- function(x){as.character(x)}


#' Read in Lowestoft VPA data
#' 
#' @export
#' @param filename The name of the file to read in
read.vpa.file <- function(filename)
{
  wk.y <- scan(filename, skip = 2, nlines = 1, quiet = TRUE)
  wk.a <- scan(filename, skip = 3, nlines = 1, quiet = TRUE)
  wk.tab <- read.delim(filename, header = FALSE, sep = "", skip = 5)
  list(y = wk.y, a = wk.a, tab = wk.tab)
}


#' Read in Lowestoft-format survey data
#' 
#' @export
#' @param filename The name of the file to read in
read.survey.file <- function(filename)
{
  wk.n <- scan(filename, skip = 1, nlines = 1, quiet = TRUE) - 100
  wk.idx <- vector("list", length = wk.n)
  wk.start <- 3
  for (wk.k in 1:wk.n)
  {
    wk.idx[[wk.k]]$name <- paste(scan(filename, skip = wk.start - 1, nlines = 1, 
                                      what = character(0), quiet = TRUE), collapse = " ")
    wk.temp <- scan(filename, skip = wk.start, nlines = 1, quiet = TRUE)
    wk.idx[[wk.k]]$y1 <- wk.temp[1]
    wk.idx[[wk.k]]$y2 <- wk.temp[2]
    wk.idx[[wk.k]]$ny <- wk.temp[2] - wk.temp[1] + 1
    wk.temp <- scan(filename, skip = wk.start + 1, nlines = 1, quiet = TRUE)
    wk.idx[[wk.k]]$rho <- 0.5 * (wk.temp[4] + wk.temp[3])
    wk.temp <- scan(filename, skip = wk.start + 2, nlines = 1, quiet = TRUE)
    wk.idx[[wk.k]]$a1 <- wk.temp[1]
    wk.idx[[wk.k]]$a2 <- wk.temp[2]
    wk.idx[[wk.k]]$na <- wk.temp[2] - wk.temp[1] + 1
    wk.idx[[wk.k]]$tab <- read.table(filename, skip = wk.start + 3, nrows = wk.idx[[wk.k]]$ny)
    wk.temp <- wk.idx[[wk.k]]$tab[,2:(wk.idx[[wk.k]]$na + 1)] 
    wk.effort <- wk.idx[[wk.k]]$tab[,1] 
    wk.idx[[wk.k]]$tab <- data.frame(wk.temp / wk.effort)
    names(wk.idx[[wk.k]]$tab) <- wk.idx[[wk.k]]$a1:wk.idx[[wk.k]]$a2
    rownames(wk.idx[[wk.k]]$tab) <- wk.idx[[wk.k]]$y1:wk.idx[[wk.k]]$y2
    wk.start <- wk.start + 4 + wk.idx[[wk.k]]$ny
  }
  list(n = wk.n, idx = wk.idx)
}


