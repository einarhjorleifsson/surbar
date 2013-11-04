#' surba_fit
#'
#' @export
#' @param dat A list that contains the ...
#' @param mqdt.control mqdt.control
#' @param fitOnly flag
#' @param boot flag
surba_fit <- function(dat,mqdt.control,fitOnly=FALSE,boot=TRUE) 
{
  if(missing(mqdt.control))
  {
    mqdt.control <- list(ftol = 0.00001, ptol = 0.00001, gtol = 0, 
                     diag = numeric(), factor = 100, maxfev = 100 * (length(dat$params0) + 1), 
                     nprint = 0, maxiter = 200)
  }
  dat$fit <- nls.lm(dat$params0, fn = surba_opt, control = mqdt.control,dat=dat)
  if(!fitOnly) {
    dat$par <- surba_par(dat)
    dat$res <- surba_res(dat)
    dat$rby <- surba_rby(dat)
  }
  if(boot) {
    dat$boot <- surba_boot(dat)
  }
  return(dat)
}

#' surba retro
#' 
#' @export
#' @param dat dat
#' @param nYears Number of retroyears
surba_retro <- function(dat,nYears=5)
  {
  x <- dat
  retro <- list()
  for (i in 1:nYears)
    {
    x <- surba_trim(x)
    x <- surba_fit(x,boot=FALSE)
    retro[[i]] <- x$rby
  }
  return(retro)
}


#' surba_opt
#' 
#' @export
#' @param wk.p wk.p
#' @param dat dat
surba_opt <- function(wk.p,dat=NULL)
{

  x <- dat$u.std
  numk <- dat$numk
  na <- length(dat$ages)
  ny <- length(dat$years)
  ref.age <- dat$ref.age
  y1 <- min(dat$years)
  qval <- dat$qval
  rho <- dat$rho
  wt <- dat$wt
  params0 <- dat$params0
  lambda <- dat$lamda
  sw <- dat$sw
  mat <- dat$mat
  a1 <- min(dat$ages)
  a2 <- max(dat$ages)
  y2 <- max(dat$years)
  
  
  wk.nk <- numk
  wk.na <- na
  wk.ny <- ny
  
  wk.s <- rep(NA, length = wk.na)
  wk.s[1:(ref.age-1)] <- wk.p[1:(ref.age-1)]
  wk.s[ref.age] <- 1.0
  wk.s[(ref.age+1):(wk.na-1)] <- wk.p[ref.age:(wk.na-2)]
  wk.s[wk.na] <- wk.s[wk.na-1]
  
  wk.f <- rep(NA, length = wk.ny)
  wk.f[1:(wk.ny-1)] <- wk.p[(wk.na-1):(wk.na+wk.ny-3)]
  wk.f[wk.ny] <- mean(wk.f[(wk.ny-3):(wk.ny-1)])
  
  wk.r <- rep(NA, length = wk.na + wk.ny - 1)
  wk.r[1:(wk.na+wk.ny-1)] <- wk.p[(wk.na+wk.ny-2):length(wk.p)]
  
  # Total mortality
  wk.z <- wk.f %o% wk.s
  
  # Abundance
  wk.n <- array(NA, dim = dim(wk.z))
  wk.n[1,] <- rev(wk.r[1:wk.na])
  wk.n[2:wk.ny,1] <- wk.r[(wk.na+1):length(wk.r)]
  
  vecs <- array(NA, dim = c(wk.na*wk.ny, 4))
  vecs[,1] <- matrix(data = wk.n, nrow = wk.na * wk.ny, ncol = 1)
  vecs[,2] <- rep(1:wk.na, each = wk.ny)
  vecs[,3] <- rep(y1:(y1 + wk.ny - 1), wk.na)
  vecs[,4] <- vecs[,3] - vecs[,2]
  
  cz.list <- tapply(wk.z, vecs[,4], cumsum)
  
  vecs.list <- lapply(levels(as.factor(vecs[,4])), function(wk){
    temp <- vecs[vecs[,4] == wk,]
    temp.rep <- dim(temp)[1]
    if (!is.null(temp.rep)) 
    {
      temp[,1] <- rep(temp[1], temp.rep)
    }
    temp
  })
  
  vecs.list <- lapply(vecs.list, function(wk){
    temp.rep <- dim(wk)[1]
    if (!is.null(temp.rep))
    {
      wk.zz <- unlist(cz.list[as.character(wk[1,4])])
      wk <- cbind(wk, wk.zz)
      wk.a <- dim(wk)[1]
      # lnN(a,y) <- lnN(a-1,y-1) - z(a-1,y-1)
      wk[2:wk.a,1] <- wk[1:(wk.a-1),1] - wk[1:(wk.a-1),5]
    }else
    {
      wk.zz <- unlist(cz.list[as.character(wk[4])])
      wk <- c(wk, wk.zz)
    }    
    wk
  })
  
  vecs.table <- do.call(rbind, vecs.list)
  vecs.table <- vecs.table[order(vecs.table[,3], vecs.table[,2]),]
  
  new.n <- matrix(vecs.table[,1], nrow = wk.ny, ncol = wk.na, byrow = TRUE)
  wk.n <- exp(new.n)
  
  # Fitted survey indices I.hat and
  # back-transformed observed survey indices I.dash.star
  i.hat <- vector("list", length = wk.nk)
  i.dash.star <- vector("list", length = wk.nk)
  for (k in 1:wk.nk)
  {
    i.hat[[k]] <- wk.n * array(unlist(qval[[k]]), dim = dim(wk.n))
    i.dash.star[[k]] <- array(unlist(x[[k]] * exp(wk.z * rho[k])), dim = dim(wk.n))
  }
  
  # Survey log residuals
  res1 <- vector("list", length = wk.nk)
  out0 <- vector("list", length = wk.nk)
  for (k in 1:wk.nk)
  {
    res1[[k]] <- sqrt(wt[[k]]) * (log(i.dash.star[[k]]) - log(i.hat[[k]]))
    out0[[k]] <- array(NA, dim = c(wk.na * wk.ny, 1))
    out0[[k]][1:(wk.na*wk.ny),1] <- unlist(res1[[k]])
    out0[[k]] <- as.vector(na.exclude(out0[[k]]))
  }
  out1 <- as.vector(unlist(out0))
  
  # Lambda smoother
  f1 <- wk.f[1:(wk.ny-2)]
  f2 <- wk.f[2:(wk.ny-1)]
  res2 <- sqrt(lambda) * (f1 - f2)
  out2 <- as.vector(res2)
  
  # Overall SSQ
  as.numeric(c(out1,out2))
}