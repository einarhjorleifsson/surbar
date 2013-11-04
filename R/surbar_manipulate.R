#' Setup for surba
#' 
#' @export
#' @param dat Input data. A list ....
surba_setup <- function(dat=NULL) 
  {
  out <- dat
  out$rho <- unlist(lapply(out$u, function(wk) {wk$rho})) # time of survey
  out$numk <- length(out$u)
  out$wt.vec <- rep(1.0, out$numk)
  
  # setup a list with length equal to the number of surveys
  x <- vector("list", length = out$numk)
  names(x) <- unlist(lapply(out$u, function(wk){wk$name}))
  
  # set up data.frame for each survey with dimensions equivalent
  #   to data.frame sw
  x <- lapply(x, function(wk)
    {
    tab <- data.frame(array(NA, dim = dim(out$sw)))
    names(tab) <- names(out$sw)
    rownames(tab) <- rownames(out$sw)
    tab
    })
  
  # populate x with survey indices
  for (k in 1:out$numk)
    {
    a <- names(out$u[[k]]$tab)    #dat$u[[k]]$a1:dat$u[[k]]$a2
    y <- rownames(out$u[[k]]$tab)  # dat$u[[k]]$y1:dat$u[[k]]$y2
    #x[[k]][ch(y),ch(a)] <- dat$u[[k]]$tab[ch(y),ch(a)]
    x[[k]][y,a] <- out$u[[k]]$tab[y,a]
  }
  
  # set these numbers also to that in the u-frames
  #for (k in 1:out$numk)
  #  {
  #  out$u[[k]]$tab <- x[[k]]
  #}
  
  # Replace zeros with minimum value
  for (k in 1:out$numk)
  {
    xx <- x[[k]]
    x.min1 <- min(xx, na.rm = TRUE)
    if (x.min1 == 0)
    {
      xx[xx == 0] <- 9999
      x.min2 <- min(xx, na.rm = TRUE)
      xx[xx == 9999.0] <- x.min2
    }
    x[[k]] <- xx
  }
  
  
  # Mean-standardise survey data
  for (k in 1:out$numk)
  {
    x.mean <- mean(colMeans(x[[k]], na.rm = TRUE), na.rm = TRUE)
    x[[k]] <- x[[k]] / x.mean
  }  
  
  # Set up standard catchability and weightings arrays
  qval <- vector("list", length = out$numk)
  wt <- vector("list", length = out$numk)
  for (k in 1:out$numk)
  {
    qval[[k]] <- data.frame(array(1.0, dim = dim(x[[k]])))
    rownames(qval[[k]]) <- rownames(x[[k]])
    names(qval[[k]]) <- names(x[[k]])
    wt[[k]] <- data.frame(array(out$wt.vec[k], dim = dim(x[[k]])))
    rownames(wt[[k]]) <- rownames(x[[k]])
    names(wt[[k]]) <- names(x[[k]])
  }
  
  out$ages <- as.numeric(names(out$sw))
  out$years <- as.numeric(rownames(out$sw))
  out$u.std <- x
  out$qval <- qval
  out$wt <- wt
  
  
  # Initital parameter setup
  #s0.a <- seq(1.0, 1.0, length = ref.age - a1 + 1)
  #s0.b <- seq(1.0, 1.0, length = a2 - ref.age)
  s0.a <- seq(1.0, 1.0, length = out$ref.age - min(out$ages) + 1)
  s0.b <- seq(1.0, 1.0, length = max(out$ages) - out$ref.age)
  s0 <- c(s0.a[1:(length(s0.a) - 1)], s0.b[2:length(s0.b)])
  if (length(na.omit(s0)) < length(s0))
  {
    stop("\nCheck reference age.\n")
  }
  
  
  f0 <- rep(1.0, length(out$years) - 1)
  
  # Initial estimates for r are taken from the log of the average (across surveys)
  # of the survey index values at the appropriate years/ages
  temp.r <- data.frame(array(NA, dim = dim(out$sw)))
  names(temp.r) <- names(out$sw)
  rownames(temp.r) <- rownames(out$sw)
  for (i in rownames(temp.r))
  {
    for (j in names(temp.r))
    {
      temp.r[i,j] <- mean(unlist(lapply(out$u.std, function(wk){wk[i, j]})), na.rm = TRUE)
    }
  }
  # Fill missing values with age-based averages
  for (j in names(temp.r))
  {
    temp.mean <- mean(temp.r[,j], na.rm = TRUE)
    temp.r[is.nan(temp.r[,j]),j] <- temp.mean
  }
  r0 <- log(as.numeric(unlist(c(rev(temp.r[1,]), temp.r[2:dim(temp.r)[1],1]))))
  out$params0 <- c(s0, f0, r0)
  names(out$params0) <- c(rep('s',length(s0)),rep('f',length(f0)),rep('r',length(r0)))
  return(out)
}



#' Get surba parameters
#' 
#' @export
#' @param dat dat
surba_par <- function(dat) {
  s <- rep(NA, length(dat$ages))
  s[1:(dat$ref.age-1)] <- dat$fit$par[1:(dat$ref.age-1)]
  s[dat$ref.age] <- 1.0
  s[(dat$ref.age+1):(length(dat$ages)-1)] <- dat$fit$par[dat$ref.age:(length(dat$ages)-2)]
  s[length(dat$ages)] <- s[length(dat$ages)-1]
  names(s) <- dat$ages
  
  f <- rep(NA, length(dat$years))
  f[1:(length(dat$years)-1)] <- dat$fit$par[(length(dat$ages)-1):(length(dat$ages) + length(dat$years) - 3)]
  f[length(dat$years)] <- mean(f[(length(dat$years)-3):(length(dat$years)-1)])
  names(f) <- dat$years
  
  r <- dat$fit$par[(length(dat$ages) + length(dat$years) - 2):length(dat$fit$par)]
  y <- c((min(dat$years)-length(dat$ages)+1):max(dat$years))
  names(r) <- y
  return(list(s=s,f=f,r=r))
}

#' Get surba residuals
#' 
#' @export
#' @param dat dat
surba_res <- function(dat)
{
  # Extract log residuals
  numk <- dat$numk
  res <- vector("list", length = numk)
  x.res <- dat$fit$fvec
  idx <- dat$u
  x <- dat$u.std
  res.start <- 1
  for (k in 1:numk)
  {
    x.temp <- x[[k]][rownames(idx[[k]]$tab),names(idx[[k]]$tab)]
    res.end <- (res.start - 1) + (dim(x.temp)[1] * dim(x.temp)[2])
    res[[k]] <- data.frame(array(x.res[res.start:res.end], dim = dim(x.temp)))
    names(res[[k]]) <- names(x.temp)
    rownames(res[[k]]) <- rownames(x.temp)
    res.start <- res.end + 1
  }
  
  smx.res <- NULL
  for (i in 1:dat$numk) {
    tmp <- res[[i]]
    tmp$year <- rownames(tmp)
    tmp$survey <- dat$u[[i]]$name
    tmp <- melt(tmp,id.vars=c('year','survey'))
    smx.res <- rbind(smx.res,tmp)
  }
  names(smx.res)[3] <- c('age')
  smx.res$year <- as.integer(smx.res$year)
  smx.res$age <-  as.integer(smx.res$age)
  res <- smx.res
  return(res)
}

#' Calculate results by year
#' 
#' @export
#' @param dat dat
surba_rby <- function(dat) {
  #if(missing(dat$Par)) dat$Par <- surba_par(dat)
  y1 <- min(dat$years)
  y2 <- max(dat$years)
  ny <- y2-y1+1
  f <- dat$par$f
  s <- dat$par$s
  r <- dat$par$r
  zbarage.1 <- dat$zbarage.1
  zbarage.2 <- dat$zbarage.2
  sw <- dat$sw
  mat <- dat$mat
  na <- length(dat$ages)
  
  stock <- data.frame(array(NA, dim = c(ny,5)))
  names(stock) <- c("year", "rec", "ssb", "tsb", "meanz")
  
  zmort <- f %o% s
  
  lnn <- array(NA, dim = dim(zmort))
  lnn[1,] <- rev(r[1:dim(lnn)[2]])
  lnn[2:dim(lnn)[1],1] <- r[(dim(lnn)[2]+1):length(r)]
  # Oldest true age
  wk.ota <- dim(zmort)[1]
  for (j in 2:dim(zmort)[2])
  {
    for (i in 2:dim(zmort)[1])
    {
      lnn[i,j] <- lnn[i-1,j-1] - zmort[i-1,j-1]
    }
  }
  n <- exp(lnn)	
  
  stock$year <- y1:y2
  stock$meanz <- apply(zmort[,zbarage.1:zbarage.2], 1, mean)
  stock$rec <- exp(r[na:length(r)])
  stock$ssb <- apply(n * sw * mat, 1, sum)	
  stock$tsb <- apply(n * sw, 1, sum)
  
  # Mean-standardise rec, ssb, tsb (if required)
  # mean.rec <- mean(stock$rec)
  # mean.ssb <- mean(stock$ssb)
  # mean.tsb <- mean(stock$tsb)
  # stock$rec <- stock$rec / mean.rec
  # stock$ssb <- stock$ssb / mean.ssb
  # stock$tsb <- stock$tsb / mean.tsb
  return(stock)
}


#' Bootstrap model fit using data simulation
#' 
#' @export
#' @param dat dat
surba_boot <- function(dat) {
  
  Fit <- dat$fit
  y1 <- min(dat$years)
  y2 <- max(dat$years)
  ny <- y2-y1+1
  na <- length(dat$ages)
  a1 <- min(dat$ages)
  a2 <- max(dat$ages)
  ref.age <- dat$ref.age
  f <- dat$par$f
  s <- dat$par$s
  r <- dat$par$r
  zbarage.1 <- dat$zbarage.1
  zbarage.2 <- dat$zbarage.2
  sw <- dat$sw
  mat <- dat$mat
  
  
  # Test for singularity
  x.eigen <- eigen(Fit$hessian, only.values = TRUE)$values
  if (length(x.eigen[x.eigen == 0]) > 0)
    {
    stop("At least one parameter cannot be estimated: \nCheck that there are data for each cohort.")
  }
  
  n.psim <- 1000
  x.psim <- mvrnorm(n = n.psim, mu = Fit$par, vcov(Fit))
  
  x.stock <- data.frame(array(NA, dim = c(ny, 5)))
  names(x.stock) <- c("year", "rec", "ssb", "tsb", "meanz")
  
  x.psim.stock <- vector("list", n.psim)
  x.psim.stock <- lapply(x.psim.stock, function(wk){wk <- x.stock})
  
  x.psim.s <- array(NA, dim = c(1000, na))
  x.psim.s[,1:(ref.age-1)] <- x.psim[,1:(ref.age-1)]
  x.psim.s[,ref.age] <- 1.0
  x.psim.s[,(ref.age+1):(na-1)] <- x.psim[,ref.age:(na-2)]
  x.psim.s[,na] <- x.psim.s[,na-1]

  x.psim.f <- array(NA, dim = c(1000, ny))
  x.psim.f[,1:(ny-1)] <- x.psim[,(na-1):(na + ny - 3)]
  x.psim.f[,ny] <- apply(x.psim.f[,(ny-3):(ny-1)], 1, mean)
  
  x.psim.r <- x.psim[,(na + ny - 2):length(Fit$par)]
  
  x.s <- rep(NA, length = na)
  x.f <- rep(NA, length = ny)
  
  for (i in 1:n.psim)
    {
    # Extract parameters
    x.s[1:(ref.age-1)] <- x.psim[i,1:(ref.age-1)]
    x.s[ref.age] <- 1.0
    x.s[(ref.age+1):(na-1)] <- x.psim[i,ref.age:(na-2)]
    x.s[na] <- x.s[na-1]
    
    x.f[1:(ny-1)] <- x.psim[i,(na-1):(na + ny - 3)]
    x.f[ny] <- mean(x.f[(ny-3):(ny-1)])
    
    x.r <- x.psim[i,(na + ny - 2):length(Fit$par)]
    
    # Generate stock summaries
    zmort <- x.f %o% x.s
    lnn <- array(NA, dim = dim(zmort))
    lnn[1,] <- rev(x.r[1:dim(lnn)[2]])
    lnn[2:dim(lnn)[1],1] <- x.r[(dim(lnn)[2]+1):length(x.r)]
    for (jj in 2:dim(zmort)[2])
      {
      for (ii in 2:dim(zmort)[1])
        {
        lnn[ii,jj] <- lnn[ii-1,jj-1] - zmort[ii-1,jj-1]
      }
    }
    
    n <- exp(lnn)	
    
    x.psim.stock[[i]]$year <- y1:y2
    x.psim.stock[[i]]$meanz <- apply(zmort[,zbarage.1:zbarage.2], 1, mean)
    x.psim.stock[[i]]$z <- zmort
    x.psim.stock[[i]]$rec <- exp(x.r[na:length(r)]) # / mean.rec
    x.psim.stock[[i]]$ssb <- apply(n * sw * mat, 1, sum) # / mean.ssb
    x.psim.stock[[i]]$tsb <- apply(n * sw, 1, sum) # / mean.tsb
  }
  rby <- do.call(rbind,x.psim.stock)[,1:5]
  rby$iter <- rep(1:n.psim,each=ny)
  z <- do.call(rbind,x.psim.stock)[,6]
  z <- as.data.frame(z)
  names(z) <- a1:a2
  z$year <- rby$year
  z$iter <- rep(1:n.psim,each=ny)
  
  rownames(x.psim.s) <- 1:n.psim
  colnames(x.psim.s) <- a1:a2
  
  rownames(x.psim.f) <- 1:n.psim
  colnames(x.psim.f) <- y1:y2
  
  rownames(x.psim.r) <- 1:n.psim
  colnames(x.psim.r) <- (y1-(a2-a1)-a1):(y2-a1)
  
  
  par <- list(s=x.psim.s,f=x.psim.f,r=x.psim.r)
  ret <- list(rby=rby,z=z,par=par)
  return(ret)
}

#' Trim data
#' 
#' @export
#' @param dat dat
surba_trim <- function(dat) 
  {
  x <- dat
  x$sw <- x$sw[1:(nrow(x$sw)-1),]
  x$mat <- x$mat[1:(nrow(x$mat)-1),]
  for (k in 1:x$numk)
  {
    x$u[[k]]$tab <- x$u[[k]]$tab[1:(nrow(x$u[[k]]$tab)-1),]
  }
  x <- surba_setup(x)
  return(x)
}
