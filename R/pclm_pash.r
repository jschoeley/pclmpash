require(Matrix)
require(splines)
require(MortalitySmooth)
require(pash)

# pclm.control ---------------------------------------------------------------

#' Control the PCLM Fitting
#'
#' Auxiliary function for controlling PCLM fitting. Use this function to set
#' control parameters of the \code{\link{pclm.fit}} and other related functions.
#'
#' @param x.div Number of sub-classes within PCLM tim/age class (default is 1).
#'   Low value of the parameter makes the PCLM computation faster. It is however
#'   recommended to set it to higher value (e.g. 10) for better \code{nax}
#'   estimates.
#' @param x.auto.trans Logical indicating if automatically multiple age
#'   intervals to remove fractions. \code{TRUE} is the recommended value. See
#'   also examples in \code{\link{pclm.fit}}.
#' @param x.max.ext Integer defining maximal multiple of an age interval. See
#'   also \code{\link{pclm.interval.multiple}}.
#' @param zero.class.add Logical indicating if additional zero count class (open
#'   interval) should be added after last age class. \code{TRUE} is the
#'   recommended value. See \code{\link{pclm.nclasses}} and
#'   \code{\link{pclm.compmat}}.
#' @param zero.class.end Positive indicating the end of zero count class =
#'   anticipated end of last (open) interval. If set to \code{NULL} and
#'   \code{zero.class.add == TRUE} then it is calculated automatically based on
#'   \code{zero.class.frac}. See \code{\link{pclm.nclasses}} and
#'   \code{\link{pclm.compmat}}.
#' @param zero.class.frac Fraction of total range of \code{x} (age/time vector)
#'   added as the last zero-count interval when \code{zero.class.end == NULL}.
#'   Used in \code{\link{pclm.compmat}}. Increase this value if the right tail
#'   of the PCLM fit is badly fitted (use \code{\link{plot.pclm}} to diagnose).
#' @param bs.use Logical indicating if use B- or P-spline basis to speed-up
#'   computations. Possible values: \code{"auto"}, \code{TRUE}, and
#'   \code{FALSE}. Used by \code{\link{pclm.compmat}} function.
#' @param bs.method Basis for B- or P-spline used by \code{\link{pclm.compmat}}
#'   function. Possible values: \itemize{ \item{\code{"MortalitySmooth"}}{ -
#'   gives "P-splines" basis based on \code{\link{MortSmooth_bbase}}
#'   \code{\{\link{MortalitySmooth}\}} (recommended)} \item{\code{"bs"}}{ -
#'   gives basic B-splines basis based on \code{\link{bs}}
#'   \code{\{\link{splines}\}}.} }
#' @param bs.df B- or P- spline degree of freedom (df, number of inner knots) or
#'   a way to its calculation used in \code{\link{pclm.compmat}} function. The
#'   value is automatically limited by the \code{bs.df.max}. It can take
#'   corresponding values: \itemize{ \item{\code{"maxprec"}}{ - df equal to the
#'   number of ungrouped (raw) age classes (recommended option).}
#'   \item{\code{"thumb"}}{ - 'rule of thumb': one knot for the B-spline basis
#'   each 4-5 observations.} \item{\code{integer}}{ - df given explicitly.} }
#' @param bs.df.max Maximal number of knots (df) for B- or P-spline basis.
#'   Defaut value is 200, but can be decreased if computations are slow. Used in
#'   \code{\link{pclm.compmat}}.
#' @param bs.deg Degree of the piecewise polynomial for B- or P-spline basis.
#'   Default and recommended value is 3. Used in \code{\link{pclm.compmat}}.
#' @param opt.method Selection criterion for \code{lambda} (smooth parameter) in
#'   \code{\link{pclm.opt}} function. Possible values are \code{"AIC"} and
#'   \code{"BIC"} (recommended).
#' @param opt.tol Tolerance for \code{\link{pclm.opt}} function that estimates
#'   smooth parameter \code{lambda}.
#' @param pclm.deg Order of differences of the components of \code{b} (PCLM
#'   coeficients, beta in the reference [1]). Default value is 2, some other
#'   values may cause algorithm to not work properly. Used by
#'   \code{\link{pclm.core}} function.
#' @param pclm.max.iter Maximal number of iterations in \code{\link{pclm.core}}
#'   function. Default is 100, but increase it when got a warning.
#' @param pclm.lsfit.tol Tolerance for \code{\link{lsfit}} function used in
#'   \code{\link{pclm.core}} function.
#' @param pclm.tol Tolerance for \code{\link{pclm.core}} function.
#' @return List with control parameters.
#' @seealso \code{\link{pclm.fit}}, \code{\link{pclm.general}},
#'   \code{\link{pclm.core}}, \code{\link{pclm.opt}},
#'   \code{\link{pclm.aggregate}}, \code{\link{pclm.compmat}},
#'   \code{\link{pclm.interval.multiple}}, \code{\link{pclm.nclasses}},
#'   \code{\link{plot.pclm}}, and \code{\link{summary.pclm}}.
#' @references \enumerate{ \item{Rizzi S, Gampe J, Eilers PHC. Efficient
#'   estimation of smooth distributions from coarsely grouped data. Am J
#'   Epidemiol. 2015;182:138?47.} }
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @export
pclm.control<-function(x.div = 1L,
                       x.auto.trans = TRUE,
                       x.max.ext = 25L,
                       zero.class.add = TRUE,
                       zero.class.end = NULL,
                       zero.class.frac = 0.2,
                       bs.use = 'auto',
                       bs.method = c('MortalitySmooth', 'bs'),
                       bs.df = c('maxprec', 'thumb'),
                       bs.df.max = 200L,
                       bs.deg = 3L,
                       opt.method = c('BIC','AIC'),
                       opt.tol = .Machine$double.eps^0.5,
                       pclm.deg  =  2L,
                       pclm.max.iter = 100L,
                       pclm.lsfit.tol = .Machine$double.eps^0.5,
                       pclm.tol = .Machine$double.eps^0.5){

  if (!(opt.method[1] %in% c('BIC','AIC'))) stop ('"AIC" or "BIC" should be used for opt.method')
  if (!(bs.df[1] %in% c('maxprec','thumb'))) if (!is.numeric(bs.df[1]))
    stop('bs.df can take "maxprec", "thumb", or any positive integer value')
  if (!(bs.method[1] %in% c('MortalitySmooth', 'bs'))) stop ('"MortalitySmooth" or "bs" can only be used for bs.method')

  list(x.max.ext = x.max.ext, x.auto.trans = x.auto.trans, x.div = x.div,
       zero.class.add = zero.class.add, zero.class.end = zero.class.end,
       zero.class.frac = zero.class.frac, bs.use = bs.use,
       bs.method = bs.method[1], bs.df = bs.df[1], bs.df.max = bs.df.max,
       bs.deg = bs.deg, opt.method = opt.method[1], opt.tol = opt.tol,
       pclm.deg = pclm.deg, pclm.max.iter = pclm.max.iter,
       pclm.lsfit.tol = pclm.lsfit.tol, pclm.tol = pclm.tol)
}

# pclm.interval.multiple ---------------------------------------------------------------

#' Calculate Smallest Factor Removing Fractions in Input Vector
#'
#' Calculates minimal multiple of the orginal age/time interval length to remove
#' fractions in \code{x}. The function (together with \code{x} and some control
#' parameters) is used to calculate the internal (raw, nonaggregated) interval
#' length during PCLM computations.
#'
#' @param x Vector with start of the interval for age/time classes.
#' @param control List with additional parameters. See
#'   \code{\link{pclm.control}}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @seealso \code{\link{pclm.general}}, \code{\link{pclm.control}}, and
#'   \code{\link{pclm.nclasses}}.
#' @export
pclm.interval.multiple <- function(x, control = list()) {
  control <- do.call("pclm.control", control)
  frac<-function(x, digits = floor(-log10(.Machine$double.eps^0.5))) round(x-floor(round(x, digits)), digits)
  for (j in 1:(control$x.max.ext)){
    y <- frac(j * x)
    if (all(y == 0)) break
  }
  j
}

# pclm.nclasses -----------------------------------------------

#' Calculate the Number of Internal PCLM Age Classes
#'
#' Calculate the number of internal PCLM age classes used in construction of the
#' composition matrix.
#'
#' @param x Vector with start of the interval for age/time classes.
#' @param control List with additional parameters. See
#'   \code{\link{pclm.control}}.
#' @examples
#' \dontrun{
#' # Use a simple data set
#' AU10 <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,
#'    nax = australia_10y$nax, nx = australia_10y$nx, last_open = TRUE)
#'
#' # Define the open interval by zero.class.frac
#' control.1 = list(x.div = 5, zero.class.frac = 0.2, zero.class.end = NULL)
#' pclm.nclasses(AU10$lt$x, control = control.1) #calculate number of raw classes
#' AU10p.1A <- pclm.fit(AU10, control = control.1)
#' length(AU10p.1A$pclm$raw$x) # the number of raw classes after fit
#' plot(AU10p.1A)
#'
#' # Define the open interval by zero.class.end
#' control.2 = list(x.div = 5, zero.class.end = 109)
#' pclm.nclasses(AU10$lt$x, control = control.2) #calculate the number of raw classes
#' AU10p.1B <- pclm.fit(AU10, control = control.2)
#' length(AU10p.1B$pclm$raw$x) # the number of raw classes after fit
#' plot(AU10p.1B)
#'
#' # **** See more examples in the help for pclm.fit() function.
#' }
#' @seealso \code{\link{pclm.fit}}, \code{\link{pclm.control}},
#'   \code{\link{pclm.interval.multiple}},
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @export
pclm.nclasses<-function(x, control = list()) {
  control <- do.call("pclm.control", control)
  if (control$zero.class.add)
    if (length(control$zero.class.end) == 0){
      drx <- diff(range(x))
      tmp <- control$x.div * drx * pclm.interval.multiple(x, control)
      tmp <- tmp * (1 + control$zero.class.frac)
      tmp <- tmp + 1
    } else {
      drx <- diff(range(x, control$zero.class.end))
      tmp <- control$x.div * drx * pclm.interval.multiple(x, control)
      tmp <- tmp + 1
    }
  tmp
}

# pclm.compmat -------------------------------------------------------------------------

#' Create Composition Matrix Object
#'
#' Construct the composition matrix object for automatically recalibrated age
#' classes.
#'
#' @param x Vector with start of the interval for age/time classes. x * x.div
#'   must be an integer. The appropriate correction for fractional intervals
#'   based on the interval multiple (\code{\link{pclm.interval.multiple}}) is
#'   performed in \code{\link{pclm.general}}.
#' @param y Vector with counts, e.g. \code{ndx}. It must have the same length as
#'   \code{x}.
#' @param exposures Vector with exposures used to calculate smoothed mortality
#'   rates (see reference [1] and \code{\link{pclm.general}}).
#' @param control List with additional parameters. See
#'   \code{\link{pclm.control}}.
#' @return List with components:
#' @return \item{\code{C}}{ Composition matrix.}
#' @return \item{\code{X}}{ B-spline base, P-spline base, or identity matrix.}
#' @return \item{\code{x}}{ Corrected age/time vector.}
#' @return \item{\code{y}}{ Corrected vector with counts.}
#' @return \item{\code{open.int.len}}{ Length of the open interval in age
#'   classes.}
#' @return \item{\code{exposures}}{ Vector with exposures if it was used to
#'   construct the composition matrix.}
#' @return \item{\code{control}}{Used control parameters, see
#'   \code{\link{pclm.control}}.}
#' @return \item{\code{warn.list}}{ List with warnings.}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @details The details of matrix construction can be found in reference [1]. if
#' \code{bs.use == TRUE} then P- or B- splines are used instead of identity
#' matrix (see \code{\link{pclm.control}}).\cr\cr The dimension of constructed
#' composition matrix can be determined before its computation. The shorter
#' dimension equals to the length of data vector + 1, whereas the longer
#' dimension is determined by the function \code{\link{pclm.nclasses}} and for
#' \code{zero.class.end == NULL} equals:\cr\cr \code{(x.div * (max(x) - min(x))
#' * m) * (1 + zero.class.frac) + 1}\cr\cr or\cr\cr \code{x.div *
#' (zero.class.end - min(x)) * m + 1}\cr\cr otherwise, where \code{m} is an
#' interval multiple calculated by \code{\link{pclm.interval.multiple}}. See
#' also \code{\link{pclm.nclasses}}.
#' @references \enumerate{ \item{Rizzi S, Gampe J, Eilers PHC. Efficient
#' estimation of smooth distributions from coarsely grouped data. Am J
#' Epidemiol. 2015;182:138?47.} \item{Camarda, C. G. (2012). MortalitySmooth: An
#' R Package for Smoothing Poisson Counts with P-Splines. Journal of Statistical
#' Software. 50, 1-24.} \item{Hastie, T. J. (1992) Generalized additive models.
#' Chapter 7 of Statistical Models in S eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.} }
#' @seealso \code{\link{pclm.general}}, \code{\link{pclm.control}},
#'   \code{\link{pclm.interval.multiple}}, and \code{\link{pclm.nclasses}}.
#' @keywords internal
pclm.compmat<-function(x, y, exposures = NULL, control = list()){
  require(Matrix)
  require(splines)
  require(MortalitySmooth)

  control <- do.call("pclm.control", control)

  if (control$bs.use == 'auto') {
    control$bs.use <- (pclm.nclasses(x, control) >= control$bs.df.max)
    message(paste('bs.use automatically set to', control$bs.use))
  } else if (!is.logical(control$bs.use)) stop('bs.use can take only "auto", TRUE or FALSE values.')
  if ((pclm.nclasses(x, control) >= control$bs.df.max) && (control$bs.use == FALSE)) warning(immediate. = TRUE, 'Big composition matrix detected. Calculations can be slow. Set "bs.use" to TRUE.')

  .MortSmooth_Bbase<-function (x, xl = min(x), xr = max(x), df = floor(length(x) / 4), deg = control$bs.deg) {
    #Modified version of MortalitySmooth:::MortSmooth_bbase
    ndx <- df-deg
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, function (x, t)  (x - t)^deg * (x > t))
    D <- diff(diag(dim(P)[2]), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
    B <- (-1)^(deg + 1) * P %*% t(D)
    B
  }

  warn.list <- list()

  if (any(abs(as.integer(x) - x) > 1e-6)) {
    WARN <- 'Fractional values in age vector. Values were rounded'
    warning(immediate. = TRUE, WARN)
    warn.list <- c(warn.list, WARN)
  }

    #x <- as.integer(round(x))
  if ((length(x) < (3 + control$bs.deg)) && (control$bs.use)) {
    WARN <- 'Not enough data classes to use B-spline basis. Exact method was used'
    warning(immediate. = TRUE, WARN)
    warn.list <- c(warn.list, WARN)
    control$bs.use <- FALSE
  }

  if (control$zero.class.add) {
    if (length(control$zero.class.end) == 0){
      mstep <- min(diff(x))
      r <- diff(range(x))
      control$zero.class.end <- ceiling(ceiling(r * control$zero.class.frac) / mstep) * mstep + max(x)
    }
    y <- c(y, 0)
    open.int.len <- control$zero.class.end - x[length(x)]
    x <- c(x, control$zero.class.end)
  } else open.int.len <- NA

  segm <- round(x * control$x.div)
  if (length(unique(segm)) != length(segm)) stop('Non-unique values in x after rounding. Try to increase x.max.ext, x.div, or round age classes manually.')
  C2 <- approx(x = segm + 1, y = 1:length(segm), method = 'constant', xout = (min(segm):(max(segm))) + 1, rule = 2)
  SparMat <- Matrix:::sparseMatrix(C2$y, C2$x)
  C. <- 1 * Matrix(SparMat, sparse = FALSE)
  if (control$bs.use) {
    if (length(exposures) > 0) stop('Exposures cannot be used together with B-Splines.')
    control$bs.method <- control$bs.method[1]
    control$bs.df <- control$bs.df[1]
    if (control$bs.df == "maxprec") control$bs.df <- min(control$bs.df.max, length(x) * control$x.div-2)
    else if (tolower(control$bs.df) == 'thumb') control$bs.df <- min(control$bs.df.max, (length(x) * control$x.div-2)/4)
    else if(is.na(as.numeric(control$bs.df))) stop(paste('Wrong bs.df parameter:', control$bs.df))
    control$bs.df <- max(length(segm)-2, min(control$bs.df, length(segm) * control$x.div-2))
    if (control$bs.method == 'bs'){
      X <- splines:::bs(C2$x, df = control$bs.df)
    } else if (tolower(control$bs.method) == 'mortalitysmooth'){
      X <- .MortSmooth_Bbase(C2$x, df = control$bs.df)
    } else stop('Unknown B-spline basis method.')
  } else{
    if (length(exposures) > 0){
      d <- length(x)-length(exposures)
      if (d > 0) exposures <- c(exposures, rep(1e-6, d))
      expo <- as.matrix(do.call('cbind', rep(list(exposures), dim(C.)[2])))
      # expo <- (rep(c(exposures), dim(C.)[2]))
      # dim(expo) <- rev(dim(C.))
      # expo <- t(expo)
      C. <- C. * expo
    } else expo <- NULL
    #control$bs.df <- NA
    X <- diag(dim(C.)[2])
  }
  list(C = as.matrix(C.),
       X = X,
       y = y,
       x = x,
       open.int.len = open.int.len,
       exposures = exposures,
       control = control,
       warn.list = warn.list)
}

# pclm.core ------------------------------------------------------------------------

#' Fit the Penalized Composite Link Model (PCLM)
#'
#' Efficient estimation of smooth distribution from coarsely grouped data based
#' on PCLM algorithm described in Rizzi et al. 2015. For further description see
#' the reference [1] and \code{\link{pclm.fit}}.
#'
#' @param CompositionMatrix Object constructed by \code{\link{pclm.compmat}}.
#' @param lambda Smoothing parameter.
#' @param control List with additional parameters. See
#'   \code{\link{pclm.control}}.
#' @return List with components:
#' @return \item{\code{gamma}}{  Ungrouped counts.}
#' @return \item{\code{dev}}{Deviance.}
#' @return \item{\code{aic}}{AIC for the fitted model.}
#' @return \item{\code{bic}}{BIC for the fitted model.}
#' @return \item{\code{warn.list}}{List with warnings.}
#' @author Silvia Rizzi (original code, see Appendix #2 of the reference [1]),
#' Maciej J. Danko (modification) <\email{danko@demogr.mpg.de}>
#' <\email{maciej.danko@gmail.com}>
#' @references \enumerate{ \item{Rizzi S, Gampe J, Eilers PHC. Efficient
#' estimation of smooth distributions from coarsely grouped data. Am J
#' Epidemiol. 2015;182:138?47.} }
#' @keywords internal
pclm.core <- function(CompositionMatrix, lambda = 1, control = list()){
  control <- do.call("pclm.control", control)
  warn.list <- list()
  C <-  CompositionMatrix$C
  X <-  CompositionMatrix$X
  y <-  CompositionMatrix$y
  y <- as.matrix(as.vector(y))
  nx <- dim(X)[2] #number of small classes
  D <- base:::diff(diag(nx), diff = control$pclm.deg)
  bstart <- log(sum(y) / nx);
  b <- rep(bstart, nx);
  was.break <- FALSE
  for (it in 1:control$pclm.max.iter) {
    b0 <- b
    eta <- X %*% b
    gam <- exp(eta)
    mu <- C %*% gam
    w <- c(1 / mu, rep(lambda, nx - control$pclm.deg))
    Gam <- gam %*% rep(1, nx)
    Q <- C %*% (Gam * X)
    z <- c(y - mu + Q %*% b, rep(0, nx - control$pclm.deg))
    ls.x <- rbind(Q, D)
    Fit <- lsfit(ls.x, z, wt = w, intercept = FALSE, tolerance = control$pclm.lsfit.tol)
    b <- Fit$coef
    db <- max(abs(b - b0))
    if (db < control$pclm.tol) {
      was.break <- TRUE
      break
    }
  }
  if (!was.break) {
    WARN <- 'Maximum iteration reached without convergence. Try to increase pclm.max.iter parameter.'
    warning(immediate. = TRUE, WARN)
    warn.list <- c(warn.list, WARN)
  }
  R <- t(Q) %*% diag(c(1 / mu)) %*% Q
  H <- solve(R + lambda * t(D) %*% D) %*% R
  fit <- list()
  .trace <- sum(diag(H))
  ok <- y > 0 & mu > 0
  fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
  fit$gamma <- gam
  fit$aic <- fit$dev + 2 * .trace
  fit$bic <- fit$dev + .trace*log(length(y))
  fit$warn.list <- warn.list
  fit
}

# pclm.opt --------------------------------------------------------------------

#' Optimize the Smoothing Parameter \code{lambda} for PCLM Method
#'
#' @param CompositionMatrix Output of \code{\link{pclm.compmat}}.
#' @param control List with additional parameters. See
#'   \code{\link{pclm.control}}.
#' @return List with components:
#' @return \item{\code{X}}{Ungrouped age/time classes.}
#' @return \item{\code{Y}}{Ungrouped counts.}
#' @return \item{\code{lambda}}{Optimal smooth parameter.}
#' @return \item{\code{fit}}{Output of the \code{\link{pclm.core}} derived for
#'   the optimal \code{lambda}.}
#' @return \item{\code{CompositionMatrix}}{Used composition matrix. See also
#'   \code{\link{pclm.compmat}}.}
#' @return \item{\code{warn.list}}{List with warnings.}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @references \enumerate{ \item{Rizzi S, Gampe J, Eilers PHC. Efficient
#' estimation of smooth distributions from coarsely grouped data. Am J
#' Epidemiol. 2015;182:138?47.} }
#' @keywords internal
pclm.opt<-function(CompositionMatrix, control = list()){
  warn.list <- list()
  control <- do.call("pclm.control", control)
  tryme<-function(G, Altern = 1e200) suppressWarnings(if (class(try(G, silent = TRUE)) == "try-error") Altern else try(G, silent = TRUE))
  if (toupper(control$opt.method) == 'AIC') opty<-function (log10lam) tryme(pclm.core(CompositionMatrix, lambda = 10^log10lam, control = control)$aic) else
    if (toupper(control$opt.method) == 'BIC') opty<-function (log10lam) tryme(pclm.core(CompositionMatrix, lambda = 10^log10lam, control = control)$bic) else
      stop('Unknown method of lambda optimization.')
  res.opt <- stats:::optimize(f = opty, interval = c(-12, 22), tol = control$opt.tol)$minimum
  if((round(res.opt) <= -11.9) || (round(res.opt) >= 21.9)) {
    WARN <- 'Lambda reached boundary values.'
    warning(immediate. = TRUE, WARN)
    warn.list <- c(warn.list, WARN)
  }
  lambda  <-  10^res.opt
  #cat(res.opt,'\n')
  fit <- pclm.core(CompositionMatrix, lambda = lambda, control = control)
  X <- seq(min(CompositionMatrix$x), max(CompositionMatrix$x), 1 / control$x.div)
  Y <- fit$gamma
  Z <- list(X = X, Y = Y, lambda = lambda, fit = fit, CompositionMatrix = CompositionMatrix,
            warn.list = warn.list)
  Z
}

# pclm.aggregate ---------------------------------------------------

#' Calculate Raw and Aggregated Life-tables
#'
#' Calculate raw and aggregated life-tables From the object returned by the
#' \code{pclm.opt()} Function.
#'
#' @param fit Object obtained by \code{\link{pclm.opt}} function.
#' @param out.step Output interval length.
#' @param count.type Type of the data, deaths(\code{"DX"})(default) or
#'   exposures(\code{"LX"}.)
#' @param exposures.used Logical indicating if exposures were used to create
#'   composition matrix.
#' @return List with components:
#' @return \item{\code{grouped}}{Life-table based on aggregated PCLM fit defined
#'   by \code{out.step}.}
#' @return \item{\code{raw}}{Life-table based on original (raw) PCLM fit.}
#' @return \item{\code{fit}}{PCLM fit used to construct life-tables.}
#' @return \item{\code{warn.list}}{List with warnings.}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @seealso \code{\link{pclm.fit}}
#' @keywords internal
pclm.aggregate<-function(fit, out.step = NULL, count.type = c('DX', 'LX'), exposures.used = FALSE){
  count.type <- count.type[1]
  p <- function(x) round(x, floor(-log10(.Machine$double.eps^0.8)))
  warn.list <- fit$warn.list
  Y <- fit$Y
  X <- fit$X
  n <- diff(X)[1]#c(diff(X), diff(X)[length(X)-1])
  if (length(out.step) == 0) x <- p(fit$CompositionMatrix$x) else {
    if (p(out.step)<p(n)) {
      WARN <- 'Output age interval length (out.step) was too small and was adjusted. It equals the smallest age class now. Re-fit PCLM with higher x.div.'
      warning(immediate. = TRUE, WARN)
      warn.list <- c(warn.list, WARN)
      out.step <- n
    }
    tmp <- round(out.step/n) * n
    if (p(tmp) != p(out.step)) {
      WARN <- 'Output age interval length (out.step) rounded to an integer multiple of the smallest age class'
      warning(immediate. = TRUE, WARN)
      warn.list <- c(warn.list, WARN)
    }
    out.step <- tmp
    x <- p(unique(c(seq(X[1], X[length(X)], out.step), X[length(X)])))
  }

  if (toupper(count.type) == 'DX'){
    if (!exposures.used){
      ax <- rep(NA, length(x)-1)
      Dx <- rep(NA, length(x)-1)
      for (j in 1:(length(x)-1)){
        ind <- (x[j] <= X) & (x[j + 1] > X)
        sDx <- sum(Y[ind])
        mX <- sum(Y[ind] * (X[ind] - X[ind][1] + n/2))
        ax[j] <- mX / sDx
        Dx[j] <- sDx
      }
      grouped <- data.frame(x = x[-length(x)], lx = sum(Dx)-c(0, cumsum(Dx)[-length(Dx)]),
                         dx = Dx, ax = ax, n = diff(x), Ax = ax / diff(x))
      raw <- data.frame(x = X, lx = sum(Y)-c(0, cumsum(Y)[-length(Y)]), dx = Y,
                     ax = c(0.5 * diff(X), 0.5*diff(X)[length(X)-1]),
                     n = n, Ax = 0.5)
    } else {
      MX = rep(NA, length(x)-1)
      for (j in 1:(length(x)-1)){
        ind <- (x[j] <= X) & (x[j + 1] > X)
        MX[j] <- sum(Y[ind], na.rm=TRUE)
      }
      grouped <- data.frame(x = x[-length(x)], mx = MX, n = out.step)
      raw <- data.frame(x = X, mx = Y, n = n)
    }
  } else if (toupper(count.type) == 'LX'){
    if (exposures.used) stop('You cannot give exposures as extra parameter if you already selected "DX" as count.type.')
    Lx <- rep(NA, length(x) - 1)
    for (j in 1:(length(x) - 1)){
      ind <- (x[j] <= X) & (x[j + 1] > X)
      sLx <- sum(Y[ind])
      Lx[j] <- sLx
    }
    grouped <- data.frame(x = x[-length(x)], Lx = Lx, n = diff(x))
    raw <- data.frame(x = X, Lx = Y, n = n)
  } else stop('Unknown life-table type')
  object <- list(grouped = grouped,
              raw = raw,
              fit = fit,
              warn.list = warn.list)
  object
}

# pclm.general -------------------------------------------------------------------------

#' PCLM De-aggregation of Life-table
#'
#' De-aggregates a life-table using the PCLM method.
#'
#' @param x Vector with start of the interval for age/time classes.
#' @param y Vector with counts, e.g. \code{ndx}. It must have the same length as
#'   \code{x}.
#' @param count.type Type of the data, deaths(\code{"DX"})(default) or
#'   exposures(\code{"LX"}.)
#' @param out.step Age interval length in output aggregated life-table. If set
#'   to \code{"auto"} then the parameter is automatically set to the length of
#'   the shortest age/time interval of \code{x}.
#' @param exposures Optional exposures to calculate smooth mortality rates. A
#'   vector of the same length as \code{x} and \code{y}. See reference [1] for
#'   further details.
#' @param control List with additional parameters. See
#'   \code{\link{pclm.control}}.
#' @details The function has four major steps: \enumerate{ \item{Calculate
#' interval multiple (\code{\link{pclm.interval.multiple}} to remove fractional
#' parts from \code{x} vector. The removal of fractional parts is necessary to
#' build composition matrix.} \item{Calculate composition matrix using
#' \code{\link{pclm.compmat}}.} \item{Fit PCLM model using
#' \code{\link{pclm.opt}}.} \item{Calculate aggregated (grouped) life-table
#' using \code{\link{pclm.aggregate}}.} } More details for PCLM algorithm can be
#' found in reference [1], but see also \code{\link{pclm.fit}} and
#' \code{\link{pclm.compmat}}.
#' @return The output is of \code{"pclm"} class with the components:
#' @return \item{\code{grouped}}{Life-table based on aggregated PCLM fit and
#'   defined by \code{out.step}.}
#' @return \item{\code{raw}}{Life-table based on original (raw) PCLM fit.}
#' @return \item{\code{fit}}{PCLM fit used to construct life-tables.}
#' @return \item{\code{m}}{Interval multiple, see
#'   \code{\link{pclm.interval.multiple}}, \code{\link{pclm.compmat}}.}
#' @return \item{\code{x.div}}{Value of \code{x.div}, see
#'   \code{\link{pclm.control}}.}
#' @return \item{\code{out.step}}{Interval length of aggregated life-table, see
#'   \code{\link{pclm.control}}.}
#' @return \item{\code{control}}{Used control parameters, see
#'   \code{\link{pclm.control}}.}
#' @return \item{\code{warn.list}}{List with warnings.}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @examples
#' # The examples with use of the \code{pash} object are presented in \link{pclm.fit}.
#' # Explicit examples of use \code{pclm.general} (especially how to use exposures)
#' # are to be written in a next package release.
#' @references \enumerate{ \item{Rizzi S, Gampe J, Eilers PHC. Efficient
#' estimation of smooth distributions from coarsely grouped data. Am J
#' Epidemiol. 2015;182:138?47.} \item{Rizzi S, Thinggaard M, Engholm G, et al.
#' Comparison of non-parametric methods for ungrouping coarsely aggregated data.
#' BMC Medical Research Methodology. 2016;16:59. doi:10.1186/s12874-016-0157-8.}
#' }
#' @seealso \code{\link{pclm.fit}}, \code{\link{pclm.compmat}},
#'   \code{\link{pclm.interval.multiple}}, and \code{\link{pclm.nclasses}}.
#' @export
pclm.general <- function(x, y, count.type = c('DX', 'LX'), out.step = 'auto', exposures = NULL, control = list()){
  control <- do.call("pclm.control", control)
  count.type <- count.type[1]
  exposures.used  <-  (length(exposures) > 0)
  if (control$bs.use == 'auto') {
    control$bs.use <- (pclm.nclasses(x, control) >= control$bs.df.max)
    message(paste('bs.use automatically set to', control$bs.use))
  } else if (!is.logical(control$bs.use)) stop('bs.use can take only "auto", TRUE or FALSE values.')
  if ((pclm.nclasses(x, control) >= control$bs.df.max) && (control$bs.use == FALSE)) warning(immediate. = TRUE, 'Big composition matrix detected. Calculations can be slow. Set "bs.use" to TRUE.')
  if (tolower(out.step) == 'auto') {
    out.step <- min(diff(x), 1) # removed: 1/pclm.interval.multiple(x, control)
    message(paste('out.step automatically set to', out.step))
  } else if (!is.numeric(out.step)) stop('Unknown command for out.step. Set out.step as "auto" or as numeric value.')
  if (!control$zero.class.add) warning(immediate. = TRUE, 'Omitting zero.class may lead to biased results.')
  if ((control$zero.class.add) && (length(control$zero.class.end)>0) && (control$zero.class.end <= max(x))) stop ("zero.class.end lower than last age class.")

  if ((toupper(count.type) == 'LX') && (length(exposures) > 0)) stop('You cannot give exposures as extra parameter if you already selected "DX" as count.type.')

  WARN <- list()
  if (all(order(x) != 1:length(x))) stop ('Age classes are not ordered')

  #Automatically change the scale to make x integer
  #Output:
  # m - magnitude of change,
  # x - transformed vector of x
  if (control$x.auto.trans) {
    m <- pclm.interval.multiple(x, control)
    if (m == control$x.max.ext){
      WARN <- 'Too small age interval found. Age classes were rounded.'
      x <- round(x * control$x.max.ext) / control$x.max.ext
      m <- pclm.interval.multiple(x)
      warning(immediate. = TRUE, WARN)
    } else {
      if (m > 1) WARN <- 'Age vector automatically transformed.'
    }
    x <- as.integer(round(x * m * control$x.div)) / control$x.div #transform
    #x <- as.integer(round(x*m)) #OLD, but x is anyway multiplied by x.div in pclm.compmat
  } else m <- 1

  # Calculate composition matrix
  CM <- pclm.compmat(x, y, exposures = exposures, control = control)
  #control <- CM$control
  #CM$control <- NULL

  # Fit PCLM model
  fit <- pclm.opt(CompositionMatrix = CM, control = control)

  # Change age/time units back to fractional (backward transformation)
  fit$CompositionMatrix$x <- fit$CompositionMatrix$x / m
  fit$X <- fit$X / m
  fit$CompositionMatrix$open.int.len <- fit$CompositionMatrix$open.int.len / m

  # Construct grouped (aggregated) and ungrouped (nonaggregated) life-tables
  GLT <- pclm.aggregate(fit, out.step = out.step, count.type = count.type[1], exposures.used = exposures.used)

  # Construct pclm object
  GLT$m <- m
  GLT$x.div <- control$x.div
  GLT$out.step <- out.step
  GLT$warn.list <- c(CM$warn.list, GLT$warn.list, WARN)
  GLT$exposures <- exposures
  GLT$control <- control
  class(GLT) <- 'pclm'
  GLT
}

# summary.pclm -------------------------------------------------------------------------

#' Summary of the Fitted PCLM Object
#'
#' @param object Fitted PCLM object.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @seealso \code{\link{pclm.fit}} \code{\link{plot.pclm}}
#' @keywords internal
#' @export
summary.pclm <- function(object){
  if (!inherits(object, 'pclm')) {
    if (inherits(object, 'pash')) pash:::summary.pash(object) else stop ('Object of class pclm needed')
  } else {
    if (inherits(object, 'pash')) {
      message('Summary of the pash object:')
      pash:::summary.pash(object)
      object <- object$pclm
      cat('\n\n')
    }
    message('Summary of the pclm object:')
    n1 <- diff(object$fit$CompositionMatrix$x)
    n1 <- c(n1, n1[length(n1)])
    z0 <- n1[1]/2
    n2 <- diff(object$fit$X)
    n2 <- c(n2, n2[length(n2)])
    if (!is.na(object$fit$CompositionMatrix$open.int.len)) ind <- -(length(n1):(length(n1) - 1)) else ind <- 1:length(n1)
    cat(paste(paste('PCLM total classes =', length(n2)),
      paste('Number of smoothing parameters for B-/P-splines =', object$fit$control$bs.df),
      paste('Original minimal interval length =', round(min(n1[ind]), 3)),
      paste('Original maximal interval length =', round(max(n1[ind]), 3)),
      paste('Open interval length =', round(object$fit$CompositionMatrix$open.int.len, 3)),
      paste('Fractional age/time class correction (multiple) =', object$m),
      paste('PCLM interval length =', round(min(n2), 3)),
      paste('PCLM class divider (x.div) =', object$x.div),
      paste('PCLM classes per original smallest interval length =', round(min(n1[ind]) / min(n2), 3)),
      paste('PCLM classes per original biggest interval length =', round(max(n1[ind]) / max(n2), 3)), sep='\n'))
    message('\nWarnings list:')
    W <- unlist(object$warn.list)
    print(W,  quote=F)
    cat('\n')
  }
  invisible()
}

# plot.pclm ---------------------------------------------------------------------------

#' Diagnostic Plot for PCLM Object
#'
#' @param object Fitted PCLM object.
#' @param type Type of PCLM plot: \itemize{ \item{\code{"aggregated"} -
#'   Aggregated PCLM fit with interval length of \code{out.step}}. See
#'   \code{\link{pclm.fit}}. \item{\code{"nonaggregated"} - Nonaggregated (raw)
#'   PCLM fit with interval of length equal to the shortest original interval
#'   length divided by \code{x.div}. See \code{\link{pclm.control}}}. }
#' @param xlab Optional label for the X-axis.
#' @param ylab Optional label for the Y-axis.
#' @param xlim Optional limits of the X-axis.
#' @param ylim Optional limits of the Y-axis.
#' @param legpos.x,legpos.y Position of the \code{\link{legend}}. If
#'   \code{legpos.x == NULL} then legend is not plotted.
#'
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
#' @seealso \code{\link{pclm.fit}} \code{\link{summary.pclm}}
#' @keywords internal
#' @export
plot.pclm<-function(object, type = c("aggregated", "nonaggregated"), xlab, ylab, xlim, ylim, legpos.x = "topleft", legpos.y = NULL){
  if (missing(xlab)) if (inherits(object, 'pash')) xlab <- attributes(object)$time_unit else xlab <- 'Age or time'
  if (!inherits(object, 'pclm')) {
    if (inherits(object, 'pash')) stop('Plot function for pash object without computed pclm is not supported.') else stop ('Object of class pclm needed')
  } else {
    if (inherits(object, 'pash')) object <- object$pclm
  }
 if (length(object$exposures) == 0){
   if (missing(ylab)) ylab <- 'Counts / interval length'

   n1 <- diff(object$fit$CompositionMatrix$x)
   n1 <- c(n1, n1[length(n1)])
   n2 <- diff(object$fit$X)
   n2 <- c(n2, n2[length(n2)])
   n3 <- diff(object$grouped$x)
   n3 <- c(n3, n3[length(n3)])

   if (missing(xlim)) xlim <- range(c(object$fit$X, object$fit$CompositionMatrix$x))
   if (missing(ylim)) {
     ylim <- range(c(object$fit$Y/n2, object$fit$CompositionMatrix$y/n1, object$grouped$dx/n3))
     ylim[1] <- 0
     ylim[2] <- ylim[2]*1.2
   }
   tmp.lwd <- par('lwd'); par(lwd = 2,xaxs = 'i', yaxs = 'i')
   barplot(width = n1, space = 0, height = object$fit$CompositionMatrix$y / n1, xlab = xlab, ylab = ylab,
           col = 'gold2', border = 'white', xlim = xlim, ylim = ylim)
   par(lwd = tmp.lwd)

   lines(object$fit$CompositionMatrix$x, object$fit$CompositionMatrix$y/n1, type = 's')
   axis(1)
   if (tolower(type[1]) == 'nonaggregated'){
     lines(object$fit$X, object$fit$Y/n2, type = 's', col = 'red', lwd = 2)
     AType <- 'PCLM nonaggregated (raw)'
   } else if (tolower(type[1]) == 'aggregated') {
     AType <- 'PCLM aggregated'
     lines(object$grouped$x, object$grouped$dx/n3, type = 's', col = 'red', lwd = 2)
   } else stop('Unknown plotting type.')
   if (length(legpos.x) > 0) legend(x = legpos.x, y = legpos.y, legend = c('Data', AType), bty = 'n', pch = c(15, NA), lty = c(NA, 1), lwd = 2, col = c('gold2', 'red'), pt.cex = 1.8)
   box()
 } else stop('Diagnostic plots for mortality smooth not available yet.')
 invisible()
}

# head.pclm ----------------------------------------------------

#' Head of PCLM Object
#'
#' @export
#' @param object PCLM object.
#' @param n A single integer. If positive, size for the resulting object: number
#'   of rows for a life-table. If negative, all but the n last/first number of
#'   elements of x.
#' @param type which life-table  should be returned. One of \code{"lt"},
#'   \code{"aggregated"} or \code{"nonaggregated"}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
head.pclm<-function(object, n = 6L, type = c("lt", "aggregated", "nonaggregated")){
  if (!inherits(object, "pash")) {
    if (type == "lt") stop('"lt" not supported for non-pash object.')
    object. <- NULL
    object.$pclm <- object
    object <- object.
  }
  type <- type[1]
  if (type == "lt")
    head(object$lt, n = n)
  else if (type == "aggregated")
    head(object$pclm$grouped, n = n)
  else if (type == "nonaggregated")
    head(object$pclm$raw, n = n)
  else stop('Unknown type')
}

# tail.pclm -----------------------------------------------------------

#' Tail of PCLM Object
#'
#' @export
#' @param object PCLM object.
#' @param n A single integer. If positive, size for the resulting object: number
#'   of rows for a life-table. If negative, all but the n last/first number of
#'   elements of x.
#' @param type which life-table  should be returned. One of \code{"lt"},
#'   \code{"aggregated"} or \code{"nonaggregated"}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
tail.pclm<-function(object, n=6L, type = c("lt", "aggregated", "nonaggregated")){
  if (!inherits(object, "pash")) {
    if (type == "lt") stop('"lt" not supported for non-pash object.')
    object. <- NULL
    object.$pclm <- object
    object <- object.
  }
  type <- type[1]
  if (type == "lt")
    tail(object$lt, n = n)
  else if (type == "aggregated")
    tail(object$pclm$grouped, n = n)
  else if (type == "nonaggregated")
    tail(object$pclm$raw, n = n)
  else stop('Unknown type')
}

# print.pclm ------------------------------------------------------------

#' Print PCLM Object
#'
#' @export
#' @param object PCLM object.
#' @param type which life-table  should be returned. One of \code{"lt"},
#'   \code{"aggregated"} or \code{"nonaggregated"}.
#' @param ... other parameters passed to \code{\link{print}}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}>
#'   <\email{maciej.danko@gmail.com}>
print.pclm<-function(object, type = c("lt", "aggregated", "nonaggregated"), ...){
  if (!inherits(object, "pash")) {
    if (type == "lt") stop('"lt" not supported for non-pash object.')
    object. <- NULL
    object.$pclm <- object
    object <- object.
  }
  type <- type[1]
  if (type == "lt")
    pash:::print.pash(object, ...)
  else if (type == "aggregated")
    print(object$pclm$grouped, ...)
  else if (type == "nonaggregated")
    print(object$pclm$raw, ...)
  else stop('Unknown type')
}
