#' Left-Truncated Gamma Distributions
#'
#' Density, distribution function, quantile function, random number generation, expected value, and variance for left-truncated gamma distributions (or right-truncated inverse gamma distributions).
#' @param x,q,p Points at which to evaluate the density (\code{dltrgamma}), distribution (\code{pltrgamma}), and quantile (\code{qltrgamma}) functions, respectively.
#' @param shape Shape parameter for the desired gamma distribution. Must be strictly positive and must exceed \code{1} and \code{2} under \code{exp_ltrgamma} and \code{var_ltrgamma}, respectively, when \code{inverse=TRUE}.
#' @param rate,scale Rate parameter for the desired gamma distribution. Must be strictly positive. Defaults to \code{1}. \code{scale} is an alternative way to specify the \code{rate} and defaults to \code{1/rate}. It is an error to supply both \code{scale} and \code{rate}.
#' @param trunc The point of truncation (corresponding to \eqn{\tau}{tau} below). Must be non-negative. The default \code{ifelse(inverse, Inf, 0)} corresponds to no left-truncation for gamma distributions and no right-truncation for inverse gamma distributions. See Details below.
#' @param inverse A logical indicating whether to use a gamma distribution left-truncated to \eqn{\left\lbrack\tau,\infty\right)} or an inverse gamma distribution right-truncated to \eqn{\left(0, \tau\right\rbrack}. Defaults to \code{FALSE}.
#' @param log,log.p Logicals indicating whether densities/probabilities should be returned on the log scale for \code{dltrgamma} (\code{log}) and \code{pltrgamma} and \code{qltrgamma} (both \code{log.p}). Both default to \code{FALSE}. See Details below.
#' @param lower.tail Logical indicating whether probabilities are \eqn{\Pr\left\lbrack X \le x\right\rbrack} when \code{TRUE} (the default) or \eqn{\Pr\left\lbrack X > x\right\rbrack} for \code{pltrgamma} and \code{qltrgamma}.
#'
#' @details The left-truncated gamma distribution has PDF:
#' \deqn{f(x|\alpha, \beta, \tau) = \frac{\beta^\alpha}{\left(\Gamma\left(\alpha\right)-\Gamma\left(\alpha, \tau\beta\right)\right)}x^{\alpha-1}e^{-x\beta}\mathbf{1}\left(x \ge \tau\right)}
#' for \eqn{0\le\tau\le x < \infty}, where \eqn{\alpha > 0}{alpha} and \eqn{\beta > 0}{beta} are the \code{shape} and \code{rate} parameters, respectively, \eqn{\tau}{tau} is the cutoff point at which left \code{trunc}ation occurs, \eqn{\mathbf{1}\left(\cdot\right)} is the usual indicator function, \eqn{\Gamma(\alpha)} is the \code{\link{gamma}} function, and \eqn{\Gamma(\alpha, \tau\beta)} is the upper incomplete gamma function.
#'
#' When \code{inverse=TRUE}, \eqn{\tau} (i.e. \code{trunc}) becomes the point of \emph{right} truncation for an inverse gamma distribution, such that \eqn{0 < x \le \tau < \infty}.
#'
#' All code is based purely on functions in \code{\link{GammaDist}}. The R code of Nadarajah and Kotz (2006) has been modified to work with \code{log=TRUE}/\code{log.p=TRUE}, &/or \code{inverse=TRUE}, and/or \code{lower.tail=FALSE}, and prevent the invalid passing of \code{log.p} and/or \code{lower.tail} arguments through the various \code{\link{GammaDist}} and \code{\link{ltrgamma}} functions called \emph{within} these functions.
#'
#' Furthermore, internal calculations for the functions \code{dltrgamma}, \code{exp_ltrgamma}, \code{var_ltrgamma}, and \code{pltrgamma} (where possible) are always \emph{initially} performed on the log scale, for numerical stability, regardless of the value(s) of \code{log} and \code{log.p}.
#'
#' Additionally, the code also enables recovery of non-truncated gamma distributions (set \code{trunc=0} and \code{inverse=FALSE}) and non-truncated inverse gamma distributions (set \code{trunc=Inf} and \code{inverse=TRUE}). In fact, the \emph{default} behaviour for \code{trunc} is already such that no truncation is performed.
#'
#' Finally, the code avails of closed-form expressions for both \code{exp_ltrgamma} and \code{var_ltrgamma}, rather than relying on numerical integration as per Nadarajah and Kotz (2006), both for \code{inverse=TRUE} and \code{inverse=FALSE}, thereby enabling proper recycling of numeric arguments throughout.
#'
#' @note \code{rltrgamma} is invoked internally for the \code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"} models to draw column shrinkage parameters for all but the first loadings column under the MGP prior when \code{truncated=TRUE} (which is \strong{not} the default) is supplied to \code{\link{mgpControl}}, at the expense of slightly longer run times. \code{exp_ltrgamma} is used within \code{\link{MGP_check}} to check the validity of the MGP hyperparameters when \code{truncated=TRUE} (which is again, \strong{not} the default). Both functions always assume \code{trunc=1} for these internal usages.
#'
#' @return \code{dltrgamma} gives the density, \code{pltrgamma} gives the distribution function, \code{qltrgamma} gives the quantile function, \code{rltrgamma} generates random deviates, \code{exp_ltrgamma} returns the expected value, and \code{var_ltrgamma} returns the variance for the left-truncated gamma distribution (or right-truncated inverse gamma distribution) with the specified \code{shape} and \code{rate}/\code{scale} parameters, and truncation point \code{trunc}.
#'
#' Values in \code{x} and \code{q} are allowed to be less than \code{trunc} (or greater than \code{trunc} when \code{inverse=TRUE}) with densities and log densities of \code{0} and \code{-Inf} respectively returned by \code{dltrgamma} and appropriate values of \code{-Inf}, \code{0}, or \code{1} returned by \code{pltrgamma} depending on the values of \code{inverse}, \code{log.p}, and \code{lower.tail}. Otherwise, invalid arguments will result in return value \code{NaN}, with a warning.
#'
#' The length of the result is determined by \code{n} for \code{rltrgamma}, and is the maximum of the lengths of the numerical arguments for the other functions as the numerical arguments other than \code{n} are recycled to the length of the result. Invalid recycling of the numeric arguments will yield a warning. Only the first elements of the logical arguments are used.
#'
#' @seealso \code{\link{GammaDist}}, \code{\link{mgpControl}}, \code{\link{MGP_check}}
#' @name ltrgamma
#' @rdname ltrgamma
#' @references Dagpunar, J. S. (1978) Sampling of variates from a truncated gamma distribution, \emph{Statistical Computation and Simulation}, 8(1): 59-64.
#'
#' Nadarajah, S. and Kotz, S. (2006) R Programs for Computing Truncated Distributions, \emph{Journal of Statistical Software, Code Snippets}, 16(2): 1–8.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @examples
#' s <- 3.1
#' r <- 2.1
#' t <- 1
#' # Generate left-truncated Ga(3.1, 2.1, 1) variates
#' rltrgamma(n=10, shape=s, rate=r, trunc=t)
#'
#' # Plot the density of a left-truncated Ga(3.1, 2.1, 1) distribution
#' curve(dltrgamma(x, s, r, t), from=0, to=5)
#' # Overlay the corresponding non-truncated Ga(3.1, 2.1) density
#' curve(dgamma(x, s, r), lty=2, add=TRUE)
#'
#' # Check validity of density
#' integrate(dltrgamma, lower=t, upper=Inf,
#'           shape=s, scale=r, trunc=t)
#'
#' # Calculate the expectation and variance of a Ga(3.1, 2.1, 1) distribution
#' exp_ltrgamma(shape=s, rate=r, trunc=t)
#' var_ltrgamma(shape=s, rate=r, trunc=t)
#'
#' # Compare expectation and variance with a similar non-truncated gamma distribution
#' exp_ltrgamma(shape=s, rate=r, trunc=0) # s/r
#' var_ltrgamma(shape=s, rate=r, trunc=0) # s/r^2
#'
#' # Repeat the above for a right-truncated inverse gamma IG(2.2, 3.3, 4.4) distribution
#' s <- 2.2
#' r <- 3.3
#' t <- 4.4
#' rltrgamma(n=10, shape=s, rate=r, trunc=t, inverse=TRUE)
#' curve(dltrgamma(x, s, r, t, inverse=TRUE), from=0, to=6.6)
#' curve(dgamma(1/x, s, r)/x^2, lty=2, add=TRUE)
#' integrate(dltrgamma, lower=0, upper=t,
#'           shape=s, rate=r, trunc=t, inverse=TRUE)
#' exp_ltrgamma(shape=s, rate=r, trunc=t, inverse=TRUE)
#' var_ltrgamma(shape=s, rate=r, trunc=t, inverse=TRUE)
#' exp_ltrgamma(shape=s, rate=r, trunc=Inf, inverse=TRUE) # r/(s - 1)
#' var_ltrgamma(shape=s, rate=r, trunc=Inf, inverse=TRUE) # (r)^2/((s - 1)^2 * (s - 2))
#'
#' # Check distribution and quantile functions
#' p <- (1:9)/10
#' t <- 1:3
#' all.equal(p, pltrgamma(qltrgamma(p, shape=2, scale=2, trunc=t),
#'                        shape=2, scale=2, trunc=t))
#' all.equal(p, pltrgamma(qltrgamma(p, shape=2, scale=2, trunc=t, inverse=TRUE),
#'                        shape=2, scale=2, trunc=t, inverse=TRUE))
NULL

#' @keywords utility
#' @rdname ltrgamma
#' @usage
#' dltrgamma(x,
#'           shape,
#'           rate = 1,
#'           trunc = ifelse(inverse, Inf, 0),
#'           inverse = FALSE,
#'           log = FALSE,
#'           scale = 1/rate)
#' @export
dltrgamma    <- function(x, shape, rate = 1, trunc = ifelse(inverse, Inf, 0), inverse = FALSE, log = FALSE, scale = 1/rate) {
  if(!is.numeric(x))                   stop("'x' must be a numeric vector",                      call.=FALSE)
  if(any(!is.numeric(shape),
         shape       <= 0))            stop("'shape' must be strictly positive and numeric",     call.=FALSE)
  if(missing(rate))   {
    rate     <- if(!missing(scale)) 1/scale else rate
  } else if(!missing(scale))           stop("Specify 'rate' or 'scale' but not both",            call.=FALSE)
  if(any(!is.numeric(rate),
         rate        <= 0))            stop("'rate' must be strictly positive and numeric",      call.=FALSE)
  if(any(!is.numeric(trunc),
         trunc        < 0))            stop("'trunc' must be strictly non-negative and numeric", call.=FALSE)
  dens       <- x
  if(isTRUE(inverse)) {
    dx       <- x <= trunc & x >= 0
    dens[dx] <- stats::dgamma(1/x[dx], shape, rate, log=TRUE) - 2 * log(x[dx]) - stats::pgamma(1/trunc, shape=shape, rate=rate, lower.tail=FALSE, log.p=TRUE)
  } else      {
    dx       <- x >= trunc
    dens[dx] <- stats::dgamma(x[dx],   shape, rate, log=TRUE) - log1p(-stats::pgamma(trunc, shape, rate))
  }
  dens[!dx]  <- if(isTRUE(log)) -Inf     else 0
  dens[dx]   <- if(isTRUE(log)) dens[dx] else exp(dens[dx])
    return(replace(dens, is.na(dens), ifelse(isTRUE(log), -Inf, 0)))
}

#' @keywords utility
#' @rdname ltrgamma
#' @usage
#' pltrgamma(q,
#'           shape,
#'           rate = 1,
#'           trunc = ifelse(inverse, Inf, 0),
#'           inverse = FALSE,
#'           log.p = FALSE,
#'           lower.tail = TRUE,
#'           scale = 1/rate)
#' @export
pltrgamma    <- function(q, shape, rate = 1, trunc = ifelse(inverse, Inf, 0), inverse = FALSE, log.p = FALSE, lower.tail = TRUE, scale = 1/rate) {
  if(!is.numeric(q))                   stop("'q' must be a numeric vector",                      call.=FALSE)
  if(any(!is.numeric(shape),
         shape       <= 0))            stop("'shape' must be strictly positive and numeric",     call.=FALSE)
  if(missing(rate))   {
    rate     <- if(!missing(scale)) 1/scale else rate
  } else if(!missing(scale))           stop("Specify 'rate' or 'scale' but not both",            call.=FALSE)
  if(any(!is.numeric(rate),
         rate        <= 0))            stop("'rate' must be strictly positive and numeric",      call.=FALSE)
  if(any(!is.numeric(trunc),
         trunc        < 0))            stop("'trunc' must be strictly non-negative and numeric", call.=FALSE)
  if(isTRUE(inverse)) {
    q        <- stats::pgamma(1/pmax(.Machine$double.eps,
                                     pmin(q, trunc)), shape=shape, rate=rate, log.p=isTRUE(lower.tail),  lower.tail=FALSE)
    G        <- stats::pgamma(1/trunc,                shape=shape, rate=rate, log.p=isTRUE(lower.tail),  lower.tail=FALSE)
  } else      {
    q        <- stats::pgamma(pmax(q, trunc),         shape=shape, rate=rate, log.p=isFALSE(lower.tail), lower.tail=isTRUE(lower.tail))
    G        <- stats::pgamma(trunc,                  shape=shape, rate=rate, log.p=isFALSE(lower.tail), lower.tail=isTRUE(lower.tail))
  }
  if(isTRUE(inverse) && isFALSE(lower.tail)) {
    q        <- if(isTRUE(log.p)) log(G - q) - log(G)    else (G - q)/G
  } else if(isFALSE(inverse)   && isTRUE(lower.tail)) {
    q        <- if(isTRUE(log.p)) log(q - G) - log1p(-G) else (q - G)/(1 - G)
  } else q   <- if(isTRUE(log.p)) q     - G              else exp(q - G)
    return(replace(q, is.na(q), ifelse(isTRUE(log.p),
                                       ifelse(isFALSE(lower.tail), 0, -Inf),
                                       isFALSE(lower.tail))))
}

#' @keywords utility
#' @rdname ltrgamma
#' @usage
#' qltrgamma(p,
#'           shape,
#'           rate = 1,
#'           trunc = ifelse(inverse, Inf, 0),
#'           inverse = FALSE,
#'           log.p = FALSE,
#'           lower.tail = TRUE,
#'           scale = 1/rate)
#' @export
qltrgamma    <- function(p, shape, rate = 1, trunc = ifelse(inverse, Inf, 0), inverse = FALSE, log.p = FALSE, lower.tail = TRUE, scale = 1/rate) {
  if(!is.numeric(p))                   stop("'p' must be a numeric vector",                      call.=FALSE)
  if(any(!is.numeric(shape),
         shape       <= 0))            stop("'shape' must be strictly positive and numeric",     call.=FALSE)
  if(missing(rate))   {
    rate     <- if(!missing(scale)) 1/scale else rate
  } else if(!missing(scale))           stop("Specify 'rate' or 'scale' but not both",            call.=FALSE)
  if(any(!is.numeric(rate),
         rate        <= 0))            stop("'rate' must be strictly positive and numeric",      call.=FALSE)
  if(any(!is.numeric(trunc),
         trunc        < 0))            stop("'trunc' must be strictly non-negative and numeric", call.=FALSE)
  if(isTRUE(inverse)) {
    G        <- stats::pgamma(1/trunc, shape, rate, lower.tail=FALSE)
      1/stats::qgamma(1   - (p * G),   shape, rate, lower.tail=isTRUE(lower.tail), log.p=isTRUE(log.p))
  } else      {
    G        <- stats::pgamma(trunc,   shape, rate, lower.tail=TRUE)
      stats::qgamma(G + p * (1 - G),   shape, rate, lower.tail=isTRUE(lower.tail), log.p=isTRUE(log.p))
  }
}

#' @keywords utility
#' @param n Number of observations to generate for \code{rltrgamma}.
#' @rdname ltrgamma
#' @usage
#' rltrgamma(n,
#'           shape,
#'           rate = 1,
#'           trunc = ifelse(inverse, Inf, 0),
#'           inverse = FALSE,
#'           scale = 1/rate)
#' @export
rltrgamma    <- function(n, shape, rate = 1, trunc = ifelse(inverse, Inf, 0), inverse = FALSE, scale = 1/rate) {
  if(!is.numeric(n)  ||
     n       <= 0    ||
     n       != floor(n))              stop("'n' must be a strictly positive integer",           call.=FALSE)
  if(any(!is.numeric(shape),
         shape       <= 0))            stop("'shape' must be strictly positive and numeric",     call.=FALSE)
  if(missing(rate))   {
    rate     <- if(!missing(scale)) 1/scale else rate
  } else if(!missing(scale))           stop("Specify 'rate' or 'scale' but not both",            call.=FALSE)
  if(any(!is.numeric(rate),
         rate        <= 0))            stop("'rate' must be strictly positive and numeric",      call.=FALSE)
  if(any(!is.numeric(trunc),
         trunc        < 0))            stop("'trunc' must be strictly non-negative and numeric", call.=FALSE)
    qltrgamma(stats::runif(n), shape=shape, rate=rate, trunc=trunc, inverse=isTRUE(inverse))
}

#' @keywords utility
#' @rdname ltrgamma
#' @usage
#' exp_ltrgamma(shape,
#'              rate = 1,
#'              trunc = ifelse(inverse, Inf, 0),
#'              inverse = FALSE,
#'              scale = 1/rate)
#' @export
exp_ltrgamma <- function(shape, rate = 1, trunc = ifelse(inverse, Inf, 0), inverse = FALSE, scale = 1/rate) {
  if(any(!is.numeric(shape),
         shape       <= 0))            stop("'shape' must be strictly positive and numeric",     call.=FALSE)
  if(missing(rate))   {
    rate     <- if(!missing(scale)) 1/scale else rate
  } else if(!missing(scale))           stop("Specify 'rate' or 'scale' but not both",            call.=FALSE)
  if(any(!is.numeric(rate),
         rate        <= 0))            stop("'rate' must be strictly positive and numeric",      call.=FALSE)
  if(any(!is.numeric(trunc),
         trunc        < 0))            stop("'trunc' must be strictly non-negative and numeric", call.=FALSE)
  if(isTRUE(inverse)) {
    if(any(shape     <= 1))            stop("'shape' must be > 1 when 'inverse'=TRUE",           call.=FALSE)
    log_exp  <- stats::pgamma(rate / trunc, shape - 1, lower.tail=FALSE, log.p=TRUE) - stats::pgamma(rate / trunc, shape, lower.tail=FALSE, log.p=TRUE) + log(rate)  - log(shape - 1)
    log_exp  <- replace(log_exp, !is.finite(trunc), log(rate) - log(shape - 1))
  } else      {
    log_exp  <- stats::pgamma(rate * trunc, shape + 1, lower.tail=FALSE, log.p=TRUE) - stats::pgamma(rate * trunc, shape, lower.tail=FALSE, log.p=TRUE) + log(shape) - log(rate)
  }
    return(exp(log_exp))
}

#' @keywords utility
#' @rdname ltrgamma
#' @usage
#' var_ltrgamma(shape,
#'              rate = 1,
#'              trunc = ifelse(inverse, Inf, 0),
#'              inverse = FALSE,
#'              scale = 1/rate)
#' @export
var_ltrgamma <- function(shape, rate = 1, trunc = ifelse(inverse, Inf, 0), inverse = FALSE, scale = 1/rate) {
  if(any(!is.numeric(shape),
         shape       <= 0))            stop("'shape' must be strictly positive and numeric",     call.=FALSE)
  if(missing(rate))   {
    rate     <- if(!missing(scale)) 1/scale else rate
  } else if(!missing(scale))           stop("Specify 'rate' or 'scale' but not both",            call.=FALSE)
  if(any(!is.numeric(rate),
         rate        <= 0))            stop("'rate' must be strictly positive and numeric",      call.=FALSE)
  if(any(!is.numeric(trunc),
         trunc        < 0))            stop("'trunc' must be strictly non-negative and numeric", call.=FALSE)
  if(isTRUE(inverse) &&
     any(shape       <= 2))            stop("'shape' must be > 2 when 'inverse'=TRUE",           call.=FALSE)
  mu         <- exp_ltrgamma(shape=shape, rate=rate, trunc=trunc, inverse=isTRUE(inverse))
  if(isTRUE(inverse)) {
    lEX2     <- stats::pgamma(rate / trunc, shape - 2, lower.tail=FALSE, log.p=TRUE) - stats::pgamma(rate / trunc, shape, lower.tail=FALSE, log.p=TRUE) + 2 * log(rate) + lgamma(shape - 2) - lgamma(shape)
  } else      {
    lEX2     <- stats::pgamma(rate * trunc, shape + 2, lower.tail=FALSE, log.p=TRUE) - stats::pgamma(rate * trunc, shape, lower.tail=FALSE, log.p=TRUE) - 2 * log(rate) + log1p(shape) + log(shape)
  }
    exp(lEX2) - mu^2
}
