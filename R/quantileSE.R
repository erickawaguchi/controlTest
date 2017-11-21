#' Find standard error for survival quantile
#'
#' @param timevar Vector of observed survival times.
#' @param censor  Vector of censoring indicators (1 = uncensored, 0 = censored).
#' @param q       Quantile of interest (Default is median).
#' @param B       Number of bootstrap samples.
#' @param alpha   Significance level for confidence interval of quantile.
#' @param seed    Seed number (for reproducibility).
#' @param plots   Logical. TRUE to show Kaplan-Meier plot
#' @return Returns quantile estimate, bootstrapped standard error, and (1 - alpha / 2) * 100% confidence interval for quantile estimate.
#' @examples
#' #Reference: Survival Analysis Techniques for Censored and Truncated Data.
#' #Klein and Moeschberger (1997) Springer.
#' #Data: Chapter 7.6 Example 7.9 (p. 211)
#' library(controlTest)
#' t1 <- c(1, 63, 105, 129, 182, 216, 250, 262, 301, 301,
#'        342, 354, 356, 358, 380, 383, 383, 338, 394, 408, 460, 489,
#'        499, 523, 524, 535, 562, 569, 675, 676, 748, 778, 786, 797,
#'        955, 968, 1000, 1245, 1271, 1420, 1551, 1694, 2363, 2754, 2950)
#' c1 <- c(rep(1, 43), 0, 0)
#' quantileSE(timevar = t1, censor = c1, q = 0.5, B = 500)
#' @export
#' @import survival
#' @importFrom graphics legend lines plot
#' @importFrom stats pnorm qnorm quantile sd
#'
quantileSE <- function(timevar, censor, q = .5, B = 1000, alpha = 0.05, seed = 1991, plots = FALSE) {

  #- Checking for silly errors
  if (q < 0 | q > 1 ) {
    stop("q should be between 0 and 1.")
  }
  if (B <= 0) {
    stop("B should be a positive integer")
  }
  if (seed <= 0) {
    stop("Seed should be a positive integer")
  }
  if (length(timevar) != length(censor)) {
    stop("Length of timevar and censor should be the same")
  }

  set.seed(seed)
  fit      <- survfit(Surv(timevar, censor) ~ 1, conf.type = 'none') # Control
  q.est    <- unname(quantile(fit, prob = 1 - q))
  if(is.na(q.est)){
    stop(paste("Estimated survival time for ", q*100, "th quantile was not found. Program Stopped."))
  }

  Qp <- function(t, c, q){
    fit1 <- survfit(Surv(t, c) ~ 1, conf.type = 'none')
    Finv <- unname(quantile(fit1, prob = 1 - q))
    if(is.na(Finv)){
      warning(paste("Estimated survival time for ", q*100, "th quantile could not be estimated for sample. Largest observed survival time was used as an estimate."))
      Finv <- max(t)
    }
    return(Finv)
  }

  quant_est <- numeric(B)
  for(i in 1:B){
    btsp         <- sample(c(1:length(timevar)), replace = TRUE)
    tmp_time     <- timevar[btsp]
    tmp_censor   <- censor[btsp]
    quant_est[i] <- Qp(tmp_time, tmp_censor, q = q) #Bootstrapped KM
  }

  if(plots == TRUE) {
    plot(fit, col = 'red', ylab = 'Estimated Survival Function', xlab = 'Time', main = 'Kaplan-Meier Estimate')
    lines(x = c(q.est, q.est), y = c(-.5, q), lty = 2)
    lines(x = c(0, q.est), y = c(q, q), lty = 2)
  }

  se        <- sd(quant_est)

  #- Output
  out <- list()
  out$quantile <- q
  out$estimate <- q.est
  out$se       <- se
  out$lower    <- round(q.est - qnorm(1 - alpha / 2) * se, 2)
  out$upper    <- round(q.est + qnorm(1 - alpha / 2) * se, 2)
  return(out)
}

