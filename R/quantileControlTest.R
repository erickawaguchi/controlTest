#' Find standard error for survival quantile
#'
#' @param timevar1 Vector of observed survival times for sample 1 (control).
#' @param censor1  Vector of censoring indicators for sample 1 (1 = uncensored, 0 = censored).
#' @param timevar2 Vector of observed survival times for sample 2 (treatment).
#' @param censor2  Vector of censoring indicators for sample 2 (1 = uncensored, 0 = censored).
#' @param q       Quantile of interest (Default is median).
#' @param B       Number of bootstrap samples.
#' @param seed    Seed number (for reproducibility).
#' @param plots   Logical. TRUE to show Kaplan-Meier plot
#' @return Returns quantile estimate, bootstrapped standard error, test statistic, and two-sided p-value.
#' @examples
#' #Reference: Survival Analysis Techniques for Censored and Truncated Data.
#' #Klein and Moeschberger (1997) Springer.
#' #Data: Chapter 7.6 Example 7.9 (p. 211)
#' library(controlTest)
#' t1 <- c(1, 63, 105, 129, 182, 216, 250, 262, 301, 301,
#'        342, 354, 356, 358, 380, 383, 383, 338, 394, 408, 460, 489,
#'        499, 523, 524, 535, 562, 569, 675, 676, 748, 778, 786, 797,
#'        955, 968, 1000, 1245, 1271, 1420, 1551, 1694, 2363, 2754, 2950)
#' t2 <- c(17, 42, 44, 48, 60, 72, 74, 95, 103, 108, 122, 144, 167, 170,
#'        183, 185, 193, 195, 197, 208, 234, 235, 254, 307, 315, 401, 445,
#'        464, 484, 528, 542, 547, 577, 580, 795, 855, 1366, 1577, 2060,
#'        2412, 2486, 2796, 2802, 2934, 2988)
#' c1 <- c(rep(1, 43), 0, 0)
#' c2 <- c(rep(1, 39), rep(0, 6))
#' quantileControlTest(t1, c1, t2, c2, q = 0.5, B = 500)
#'
#'@details It is important to note the possiblilty that the estimated quantile may not be estimable in our bootstrap samples. In such cases
#' the largest observed survival time will be considered as an estimate for the quantile.
#'
#' @references
#' Li, G., Tiwari, R.C., and Wells, M. (1996). "Quantile Comparison Functions in Two-Sample Problems: With Applications to Comparisons of Diagnostic Markers." Journal of the American Statistical Association, 91, 689-698.
#'
#' Chakraborti, S., and Mukerjee, R. (1989), "A Confidence Interval for a Measure Associated With the Comparison of a Treatment With a Control," South African Statistical Journal, 23, 219-230.
#'
#' Gastwirth, J. L., and Wang, J. L. (1988), "Control Percentile Test for Censored Data," Journal of Statistical Planning and Inference, 18, 267-276.
#' @export
#' @import survival
#' @importFrom graphics legend lines plot
#' @importFrom stats pnorm qnorm quantile sd
#'
quantileControlTest <- function(timevar1, censor1, timevar2, censor2, q = 0.5, B = 1000, seed = 1234, plots = FALSE) {

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
  if (length(timevar1) != length(censor1)) {
    stop("Length of timevar1 and censor1 should be the same")
  }
  if (length(timevar2) != length(censor2)) {
    stop("Length of timevar2 and censor2 should be the same")
  }

  set.seed(seed)
  fit1 <- survfit(Surv(timevar1, censor1) ~ 1, conf.type = "none")
  fit2 <- survfit(Surv(timevar2, censor2) ~ 1, conf.type = "none")
  F.inv <- unname(quantile(fit1, prob = 1 - q))
  G.inv <- unname(quantile(fit2, prob = 1 - q))

  if(is.na(F.inv)) {
    stop(paste0("Estimated survival time for ", q*100, "-th quantile for sample 1 (control) not found. Program Stopped."))
  } else if(is.na(G.inv)) {
    stop(paste0("Estimated survival time for ", q*100, "-th quantile for sample 2 (treatment) not found. Program Stopped."))
  }


  Qp <- function(t1, c1, t2, c2) {
    fit1 <- survfit(Surv(t1, c1) ~ 1, conf.type = "none")
    fit2 <- survfit(Surv(t2, c2) ~ 1, conf.type = "none")
    F.inv <- unname(quantile(fit1, prob = 1 - q))
    if (is.na(F.inv)) {
      warning(paste0("Estimated survival time for ", q*100, "-th quantile could not be estimated for bootstrap sample. Largest observed survival time was used instead."))
      F.inv <- max(t1)
    }
    t0 <- which.max(fit2$time[fit2$time <= F.inv])
    Q  <- fit2$surv[t0]
    #G <- stepfun(fit2$time, c(1, fit2$surv))
    #Q <- G(F.inv)
    return(Q)
  }

  Q <- Qp(timevar1, censor1, timevar2, censor2)

  # Bootstrap
  b.est <- numeric(B)
  for (i in 1:B) {
    boot1    <- sample(1:length(timevar1), replace = TRUE)
    t1.boot  <- timevar1[boot1]
    c1.boot  <- censor1[boot1]
    boot2    <- sample(1:length(timevar2), replace = TRUE)
    t2.boot  <- timevar2[boot2]
    c2.boot  <- censor2[boot2]
    b.est[i] <- Qp(t1.boot, c1.boot, t2.boot, c2.boot)
  }

  se   <- sd(b.est)
  Z    <- (Q - (1 - q)) / se
  pval <- 2 * (1 - pnorm(abs(Z)))

  #- Plots
  if (plots == TRUE) {
    plot(fit1, col = "red", ylab = "Estimated Survival Function",
         xlab = "Time", main = "Kaplan-Meier Estimates")
    lines(fit2, lty = 2, col = "blue")
    legend("topright", c("KM-Estimate for Control Group", "KM-Estimate for Trt. Group"),
           lty = c(1, 2), col = c("red", "blue"), bty = "n", cex = 0.8)
  }

  out <- list()
  out$quantile <- q
  out$sample1  <- F.inv
  out$sample2  <- G.inv
  out$Z        <- Z
  out$se       <- se
  out$pval     <- round(pval, 3)
  out$B        <- B
  return(out)
}
