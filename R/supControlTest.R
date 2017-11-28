#' Supremum-type test for two-sample comparison of survival quantiles.
#'
#' @param timevar1 Vector of observed survival times for sample 1 (control).
#' @param censor1  Vector of censoring indicators for sample 1 (1 = uncensored, 0 = censored).
#' @param timevar2 Vector of observed survival times for sample 2 (treatment).
#' @param censor2  Vector of censoring indicators for sample 2 (1 = uncensored, 0 = censored).
#' @param q.min    Smallest quantile (in terms of CDF) to test. Default is the time to earliest event for sample 1.
#' @param q.max    Largest quantile (in terms of CDF) to test.
#' @param gridpts Number of grid points between q.min and q.max to test.
#' @param B       Number of bootstrap samples.
#' @param seed    Seed number (for reproducibility).
#' @param plots   Logical. TRUE to show plot of cumulative distribution functions.
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
#' supControlTest(t1, c1, t2, c2, q.max = 0.5, B = 500)
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
#' @importFrom stats pnorm qnorm quantile sd stepfun
#'
supControlTest <- function(timevar1, censor1, timevar2, censor2, q.min = NULL, q.max = 0.5, gridpts = 50, B = 500, seed = 1234, plots = FALSE) {

  #- Checking for silly errors
  if (q.max < 0 | q.max > 1 ) {
    stop("q.max should be between 0 and 1.")
  }
  if (gridpts < 1) {
    stop("gridpts should be an integer larger than 0.")
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
  fit1  <- survfit(Surv(timevar1, censor1) ~ 1, conf.type = "none")
  fit2  <- survfit(Surv(timevar2, censor2) ~ 1, conf.type = "none")
  q.min <- min(1 - fit1$surv[fit1$n.censor == 0]) #Get first observed survival time

  if (q.max < q.min) {
    stop("Choose larger value for q.max.")
  }

  q.seq <- sort(seq(q.min, q.max, length = gridpts), decreasing = TRUE) #Discretize the quantiles

  test.stat <- numeric(length(q.seq))
  test.pval <- numeric(length(q.seq))

  for(j in 1:length(q.seq)) {
    fit <- quantileControlTest(timevar1, censor1, timevar2, censor2,
                               q = q.seq[j], B = B, plot = FALSE)
    test.stat[j] <- abs(fit$Z)
    test.pval[j] <- fit$pval
  }

  #- Plots
  if (plots == TRUE) {
    plot( (1 - fit1$surv) ~ fit1$time, col = "red", type = "s", ylab = "F(x)",
         xlab = "Time", main = "Estimated CDF for Control and Trt. Group")
    lines((1 - fit2$surv) ~ fit2$time, type = "s", lty = 2, col = "blue")
    legend("bottomright", c("CDF Estimate for Control Group", "CDF Estimate for Trt. Group"),
           lty = c(1, 2), col = c("red", "blue"), bty = "n", cex = 0.8)
  }

  out <- list()
  dat <- data.frame(cdf_quantile = round(q.seq, 3),
                    Z = round(abs(test.stat), 3),
                    pval = round(test.pval, 4))
  out$Z        <- max(dat$Z)
  out$pval     <- dat$pval[which.max(dat$Z)]
  out$results  <- dat
  return(out)
}
