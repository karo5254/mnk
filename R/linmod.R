

#' @title Projekt MNK
#' @author Karolina Augustyniak, Aleksandra Chwalek
#' @description Estymacja MNK
#' @param x matrix
#' @param ... additional arguments to be passed in the future update of this package
#' @return Funkcja zwaraca w wyniku parametry MNK
#' @rdname linmod
#' @export

linmod <- function(x, ...) UseMethod("linmod")

#' @title Estymacja MNK
#' @author Karolina Augustyniak, Aleksandra Chwalek
#' @description Estymacja MNK
#' @param x matrix
#' @param y matrix
#' @return Funkcja zwaraca w wyniku parametry MNK
#' @rdname linmodEst
#' @export

linmodEst <- function(x, y)
{
  ## compute QR-decomposition of x
  qx <- qr(x)
  ## compute (x'x)^(-1) x'y
  coef <- solve.qr(qx, y)
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)
  sigma2 <- sum((y - x%*%coef)^2)/df
  ## compute sigma^2 * (x'x)^-1
  vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       df = df)
}

#' @title Estymacja MNK
#' @author Karolina Augustyniak, Aleksandra Chwalek
#' @description Estymacja MNK
#' @param x matrix
#' @param y numeric
#' @param ... additional arguments to be passed in the future update of this package
#' @return Funkcja zwaraca w wyniku parametry MNK
#' @rdname linmod.default
#' @export

linmod.default <- function(x, y, ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  est <- linmodEst(x, y)
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "linmod"
  est
}

#' @title Wydruk MNK
#' @author Karolina Augustyniak, Aleksandra Chwalek
#' @description Estymacja MNK
#' @param x matrix
#' @param ... additional arguments to be passed in the future update of this package
#' @return Funkcja zwaraca w wyniku parametry MNK
#' @rdname print.linmod
#' @export

print.linmod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#' @title Podsumowanie MNK
#' @author Karolina Augustyniak, Aleksandra Chwalek
#' @description Estymacja MNK
#' @param object matrix
#' @param ... additional arguments to be passed in the future update of this package
#' @return Funkcja zwaraca w wyniku podsumowanie MNK
#' @rdname summary.linmod
#' @export

summary.linmod <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.linmod"
  res
}

#' @title Wydruk podsumowania MNK
#' @author Karolina Augustyniak, Aleksandra Chwalek
#' @description Estymacja MNK
#' @param x matrix
#' @param ... additional arguments to be passed in the future update of this package
#' @return Funkcja zwaraca w wyniku podsumowanie MNK
#' @rdname print.summary.linmod
#' @export

print.summary.linmod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}


