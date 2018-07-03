#' Smoothing kernels used in the package
#'
#' @keywords internal
#' @noRd

K <- function(x)
{
  (abs(x)<1)*0.75*(1-(x)^2)
}

sigk1 <- sqrt(0.2)
Kq <- function(x,q)
{
  2/(q+1)*K(2/(q+1)*(x-(q-1)/2))*(1+((q-1)/(q+1)/sigk1)^2+2/sigk1^2*(1-q)/(1+q)^2*x)
}
