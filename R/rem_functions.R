#'Internal helper function of rainforest(): Summary effect (REM)
#'
#'Computes the meta-analytic summary effect(s) for the random effects model (for method = "REM")
#'@export
#'@keywords internal
rem_effect <- function(es, se) {
  summary_es_FEM <- sum((1/se^2)*es)/sum(1/se^2)
  ni <- length(es)
  if(ni == 1) {
    warning("Some subgroups consist only of one study.
            Between study variance in these subgroups cannot be estimated and is set to zero")
    t2 <- 0
  } else {
    Q <- sum((1 / se^2) * (es - summary_es_FEM)^2)
    t2 <- max(c(0, (Q - (ni - 1)) / (sum(1 / se^2) - sum((1 / se^2)^2) / sum(1/se^2))))
  }
  w <- 1/(se^2 + t2)
  sum(w*es)/sum(w)
}
#'Internal helper function of rainforest(): Summary effect (REM)
#'
#'Computes the standard error of the meta-analytic summary effect(s) for the random effects model (for method = "REM")
#'@export
#'@keywords internal
rem_err <- function(es, se) {
  summary_es_FEM <- sum((1/se^2)*es)/sum(1/se^2)
  ni <- length(es)
  if(ni == 1) {
    t2 <- 0
  } else {
    Q <- sum((1 / se^2) * (es - summary_es_FEM)^2)
    t2 <- max(c(0, (Q - (ni - 1)) / (sum(1 / se^2) - sum((1 / se^2)^2) / sum(1/se^2))))
  }
  w <- 1/(se^2 + t2)
  sqrt(1/sum(w))
}
