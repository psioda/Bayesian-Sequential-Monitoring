integrate_debug <- function(fun, xmin, xmax, ymin, ymax){
  tryCatch(
    tryCatch(
      tryCatch(
        integral2(fun, xmin, xmax, ymin, ymax)$Q,
        error = function(e) integral2(fun, xmin, xmax, ymin, ymax, singular = T)$Q),
      error = function(e) integral2(fun, xmin, xmax, ymin, ymax, abstol = 1E-6)$Q),
    error = function(e) integral2(fun, xmin, xmax, ymin, ymax, abstol = 1E-4)$Q)
}