is.scalar <- function(x) is.atomic(x) && length(x) == 1L
is.numeric.scalar <- function(x) is.scalar(x) && is.numeric(x)