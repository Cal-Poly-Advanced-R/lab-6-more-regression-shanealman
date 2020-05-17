#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#' @param explanatory The name of the explanatory variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#' @import data.table
#'
#' @export

slr_gd <- function(dat, response, explanatory){

  y = dat %>% pull({{response}})
  x = dat %>% pull({{explanatory}})

  sdy = sd(y)
  sdx = sd(x)
  meanx = mean(x)
  meany = mean(y)

  x = as.numeric(scale(x))
  y = as.numeric(scale(y))

  psi = 0.00001
  m = runif(1,0); b = runif(1,0); n = length(x)
  yhat = m*x + b

  for(i in 1:1000000){

    m = m - psi * (2/n) * sum((yhat - y) * x)
    b = b - psi * (2/n) * sum(yhat - y)

    yhat =  m * x + b

  }

  results = data.table(m, 'Intercept' = b)
  names(results)[1] = mtcars %>% select({{explanatory}}) %>% names()


  results[,1] = results[,1] * sdy / sdx
  results[,2] = results[,2]*sdy + meany - results[,1] * meanx

  return(results)

}



#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#' @impor data.table
#'
#' @export

mlr_gd <- function(dat, response) {

  dt = data.table(dat)
  y = dt %>%  pull({{response}})
  x = dt %>% select(-{{response}})

  sdy = sd(y)
  sdx = sapply(x, sd)
  meany = mean(y)
  meanx = sapply(x, mean)

  y = as.numeric(scale(y))
  x = scale(x)

  x = cbind(1, x)
  x = as.matrix(x)

  betas = as.matrix(runif(dim(x)[2],0))
  n = length(y)
  psi = 0.0001

  for(i in 1:1000000){

    betas = betas - psi * (1/n) * t(x) %*% (x %*% betas - y)

  }

  rownames(betas)[1] = 'Intercept'
  results = data.table(t(betas))

  results[,2:ncol(x)] = results[,2:ncol(x)] * sdy / sdx
  results[,1] = results[,1] * sdy + meany - sum(results[,2:ncol(x)] * meanx)

  return(results)
}



#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#' @export

mlr_qr <- function(dat, response) {

  dt <- data.table(dat)
  y <- dt %>%  pull({{response}})
  x <- dt %>% select(-{{response}})

  sdy = sd(y)
  sdx = sapply(x, sd)
  meany = mean(y)
  meanx = sapply(x, mean)

  y = as.numeric(scale(y))
  x = scale(x)

  x <- cbind(1, x)
  x <- as.matrix(x)

  QR <- qr(x)
  betas <- solve.qr(QR, y)

  results <- data.table(t(betas))
  names(results)[1] <- 'Intercept'

  results[,2:ncol(x)] = results[,2:ncol(x)] * sdy / sdx
  results[,1] = results[,1] * sdy + meany - sum(results[,2:ncol(x)] * meanx)

  return(results)

}


