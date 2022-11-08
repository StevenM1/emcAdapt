
ll_func <- function(x, data, sample = FALSE, updateOnly=TRUE) {
  #       data: column-vector of outcomes
  #       0<lambda<1, volatility learning rate
  #       v0>0, initial volatility
  #       sigma2>0, outcome noise
  #       sigmar = simulation response noise (sd)
  #       mu0 = initial value of mu (assumed to be 0 in code below)
  # also see https://github.com/payampiray/VKF/blob/master/vkf.m

  if (is.null(names(x))) {
    names(x) <- c("lambda", "v0", "sigma2", "sigmar", "mu0")
  }

  # extract parameters
  lambda <- x[["lambda"]]
  v0 <- x[["v0"]]
  sigma2 <- x[["sigma2"]]

  # number of trials
  nt <- nrow(data)

  # initial values
  m <- 300 * x[["mu0"]]
  w0 <- x[["sigma2"]]
  w <- w0
  v <- x[["v0"]]

  predictions <- numeric(nt)
  learning_rate <- numeric(nt)
  volatility <- numeric(nt)
  prediction_error <- numeric(nt)
  volatility_error <- numeric(nt)
  uncertainty <- numeric(nt)

  for (t in seq_len(nt)) {

    o <- data$o[t]
    predictions[t] <- m
    volatility[t] <- v
    uncertainty[t] <- w

    mpre <- m
    wpre <- w

    delta_m <- o - m
    k       <- (w + v) / (w + v + sigma2)                         # Eq 9
    m       <- m + k * delta_m                                    # Eq 10
    w       <- (1 - k) * (w + v)                                  # Eq 11

    wcov    <-  (1 - k) * wpre                                    # Eq 12
    delta_v <-  (m - mpre)^2 + w + wpre - 2 * wcov - v
    v       <-  v + lambda * delta_v                              # Eq 13

    learning_rate[t] <- k
    prediction_error[t] <- delta_m
    volatility_error[t] <- delta_v

  }

  if(updateOnly) {
    return(list(predictions = predictions,
         volatility = volatility,
         learning_rate = learning_rate,
         prediction_error = prediction_error,
         volatility_error = volatility_error,
         uncertainty = uncertainty))
  }
  # response model
  if (sample) {

    # generate responses
    data$r <- rnorm(nt, mean = predictions, sd = x[["sigmar"]])
    attr(data, "latent_state_pars") <- list(predictions = predictions,
                                            volatility = volatility,
                                            learning_rate = learning_rate,
                                            prediction_error = prediction_error,
                                            volatility_error = volatility_error,
                                            uncertainty = uncertainty)
    return(data)


  } else {

    out <- sum(dnorm(data$r,
                     mean = predictions,
                     sd = x[["sigmar"]],
                     log = TRUE))
    attr(data, "latent_state_pars") <- list(predictions = predictions,
                                            volatility = volatility,
                                            learning_rate = learning_rate,
                                            prediction_error = prediction_error,
                                            volatility_error = volatility_error,
                                            uncertainty = uncertainty)
    return(out)

  }
}

## testing

x <- c("lambda"=.5, "v0"=.1, "sigma2"=1, "sigmar"=1, "mu0"=2)
data <- data.frame('o'=rnorm(50))
output <- ll_func(x, data, sample=TRUE)
attr(output, 'latent_state_pars')
