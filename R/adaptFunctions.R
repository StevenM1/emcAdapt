adapt.c.emc <- function(feedback, arguments, learningRule='delta') {
  nTrials <- nrow(feedback)
  nAdapt <- ncol(feedback)

  if(learningRule == 'delta') {
    # Simple Delta rule
    startValues <- arguments$startValues
    learningRates <- arguments$learningRates

    ## empty output arrays
    adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)
    out = .C('adaptDelta',
             nTrials=nTrials,
             nChoices=nAdapt,
             values=as.double(startValues),
             adaptedValues=as.double(adaptedValues),
             predictionErrors=as.double(predictionErrors),
             outcomes=as.double(feedback),
             learningRates=as.double(learningRates),
             NAOK=TRUE)

    adaptedValues <- matrix(out$adaptedValues, nrow=nTrials, ncol=nAdapt)
    predictionErrors <- matrix(out$predictionErrors, nrow=nTrials, ncol=nAdapt)
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors))
  } else if(learningRule == 'vkf' | learningRule == 'vkfbinary') {
    # Variational Kalman filter
    predictionsStartValues <- arguments$predictionsStartValues
    volatilitiesStartValues <- arguments$volatilitiesStartValues
    volatilityLearningRates <- arguments$volatilityLearningRates
    uncertaintiesStartValues <- arguments$uncertaintiesStartValues

    if(learningRule == 'vkfbinary') {
      ccall <- 'adaptVKFbinary'
    } else {
      ccall <- 'adaptVKF'
    }

    ## empty output arrays
    predictions <- predictionErrors <- learningRates <- volatilities <- volatilityPredictionErrors <- uncertainties <- matrix(nrow=nTrials, ncol=nAdapt)

    out = .C(ccall,
             nTrials=nTrials,
             nChoices=nAdapt,
             # Predictions
             predictions=as.double(predictionsStartValues),
             adaptedPredictions=as.double(predictions),
             predictionErrors=as.double(predictionErrors),
             learningRates=as.double(learningRates),
             # volatility
             volatilities=as.double(volatilitiesStartValues),
             adaptedVolatilities=as.double(volatilities),
             volatilityPredictionErrors=as.double(volatilityPredictionErrors),
             volatilityLearningRates=as.double(volatilityLearningRates),
             # uncertainty
             uncertainties=as.double(uncertaintiesStartValues),
             adaptedUncertainties=as.double(uncertainties),
             # feedback
             outcomes=as.double(feedback),
             NAOK=TRUE)

    adaptedPredictions <- matrix(out$adaptedPredictions, nrow=nTrials, ncol=nAdapt)
    predictionErrors <- matrix(out$predictionErrors, nrow=nTrials, ncol=nAdapt)

    adaptedVolatilities <- matrix(out$adaptedVolatilities, nrow=nTrials, ncol=nAdapt)
    volatilityPredictionErrors <- matrix(out$volatilityPredictionErrors, nrow=nTrials, ncol=nAdapt)

    learningRates <- matrix(out$learningRates, nrow=nTrials, ncol=nAdapt)
    adaptedUncertainties <- matrix(out$adaptedUncertainties, nrow=nTrials, ncol=nAdapt)

    return(list(adaptedPredictions=adaptedPredictions, predictionErrors=predictionErrors, learningRates=learningRates,
                adaptedVolatilities=adaptedVolatilities, volatilityPredictionErrors=volatilityPredictionErrors, adaptedUncertainties=adaptedUncertainties))
  }
}

adapt.r.test <- function(startValues, learningRates, feedback, learningRule='delta',
                         learningRatesNeg=NULL, riskStartValues=NULL, riskLearningRates=NULL) {
  nTrials <- nrow(feedback)
  nAdapt <- ncol(feedback)

  # declare output array
  adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)

  if(learningRule == 'SARSARisk') {
    adaptedRiskValues <- riskPredictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)
    riskValues <- riskStartValues
  }

  values <- startValues
  # Update all values (can we vectorize this? Hopefully)
  for(trial in 1:nTrials) {
    c_ = which(!is.na(feedback[trial,])) #choice[trial] # c_ = choice
    o_ = feedback[trial,c_] # o_ = outcome

    # Keep track of current value, used later for drift rates
    adaptedValues[trial, c_] = values[c_]
    adaptedValues[trial, -c_] = values[-c_]

    if(learningRule == 'delta' | learningRule == 'SARSAvarLR' | learningRule == 'SARSARisk') {
      predicted <- values[c_]
      if(learningRule == 'SARSARisk') {
        adaptedRiskValues[trial, c_] = riskValues[c_]
        adaptedRiskValues[trial, -c_] = riskValues[-c_]
      }
    } else if(learningRule == 'Qlearning')  {
      predicted <- max(values)
    }
    # calculate PE
    dv = o_-predicted  # prediction error = outcome (reward) - predicted value
    predictionErrors[trial,c_] <- dv  # keep track of this

    if(!is.null(learningRatesNeg)) {
      LR = ifelse(dv>0, learningRates[trial,c_], learningRatesNeg[trial,c_])
    } else if(learningRule == 'SARSAvarLR') {
      LR = learningRates[trial,c_] * (predicted*(1-predicted))
    } else if (learningRule == 'delta') {
      LR = learningRates[trial,c_]
    } else if(learningRule == 'SARSARisk') {
      LR = learningRates[trial,c_]
      riskLR = riskLearningRates[trial, c_]
      # update risk
      riskPE <- dv^2 - riskValues[c_]
      riskValues[c_] = riskValues[c_]+riskLR*riskPE
      riskPredictionErrors[trial,c_] = riskPE
      dv <- dv /sqrt(riskValues[c_])  # scale PE
    }
    values[c_] = values[c_] + LR*dv

    # # update values
    # if(!is.null(learningRatesNeg)) {
    #   values[c_] = values[c_] + ifelse(dv>0, learningRates[trial,c_], learningRatesNeg[trial,c_])*dv
    # } else {
    #   if(learningRule == 'SARSAvarLR') {
    #     LR
    #   }
    #   values[c_] = values[c_] + learningRates[trial,c_]*dv
    # }
  }
  if(learningRule == 'SARSARisk') {
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors,
           adaptedRiskValues=adaptedRiskValues, riskPredictionErrors=riskPredictionErrors))
  } else {
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors))
  }
}

vkf.r <- function(x, data, sample = FALSE, updateOnly=TRUE) {
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

vkf_binary.r <- function(x, data, sample = FALSE, updateOnly=TRUE) {
  #       data: column-vector of outcomes
  #       0<lambda<1, volatility learning rate
  #       v0>0, initial volatility
  #       sigma2>0, outcome noise
  #       sigmar = simulation response noise (sd)
  #       mu0 = initial value of mu (assumed to be 0 in code below)
  # also see https://github.com/payampiray/VKF/blob/master/vkf.m

  if (is.null(names(x))) {
    names(x) <- c("lambda", "v0", "omega", "sigmar", "mu0")
  }

  # extract parameters
  lambda <- x[["lambda"]]
  v0 <- x[["v0"]]
  omega <- x[["omega"]]

  # number of trials
  nt <- nrow(data)

  # initial values
  m <- 300 * x[["mu0"]]
  w0 <- x[["omega"]]
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

    k       <- (w + v) / (w + v + omega)                          # Eq 14
    alpha   <- sqrt(w + v)                                        # Eq 15
    delta_m <- o - (1/(1+exp(-m)))                                # Eq 16 (sigmoid)
    m       <- m + alpha * delta_m                                # Eq 16
    w       <- (1 - k) * (w + v)                                  # Eq 17

    wcov    <-  (1 - k) * wpre                                    # Eq 18
    delta_v <-  (m - mpre)^2 + w + wpre - 2 * wcov - v            # Eq 19.1
    v       <-  v + lambda * delta_v                              # Eq 19.2

    learning_rate[t] <- alpha                                     # learning rate is now alpha!!
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


# Testing VKF -------------------------------------------------------------
# feedback = matrix(rnorm(20), ncol=2)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# out = adapt.c.emc(feedback=feedback, learningRule='vkf',
#                   arguments=list(predictionsStartValues=predictionsStartValues,
#                                  volatilitiesStartValues=volatilitiesStartValues,
#                                  uncertaintiesStartValues=uncertaintiesStartValues,
#                                  volatilityLearningRates=volatilityLearningRates))
#
# ## compare with Q's implementation
# column <- 2
# x <- c("lambda"=volatilityLearningRates[1,column],    # volatility learning rate
#        "v0"=volatilitiesStartValues[column],        # volatility start value
#        "sigma2"=uncertaintiesStartValues[column],    # uncertainty start value
#        "sigmar"=uncertaintiesStartValues[column],    # simulation noise, unused for updating
#        "mu0"=predictionsStartValues[column]      # prediction start value
#        )
# out2 <- vkf.r(x, data=data.frame(o=feedback[,column]), updateOnly=TRUE)
#
# all(out2$predictions==out$adaptedPredictions[,column])
# all(out2$volatility==out$adaptedVolatilities[,column])
# all(out2$uncertainty==out$adaptedUncertainties[,column])
#

## speed comparison
# feedback = matrix(rnorm(2000), ncol=1)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# library(microbenchmark)
# microbenchmark(adapt.c.emc(feedback=feedback, learningRule='vkf',
#                            arguments=list(predictionsStartValues=predictionsStartValues,
#                                           volatilitiesStartValues=volatilitiesStartValues,
#                                           uncertaintiesStartValues=uncertaintiesStartValues,
#                                           volatilityLearningRates=volatilityLearningRates)),
#                vkf.r(x, data=data.frame(o=feedback), updateOnly=TRUE))
# Unit: microseconds
#       min       lq      mean   median        uq      max neval
#   67.978   83.517  107.8521   92.824  103.9555 1274.075   100
# 1183.014 1215.014 1349.7270 1252.796 1352.3645 4342.310   100


# Testing BINARY VKF -------------------------------------------------------------
# feedback = matrix(rnorm(20), ncol=2)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# out = adapt.c.emc(feedback=feedback, learningRule='vkfbinary',
#                   arguments=list(predictionsStartValues=predictionsStartValues,
#                                  volatilitiesStartValues=volatilitiesStartValues,
#                                  uncertaintiesStartValues=uncertaintiesStartValues,
#                                  volatilityLearningRates=volatilityLearningRates))
#
# ## compare with Q's implementation
# column <- 2
# x <- c("lambda"=volatilityLearningRates[1,column],    # volatility learning rate
#        "v0"=volatilitiesStartValues[column],        # volatility start value
#        "omega"=uncertaintiesStartValues[column],    # uncertainty start value
#        "sigmar"=uncertaintiesStartValues[column],    # simulation noise, unused for updating
#        "mu0"=predictionsStartValues[column]      # prediction start value
#        )
# out2 <- vkf_binary.r(x, data=data.frame(o=feedback[,column]), updateOnly=TRUE)
#
# all(out2$predictions==out$adaptedPredictions[,column])
# all(out2$volatility==out$adaptedVolatilities[,column])
# all(out2$uncertainty==out$adaptedUncertainties[,column])
#
#
# # speed comparison
# feedback = matrix(rnorm(2000), ncol=1)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# library(microbenchmark)
# microbenchmark(adapt.c.emc(feedback=feedback, learningRule='vkfbinary',
#                            arguments=list(predictionsStartValues=predictionsStartValues,
#                                           volatilitiesStartValues=volatilitiesStartValues,
#                                           uncertaintiesStartValues=uncertaintiesStartValues,
#                                           volatilityLearningRates=volatilityLearningRates)),
#                vkf_binary.r(x, data=data.frame(o=feedback), updateOnly=TRUE))
# Unit: microseconds
# min        lq       mean   median        uq      max neval
# 75.973   89.7285   98.45125   97.088  103.3405  168.428   100
# 1346.276 1365.3820 1481.72811 1382.028 1410.5435 3725.588   100



# Testing delta rule ------------------------------------------------------
# #
# nTrials <- 100
# notChosen <- sample(c(1, 2), size=nTrials, replace=TRUE)
# eta1 <- .3
# eta2 <- .4
# startValues <- c(.5, .5)
# feedback <- matrix(rnorm(nTrials*2, 5, 3), nrow=nTrials)
# feedback[cbind(1:nTrials, as.numeric(notChosen))] <- NA
# choice <- ifelse(notChosen==1, 2, 1)
# #
# tmpSarsaC = adapt.c.emc(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='delta')
# tmpSarsaR = adapt.r.test(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='delta')
#
# tmpSarsaC$adaptedValues == tmpSarsaR$adaptedValues
#
# #
# ## negative learning rates?
# posLR <- matrix(eta1, nrow=nTrials, ncol=2)
# negLR <- matrix(eta2, nrow=nTrials, ncol=2)
# tmpSarsaC = adapt.c.emc(feedback=feedback, startValues=startValues, posLR, learningRule='SARSA', negLR)
# tmpSarsaR = adapt.r.test(feedback=feedback, startValues=startValues, posLR, learningRule='SARSA', negLR)
#
#
# ## variable learning rates?
# tmpSarsaC = adapt.c.emc(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='SARSAvarLR')
# tmpSarsaR = adapt.r.test(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='SARSAvarLR')
#
# tmpSarsaC$adaptedValues == tmpSarsaR$adaptedValues
#
#
# ## Risk learning?
# riskStartValues <- c(1,1)
# tmpSarsaRiskC = adapt.c.emc(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2),
#                             learningRule='SARSARisk', riskStartValues=riskStartValues, riskLearningRates = matrix(eta1, nrow=nTrials, ncol=2))
# tmpSarsaRiskR = adapt.r.test(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2),
#                              learningRule='SARSARisk', riskStartValues=riskStartValues, riskLearningRates = matrix(eta1, nrow=nTrials, ncol=2))
# tmpSarsaRiskC$adaptedValues == tmpSarsaRiskR$adaptedValues
# tmpSarsaRiskC$adaptedRiskValues == tmpSarsaRiskR$adaptedRiskValues
# tmpSarsaRiskC$riskPredictionErrors == tmpSarsaRiskR$riskPredictionErrors
