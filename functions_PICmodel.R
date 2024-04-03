# Model functions March 2024
library(Rlab)
library(numDeriv)
library(dplyr)

model.fit <- function(l1_x, l2_x, pi_x, data, epsilon=1e-08, short.epsilon=1e-1, short.iter=10, short.runs=20,
                      silent=T, include.h=T, gen.inits=T, init=NULL, init.h.only=F, two.step.h=F){
  sapply(c("Rlab", "dplyr", "numDeriv"), require, character.only = TRUE)
  # g(): the competing risks function used in the 'incidence' part of the mixture model
  # Input:
  #   - l1: the rate parameter for progression to CIN2+
  #   - l2: the rate parameter for viral clearance
  #   - t: time
  # Output:
  #   - the cumulative risk of CIN2+ among HPV-positive women wihtout CIN2+ at baseline
  g <- function(l1, l2, t){
    return(l1 / (l1 + l2) * (1 - exp(-(l1 + l2)* t)))
  }

  # create.par(): function to create the parameter names which will be used in the mstep expressions
  # Input:
  #   - n1: number of covariates used for the progression rate parameter + 1
  #   - n2: number of covariates used for the clearance rate parameter + 1
  #   - n3: number of covariates used in the parameter for the probability of prevalent disease + 1
  # Output:
  #   - vector of names of the parameters with (1) g0, g1... for l1, (2) w0, w1... for l2 and (3) p0, p1... for pi
  create.par <- function(n1, n2, n3) {
    par1 <- paste0(c("g"), seq(0, n1 - 1))
    par2 <- paste0(c("w"), seq(0, n2 - 1))
    par3 <- paste0(c("p"), seq(0, n3 - 1))
    return(c(par1, par2, par3))
  }

  # create.covariate.data(): creates matrices for each parameter of observed independent variables (like the 'X' matrices in a simple regression model)
  # Input:
  #   - data: first three columns must be (1) left interval, (2) right interval and (3) z indicator for prevalent/incident disease,
  #           following columns must match the column names given in l1_x, l2_x and pi_x
  # Output:
  #   - list of 3 matrices corresponding to relevant covariates for each of the parameters lambda_1, lambda_2 and pi
  create.covariate.data <- function(data){
    data1 <- matrix(c(rep(1, dim(data)[1]*n1)), ncol=n1) # create empty matrix
    if (n1 != 1){
      for (i in 1:(n1-1)){
        data1[,(i+1)] <- data[[l1_x[i]]] # use relevant column from the main data
      }
    }
    data2 <- matrix(c(rep(1, dim(data)[1]*n2)), ncol=n2) # create empty matrix
    if (n2 != 1){
      for (i in 1:(n2-1)){
        data2[,(i+1)] <- data[[l2_x[i]]] # use relevant column from the main data
      }
    }
    data3 <- matrix(c(rep(1, dim(data)[1]*n3)), ncol=n3) # create empty matrix
    if (n3 != 1){
      for (i in 1:(n3-1)){
        data3[,(i+1)] <- data[[pi_x[i]]] # use relevant column from the main data
      }
    }
    return(list(data1, data2, data3))
  }

  # log.likelihood(): function for the log-likelihood
  # Input:
  #   - current_par: the current value of the parameters
  #   - data: data frame with columns (1) left, (2) right, and (3) z, representing the left and right interval
  #               and the prevalent/incident indicator variable
  # Output:
  #   - value of the log-likelihood for the current parameter values
  log.likelihood.h <- function(current_par, data){
    z <- data[[3]] # indicator variable for prevalent or incident disease
    if (include.h){
      h <- exp(current_par[1])
      current_par <- current_par[-1]
    }else {
      h <-0
    }
    l1 <- exp(data1 %*% current_par[1:n1]) # create current lambda 1 vector for each value in the data
    l2 <- exp(data2 %*% current_par[(n1+1):(n1+n2)]) # create current lambda 2 vector for each value in the data
    if (n3 > 1) {
      # if there are covariates for pi then calculate pi for each data value using transformation from logit
      p <- exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)]))
    }else if (n3 == 1) {
      # if there are no covariates for pi (n==1 since only intercept), then calculate pi for each data value
      p <- data3 %*% current_par[n1+n2+n3]
    }

    f <- function(t) return(ifelse(t==Inf, 1, (1-exp(-h*t)) + exp(-h*t)*(l1/(l1 + l2))*(1-exp(-(l1 + l2)*t))))

    llk <- rep(0, length(right)) # create empty vector to store log-likelihood values
    llk[which(z==1)] <- I(log(p))[which(z==1)]
    llk[which(z==0 & right<Inf)] <- I(log((1-p)*(f(right) - f(left))))[which(z==0 & right <Inf)]
    llk[which(z==0 & right==Inf)] <- I(log(1-p)*(1- f(left)))[which(z==0 & right==Inf)]
    llk[which(is.na(z) & right<Inf)] <- I(log(p + (1-p)*f(right)))[which(is.na(z) & right<Inf)]
    llk[which(is.na(z) & right==Inf)] <- 0
    return(sum(llk))
  }

  estep.h <- function(current_par, data){
    z <- data[[3]] # indicator variable for prevalent or incident disease
    if (include.h){
      h <- exp(current_par[1])
      current_par <- current_par[-1]
    }else h <- 0

    # multiply data by corresponding parameters to calculate l1, l2 and pi
    l1 <- exp(data1 %*% current_par[1:(n1)])
    l2 <- exp(data2 %*% current_par[(n1+1):(n1+n2)])
    if (n3 == 1){
      # if pi has no covariates then the value of p is the same for all women
      p <- data3 %*% current_par[n1+n2+n3]
    }else if (n3 > 1){
      p <- exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)]))
    }
    # the expected value of z only depends on the right interval since the left is zero when z is unknown
    est <- p/(p + (1 - p)*(1-exp(-h*right) + exp(-h*right)*g(l1, l2, right)))
    est[which(right==Inf)] <- p[which(right==Inf)]
    # the expected value of the latent variable is the estimate above if z is unknown otherwise it is the true value of z
    z[is.na(z)] <-  est[is.na(z)]
    return(z)
  }

  build.expr.h <- function(model_par){
    # function to 'build' the expected log-likelihood function as an expression which
    # allows it to be differentiated w.r.t. different parameters in the mstep() function

    # we loop through the names of covariate parameters (e.g., g0, g1, w0, p0, p1, p2) for each of l1, l2 and pi
    # and paste them together with the corresponding data matrix, this string is then used to create the
    # expression which will be differentiated w.r.t each parameter in the mstep() function
    if(include.h) model_par <- model_par[-1]
    # expression for lambda_1 (progression rate parameter)
    expr1 <- model_par[1]
    if (n1 > 1){
      for (i in 1:(n1-1)){
        expr1 <- paste0(expr1, "+", model_par[1+i], "*", "data1", i)
      }}

    # expression for lambda_2 (clearance rate parameter)
    expr2 <- model_par[n1+1]
    if (n2>1){
      for (i in 1:(n2-1)){
        expr2 <- paste0(expr2, "+", model_par[n1+1+i], "*", "data2", i)
      }}

    # expression for pi (prevalent probability parameter)
    expr3 <- model_par[n1+n2+1]
    if (n3 > 1){
      for (i in 1:(n3-1)){
        expr3 <- paste0(expr3, "+", model_par[n1+n2+1+i], "*", "data3", i)
      }}
    p.expr <- ifelse(n3==1, str2expression(paste0("z*log(", expr3, ") + (1-z)*log(1-", expr3, ")")),
                     str2expression(paste0("z*(", expr3, ") - log(1 + exp(",expr3, "))")))
    f.expr <- str2expression(paste0("(1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",expr1 ,")/(exp(", expr1, ") + exp(",expr2,")))*(1-exp(-(exp(",expr1,") + exp(", expr2,"))*t))"))
    return(list(p.expr, f.expr))
  }


  mstep.h <- function(current_par, expected_z, data){
    z <- expected_z # expected value of z = output from the estep.h() function
    m1 <- n1+n2
    if (include.h) {m1 <- n1+n2 + 1} else {h <- -Inf}

    # assign the values of current_par to vectors with names from the 'pars' list (this will or won't have h and the number of parameters will correspond)
    for (i in 1:length(pars)){
      assign(pars[i], current_par[i], envir=environment())
    }

    f <- function(t) return(ifelse(t==Inf, 1, eval(f.expr)))     # incident function
    deriv.f <- function(t, param, order=1) { # derivative of incident function
      if (order ==1){
        return(ifelse(t==Inf, 0, eval(D(f.expr, param))))
      }else if (order==2){
        return(ifelse(t==Inf, 0, eval(D(D(f.expr, param[1]), param[2]))))
      }
    }
    exprs <- build.expr.h(pars)
    p.expr <- exprs[[1]]
    f.expr <- exprs[[2]]

    hess <- matrix(rep(0, (length(current_par))^2), nrow=length(current_par)) # empty hessian matrix
    grad <- rep(0, length(current_par)) # empty gradient matrix

    for (i in 1:m1){
      # we loop through the parameter list and calculate the first derivatives for gradient vector for the l1 and l2 parameters
      grad[i] <- sum(((1-z)*(deriv.f(right, pars[i]) - deriv.f(left, pars[i]))/(f(right) - f(left)))[z!=1])
      for (j in 1:i){
        # calculate 2nd derivatives for the Hessian matrix for h, l1 and l2 parameters; the hessian is symmetric so we can fill both hess[i,j] and hess[j,i] together to save time
        hess[i,j] <- hess[j,i] <- sum(((1-z)*(f(right) - f(left))^(-2)*
                                         ((deriv.f(right, c(pars[i], pars[j]), 2) - deriv.f(left, c(pars[i], pars[j]),2))*(f(right) - f(left)) -
                                            (deriv.f(right, pars[i], 1) - deriv.f(left, pars[i],1))*(deriv.f(right, pars[j],1) - deriv.f(left, pars[j], 1))))[z!=1])
      }
    }
    # parameters for pi (for example p0, p1, p2, p3, p4) do not depend on left/right so we can calculate it separately
    for (i in (m1+1):length(current_par)){
      # calculate first derivatives for gradient vector for the pi parameters
      grad[i] <- sum(eval(D(p.expr, pars[i])))
      for (j in (m1+1):i){
        # calculate 2nd derivatives for Hessian matrix for the pi parameters; the hessian is symmetric so we can fill both hess[i,j] and hess[j,i] together to save time
        hess[i,j] <- hess[j,i] <- sum(eval(D(D(p.expr, pars[i]), pars[j])))
      }
    }
    new_par <-  current_par - solve(hess)%*%grad # single Newton step to calculate updated parameters
    return(list(as.vector(new_par), hess, grad)) # return the new parameters
  }


  # em.function.h(): combining the E- and M-step and repeating until the parameter estimates converge
  em.function.h <- function(init, data){
    new_theta <- init
    old_llk <- 0
    new_llk <- 100 # make big to start with
    iter <- 1

    while ((abs(new_llk - old_llk) > epsilon)){
      current_theta <- new_theta
      old_llk <- new_llk
      next_em_step <- mstep.h(current_theta, estep.h(current_theta, data), data)
      new_theta <- as.vector(next_em_step[[1]])
      new_llk <- log.likelihood.h(new_theta, data)
      if (!silent & iter %/% 5 ==0) print(round(c(new_theta, new_llk), 4))
      iter <- iter + 1
    }
    hess <- next_em_step[[2]]
    names(new_theta) <- pars
    std.dev <- round(sqrt(-diag(solve(hess))),4)
    summary <- round(data.frame(theta.hat = new_theta, std.dev, lower = new_theta - 1.96*std.dev, upper = new_theta + 1.96*std.dev), 4)
    rownames(summary) <- names(new_theta)
    return(list(initial.values = init, theta.hat = new_theta, num.iterations=iter, log.likelihood = new_llk, hess=hess, grad = next_em_step[[3]], std.dev=std.dev, summary=summary))
  }

  # short.em():  lots of short runs of the em function to generate appropriate starting values
  short.em <- function(data){
    if(!init.h.only){
      new_theta <- log(runif(n1+n2+n3))
      if (n3==1) new_theta[n1+n2+n3] <- exp(new_theta[n1+n2+n3])
      if(include.h) {
        new_theta <- c(log(runif(1, 0, 0.2)), new_theta)
      }
    }else{
      new_theta <- c(log(runif(1, 0, 0.008)), init)
    }
    new_llk <- 100
    old_llk <- 0
    iter <- 1

    while ((abs(new_llk - old_llk) > short.epsilon) & iter < short.iter){
      current_theta <- new_theta
      old_llk <- new_llk
      iter <- iter + 1
      new_theta <- try(as.vector(mstep.h(current_theta, estep.h(current_theta, data), data)[[1]]), silent=T)
      if (class(new_theta) == 'try-error') return(short.em(data))
      new_llk <- log.likelihood.h(new_theta, data)
      if (abs(new_llk) == Inf | is.nan(new_llk) | is.na(new_llk)) return(short.em(data))
    }
    return(list(theta.hat = new_theta, log.likelihood = new_llk))
  }

  # init.generator(): performs short runs of the em function using short.em and then returns the most common set of starting values as the initial values
  init.generator <- function(data){
    if (!silent) pb <- txtProgressBar(min = 0, max = short.runs, style = 3, width = 50, char = "=")
    short.inits <- list()
    for(k in 1:short.runs) {
      short.inits[[k]] <- short.em(data)[c("theta.hat", "log.likelihood")]
      if (!silent) setTxtProgressBar(pb, k)
    }
    if (!silent) close(pb)
    short.inits.mat <- matrix(unlist(short.inits), nrow=length(short.inits), byrow=T)

    ncols <- ncol(short.inits.mat)
    # find the set of parameter values that results in the maximum likelihood
    init.counts <- cbind(round(short.inits.mat[rev(order(short.inits.mat[,ncols])),1:(ncols-1)],1), round(short.inits.mat[rev(order(short.inits.mat[,ncols])),ncols])) %>%
      data.frame()%>% group_by_all() %>% dplyr::count() %>% as.data.frame()

    if (all(init.counts$n==1)){
      # if there are no repeated starting values, we take the one with the highest likelihood
      inits <- short.inits.mat[which.max(short.inits.mat[,dim(short.inits.mat)[2]]),1:(dim(short.inits.mat)[2]-1)]
    }else{
      # of the starting values which are repeated more than once, take the one with the highest likelihood
      inits <- init.counts %>% arrange(., n) %>% filter(n > 1) %>% slice_max(n=1, order_by = get(noquote(paste0("X", ncols)))) %>% head(1) %>% as.numeric() %>% head(-2)
    }
    if (any(inits ==0)){
      inits[inits==0] <- 0.01
    }
    return(inits)
  }


  ##########################
  # run the function
  ##########################

  if (two.step.h){
    if(!silent) print(noquote("Generating inital values without background risk."))
    init.without.h <- model.fit(l1_x, l2_x, pi_x, data, epsilon, short.epsilon, short.iter, short.runs, silent, include.h=F, two.step.h = F)$theta.hat
    if(!silent) print(noquote("Generating inital values with background risk."))
    return(model.fit(l1_x, l2_x, pi_x, data, epsilon, short.epsilon, short.iter, short.runs, silent, include.h=T, init.h.only = T, init=init.without.h, two.step.h = F))
  }

  left <- data[[1]]
  right <- data[[2]] # right intervals from the data

  n1 <- length(l1_x) + 1 # number of parameters for lambda_1 (add 1 for intercept)
  n2 <- length(l2_x) + 1 # number of parameters for lambda_2 (add 1 for intercept)
  n3 <- length(pi_x) + 1 # number of parameters for pi (add 1 for intercept)

  covariate_data <- create.covariate.data(data)
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  if (n1 !=1){
    for (i in 1:(n1-1)){
      assign(paste0("data1", i), data1[,i+1], envir=environment())
    }}
  if (n2 !=1){
    for (i in 1:(n2-1)){
      assign(paste0("data2", i), data2[,i+1], envir=environment())
    }}
  if (n3!=1){
    for (i in 1:(n3-1)){
      assign(paste0("data3", i), data3[,i+1], envir=environment())
    }}

  pars <- create.par(n1, n2, n3)
  if(include.h) {
    pars <- c('h', pars)  # creates a vector of the names of parameters that need to be estimated
  }
  if (gen.inits) init <- init.generator(data)
  if(!silent) print(noquote("Running EM algorithm."))
  final.res <- em.function.h(init, data)
  return(final.res)
}

screening.simulator <- function(n, l1_x, l2_x, pi_x, params, show_prob = 0.9, i=5, include.h=T){
  d <- 1
  screening_times <- data.frame(x1 = ifelse(rbern(n, show_prob), 0, NA),
                                x2 = ifelse(rbern(n, show_prob), rnorm(n, i*1, sqrt(d)), NA),
                                x3 = ifelse(rbern(n, show_prob), rnorm(n, i*2, sqrt(d)), NA),
                                x4 = ifelse(rbern(n, show_prob), rnorm(n, i*3, sqrt(d)), NA),
                                x5 = ifelse(rbern(n, show_prob), rnorm(n, i*4, sqrt(d)), NA))
  age <- rbern(n, 0.5) # for older or younger than 40

  # Cytology Results: this is an indicator variable so 1 means abnormal cytology and 0 means not abnormal (or unknown for z=NA)
  # if they did not show up for screening at time 0 then their cytology result is 0 because it is unknown
  cytology <- rbern(n, 0.4) #ifelse(is.na(screening_times$x1), 0, rbern(n, 0.4))

  # HPV genotype (HPV 16 or other) - this is an indicator variable so 1 means they have HPV16 and 0 means other HPV type
  hpv <- rbern(n, 1/3)
  create.covariate.data <- function(data){
    data1 <- matrix(c(rep(1, dim(data)[1]*n1)), ncol=n1) # create empty matrix
    if (n1 != 1){
      for (i in 1:(n1-1)){
        data1[,(i+1)] <- data[[l1_x[i]]] # use relevant column from the main data
      }
    }
    data2 <- matrix(c(rep(1, dim(data)[1]*n2)), ncol=n2) # create empty matrix
    if (n2 != 1){
      for (i in 1:(n2-1)){
        data2[,(i+1)] <- data[[l2_x[i]]] # use relevant column from the main data
      }
    }
    data3 <- matrix(c(rep(1, dim(data)[1]*n3)), ncol=n3) # create empty matrix
    if (n3 != 1){
      for (i in 1:(n3-1)){
        data3[,(i+1)] <- data[[pi_x[i]]] # use relevant column from the main data
      }
    }
    return(list(data1, data2, data3))
  }

  create.covariate.data <- function(data){
    data1 <- matrix(c(rep(1, dim(data)[1]*n1)), ncol=n1) # create empty matrix
    if (n1 != 1){
      for (i in 1:(n1-1)){
        data1[,(i+1)] <- data[[l1_x[i]]] # use relevant column from the main data
      }
    }
    data2 <- matrix(c(rep(1, dim(data)[1]*n2)), ncol=n2) # create empty matrix
    if (n2 != 1){
      for (i in 1:(n2-1)){
        data2[,(i+1)] <- data[[l2_x[i]]] # use relevant column from the main data
      }
    }
    data3 <- matrix(c(rep(1, dim(data)[1]*n3)), ncol=n3) # create empty matrix
    if (n3 != 1){
      for (i in 1:(n3-1)){
        data3[,(i+1)] <- data[[pi_x[i]]] # use relevant column from the main data
      }
    }
    return(list(data1, data2, data3))
  }


  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  covariate_data <- create.covariate.data(data=data.frame(age=age, hpv=hpv, cyt=cytology))
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  if (include.h){
    h <- params[1]
    params <- params[-1]
  }

  l1 <- exp(data1 %*% params[1:(n1)])
  l2 <- exp(data2 %*% params[(n1+1):(n1+n2)])
  if (n3 == 1){
    # if pi has no covariates then the value of p is the same for all women
    p <- params[n1+n2+n3]
  }else if (n3 > 1){
    p <- exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)]))
  }

  # disease process
  t1 <- ifelse(rbern(n, p)==1, 0, Inf) # prevalent CIN2/3
  t2 <- ifelse(rbern(n, l1/(l1+l2))==1, rexp(n, rate=(l1+l2)), Inf) # due to the HPV infection at baseline
  if (include.h) {
    t3 <- rexp(n, rate=exp(h)) # due to the background risk (simulate value for everyone)
    t <- pmin(t1, t2, t3) # keep the minimum value as the actual event time
    cause <- apply(data.frame(t1, t2, t3), 1, FUN = which.min)
  }else{
    t <- pmin(t1, t2) # keep the minimum value as the actual event time
    cause <- apply(data.frame(t1, t2), 1, FUN = which.min)
  }

  screening_times$actual <- t
  screening_times$cause <- cause
  screening_times$cause <- NULL
  # create a list of the screening times with NA values removed for each subject
  # the code makes the following change: 0   NA    t2    t3    NA  ---> 0    t2    t3
  # this allows us to find the interval where CIN2/3 was detected
  if (show_prob ==1) {
    screens <- lapply(seq(1, n), function(x) unlist((screening_times)[x,] ))
  }
  else {
    screens <- apply(screening_times, 1, function(x) c(na.omit(unlist(x, use.names=FALSE))))
  }


  z <- rep(0, n) # create the indicator variable Z

  # create left intervals by finding the last value in the list of screens that is smaller than the actual event time ß
  left <- vapply(screens, function(x) x[Position(function(z) z <= x[length(x)], x[c(1:(length(x))-1)], right=TRUE)], 1)

  # if left interval is NA then disease was unknown at baseline because it was not checked
  z[is.na(left)] <- NA

  left[is.na(left)] <- 0 # set unknown left intervals to 0 because CIN2/3 could have been there at baseline

  # create a list of right intervals by finding the first value in the
  # list of screen times that is greater than the actual event time
  right <- vapply(screens, function(x) x[Position(function(z) z > x[length(x)], x[c(1:(length(x))-1)])], 1)

  # if the actual event time t=0 and left interval l=0 and the indicator is not unknown
  # (meaning disease was checked at baseline), then the right interval is also zero
  right[left==0 & t==0 & !is.na(z)] <-  0

  z[which(right==0)] <- 1 # right is only zero when disease is prevalent (defined above)

  # if the actual time of CIN2/3 development is after the last screening time, then the set the time to Inf
  last_screening_time <- vapply(screens, function(x) tail(x, 2)[1], 1)
  right[screening_times$actual > last_screening_time] <- Inf

  # if the right interval is NA then set it to infinity - this happens if all screening
  # rounds were NA (very rare, this is just to avoid errors in case it happens)
  right[is.na(right)] <- Inf

  return(data.frame(left, right, z = z, age = age, hpv = hpv, cyt=cytology, cause=cause, actual=t))
}

model.predict <- function(l1_x, l2_x, pi_x, data, time.points, fit, calc.CI=F, include.h=T){
  g <- function(l1, l2, t){
    return(l1 / (l1 + l2) * (1 - exp(-(l1 + l2)* t)))
  }

  requireNamespace("dlpyr", quietly = TRUE)
  theta.hat <- fit$theta.hat
  if (include.h) {
    h <- exp(theta.hat[1])
    theta.hat <- theta.hat[-1]
  }else{
    h <- 0
  }

  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  preds <- list()
  data <- unique(data)
  rownames(data) <- seq(1:nrow(data))
  # calculate the cumulative risk for each unique combination of covariates in the input prediction data and store in a list
  for (row in 1:nrow(data)){
    sub.data <- subset(data,  subset=rownames(data) == row)

    create.covariate.data <- function(data){
      data1 <- matrix(c(rep(1, dim(data)[1]*n1)), ncol=n1) # create empty matrix
      if (n1 != 1){
        for (i in 1:(n1-1)){
          data1[,(i+1)] <- data[[l1_x[i]]] # use relevant column from the main data
        }
      }
      data2 <- matrix(c(rep(1, dim(data)[1]*n2)), ncol=n2) # create empty matrix
      if (n2 != 1){
        for (i in 1:(n2-1)){
          data2[,(i+1)] <- data[[l2_x[i]]] # use relevant column from the main data
        }
      }
      data3 <- matrix(c(rep(1, dim(data)[1]*n3)), ncol=n3) # create empty matrix
      if (n3 != 1){
        for (i in 1:(n3-1)){
          data3[,(i+1)] <- data[[pi_x[i]]] # use relevant column from the main data
        }
      }
      return(list(data1, data2, data3))
    }

    covariate.data <- create.covariate.data(sub.data)
    data1 <- covariate.data[[1]]
    data2 <- covariate.data[[2]]
    data3 <- covariate.data[[3]]

    l1 <- as.numeric(exp(data1 %*% theta.hat[1:(n1)]))
    l2 <- as.numeric(exp(data2 %*% theta.hat[(n1+1):(n1+n2)]))
    if (n3 == 1){
      #if pi has no covariates then the value of p is the same for all women
      p <- theta.hat[n1+n2+n3]
    }else if (n3 > 1){
      p <- as.numeric(exp(data3 %*% theta.hat[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% theta.hat[(n1+n2+1):(n1+n2+n3)])))
    }
    prob <- p + (1-p)*(1 - exp(-h*time.points) + exp(-h*time.points)*g(l1, l2, time.points))

    ## calculate confidence errors for cumulative risk using the complementary log-log cumulative risk

    if (calc.CI){
      cr.est <- log(-log(prob))
      cr.se <- sapply(time.points,  se.cumrisk, data1, data2, data3, n1, n2, n3, fit, include.h)
      cr.lci <- exp(-exp(cr.est + 1.96*cr.se))
      cr.uci <- exp(-exp(cr.est - 1.96*cr.se))
      preds[[row]] <- data.frame(Time=time.points, CR = prob, CR.se = cr.se, CR.lower95 = cr.lci, CR.upper95 = cr.uci)
    }else{
      preds[[row]] <- data.frame(Time=time.points, CR = prob)
    }

  }
  return(preds)
}

se.cumrisk <- function(t, data1, data2, data3, n1, n2, n3, fit, include.h){
  create.par <- function(n1, n2, n3) {
    par1 <- paste0(c("g"), seq(0, n1 - 1))
    par2 <- paste0(c("w"), seq(0, n2 - 1))
    par3 <- paste0(c("p"), seq(0, n3 - 1))
    return(c(par1, par2, par3))
  }

  if (include.h){
    pars <- c("h", create.par(n1, n2, n3))
  }else{
    pars <- create.par(n1, n2, n3)
    h <- 0
  }

  if (n1 !=1){
    for (i in 1:(n1-1)){
      assign(paste0("data1", i), data1[,i+1], envir=environment())
    }}
  if (n2 !=1){
    for (i in 1:(n2-1)){
      assign(paste0("data2", i), data2[,i+1], envir=environment())
    }}
  if (n3!=1){
    for (i in 1:(n3-1)){
      assign(paste0("data3", i), data3[,i+1], envir=environment())
    }}


  # assign the values of current_par to vectors with names from the 'pars' list
  for (i in 1:length(pars)){
    assign(pars[i], fit[["theta.hat"]][i], envir=environment())
  }


  # create the expression for the complementary log log cumulative risk, it is coded this way to allow for user-defined covariates
  # and then we take the derivative so it needs to be in expression form
  create.cllcr.expr <- function(n1, n2, n3, pars, data1, data2, data3){
    pars <- pars[-1]
    # expression for lambda_1 (progression rate parameter)
    expr1 <- pars[1]
    if (n1 != 1){
      for (i in 1:(n1-1)){
        expr1 <- paste0(expr1, "+", pars[1+i], "*", "data1", i)
      }}

    # expression for lambda_2 (clearance rate parameter)
    expr2 <- pars[n1+1]
    if (n2!=1){
      for (i in 1:(n2-1)){
        expr2 <- paste0(expr2, "+", pars[n1+1+i], "*", "data2", i)
      }}

    # expression for pi (prevalent probability parameter)
    expr3 <- pars[n1+n2+1]
    if (n3 > 1){
      for (i in 1:(n3-1)){
        expr3 <- paste0(expr3, "+", pars[n1+n2+1+i], "*", "data3", i)
      }

      #  expected complete log-likelihood for Ri!=Inf
      llcr <- paste0("log(- log( (exp(", expr3, ")/(1+ exp(", expr3, "))) + (1/(1+exp(", expr3, ")))*((1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",
                     expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,"))*t)))))" )


    }else if (n3 == 1){ # if there are no covariates for pi we use a different expression (don't need to use logit function)
      #  expected complete log-likelihood for Ri!=Inf for when there is only intercept for pi parameter
      # llcr <- paste0("log(- log( log(", expr3, ") + (log(1-", expr3, ") + log((1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",
      #                expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,"))*t))))))" )
      llcr <- paste0("log(- log( ", expr3, " + (1 - ", expr3, ")*((1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",
                     expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,"))*t)))))" )
    }

    # convert from string to expression to be used in the differentiate function in the gradient calculation
    return(str2expression(llcr))
  }

  # function to calculate the gradient of the log log cumulative risk using the built in R derivative function D()
  grad.loglog <- function(par){
    return(eval(D(create.cllcr.expr(n1, n2, n3, pars, data1, data2, data3), par)))
  }

  return(sqrt(as.numeric(t(sapply(pars, grad.loglog)) %*% (solve(-fit[["hess"]])) %*% (sapply(pars, grad.loglog)))))
}

