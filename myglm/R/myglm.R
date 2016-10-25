myglm = function(formula, data = list(), family, ...){
  # Extract model matrix & responses
  mf = model.frame(formula = formula, data = data)
  X  = model.matrix(attr(mf, "terms"), data = mf)
  y  = model.response(mf)
  terms = attr(mf, "terms")

  if(family == "gaussian"){
  if (attr(terms,"intercept") == 1){
    rss_null <- t(y-mean(y))%*%(y-mean(y))
    p = ncol(X) - 1
  } else{
    rss_null <- t(y)%*%y
    p = ncol(X)
  }

  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y;      residuals = y - X%*%beta_hat;       rss_beta = t(y - X%*%beta_hat)%*%(y - X%*%beta_hat)
  n = nrow(X);                                sigma2 = as.numeric(rss_beta/n);
  beta_hat_cov = solve((t(X)%*%X)) * sigma2;  fitted = as.vector(X %*% beta_hat)
  R_squared = 1 - rss_beta/rss_null;          adj_R_squared = R_squared - (1 - R_squared) * (p)/(n - p)

  est = list(terms = terms, model = mf, coefficients = beta_hat, beta_cov = beta_hat_cov,
             residuals = residuals, n=n, p = p, sigma2=sigma2,
             rss_beta = rss_beta, rss_null = rss_null, fitted = fitted, y = y,
             R_squared = R_squared, adj_R_squared = adj_R_squared, x = X)
  # Store call and formula used
  est$call = match.call()
  est$formula = formula
  # Set class name. This is very important!
  class(est) = 'myglm'
  # Return the object with all results
  return(est)
  }
  else if(family == "poisson"){
    offset = model.offset(mf)
    loglik = function(offset, beta, X, y) {
      return (sum(dpois(y, lambda = exp(offset + X%*%beta),  log=TRUE)))
    }
    loglik1 = function(offset, beta, X, y) {
      return (sum(dpois(y, lambda = offset*exp(beta),  log=TRUE)))
    }
    beta = optim(par = numeric(dim(t(X))[1]), fn = loglik, method="BFGS", hessian = TRUE, control = list(fnscale=-1), X = X, offset = offset, y = y)
    betanull = optim(par = numeric(1), fn = loglik1, method="BFGS", hessian = TRUE, control = list(fnscale=-1), X = NULL, offset = offset, y = y)
    logsum <- function(y){
      S <- 0
      for(i in 1:y){
        S <- S + log(i)
      }
      return(S)
    }
    logsumy <- c()
    for(i in 1:length(y)){
      logsumy <- c(logsumy, logsum(y[i]))
    }
    lf = sum(y * log(y) - y - logsumy)
    ln = sum(y * log(offset*exp(betanull$par)) - offset*exp(betanull$par) - logsumy)
    print(2*(lf - ln))


    residual_df = length(dim(X)[1]) - length(beta$par)
    est = list(terms = terms, y = y, x = X, model = mf, offset = offset,
               coefficients = matrix(c(attr(X,"dimnames")[[2]], beta$par), ncol = length(beta$par),nrow = 2, byrow=TRUE),
               beta_cov = solve(-beta$hessian), residual_df = residual_df)


    est$call = match.call()
    est$formula = formula
    class(est) = 'myglm'
    return(est)
  }
}

print.myglm = function(x, ...){
  #cat('Call: \n')
  #print(x$call)
  #cat('\nCoefficients: \n')
  #print(x$coefficients[,])
}

summary.myglm = function(object, ...){
  coeff = object$coefficients; resid = object$residuals; p = object$p; n = object$n
  std_error = z_value = z_score = numeric(length(coeff))
  for (i in 1:length(coeff)){
    z_value[i]   = coeff[i] / sqrt(object$beta_cov[i,i])
    std_error[i] = sqrt(object$beta_cov[i,i])
  }
  z_score   = 1 - pnorm(abs(z_value), mean = 0, sd = 1)
  chi_stat = ((object$rss_null - object$rss_beta)/p) / (object$rss_beta/(n - p))
  chi_p_value = 1 - pchisq(chi_stat, p)
  chi_p_value[chi_p_value < 2e-16] <- '< 2e-16';   z_score[z_score < 2e-16] <- '< 2e-16'


  options(scipen=0)
  cat("Call: \n")
  print(object$call)
  cat('\n Residuals: \n')
  frame = data.frame(min(resid), quantile(resid, 0.25), median(resid),
                     quantile(resid, 0.75), max(resid), row.names = '')
  colnames(frame) <- c('Min','1stQ','Median','3rdQ','Max')
  print(frame)
  cat("\n Coefficients: \n")
  frame = data.frame(round(coeff,3), round(std_error,3), round(z_value,3), z_score)
  colnames(frame) <- c('Estimate','Std. Error','Z-value','Z-score')
  print(frame)
  cat('--- \n')
  cat("Residual standard error: ", sqrt(object$sigma2), " on ", n-p, " degrees of freedom \n")
  cat("Multiple R-squared: ", round(object$R_squared,5), "\t Adjusted R-squared:  ", round(object$adj_R_squared,5), "\n")
  cat("'Chi2'-statistic: ", round(chi_stat,3), " on ", p, "DF,  p-value: ", chi_p_value)
}

plot.myglm = function(x, ...){
  plot(x$fitted, x$y, xlab='Fitted', ylab='Observed', main='Observed vs fitted values')
}

anova.myglm = function(object, ...){
  comp = attr(object$terms, "term.labels")
  dof = Sum_Sq = numeric(length(comp)+1)
  # Name of response
  response = deparse(object$terms[[2]])
  # Fit the sequence of models
  txtFormula = paste(response, "~", sep = "")
  model = list()
  extra <- 0
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula = paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula = paste(txtFormula, comp[numComp], sep = "+")
    }
    formula = formula(txtFormula)
    model[[numComp]] = myglm(formula=formula, data = object$model)

    intercept_bool <- attr(object$terms, "intercept")
    num_levels <- length(levels(model[[numComp]]$model[[comp[numComp]]]))

    if (num_levels == 0){
      dof[numComp] = 1
    } else {
      dof[numComp] = num_levels - 1

      if (extra == 0 && intercept_bool == 0){
        extra = numComp
      }

    }

    if (numComp == 1){
      Sum_Sq[numComp] = object$rss_null - model[[numComp]]$rss_beta
    } else {
      Sum_Sq[numComp] = model[[numComp - 1]]$rss_beta - model[[numComp]]$rss_beta
    }

    factors_combi <- attr(model[[numComp]]$terms,"factors")[,numComp]
    if(sum(factors_combi) > 1){
      factors_combi_bool = as.logical(attr(model[[numComp]]$terms,"factors")[-1,numComp])
      multi_dof <- 1
      for (i in 1:length(factors_combi_bool)){
        if (factors_combi_bool[i] != 0){
          multi_dof = multi_dof * dof[i]
        }
      }
      dof[numComp] = multi_dof
    }
  }

  if (attr(object$terms, "intercept") == 0){
    dof[extra] = dof[extra] + 1
  }
  Sum_Sq[length(Sum_Sq)] = object$rss_null - sum(Sum_Sq);           dof[length(dof)] = nrow(object$x) - ncol(object$x);
  Mean_Sq = Sum_Sq/dof;                                             Chi2_value = Mean_Sq[1:length(Mean_Sq)-1] / Mean_Sq[length(Mean_Sq)];
  Chi2_score = 1 - pchisq(Chi2_value,dof[length(dof) - 1])

  Chi2_score[Chi2_score < 2e-16] <- '< 2e-16'

  frame = data.frame(dof, Sum_Sq, Mean_Sq, c(round(Chi2_value,5), ''), c(Chi2_score,''))
  colnames(frame) <- c('Df','Sum Sq','Mean Sq','Chi2-value','Chi2-score')
  rownames(frame) <- c(comp, "Residuals")
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep='')
  print(frame)
  cat('--- \n')
  #return(model)
}
