#load("~/Documents/GLM/myglm/Workspace.RData")
credit_data = read.table("https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment3/credit.txt",header=TRUE)


cross_validation = function(nr_folds, price_missed_bad, price_missed_good, intercept, brute){
  acceptance_prob = price_missed_good/(price_missed_bad + price_missed_good)
  size_fold = ceiling(dim(credit_data)[1]/nr_folds)
  comp = c(colnames(credit_data)[1:5],
           paste(colnames(credit_data)[6:11], collapse = "+"),
           paste(colnames(credit_data)[12:17], collapse = "+"),
           paste(colnames(credit_data)[18:23], collapse = "+"))

  # Name of response
  response = credit_data$Y
  # Fit the sequence of models
  n = 0
  model_error_rates = numeric()
  if (brute){
    for (m in 1:length(comp)){
      combinations = combn(1:length(comp), m)
      for (i in 1:dim(combinations)[2]){
        if (intercept){
          txtFormula = paste("Y", "~", sep = " ")
        } else {
          txtFormula = paste("Y", "~ - 1 +", sep = " ")
        }
        txtFormula = paste(txtFormula, paste(comp[combinations[,i]], collapse = "+"), sep = "")
        n = n + 1
        cat(txtFormula, 'Model nr: ', n, '/', 2^(length(comp)) - 1, '\n')
        
        formula = formula(txtFormula)
        mf = model.frame(formula = formula, data = credit_data)
        X = model.matrix(formula, data = mf) 
        y  = model.response(mf)
        model_error_rates[[txtFormula]] = 0
        for (fold in 1:nr_folds){
          #Select test data
          test_data <- ((fold-1)*size_fold + 1) : (fold*size_fold)
          #Select training data
          training_data <- 1:dim(credit_data)[1]
          training_data[test_data] = 0
          #Run estimate to achieve parameter
          beta <- IRLS(numeric(dim(credit_data)[1] - size_fold) + 1, y[training_data], X[training_data,])
          #Get matrix of X to calculate probabilities of test set
          mu <- X[test_data,] * beta$par
          if (length(beta$par) > 1) {
            mu = X[test_data,] %*% beta$par
          }
          probabilities = exp(mu) / (1 + exp(mu))
          #we achieve an estimate of the probability that the customer will default a loan
          predictions = probabilities > acceptance_prob
          model_error_rates[txtFormula] <- as.numeric(model_error_rates[txtFormula]) + size_fold - sum(y[test_data] == predictions)
        }
        
      }
    }
  
  }
  else{
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula = paste(txtFormula, comp[numComp])
      }
    else{
      txtFormula = paste(txtFormula, comp[numComp], sep = "+")
    }
    formula = formula(txtFormula)
    mf = model.frame(formula = formula, data = credit_data)
    X = model.matrix(attr(mf, "terms")[1:numComp], data = mf)
    y  = model.response(mf)
    model_error_rates[[txtFormula]] = 0
      for (fold in 1:nr_folds){
        #Select test data
        test_data <- ((fold-1)*size_fold + 1) : (fold*size_fold)
        #Select training data
        training_data <- 1:dim(credit_data)[1]
        training_data[test_data] = 0
        #Run estimate to achieve parameter
        beta <- IRLS(numeric(dim(credit_data)[1] - size_fold) + 1, y[training_data], X[training_data,])
        #Get matrix of X to calculate probabilities of test set
        mu <- X[test_data,] * beta$par
        if (length(beta$par) > 1) {
          mu = X[test_data,] %*% beta$par
        }
        probabilities = exp(mu) / (1 + exp(mu))
        #we achieve an estimate of the probability that the customer will default a loan
        predictions = probabilities > acceptance_prob
        model_error_rates[txtFormula] <- as.numeric(model_error_rates[txtFormula]) + size_fold - sum(y[test_data] == predictions)
      }
    }  
  }
  return (model_error_rates)
}
  
IRLS = function(n, y, X){
  z <- log(y + 0.5) - log(n - y + 0.5)
  beta <- solve(t(X) %*% X) %*% t(X) %*% z
  beta_prev = 10000000
  epsilon = 1e-5
  while(sum(abs(beta - beta_prev)) > epsilon) {
    eta <- X %*% beta;
    mu <- (exp(eta) / (1 + exp(eta)))*n;
    W <- as.numeric(mu*(n - mu)/n);   
    z <- eta + n*(y - mu)/(mu*(n - mu));
    beta_prev <- beta;              beta <- solve(t(X*W) %*% X) %*% t(X*W) %*% z;
  }
  beta_covmatrix = solve(t(X*W) %*% X)
  return (list(par = as.numeric(beta), covmatrix = beta_covmatrix))
}