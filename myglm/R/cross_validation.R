credit_data = read.table("https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment3/credit.txt",header=TRUE)


cross_validation = function(nr_folds){
  
  size_fold = ceiling(dim(credit_data)[1]/nr_folds)
  comp = c(colnames(credit_data)[1:5],
           paste(colnames(credit_data)[6:11], collapse = "+"),
           paste(colnames(credit_data)[12:17], collapse = "+"),
           paste(colnames(credit_data)[17:23], collapse = "+"))
  boolean = numeric(dim(credit_data)[1]) + 1

  # Name of response
  response = credit_data$Y
  # Fit the sequence of models
  txtFormula = paste("Y", "~", sep = "")
  model_error_rate=list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula = paste(txtFormula, comp[numComp])
      }
    else{
      txtFormula = paste(txtFormula, comp[numComp], sep = "+")
    }
    formula = formula(txtFormula)
    mf = model.frame(formula = formula, data = credit_data)
      for (fold in 2:nr_folds){
        #Select test data
        test_data = ((fold-1)*size_fold + 1) : (fold*size_fold)
        #Select training data
        training_data = 1:dim(credit_data)[1]
        training_data[test_data] = 0
        #Run estimate to achieve parameter
        current = myglm(formula=formula, data = na.omit(credit_data[training_data,]), family="geometric")
        #Get matrix of X to calculate probabilities of test set
        X  = model.matrix(attr(mf, "terms")[1:numComp], data = mf[test_data,])
        predictions = exp(X %*% current$par) / (1 + exp(X %*% current$par))
      }
  }
}
  
  
  #kj??rer over alle parameterene, her bruker vi samme framgangsm??te som anova funksjonen i myglm
  #"popper" dim(datasett)[1]/10
  #-> estimerer p-sannsynligheten, tilpasser til dataene i testfold. klassifiserer de etter hva som minimerer forventet tap
  #-> (alts?? ta h??yde for at en som defaulter blir merket ok er 10x dyrere enn motsatt, velg en p etter dette)
  #-> sjekker dette opp mot faktiske klassifikasjoner i folden og beregner hvor mye vi bommet
  #-> gj??r dette for alle foldene
  #-> itererer videre til en ekstra parameter
  #Det settet med parametere som

  #Mulig feilkilde: ?? avgj??re klassifikasjonen etter en tapsfunksjon blir litt underlig, da vi egentlig gj??r regresjon over
  #hva som er sannsynlig status. Gj??r vi klassifikasjonen over "midten" p = 0.5, s?? vil vi f?? klassifisert hva de mest sannsynlig er.
  #Gj??r vi det via en tapsfunksjon vil vi klassifisere slik det er lurt av bedriften ?? klassifisere dem.
