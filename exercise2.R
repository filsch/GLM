data.file = "https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment2/PremierLeague2015.txt"
d = read.table(data.file,
               col.names = c("home", "away", "x", "y"),
               colClasses = c("character", "character", "numeric","numeric"))

#testing independence
names = c('0','1','2','3','4+')
data = matrix(data <- numeric(25), nrow = 5,ncol = 5, dimnames = list(names,names))

for(i in 1:(dim(d)[1])){
  if(d[i,3] < 4 & d[i,4] < 4){
    data[d[i,4]+1, d[i,3]+1] = data[d[i,4]+1, d[i,3]+1] + 1
  }else if(d[i,3] < 4 & d[i,4] >= 4){
    data[5, d[i,3]+1] = data[5, d[i,3]+1] + 1
  }else if(d[i,3] >= 4 & d[i,4] < 4){
    data[d[i,4]+1, 5] = data[d[i,4]+1, 5] + 1
  }else{
    data[5,5] = data[5,5] + 1
  }
}
print(data)
O <- data
O = matrix(O, nrow = 5)
O = t(O)
print(O)
E = numeric(25)
print(E)
E = matrix(E, nrow = 5)
for(i in 1:5){
  for(j in 1:5){
    E[i,j] = sum(O[i,])*sum(O[,j])/sum(O)
  }
}
print(O)
print(E)
X = 0
for(i in 1:5){
  for(j in 1:5){
    X = X + (O[i,j] - E[i,j])^2/E[i,j]
  }
}
print(X)

#creating full ranking
createTable <- function(){
  data.file = "https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment2/PremierLeague2015.txt"
  d = read.table(data.file,
                 col.names = c("home", "away", "x", "y"),
                 colClasses = c("character", "character", "numeric","numeric"))
  teamnames <- c()
  points <- numeric(20)
  for(i in 1:(dim(d)[1])){
    if(is.element(d[i,1],teamnames) == FALSE){
      teamnames <- c(teamnames, d[i,1])
    }
  }
  for(i in 1:(dim(d)[1])){
    if(d[i,3] == d[i,4]){
      points[grep(d[i,1],teamnames)] = points[grep(d[i,1],teamnames)] + 1
      points[grep(d[i,2],teamnames)] = points[grep(d[i,2],teamnames)] + 1
    }
    else if(d[i,3] > d[i,4]){
      points[grep(d[i,1],teamnames)] = points[grep(d[i,1],teamnames)] + 3
    }
    else if(d[i,3] < d[i,4]){
      points[grep(d[i,2],teamnames)] = points[grep(d[i,2],teamnames)] + 3
    }
  }
  names(points) <- teamnames
  points = sort(points,decreasing = TRUE)
  return(points)
}

#creating X-matrix

createY <- function(){
  data.file = "https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment2/PremierLeague2015.txt"
  d = read.table(data.file,
                 col.names = c("home", "away", "x", "y"),
                 colClasses = c("character", "character", "numeric","numeric"))
  for(i in 1:dim(d)[1]){
    y <- c(y, d[i,3], d[i,4])
  }
  return(y)
}

createX <- function(){
  teamnames <- c()
  for(i in 1:(dim(d)[1])){
    if(is.element(d[i,1],teamnames) == FALSE){
      teamnames <- c(teamnames, d[i,1])
    }
  }
  colnames <- c("HomeAdvantage", teamnames)
  X <- numeric(2*dim(d)[1]*(length(teamnames) + 1))
  X <- matrix(X, nrow = 2*dim(d)[1], ncol = length(teamnames) + 1)
  colnames(X) <- colnames
  for(i in 1:(dim(X)[1])){
    if(i %% 2 == 1){
      X[i,1] = 1
      X[i+1,1] = 0
      X[i,grep(d[i,1],colnames)] = 1
      X[i+1,grep(d[i,2],colnames)] = 1
    }
  }
  return(X)
}
