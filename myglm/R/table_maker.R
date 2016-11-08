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
  points = data.frame(sort(points,decreasing = TRUE))
  colnames(points) <- "Ranking"
  return(points)
}
