createX <- function(simple){
  data.file = "https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment2/PremierLeague2015.txt"
  d = read.table(data.file,
                 col.names = c("home", "away", "x", "y"),
                 colClasses = c("character", "character", "numeric","numeric"))
  teamnames <- c()
  parameters <- c()
  for(i in 1:(dim(d)[1])){
    if(is.element(d[i,1],teamnames) == FALSE){
      teamnames <- c(teamnames, d[i,1])
      parameters <- c(parameters, paste(d[i,1],"Attack", sep=""), paste(d[i,1],"Defense", sep=""))
    }
  }
  if (simple){
  colnames <- c("HomeAdvantage", teamnames)
  X <- matrix(numeric(2*dim(d)[1]*(length(teamnames) + 1)), nrow = 2*dim(d)[1], ncol = length(teamnames) + 1)
  colnames(X) <- colnames
  y <- c()
  for(i in 1:dim(d)[1]){
    y <- c(y, d[i,3], d[i,4])
  }
  for(i in 1:(dim(d)[1])){
    X[2*i-1,1] = 1
    X[2*i,1] = 0
    X[2*i-1,grep(d[i,1],colnames)] = 1
    X[2*i,grep(d[i,2],colnames)] = 1
  }
  } else {
    X <- matrix(numeric(2*dim(d)[1]*(length(parameters) + 1)), nrow = 2*dim(d)[1], ncol = length(parameters) + 1)
    colnames(X) <- c("HomeAdvantage", parameters)
    y <- c()
    for(i in 1:dim(d)[1]){
      y <- c(y, d[i,3], d[i,4])
    }
    for(i in 1:(dim(d)[1])){
      X[2*i-1,1] = 1
      X[2*i,1] = 0
      X[2*i-1,grep(paste(d[i,1],"Attack", sep=""),colnames(X))] = 1
      X[2*i-1,grep(paste(d[i,2],"Defense", sep=""),colnames(X))] = 1
      X[2*i,grep(paste(d[i,2],"Attack", sep=""),colnames(X))] = 1
      X[2*i,grep(paste(d[i,1],"Defense", sep=""),colnames(X))] = 1
    }
  }
  return(list(X = X,Y = y, teamnames = teamnames))
}


season_simulation <- function(n_sim, simple){
   loglik = function(beta, Y, X) {
    return (sum(dpois(Y, lambda = exp(X%*%beta),  log=TRUE)))
   }
    if (simple){
    team_matrix = createX(simple=TRUE)
    teamnames = team_matrix$teamnames
    par = numeric(21);    names(par) <- c("Home advantage", teamnames)

    strength = optim(par, fn = loglik, method="BFGS", hessian = TRUE,
                      control = list(fnscale=-1), X = team_matrix$X, Y = team_matrix$Y)
    home_advantage = strength$par[1]

    simulated_positions <- matrix(numeric(n_sim*20), nrow = n_sim, ncol=20)
    points <- numeric(20)
    team_strength <- strength$par[2:21]

    names(points) = colnames(simulated_positions) <- teamnames
    reset <- points

    for (i in 1:n_sim){
        for (j in teamnames){
          for (k in teamnames){
            if (j != k){
            goals_home = rpois(1, exp(home_advantage + team_strength[[j]] - team_strength[[k]]))
            goals_away = rpois(1, exp(team_strength[[k]] - team_strength[[j]]))
            #cat('Goals, home team:', j, exp(home_advantage + team_strength[[j]] - team_strength[[k]]), '   Goals scored', goals_home, '\n')
            #cat('Goals, away team:', k, exp(team_strength[[k]] - team_strength[[j]]), '     Goals scored', goals_away, '\n')
              if(goals_home == goals_away){
                points[[j]] = points[[j]] + 1
                points[[k]] = points[[k]] + 1
              } else if (goals_home > goals_away){
                points[[j]] = points[[j]] + 3
              } else {
                points[[k]] = points[[k]] + 3
              }
            #cat('Points hometeam:    ', points[[j]])
            #readline("Continue ")
            }
          }
        }
      points = sort(points,decreasing = TRUE)
      place = 1
      for (l in names(points)){
        simulated_positions[[i,l]] = place
        place = place + 1
      }
      points = reset
    }
    } else {
      team_matrix = createX(simple=FALSE)
      teamnames = team_matrix$teamnames
      par = numeric(41);    names(par) <- c("Home advantage", attr(team_matrix$X,"dimnames")[[2]][2:41])

      strength = optim(par, fn = loglik, method="BFGS", hessian = TRUE,
                       control = list(fnscale=-1), X = team_matrix$X, Y = team_matrix$Y)
      home_advantage = strength$par[1]

      simulated_positions <- matrix(numeric(n_sim*20), nrow = n_sim, ncol=20)

      colnames(simulated_positions) <- teamnames
      points <- numeric(20);  names(points) <- teamnames
      reset = points

      team_strength <- strength$par[2:41]
      for (i in 1:n_sim){
        for (j in teamnames){
          for (k in teamnames){
            if (j != k){
              goals_home = rpois(1, exp(home_advantage + team_strength[[paste(j, "Attack",sep="")]] + team_strength[[paste(k, "Defense",sep="")]]))
              goals_away = rpois(1, exp(team_strength[[paste(k, "Attack",sep="")]] + team_strength[[paste(j, "Defense",sep="")]]))
              if(goals_home == goals_away){
                points[[j]] = points[[j]] + 1
                points[[k]] = points[[k]] + 1
              } else if (goals_home > goals_away){
                points[[j]] = points[[j]] + 3
              } else {
                points[[k]] = points[[k]] + 3
              }
            }
          }
        }
      points = sort(points,decreasing = TRUE)
      place = 1
      for (l in names(points)){
        simulated_positions[[i,l]] = place
        place = place + 1
      }
      points = reset
      }
    }
    total = data.frame(colMeans(simulated_positions,dims = 1),
                      100*colMeans(simulated_positions == 1,dims=1),
                      sqrt(colMeans(simulated_positions^2,dims = 1) - (colMeans(simulated_positions,dims = 1)^2)))
    rownames(total) = colnames(simulated_positions)
    total= total[order(total[,1]),]
    colnames(total) = c('Average position   ','% top position   ','Est. std. dev. on position   ')
    cat('Number of seasons simulations performed: ', n_sim , '\n \n')
    print(total,digits=3, right = FALSE)
}
