####
# Covariances' models
####

linear_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  ifelse(h <= range, sill * (1 - h / range), 0)
}
spherical_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  ifelse(h <= range, sill * (1 - (1.5 * (h / range)) + (0.5 * (h / range)^3)), 0)
}
exponential_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill * exp(-h/range))
}
mariah_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill * ((log(h/range + 1))/(h/range)))
}

####
# Corresponding models for the variograms
####

linear_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  ifelse(h <= range, sill - sill * (1 - h / range), sill)
}
spherical_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  ifelse(h <= range, sill - sill * (1 - (1.5 * (h / range)) + (0.5 * (h / range)^3)), sill)
}
exponential_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill - sill * exp(-h/range))
}
mariah_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill - sill * ((log(h/range + 1))/(h/range)))
}

####
# Simulation functions
####

build_covariances <- function(sill, range, paths, paths2 = NULL, fun)
{
  B <- length(paths)
  covs1 <- matrix(0, B, B)
  diag(covs1) <- sill
  
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      sums_1 <- sums_2 <- 0
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        l <- dim(mat)[1]
        for (k in 1:l)
        {
          sums_1 <- sums_1 + mat[k,2]*fun(c(sill, range), mat[k,1])
        }
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        l <- dim(mat)[1]
        for (k in 1:l)
        {
          sums_2 <- sums_2 + mat[k,2]*fun(c(sill, range), mat[k,1])
        }
      }
      covs1[i,j] <- covs1[j,i] <- sums_1 + sums_2
    }
  }
  
  if (!is.null(paths2))
  {
    covs2 <- matrix(0, B, B)
    diag(covs2) <- sill
    
    for (i in 1:(B-1))
    {
      for (j in (i+1):B)
      {
        sums_1 <- sums_2 <- 0
        if (!is.null(paths2[[i]]$lengths[[j]]))
        {
          mat <- paths2[[i]]$lengths[[j]]
          l <- dim(mat)[1]
          for (k in 1:l)
          {
            sums_1 <- sums_1 + mat[k,2]*fun(c(sill, range), mat[k,1])
          }
        }
        if (!is.null(paths2[[j]]$lengths[[i]]))
        {
          mat <- paths2[[j]]$lengths[[i]]
          l <- dim(mat)[1]
          for (k in 1:l)
          {
            sums_2 <- sums_2 + mat[k,2]*fun(c(sill, range), mat[k,1])
          }
        }
        covs2[i,j] <- covs2[j,i] <- sums_1 + sums_2
      }
    }
  }
  if (!is.null(paths2))
    return(list(covs1, covs2))
  else
    return(covs1)
}

create_process_covariance <- function(seed, cov)
{
  B <- dim(cov)[1]
  set.seed(seed)
  wn <- rnorm(B, 0, 1)
  eigV <- eigen(cov)
  S <- eigV$vectors %*% diag(sqrt(eigV$values)) %*% t(eigV$vectors)
  return(S %*% wn)
}

####
# Geostatistical anaysis
####

lengths_distribution <- function(paths)
{
  B <- length(paths)
  dists <- list()
  inx <- 1
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        dists[[inx]] <- mat[, 1]
        inx <- inx + 1
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        dists[[inx]] <- mat[, 1]
        inx <- inx + 1
      }
    }
  }
  return(dists)
}

evaluate_variogram_fcwa <- function(ma, paths, l, cutoff = NULL, return_fuv = F)
{
  B <- length(paths)
  if (B > 50)
  {
    set.seed(1)
    inx <- sample(1:B, 50)
  }
  else
  {
    inx <- 1:B
  }
  D <- length(inx)
  fuv <- NULL
  maxs <- 0
  for (i in 1:(D-1))
  {
    for (j in (i+1):D)
    {
      if (is.null(paths[[i]]$lengths[[j]]) & is.null(paths[[j]]$lengths[[i]]))
      {
        fuv <- c(fuv, (ma[i] - ma[j])^2)
      }
      else
      {
        if(is.null(paths[[i]]$lengths[[j]]))
          maxs <- max(maxs,max(paths[[j]]$lengths[[i]][, 1]))
        if(is.null(paths[[j]]$lengths[[i]]))
          maxs <- max(maxs,max(paths[[i]]$lengths[[j]][, 1]))
      }
    }
  }
  fuv <- mean(fuv) / 2
  lags <- seq(0, ifelse(is.null(cutoff), maxs, cutoff), length = l)
  cov <- vector("list", l)
  dists <- vector("list", l)
  np <- rep(0, l)
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        interval <- findInterval(mat[, 1], lags)
        estimation <- (ma[i] - ma[j]) ^ 2
        total <- sum(mat[, 2])
        tab <- table(interval)
        index <- as.numeric(names(tab)[which.max(tab)])
        indices <- which(interval == index)
        estimation <- (2 * fuv - estimation) / (2*sum(mat[indices, 2]))
        estimation <- fuv - estimation
        np[index] <- np[index]+1
        cov[[index]] <- rbind(cov[[index]], (estimation))
        dists[[index]] <- rbind(dists[[index]], mean(mat[indices,1]))
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        interval <- findInterval(mat[, 1], lags)
        estimation <- (ma[i] - ma[j]) ^ 2
        total <- sum(mat[, 2])
        tab <- table(interval)
        index <- as.numeric(names(tab)[which.max(tab)])
        indices <- which(interval == index)
        estimation <- (2 * fuv - estimation) / (2*sum(mat[indices, 2]))
        estimation <- fuv - estimation
        np[index] <- np[index]+1
        cov[[index]] <- rbind(cov[[index]], (estimation))
        dists[[index]] <- rbind(dists[[index]], mean(mat[indices,1]))
      }
    }
  }
  variogram <- unlist(lapply(cov, function(x) mean(x[which(x>=0 & x<=fuv)])))
  weights <- np/unlist(lapply(dists, mean, na.rm=T))
  if (any(is.na(weights)) || any(is.na(variogram)))
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)[-unique(c(which(is.na(weights)), 
                                                    which(is.na(variogram)))) ,]
  }
  else
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)
  }
  variogram <- rbind(c(0, 0, max(variogram$np)), variogram)
  if(return_fuv)
  {
    return(list(fuv, variogram))
  }
  else
    return(variogram)
}

evaluate_variogram_unadjusted <- function(ma, paths, l, cutoff)
{
  B <- length(paths)
  maxs <-0
  fuv <- NULL
  lags <- seq(0, cutoff, length = l)
  vars <- vector("list", l)
  dists <- vector("list", l)
  np <- rep(0, l)
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        estimation <- (ma[i] - ma[j]) ^ 2
        total <- sum(mat[, 2])
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        estimation <- (ma[i] - ma[j]) ^ 2
        total <- sum(mat[, 2])
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
    }
  }
  variogram <- unlist(lapply(vars, mean))/2
  weights <- np/unlist(lapply(dists, mean, na.rm=T))
  if (any(is.na(weights)) || any(is.na(variogram)))
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)[-unique(c(which(is.na(weights)), 
                                                    which(is.na(variogram)))) ,]
  }
  else
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)
  }
  # variogram <- rbind(c(0, 0, max(variogram$np)), variogram)
  return(variogram)
}

evaluate_variogram_penalization <- function(ma, paths, l, cutoff = NULL, return_fuv = F, lambda)
{
  B <- length(paths)
  if (B > 50)
  {
    set.seed(1)
    inx <- sample(1:B, 50)
  }
  else
  {
    inx <- 1:B
  }
  D <- length(inx)
  fuv <- NULL
  maxs <- 0
  for (i in 1:(D-1))
  {
    for (j in (i+1):D)
    {
      if (is.null(paths[[i]]$lengths[[j]]) & is.null(paths[[j]]$lengths[[i]]))
      {
        fuv <- c(fuv, (ma[i] - ma[j])^2)
      }
      else
      {
        if(is.null(paths[[i]]$lengths[[j]]))
          maxs <- max(maxs,max(paths[[j]]$lengths[[i]][, 1]))
        if(is.null(paths[[j]]$lengths[[i]]))
          maxs <- max(maxs,max(paths[[i]]$lengths[[j]][, 1]))
      }
    }
  }
  fuv <- mean(fuv) / 2
  lags <- seq(0, maxs, length = l)
  vars <- vector("list", l)
  omega <- vector("list", l)
  dists <- vector("list", l)
  np <- rep(0, l)
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        estimation <- (ma[i] - ma[j]) ^ 2
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        total <- sum(mat[, 2])
        omega[[interval]] <- rbind(omega[[interval]], total)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation/2)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        estimation <- (ma[i] - ma[j]) ^ 2
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        total <- sum(mat[, 2])
        omega[[interval]] <- rbind(omega[[interval]], total)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation/2)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
    }
  }
  weights <- unlist(lapply(omega, mean))
  variogram <- fuv - (fuv - unlist(lapply(vars, mean)))*(weights/(lambda + weights^2))
  weights <- np/unlist(lapply(dists, mean, na.rm=T))
  if (any(is.na(weights)) || any(is.na(variogram)))
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)[-unique(c(which(is.na(weights)), 
                                                    which(is.na(variogram)))) ,]
  }
  else
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)
  }
  variogram <- rbind(c(0, 0, max(variogram$np)), variogram)
  if(return_fuv)
  {
    return(list(fuv, variogram))
  }
  else
    return(variogram)
}

evaluate_variogram_penalization_best_lambda <- function(ma, paths, l, 
                                                        cutoff = NULL, return_fuv = F)
{
  B <- length(paths)
  if (B > 50)
  {
    set.seed(1)
    inx <- sample(1:B, 50)
  }
  else
  {
    inx <- 1:B
  }
  D <- length(inx)
  fuv <- NULL
  maxs <- 0
  for (i in 1:(D-1))
  {
    for (j in (i+1):D)
    {
      if (is.null(paths[[i]]$lengths[[j]]) & is.null(paths[[j]]$lengths[[i]]))
      {
        fuv <- c(fuv, (ma[i] - ma[j])^2)
      }
      else
      {
        if(is.null(paths[[i]]$lengths[[j]]))
          maxs <- max(maxs,max(paths[[j]]$lengths[[i]][, 1]))
        if(is.null(paths[[j]]$lengths[[i]]))
          maxs <- max(maxs,max(paths[[i]]$lengths[[j]][, 1]))
      }
    }
  }
  fuv <- mean(fuv) / 2
  lags <- seq(0, maxs, length = l)
  vars <- vector("list", l)
  omega <- vector("list", l)
  dists <- vector("list", l)
  np <- rep(0, l)
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        estimation <- (ma[i] - ma[j]) ^ 2
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        total <- sum(mat[, 2])
        omega[[interval]] <- rbind(omega[[interval]], total)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation/2)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        estimation <- (ma[i] - ma[j]) ^ 2
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        total <- sum(mat[, 2])
        omega[[interval]] <- rbind(omega[[interval]], total)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation/2)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
    }
  }
  weights <- unlist(lapply(omega, mean))
  variogram <- rep(NA, l)
  lambda <- NULL
  gamma_hat <- unlist(lapply(vars, mean))
  for (k in 1:l)
  {
    if(!is.na(gamma_hat[k]))
    {
      lam <- 0
      if (gamma_hat[k] <= fuv)
      {
        lam <- (weights[k] * (fuv - gamma_hat[k])/(fuv)) - (weights[k]^2)
        if (lam < 0)
        {
          lam <- 0
        }
      }
      else
        lam <- weights[k] - (weights[k]^2)
      lambda <- c(lambda, lam)
    }
  }
  lambda <- max(lambda)
  variogram <- fuv - (fuv - gamma_hat)*((weights)/(lambda + (weights^2)))

  weights <- np/unlist(lapply(dists, mean, na.rm=T))
  if (any(is.na(weights)) || any(is.na(variogram)))
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)[-unique(c(which(is.na(weights)), 
                                                    which(is.na(variogram)))) ,]
  }
  else
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)
  }
  # variogram <- rbind(c(0, 0, max(variogram$np)), variogram)
  if(return_fuv)
  {
    return(list(fuv, variogram))
  }
  else
    return(variogram)
}

evaluate_variogram_penalization_train <- function(ma, paths, train, l, cutoff = NULL, return_fuv = F, lambda)
{
    B <- length(train)
    if (B > 50)
    {
      set.seed(1)
      inx <- sample(1:B, 50)
    }
    else
    {
      inx <- 1:B
    }
    D <- length(inx)
    fuv <- NULL
    maxs <- 0
    for (in1 in 1:(D-1))
    {
      for (in2 in (in1+1):D)
      {
        i <- train[in1]
        j <- train[in2]
        if (is.null(paths[[i]]$lengths[[j]]) & is.null(paths[[j]]$lengths[[i]]))
        {
          fuv <- c(fuv, (ma[i] - ma[j])^2)
        }
        else
        {
          if(is.null(paths[[i]]$lengths[[j]]))
            maxs <- max(maxs,max(paths[[j]]$lengths[[i]][, 1]))
          if(is.null(paths[[j]]$lengths[[i]]))
            maxs <- max(maxs,max(paths[[i]]$lengths[[j]][, 1]))
        }
      }
    }
    fuv <- mean(fuv) / 2
    mean <- mean(ma)
    lags <- seq(0, maxs, length = l)
    vars <- vector("list", l)
    omega <- vector("list", l)
    dists <- vector("list", l)
    np <- rep(0, l)
    for (in1 in 1:(B-1))
    {
      for (in2 in (in1+1):B)
      {
        i <- train[in1]
        j <- train[in2]
        if (!is.null(paths[[i]]$lengths[[j]]))
        {
          mat <- paths[[i]]$lengths[[j]]
          estimation <- (ma[i] - ma[j]) ^ 2
          dist <- mean(mat[, 1])
          interval <- findInterval(dist, lags)
          total <- sum(mat[, 2])
          omega[[interval]] <- rbind(omega[[interval]], total)
          np[interval] <- np[interval]+1
          vars[[interval]] <- rbind(vars[[interval]], estimation/2)
          dists[[interval]] <- rbind(dists[[interval]], dist)
        }
        if (!is.null(paths[[j]]$lengths[[i]]))
        {
          mat <- paths[[j]]$lengths[[i]]
          estimation <- (ma[i] - ma[j]) ^ 2
          dist <- mean(mat[, 1])
          interval <- findInterval(dist, lags)
          total <- sum(mat[, 2])
          omega[[interval]] <- rbind(omega[[interval]], total)
          np[interval] <- np[interval]+1
          vars[[interval]] <- rbind(vars[[interval]], estimation/2)
          dists[[interval]] <- rbind(dists[[interval]], dist)
        }
      }
    }
    weights <- unlist(lapply(omega, mean))
    variogram <- fuv - (fuv - unlist(lapply(vars, mean)))*(weights/(lambda + weights^2))
    weights <- np/unlist(lapply(dists, mean, na.rm=T))
    if (any(is.na(weights)) || any(is.na(variogram)))
    {
      variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                              squared_diff = variogram, 
                              np = weights)[-unique(c(which(is.na(weights)), 
                                                      which(is.na(variogram)))) ,]
    }
    else
    {
      variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                              squared_diff = variogram, 
                              np = weights)
    }
    variogram <- rbind(c(0, 0, max(variogram$np)), variogram)
    if(return_fuv)
    {
      return(list(fuv, variogram))
    }
    else
      return(variogram)
}

evaluate_variogram_unadjusted_train <- function(ma, train, paths, l, cutoff)
{
  B <- length(train)
  maxs <- 0
  fuv <- NULL
  lags <- seq(0, cutoff, length = l)
  vars <- vector("list", l)
  dists <- vector("list", l)
  np <- rep(0, l)
  for (in1 in 1:(B-1))
  {
    for (in2 in (in1+1):B)
    {
      i <- in1
      j <- in2
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        estimation <- (ma[i] - ma[j]) ^ 2
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        total <- sum(mat[, 2])
        omega[[interval]] <- rbind(omega[[interval]], total)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation/2)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        estimation <- (ma[i] - ma[j]) ^ 2
        dist <- mean(mat[, 1])
        interval <- findInterval(dist, lags)
        total <- sum(mat[, 2])
        omega[[interval]] <- rbind(omega[[interval]], total)
        np[interval] <- np[interval]+1
        vars[[interval]] <- rbind(vars[[interval]], estimation/2)
        dists[[interval]] <- rbind(dists[[interval]], dist)
      }
    }
  }
  variogram <- unlist(lapply(vars, mean))/2
  weights <- np/unlist(lapply(dists, mean, na.rm=T))
  if (any(is.na(weights)) || any(is.na(variogram)))
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)[-unique(c(which(is.na(weights)), 
                                                    which(is.na(variogram)))) ,]
  }
  else
  {
    variogram <- data.frame(dist = unlist(lapply(dists, mean, na.rm=T)), 
                            squared_diff = variogram, 
                            np = weights)
  }
  variogram <- rbind(c(0, 0, max(variogram$np)), variogram)
  return(variogram)
}

fit_variogram <- function(var, kernel, initials)
{
  objective_function <- function(params, h, observed_variogram, np) {
    sill <- params[1]
    range <- params[2]
    predicted_variogram <- kernel(params, h)
    weights <- sqrt(np)
    return(sum(weights * (observed_variogram - predicted_variogram)^2))
  }
  
  optim_result <- optim(par = initials, fn = objective_function, 
                        h = var$dist, 
                        observed_variogram = var$squared_diff,
                        np = var$np)
  
  return(c(optim_result$par[1], optim_result$par[2]))
}

loo_predictions <- function(ma, covs, covs2 = NULL, eucl = F)
{
  D <- length(ma)
  ones <- rep(1, D-1)
  preds1 <- preds2 <- matrix(0, D, 1)
  vars1 <- vars2 <- matrix(0, D, 1)
  for (i in 1:D)
  {
    sigma <- covs[-i,-i]
    cc <- covs[i, -i]
    lam <- solve(sigma, cbind(cc, ones))
    ViX = lam[,-1]
    skwts = lam[,1]
    beta = solve(t(ones) %*% ViX, t(ViX) %*% ma[-i])
    preds1[i] <- 1 %*% beta + t(skwts) %*% (ma[-i] - c(1 %*% beta))
    
    Q = t(1) - t(ViX) %*% cc
    vars1[i] <- covs[1,1] - apply(cbind(cc*skwts), 2, sum) + 
      apply(Q * solve(t(ones) %*% ViX, Q), 2, sum) 
    
    if (!is.null(covs2))
    {
      sigma <- covs2[-i,-i]
      cc <- covs2[i, -i]
      lam <- solve(sigma, cbind(cc, ones))
      ViX = lam[,-1]
      skwts = lam[,1]
      beta = solve(t(ones) %*% ViX, t(ViX) %*% ma[-i])
      preds2[i] <- 1 %*% beta + t(skwts) %*% (ma[-i] - c(1 %*% beta))
      
      Q = t(1) - t(ViX) %*% cc
      vars2[i] <- covs2[1,1] - apply(cbind(cc*skwts), 2, sum) + 
        apply(Q * solve(t(ones) %*% ViX, Q), 2, sum) 
    }
  }
  mean <- mean(ma)
  pred_eucl <- data.frame(cbind(rep(mean, D), rep(var(ma), D)))
  names(pred_eucl) <- c("preds", "vars")
  pred_1 <- data.frame(cbind(preds1, vars1))
  names(pred_1) <- c("preds", "vars")
  if(!is.null(covs2))
  {
    pred_2 <- data.frame(cbind(preds2, vars2))
    names(pred_2) <- c("preds", "vars")
    if (eucl)
      return(list(pred_1, pred_2, pred_eucl))
    else
      return(return(list(pred_1, pred_2)))
  }
  else
  {
    if (eucl)
      return(list(pred_1, pred_eucl))
    else
      return(pred_1)
  }
}

compute_cv_indices <- function(ma, covs, covs2 = NULL, eucl = F)
{
  D <- length(ma)
  ones <- rep(1, D-1)
  err1 <- err2 <- 0
  preds1 <- preds2 <- matrix(0, D, 1)
  for (i in 1:D)
  {
    sigma <- covs[-i,-i]
    cc <- covs[i, -i]
    lam <- solve(sigma, cbind(cc, ones))
    ViX = lam[,-1]
    skwts = lam[,1]
    beta = solve(t(ones) %*% ViX, t(ViX) %*% ma[-i])
    preds1[i] <- 1 %*% beta + t(skwts) %*% (ma[-i] - c(1 %*% beta))
    
    if (!is.null(covs2))
    {
      sigma <- covs2[-i,-i]
      cc <- covs2[i, -i]
      lam <- solve(sigma, cbind(cc, ones))
      ViX = lam[,-1]
      skwts = lam[,1]
      beta = solve(t(ones) %*% ViX, t(ViX) %*% ma[-i])
      preds2[i] <- 1 %*% beta + t(skwts) %*% (ma[-i] - c(1 %*% beta))
    }
  }
  mean <- mean(ma)
  err <- mean((ma-mean)^2)
  err1 <- mean((ma - preds1)^2)
  if(!is.null(covs2))
  {
    err2 <- mean((ma - preds2)^2)
    if (eucl)
      return(c(err1, err2, err))
    else
      return(return(c(err1, err2)))
  }
  else
  {
    if (eucl)
      return(c(err1, err))
    else
      return(err1)
  }
}

evaluate_test <- function(values, covs, train, test)
{
  sigma <- covs[train,train]
  D <- length(test)
  B <- length(train)
  ones <- rep(1, B)
  predictions <- matrix(0, 1, D)
  vars <- matrix(0, 1, D)
  for (i in 1:D)
  {
    cc <- covs[test[i], train]
    lam <- solve(sigma, cbind(cc, ones))
    ViX = lam[,-1]
    skwts = lam[,1]
    beta = solve(t(ones) %*% ViX, t(ViX) %*% values[train])
    predictions[i] <- 1 %*% beta + t(skwts) %*% (values[train] - c(1 %*% beta))
    
    Q = t(1) - t(ViX) %*% cc
    vars[i] <- covs[1,1] - apply(cbind(cc*skwts), 2, sum) + 
      apply(Q * solve(t(ones) %*% ViX, Q), 2, sum)  
  }
  
  return(data.frame(preds = t(predictions), vars = t(vars)))
}

