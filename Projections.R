## Libraries ##

library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)
library(progress)
library(tidyverse)
library(plotly)
library(reshape2)

## Functions ##

source("functions_network.R")
source("functions_simulation.R")
col1 <- "cyan"
col2 <- "darkblue"
third <- rgb(63/255, 128/255, 97/255)
fourth <- "darkred"
grey <- "grey80"
custom_colorscale <- list(
  c(0, 1),
  c(col2, col1)
)

## Introductive analysis ##

data <- read.csv('Data/total_august.csv')
data_real <- read.csv('Data/august_LigurianSea.csv')

df <- data[which(data$year==2006 & data$RCP == 'rcp45'), 1:5]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
       aes(x=lon, y=lat, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
       aes(x = lon, y = lat, fill=value)) +
  geom_raster() + 
  geom_segment(data = df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + nord),
               arrow = arrow(length = unit(0.2, "cm")), color = third, alpha=0.8) + 
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = col2, high = col1, name = "Temperature")

inx <- which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value))

results <- compute_matrices(df = df)
lines <- results[[1]]
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]

g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
has_cycles <- !is_acyclic(g)
print("The built network has a cycle: ")
print(has_cycles)

df_real <- data_real[which(data_real$year==2006), 1:5]

ggplot(df_real[which(!is.na(df_real$temperature)) ,], 
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1) +
  theme_minimal()

ggplot(df_real[which(!is.na(df_real$east) & 
                       !is.na(df_real$nord) & !is.na(df_real$temperature)) ,], 
       aes(x = lon, y = lat, fill=temperature)) +
  geom_raster() + 
  geom_segment(data = df_real[which(!is.na(df_real$temperature)) ,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + nord),
               arrow = arrow(length = unit(0.2, "cm")), color = third, alpha=0.8) + 
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = col2, high = col1, name = "Temperature") + 
  theme_minimal()

results <- compute_matrices(df = df_real)
lines <- results[[1]]
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]

g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
has_cycles <- !is_acyclic(g)
print("The built network has a cycle: ")
print(has_cycles)

plot_lin_net(lines, df_real[, 1:2])

#### Initial analysis on the trend for the temperature ####

years <- unique(data$year)
rcp45 <- c()
rcp85 <- c()
for (y in years)
{
  df <- data[which(data$year==y & data$RCP == 'rcp45'), 1:5]
  rcp45 <- c(rcp45, mean(df$value, na.rm = T))
  df <- data[which(data$year==y & data$RCP == 'rcp85'), 1:5]
  rcp85 <- c(rcp85, mean(df$value, na.rm = T))
}

# x11()
# plot(years, rcp45, type = 'l', ylim = c(min(c(rcp45, rcp85))-2, max(c(rcp45, rcp85))+2),
#      col = third, ylab = 'Mean temperature', lwd = 5)
# lines(years, rcp85, col = fourth, lwd = 5)
# legend('topleft', legend = c('rcp 4.5', 'rcp 8.5'), fill = c(third, fourth))

#### Trying to extract the residuals ####

# Collect the real values, analyze the difference with the closest point in the simulation.
# Then evaluate the variogram for these residuals and analyze which is the best model.

#### rcp45 ####

#### Bistochastic ####

years <- unique(data_real$years)
len <- length(years)
means <- rep(0, len)
sill <- rep(0, len)
range <- rep(0, len)

variograms_45_bistochastic <- NULL

for(y in 1:len)
{
  df <- data[which(data$year==years[y] & data$RCP == 'rcp45'), 1:5]
  df_real <- data_real[which(data_real$year==years[y]), 1:5]
  results <- compute_matrices(df = df_real)
  dist <- results[[2]]
  PI_in <- results[[3]]
  PI_out <- results[[4]]
  g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
  has_cycles <- !is_acyclic(g)
  if(has_cycles)
    print(c("Iteration of year ", years[y], "has a cycle."))
  else
  {
    print(c("Comuputing iteration for year ", years[y]))
    inx <- which(!is.na(df_real$temperature))
    df_real <- df_real[inx,]
    dist <- dist[inx, inx]
    PI_in <- PI_in[inx, inx]
    PI_out <- PI_out[inx, inx]
    
    coords <- df_real[, 1:2]
    # I want a vector that for each location, find the nearest location based on data's locations
    
    coords_data <- expand.grid(unique(df$lon), unique(df$lat))
    n <- dim(coords_data)[1]
    
    nearest <- apply(coords, MARGIN = 1, function(x) 
      which.min(apply((coords_data - 
                         matrix(as.numeric(rep(x, times = n)), nrow = n, byrow = TRUE))^2, 
                      MARGIN = 1, sum)))
    
    inx2 <- which(!is.na(df$value[nearest]))
    
    res <- df_real$temperature[inx2] - df$value[nearest][inx2]
    means[y] <- mean(res)
    dist <- dist[inx2, inx2]
    PI_in <- PI_in[inx2, inx2]
    PI_out <- PI_out[inx2, inx2]
    
    P <- build_bistochastic(PI_out, 250, 2405)
    B <- dim(PI_out)[1]
    for(i in 1:B)
    {
      for(j in 1:B)
      {
        if(i!=j & !is.na(P[i,j]))
        {
          P[i,j] <- P[i,j]/sqrt(1-P[j,j])
        }
      }
    }
    P[which(is.na(P))] <- 0
    bistochastic <- initialize(dist, P)
    updated <- TRUE
    while (updated)
    {
      updated <- FALSE
      bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
    }
    
    variogram <- evaluate_variogram_penalization(res, bistochastic, l = 15, 
                                                 lambda = 0.2, return_fuv = T)
    fuv <- variogram[[1]]
    variogram <- variogram[[2]]
    variograms_45_bistochastic <- rbind(variograms_45_bistochastic, 
             c(variogram$squared_diff, rep(fuv, 15-length(variogram$squared_diff))))
    
    initial_params <- c(1, 3e5)
    params <- fit_variogram(variogram, spherical_kernel, initial_params)
    sill[y] <- params[1]
    range[y] <- params[2]
  }
}

inx <- (1:len)[-c(7,11)]
plot(years[inx], means[inx], pch = 16, col = third)
plot(years[inx], sill[inx], pch = 16, col = third)
plot(years[inx], range[inx], pch = 16, col = third)

mean_variogram_45_bistochastic <- data.frame(dist = variogram$dist, 
                             squared_diff = colMeans(variograms_45_bistochastic)[-15], 
                             np = variogram$np)

variograms_45_bistochastic <- melt(variograms_45_bistochastic)
colnames(variograms_45_bistochastic) <- c("Row", "Column", "Value")
variograms_45_bistochastic$Distance <- 
  round(variogram$dist[variograms_45_bistochastic$Column])

ggplot(variograms_45_bistochastic, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = col2) +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

ggplot(variograms_45_bistochastic, aes(x = factor(Distance), y = Value)) +
  geom_boxplot(fill = third) +
  geom_line(data = mean_variogram_45_bistochastic, 
            aes(x = factor(round(dist)), y = squared_diff, group = 1), col = fourth, 
            linewidth = 2) +
  labs(x = "Distance", y = "Values") +
  theme_minimal()

#### New Matrix ####

years <- unique(data_real$years)
len <- length(years)
means <- rep(0, len)
sill <- rep(0, len)
range <- rep(0, len)

variograms_45_newmatrix <- NULL

for(y in 1:len)
{
  df <- data[which(data$year==years[y] & data$RCP == 'rcp45'), 1:5]
  df_real <- data_real[which(data_real$year==years[y]), 1:5]
  results <- compute_matrices(df = df_real)
  dist <- results[[2]]
  PI_in <- results[[3]]
  PI_out <- results[[4]]
  g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
  has_cycles <- !is_acyclic(g)
  if(has_cycles)
    print(c("Iteration of year ", years[y], "has a cycle."))
  else
  {
    print(c("Comuputing iteration for year ", years[y]))
    inx <- which(!is.na(df_real$temperature))
    df_real <- df_real[inx,]
    dist <- dist[inx, inx]
    PI_in <- PI_in[inx, inx]
    PI_out <- PI_out[inx, inx]
    
    coords <- df_real[, 1:2]
    # I want a vector that for each location, find the nearest location based on data's locations
    
    coords_data <- expand.grid(unique(df$lon), unique(df$lat))
    n <- dim(coords_data)[1]
    
    nearest <- apply(coords, MARGIN = 1, function(x) 
      which.min(apply((coords_data - 
                         matrix(as.numeric(rep(x, times = n)), nrow = n, byrow = TRUE))^2, 
                      MARGIN = 1, sum)))
    
    inx2 <- which(!is.na(df$value[nearest]))
    
    res <- df_real$temperature[inx2] - df$value[nearest][inx2]
    means[y] <- mean(res)
    dist <- dist[inx2, inx2]
    PI_in <- PI_in[inx2, inx2]
    PI_out <- PI_out[inx2, inx2]
    
    P <- PI_in
    B <- dim(P)[1]
    for(i in 1:B)
    {
      for(j in 1:B)
      {
        if(i!=j & !is.na(P[i,j]))
        {
          P[i,j] <- sqrt(PI_in[i,j]*PI_out[i,j])
        }
      }
    }
    distances <- initialize(dist, P)
    updated <- TRUE
    while (updated)
    {
      updated <- FALSE
      distances <- lapply(distances, function(p) update(p,dist, P))
    }
    
    variogram <- evaluate_variogram_penalization(res, 
              distances, l = 15, lambda = 0.2, return_fuv = T)
    fuv <- variogram[[1]]
    variogram <- variogram[[2]]
    variograms_45_newmatrix <- rbind(variograms_45_newmatrix, 
               c(variogram$squared_diff, rep(fuv, 15-length(variogram$squared_diff))))
    
    initial_params <- c(1, 3e5)
    params <- fit_variogram(variogram, spherical_kernel, initial_params)
    sill[y] <- params[1]
    range[y] <- params[2]
  }
}

inx <- (1:len)[-c(7,11)]
plot(years[inx], means[inx], pch = 16, col = third)
plot(years[inx], sill[inx], pch = 16, col = third)
plot(years[inx], range[inx], pch = 16, col = third)

mean_variogram_45_newmatrix <- data.frame(dist = variogram$dist, 
                             squared_diff = colMeans(variograms_45_newmatrix)[-15], 
                             np = variogram$np)

variograms_45_newmatrix <- melt(variograms_45_newmatrix)
colnames(variograms_45_newmatrix) <- c("Row", "Column", "Value")
variograms_45_newmatrix$Distance <- 
  round(variogram$dist[variograms_45_newmatrix$Column])

ggplot(variograms_45_newmatrix, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = col2) +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

ggplot(variograms_45_newmatrix, aes(x = factor(Distance), y = Value)) +
  geom_boxplot(fill = third) +
  geom_line(data = mean_variogram_45_newmatrix, 
            aes(x = factor(round(dist)), y = squared_diff, group = 1), col = fourth, 
            linewidth = 2) +
  labs(x = "Distance", y = "Values") +
  theme_minimal()

#### Fit of the variogram ####

# The choice is to use the bistochastic model, because of the better shape and the higher
# sill

res <- fit_variogram(mean_variogram_45_bistochastic, 
                     spherical_kernel, initials = c(0.5, 2.5e5))

sill45 <- res[1]
range45 <- res[2]

#### rcp85 ####

#### Bistochastic ####

years <- unique(data_real$years)
len <- length(years)
means <- rep(0, len)
sill <- rep(0, len)
range <- rep(0, len)

variograms_85_bistochastic <- NULL

for(y in 1:len)
{
  df <- data[which(data$year==years[y] & data$RCP == 'rcp85'), 1:5]
  df_real <- data_real[which(data_real$year==years[y]), 1:5]
  results <- compute_matrices(df = df_real)
  dist <- results[[2]]
  PI_in <- results[[3]]
  PI_out <- results[[4]]
  g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
  has_cycles <- !is_acyclic(g)
  if(has_cycles)
    print(c("Iteration of year ", years[y], "has a cycle."))
  else
  {
    print(c("Comuputing iteration for year ", years[y]))
    inx <- which(!is.na(df_real$temperature))
    df_real <- df_real[inx,]
    dist <- dist[inx, inx]
    PI_in <- PI_in[inx, inx]
    PI_out <- PI_out[inx, inx]
    
    coords <- df_real[, 1:2]
    # I want a vector that for each location, find the nearest location based on data's locations
    
    coords_data <- expand.grid(unique(df$lon), unique(df$lat))
    n <- dim(coords_data)[1]
    
    nearest <- apply(coords, MARGIN = 1, function(x) 
      which.min(apply((coords_data - 
                         matrix(as.numeric(rep(x, times = n)), nrow = n, byrow = TRUE))^2, 
                      MARGIN = 1, sum)))
    
    inx2 <- which(!is.na(df$value[nearest]))
    
    res <- df_real$temperature[inx2] - df$value[nearest][inx2]
    means[y] <- mean(res)
    dist <- dist[inx2, inx2]
    PI_in <- PI_in[inx2, inx2]
    PI_out <- PI_out[inx2, inx2]
    
    P <- build_bistochastic(PI_out, 250, 2405)
    B <- dim(PI_out)[1]
    for(i in 1:B)
    {
      for(j in 1:B)
      {
        if(i!=j & !is.na(P[i,j]))
        {
          P[i,j] <- P[i,j]/sqrt(1-P[j,j])
        }
      }
    }
    P[which(is.na(P))] <- 0
    bistochastic <- initialize(dist, P)
    updated <- TRUE
    while (updated)
    {
      updated <- FALSE
      bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
    }
    
    variogram <- evaluate_variogram_penalization(res, 
                        bistochastic, l = 15, lambda = 0.2, return_fuv = T)
    fuv <- variogram[[1]]
    variogram <- variogram[[2]]
    variograms_85_bistochastic <- rbind(variograms_85_bistochastic, 
                     c(variogram$squared_diff, rep(fuv, 15-length(variogram$squared_diff))))
    
    initial_params <- c(1, 3e5)
    params <- fit_variogram(variogram, spherical_kernel, initial_params)
    sill[y] <- params[1]
    range[y] <- params[2]
  }
}

inx <- (1:len)[-c(7,11)]
plot(years[inx], means[inx], pch = 16, col = third)
plot(years[inx], sill[inx], pch = 16, col = third)
plot(years[inx], range[inx], pch = 16, col = third)

mean_variogram_85_bistochastic <- data.frame(dist = variogram$dist, 
             squared_diff = colMeans(variograms_85_bistochastic)[-15], 
                                             np = variogram$np)

variograms_85_bistochastic <- melt(variograms_85_bistochastic)
colnames(variograms_85_bistochastic) <- c("Row", "Column", "Value")
variograms_85_bistochastic$Distance <- 
  round(variogram$dist[variograms_85_bistochastic$Column])

ggplot(variograms_85_bistochastic, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = col2) +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

ggplot(variograms_85_bistochastic, aes(x = factor(Distance), y = Value)) +
  geom_boxplot(fill = third) +
  geom_line(data = mean_variogram_85_bistochastic, 
            aes(x = factor(round(dist)), y = squared_diff, group = 1), col = fourth, 
            linewidth = 2) +
  labs(x = "Distance", y = "Values") +
  theme_minimal()

#### New Matrix ####

years <- unique(data_real$years)
len <- length(years)
means <- rep(0, len)
sill <- rep(0, len)
range <- rep(0, len)

variograms_85_newmatrix <- NULL

for(y in 1:len)
{
  df <- data[which(data$year==years[y] & data$RCP == 'rcp85'), 1:5]
  df_real <- data_real[which(data_real$year==years[y]), 1:5]
  results <- compute_matrices(df = df_real)
  dist <- results[[2]]
  PI_in <- results[[3]]
  PI_out <- results[[4]]
  g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
  has_cycles <- !is_acyclic(g)
  if(has_cycles)
    print(c("Iteration of year ", years[y], "has a cycle."))
  else
  {
    print(c("Comuputing iteration for year ", years[y]))
    inx <- which(!is.na(df_real$temperature))
    df_real <- df_real[inx,]
    dist <- dist[inx, inx]
    PI_in <- PI_in[inx, inx]
    PI_out <- PI_out[inx, inx]
    
    coords <- df_real[, 1:2]
    # I want a vector that for each location, find the nearest location based on data's locations
    
    coords_data <- expand.grid(unique(df$lon), unique(df$lat))
    n <- dim(coords_data)[1]
    
    nearest <- apply(coords, MARGIN = 1, function(x) 
      which.min(apply((coords_data - 
                         matrix(as.numeric(rep(x, times = n)), nrow = n, byrow = TRUE))^2, 
                      MARGIN = 1, sum)))
    
    inx2 <- which(!is.na(df$value[nearest]))
    
    res <- df_real$temperature[inx2] - df$value[nearest][inx2]
    means[y] <- mean(res)
    dist <- dist[inx2, inx2]
    PI_in <- PI_in[inx2, inx2]
    PI_out <- PI_out[inx2, inx2]
    
    P <- PI_in
    B <- dim(P)[1]
    for(i in 1:B)
    {
      for(j in 1:B)
      {
        if(i!=j & !is.na(P[i,j]))
        {
          P[i,j] <- sqrt(PI_in[i,j]*PI_out[i,j])
        }
      }
    }
    distances <- initialize(dist, P)
    updated <- TRUE
    while (updated)
    {
      updated <- FALSE
      distances <- lapply(distances, function(p) update(p,dist, P))
    }
    
    variogram <- evaluate_variogram_penalization(res, 
                              distances, l = 15, lambda = 0.2, return_fuv = T)
    fuv <- variogram[[1]]
    variogram <- variogram[[2]]
    variograms_85_newmatrix <- rbind(variograms_85_newmatrix, 
                  c(variogram$squared_diff, rep(fuv, 15-length(variogram$squared_diff))))
    
    initial_params <- c(1, 3e5)
    params <- fit_variogram(variogram, spherical_kernel, initial_params)
    sill[y] <- params[1]
    range[y] <- params[2]
  }
}

inx <- (1:len)[-c(7,11)]
plot(years[inx], means[inx], pch = 16, col = third)
plot(years[inx], sill[inx], pch = 16, col = third)
plot(years[inx], range[inx], pch = 16, col = third)

mean_variogram_85_newmatrix <- data.frame(dist = variogram$dist, 
          squared_diff = colMeans(variograms_85_newmatrix)[-15], np = variogram$np)

variograms_85_newmatrix <- melt(variograms_85_newmatrix)
colnames(variograms_85_newmatrix) <- c("Row", "Column", "Value")
variograms_85_newmatrix$Distance <- 
  round(variogram$dist[variograms_85_newmatrix$Column])

ggplot(variograms_85_newmatrix, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = col2) +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

ggplot(variograms_85_newmatrix, aes(x = factor(Distance), y = Value)) +
  geom_boxplot(fill = third) +
  geom_line(data = mean_variogram_85_newmatrix, 
            aes(x = factor(round(dist)), y = squared_diff, group = 1), col = fourth, 
            linewidth = 2) +
  labs(x = "Distance", y = "Values") +
  theme_minimal()

#### Fit of the variogram ####

# The choice is to use the bistochastic model, because of the better shape and the higher
# sill

res <- fit_variogram(mean_variogram_85_bistochastic, 
                     spherical_kernel, initials = c(0.5, 2.5e5))

sill85 <- res[1]
range85 <- res[2]

#### Projections with the identified kernel ####

compute_matrices1 <- function(df) #This function is ok if the grid is regular - COPERNICUS
{
  coord_sub <- df[, 1:2]
  coord_utm <- coord_sub
  coordinates(coord_utm) <- c("lon","lat")
  proj4string(coord_utm) <- CRS("+proj=longlat +datum=WGS84")
  coord_utm <- spTransform(coord_utm, CRS("+proj=utm +zone=33 +datum=WGS84"))
  dist_eucl <- as.matrix(dist(coord_utm@coords, method = 'euclidean'))
  # dist_eucl <- as.matrix(dist(coord_sub, method = 'euclidean'))
  B <- dim(df)[1]
  len_lon <- length(unique(coord_sub[,1]))
  dist_mat <- matrix(NA, B,B)
  mat <- matrix(F, nrow=B, ncol=B)
  PI_in <- matrix(NA, B, B)
  PI_out <- matrix(NA, B, B)
  PI <- matrix(NA, B, B)
  lines <- list()
  
  positive_east <- !is.na(df$east) & df$east>0 & df$lon!= max(df$lon)
  positive_nord <- !is.na(df$nord) & df$nord>0 & df$lat!= max(df$lat)
  negative_east <- !is.na(df$east) & df$east<0 & df$lon!= min(df$lon)
  negative_nord <- !is.na(df$nord) & df$nord<0 & df$lat!= min(df$lat)
  
  lines <- list()
  for (i in 1:B) {
    line_id <- generate_line_id(coord_sub[i, ])
    lines[[i]] <- sp::Lines(slinelist = list(), ID = line_id)
  }
  
  # La distanza dalla riga alla colonna
  trues <- which(positive_east)
  for (i in trues)
  {
    temp <- list()
    line <- NULL
    line2 <- NULL
    if (df$east[i] > df$nord[i] || is.na(df$nord[i]))
    {
      dist_mat[i, i+1] <- dist_eucl[i, i+1]
      PI[i, i+1] <- df$east[i]
      line <- sp::Line(rbind(coord_sub[i,], coord_sub[i+1 ,]))
    }
    if (positive_nord[i])
    {
      dist_mat[i, (i+len_lon+1)] <- dist_eucl[i, (i+len_lon+1)]
      PI[i, (i+len_lon+1)] <- (sqrt(2)/2)*(abs(df$east[i])+df$nord[i])
      line2 <- sp::Line(rbind(coord_sub[i,], coord_sub[(i+len_lon+1) ,]))
    }
    temp <- list(line, line2)
    id <- paste(coord_sub[i,1], coord_sub[i,2])
    inx <- which(sapply(lines, function(line) attr(line, "ID") == id))
    if ((!is.null(temp[[1]]) || !is.null(temp[[2]])) & !length(lines[[inx]]@Lines))
    {
      lines[[inx]] <- sp::Lines(temp[!sapply(temp, is.null)], ID = id)
    }
    if ((!is.null(temp[[1]]) || !is.null(temp[[2]])) & length(lines[[inx]]@Lines))
    {
      lines[[inx]]@Lines <- c(lines[[inx]]@Lines, temp[!sapply(temp, is.null)])
    }
  }
  trues <- which(positive_nord)
  for (i in trues)
  {
    line <- NULL
    if (df$nord[i] > df$east[i] || is.na(df$east[i]))
    {
      dist_mat[i, (i+len_lon)] <- dist_eucl[i, (i+len_lon)]
      PI[i, (i+len_lon)] <- df$nord[i]
      line <- sp::Line(rbind(coord_sub[i,], coord_sub[(i+len_lon) ,]))
    }
    id <- paste(coord_sub[i,1], coord_sub[i,2])
    inx <- which(sapply(lines, function(line) attr(line, "ID") == id))
    if ((!is.null(line)) & !length(lines[[inx]]@Lines))
    {
      lines[[inx]] <- sp::Lines(line, ID = id)
    }
    if ((!is.null(line)) & length(lines[[inx]]@Lines))
    {
      lines[[inx]]@Lines <- c(lines[[inx]]@Lines, line)
    }
  }
  
  PI_out <- t(apply(PI, MARGIN =  1, normalize_column))
  PI_in <- apply(PI, 2, normalize_column)
  
  for (i in B:1)
  {
    if (!length(lines[[i]]@Lines))
      lines[[i]] <- NULL
  }
  
  return(list(lines, dist_mat, PI_in, PI_out))
}
compute_matrices2 <- function(df) #This function is ok if the grid is regular - COPERNICUS
{
  coord_sub <- df[, 1:2]
  coord_utm <- coord_sub
  coordinates(coord_utm) <- c("lon","lat")
  proj4string(coord_utm) <- CRS("+proj=longlat +datum=WGS84")
  coord_utm <- spTransform(coord_utm, CRS("+proj=utm +zone=33 +datum=WGS84"))
  dist_eucl <- as.matrix(dist(coord_utm@coords, method = 'euclidean'))
  # dist_eucl <- as.matrix(dist(coord_sub, method = 'euclidean'))
  B <- dim(df)[1]
  len_lon <- length(unique(coord_sub[,1]))
  dist_mat <- matrix(NA, B,B)
  mat <- matrix(F, nrow=B, ncol=B)
  PI_in <- matrix(NA, B, B)
  PI_out <- matrix(NA, B, B)
  PI <- matrix(NA, B, B)
  lines <- list()
  
  positive_east <- !is.na(df$east) & df$east>0 & df$lon!= max(df$lon)
  positive_nord <- !is.na(df$nord) & df$nord>0 & df$lat!= max(df$lat)
  negative_east <- !is.na(df$east) & df$east<0 & df$lon!= min(df$lon)
  negative_nord <- !is.na(df$nord) & df$nord<0 & df$lat!= min(df$lat)
  
  lines <- list()
  for (i in 1:B) {
    line_id <- generate_line_id(coord_sub[i, ])
    lines[[i]] <- sp::Lines(slinelist = list(), ID = line_id)
  }
  
  # La distanza dalla riga alla colonna
  trues <- which(negative_east)
  for (i in trues)
  {
    temp <- list()
    line <- NULL
    line2 <- NULL
    if (abs(df$east[i]) > (-df$nord[i]) || is.na(df$nord[i]))
    {
      dist_mat[i, i-1] <- dist_eucl[i, i-1]
      PI[i, i-1] <- abs(df$east[i])
      line <- sp::Line(rbind(coord_sub[i,], coord_sub[i-1 ,]))
    }
    if (negative_nord[i])
    {
      dist_mat[i, (i-len_lon-1)] <- dist_eucl[i, (i-len_lon-1)]
      PI[i, (i-len_lon-1)] <- (sqrt(2)/2)*(abs(df$east[i])+abs(df$nord[i]))
      line2 <- sp::Line(rbind(coord_sub[i,], coord_sub[(i-len_lon-1) ,]))
    }
    temp <- list(line, line2)
    id <- paste(coord_sub[i,1], coord_sub[i,2])
    inx <- which(sapply(lines, function(line) attr(line, "ID") == id))
    if ((!is.null(temp[[1]]) || !is.null(temp[[2]])) & !length(lines[[inx]]@Lines))
    {
      lines[[inx]] <- sp::Lines(temp[!sapply(temp, is.null)], ID = id)
    }
    if ((!is.null(temp[[1]]) || !is.null(temp[[2]])) & length(lines[[inx]]@Lines))
    {
      lines[[inx]]@Lines <- c(lines[[inx]]@Lines, temp[!sapply(temp, is.null)])
    }
  }
  trues <- which(negative_nord)
  for (i in trues)
  {
    line <- NULL
    if (abs(df$nord[i]) > (-df$east[i]) || is.na(df$east[i]))
    {
      dist_mat[i, (i-len_lon)]  <- dist_eucl[i, (i-len_lon)]
      PI[i, (i-len_lon)] <- abs(df$nord[i])
      line <- sp::Line(rbind(coord_sub[i,], coord_sub[(i-len_lon) ,]))
    }
    id <- paste(coord_sub[i,1], coord_sub[i,2])
    inx <- which(sapply(lines, function(line) attr(line, "ID") == id))
    if ((!is.null(line)) & !length(lines[[inx]]@Lines))
    {
      lines[[inx]] <- sp::Lines(line, ID = id)
    }
    if ((!is.null(line)) & length(lines[[inx]]@Lines))
    {
      lines[[inx]]@Lines <- c(lines[[inx]]@Lines, line)
    }
  }
  
  PI_out <- t(apply(PI, MARGIN =  1, normalize_column))
  PI_in <- apply(PI, 2, normalize_column)
  
  for (i in B:1)
  {
    if (!length(lines[[i]]@Lines))
      lines[[i]] <- NULL
  }
  
  return(list(lines, dist_mat, PI_in, PI_out))
}

#### Year 2050 ####

#### RCP 4.5 ####

df <- data[which(data$year==2050 & data$RCP == 'rcp45'), 1:5]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
       aes(x=lon, y=lat, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

# Constructing the bistochastic object

inx <- which(!is.na(df$value))

results <- compute_matrices1(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 500, 2405, control_bar = T)
B <- dim(P)[1]
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov1 <- build_covariances(sill45/2, range45, bistochastic, fun = spherical_covariance)
results <- compute_matrices2(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 500, 2405,control_bar = T)
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov2 <- build_covariances(sill45/2, range45, bistochastic, fun = spherical_covariance)
cov <- cov1 + cov2

simulation <- NULL
K <- 500
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (i in 1:500)
{
  simulation <- rbind(simulation,t(create_process_covariance(i, cov)))
  pb$tick()
}

simulation <- simulation + 
    matrix(rep(df$value[inx], times = 500), nrow = 500, ncol = length(inx), byrow = T)

means <- colMeans(simulation)
up <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.95))
low <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.05))

df1 <- df
df1$value[inx] <- means
df1$up <- NA
df1$down <- NA
df1$up[inx] <- up
df1$down[inx] <- low
df_filtered <- df1 %>% filter(!is.na(value))
lon_unique <- sort(unique(df_filtered$lon))
lat_unique <- sort(unique(df_filtered$lat))
values_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
upper_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
lower_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
for (i in 1:nrow(df_filtered)) {
  x_index <- which(lon_unique == df_filtered$lon[i])
  y_index <- which(lat_unique == df_filtered$lat[i])
  if (length(x_index) > 0 && length(y_index) > 0) {
    values_matrix[y_index, x_index] <- df_filtered$value[i]
    upper_matrix[y_index, x_index] <- df_filtered$up[i]
    lower_matrix[y_index, x_index] <- df_filtered$down[i]
  }
}

plot_ly() %>%
  # add_surface(z = values_matrix, x = lon_unique, y = lat_unique,
  #             colorscale = custom_colorscale,
  #             cmin = 22, cmax = 32,
  #             name = 'Values', showlegend = TRUE, showscale = FALSE, opacity = 0.5) %>%
  add_surface(z = upper_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Upper Bound', 
              showlegend = FALSE, showscale = TRUE) %>%
  add_surface(z = lower_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Lower Bound', 
              showlegend = FALSE, showscale = FALSE) %>%
  
  add_surface(z = matrix(27, nrow = nrow(values_matrix), ncol = ncol(values_matrix)), 
              x = lon_unique, y = lat_unique, 
              colorscale = list(c(0, 1), c("grey", "grey")), name = 'Temperature = 27', 
              opacity = 0.5, showlegend = FALSE, showscale = FALSE) %>%
  
  layout(scene = list(xaxis = list(title = 'Longitude'),
                      yaxis = list(title = 'Latitude'),
                      zaxis = list(title = 'Temperature', range = c(22, 32)),
                      camera = list(eye = list(x=-2,y=-2,z=1))),
          title = "Temperature - North Tyrrhenian - year 2050 - RCP 4.5")

#### RCP 8.5 ####

df <- data[which(data$year==2050 & data$RCP == 'rcp85'), 1:5]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
       aes(x=lon, y=lat, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

# Constructing the bistochastic object

inx <- which(!is.na(df$value))

results <- compute_matrices1(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 250, 2405)
B <- dim(P)[1]
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov1 <- build_covariances(sill85/2, range85, bistochastic, fun = spherical_covariance)
results <- compute_matrices2(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 250, 2405)
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov2 <- build_covariances(sill85/2, range85, bistochastic, fun = spherical_covariance)
cov <- cov1 + cov2

simulation <- NULL
K <- 500
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (i in 1:500)
{
  simulation <- rbind(simulation,t(create_process_covariance(i, cov)))
  pb$tick()
}

simulation <- simulation + 
  matrix(rep(df$value[inx], times = 500), nrow = 500, ncol = length(inx), byrow = T)

means <- colMeans(simulation)
up <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.95))
low <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.05))

df1 <- df
df1$value[inx] <- means
df1$up <- NA
df1$down <- NA
df1$up[inx] <- up
df1$down[inx] <- low
df_filtered <- df1 %>% filter(!is.na(value))
lon_unique <- sort(unique(df_filtered$lon))
lat_unique <- sort(unique(df_filtered$lat))
values_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
upper_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
lower_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
for (i in 1:nrow(df_filtered)) {
  x_index <- which(lon_unique == df_filtered$lon[i])
  y_index <- which(lat_unique == df_filtered$lat[i])
  if (length(x_index) > 0 && length(y_index) > 0) {
    values_matrix[y_index, x_index] <- df_filtered$value[i]
    upper_matrix[y_index, x_index] <- df_filtered$up[i]
    lower_matrix[y_index, x_index] <- df_filtered$down[i]
  }
}

plot_ly() %>%
  # add_surface(z = values_matrix, x = lon_unique, y = lat_unique,
  #             colorscale = custom_colorscale,
  #             cmin = 22, cmax = 32,
  #             name = 'Values', showlegend = TRUE, showscale = FALSE, opacity = 0.5) %>%
  add_surface(z = upper_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Upper Bound', 
              showlegend = FALSE, showscale = TRUE) %>%
  add_surface(z = lower_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Lower Bound', 
              showlegend = FALSE, showscale = FALSE) %>%
  
  add_surface(z = matrix(27, nrow = nrow(values_matrix), ncol = ncol(values_matrix)), 
              x = lon_unique, y = lat_unique, 
              colorscale = list(c(0, 1), c("grey", "grey")), name = 'Temperature = 27', 
              opacity = 0.5, showlegend = FALSE, showscale = FALSE) %>%
  
  layout(scene = list(xaxis = list(title = 'Longitude'),
                      yaxis = list(title = 'Latitude'),
                      zaxis = list(title = 'Temperature', range = c(22, 32)),
                      camera = list(eye = list(x=-2,y=-2,z=1))),
         title = "Temperature - North Tyrrhenian - year 2050 - RCP 8.5")

#### Year 2099 ####

#### RCP 4.5 ####

df <- data[which(data$year==2099 & data$RCP == 'rcp45'), 1:5]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
       aes(x=lon, y=lat, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

# Constructing the bistochastic object

inx <- which(!is.na(df$value))

results <- compute_matrices1(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 250, 2405)
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov1 <- build_covariances(sill45/2, range45, bistochastic, fun = spherical_covariance)
results <- compute_matrices2(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 250, 2405)
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov2 <- build_covariances(sill45/2, range45, bistochastic, fun = spherical_covariance)
cov <- cov1 + cov2

simulation <- NULL
K <- 500
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (i in 1:500)
{
  simulation <- rbind(simulation,t(create_process_covariance(i, cov)))
  pb$tick()
}

simulation <- simulation + 
  matrix(rep(df$value[inx], times = 500), nrow = 500, ncol = length(inx), byrow = T)

means <- colMeans(simulation)
up <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.95))
low <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.05))

df1 <- df
df1$value[inx] <- means
df1$up <- NA
df1$down <- NA
df1$up[inx] <- up
df1$down[inx] <- low
df_filtered <- df1 %>% filter(!is.na(value))
lon_unique <- sort(unique(df_filtered$lon))
lat_unique <- sort(unique(df_filtered$lat))
values_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
upper_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
lower_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
for (i in 1:nrow(df_filtered)) {
  x_index <- which(lon_unique == df_filtered$lon[i])
  y_index <- which(lat_unique == df_filtered$lat[i])
  if (length(x_index) > 0 && length(y_index) > 0) {
    values_matrix[y_index, x_index] <- df_filtered$value[i]
    upper_matrix[y_index, x_index] <- df_filtered$up[i]
    lower_matrix[y_index, x_index] <- df_filtered$down[i]
  }
}

plot_ly() %>%
  # add_surface(z = values_matrix, x = lon_unique, y = lat_unique,
  #             colorscale = custom_colorscale,
  #             cmin = 22, cmax = 32,
  #             name = 'Values', showlegend = TRUE, showscale = FALSE, opacity = 0.5) %>%
  add_surface(z = upper_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Upper Bound', 
              showlegend = FALSE, showscale = TRUE) %>%
  add_surface(z = lower_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Lower Bound', 
              showlegend = FALSE, showscale = FALSE) %>%
  
  add_surface(z = matrix(27, nrow = nrow(values_matrix), ncol = ncol(values_matrix)), 
              x = lon_unique, y = lat_unique, 
              colorscale = list(c(0, 1), c("grey", "grey")), name = 'Temperature = 27', 
              opacity = 0.5, showlegend = FALSE, showscale = FALSE) %>%
  
  layout(scene = list(xaxis = list(title = 'Longitude'),
                      yaxis = list(title = 'Latitude'),
                      zaxis = list(title = 'Temperature', range = c(22, 32)),
                        camera = list(eye = list(x=-2,y=-2,z=1))),
    title = "Temperature - North Tyrrhenian - year 2099 - RCP 4.5")

#### RCP 8.5 ####

df <- data[which(data$year==2099 & data$RCP == 'rcp85'), 1:5]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$value)) ,], 
       aes(x=lon, y=lat, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

# Constructing the bistochastic object

inx <- which(!is.na(df$value))

results <- compute_matrices1(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 250, 2405)
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov1 <- build_covariances(sill85/2, range85, bistochastic, fun = spherical_covariance)
results <- compute_matrices2(df = df)
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]
PI_in <- PI_in[inx, inx]
PI_out <- PI_out[inx, inx]
dist <- dist[inx, inx]
P <- build_bistochastic(PI_out, 250, 2405)
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(1-P[j,j])
    }
  }
}
P[which(is.na(P))] <- 0
bistochastic <- initialize(dist, P)
updated <- TRUE
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  bistochastic <- lapply(bistochastic, function(p) update(p,dist, P))
  order <- order + 1
}
cov2 <- build_covariances(sill85/2, range85, bistochastic, fun = spherical_covariance)
cov <- cov1 + cov2

simulation <- NULL
K <- 500
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (i in 1:500)
{
  simulation <- rbind(simulation,t(create_process_covariance(i, cov)))
  pb$tick()
}

simulation <- simulation + 
  matrix(rep(df$value[inx], times = 500), nrow = 500, ncol = length(inx), byrow = T)

means <- colMeans(simulation)
up <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.95))
low <- apply(simulation, MARGIN = 2, function(x) quantile(x, prob=.05))

df1 <- df
df1$value[inx] <- means
df1$up <- NA
df1$down <- NA
df1$up[inx] <- up
df1$down[inx] <- low
df_filtered <- df1 %>% filter(!is.na(value))
lon_unique <- sort(unique(df_filtered$lon))
lat_unique <- sort(unique(df_filtered$lat))
values_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
upper_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
lower_matrix <- matrix(NA, nrow = length(lat_unique), ncol = length(lon_unique))
for (i in 1:nrow(df_filtered)) {
  x_index <- which(lon_unique == df_filtered$lon[i])
  y_index <- which(lat_unique == df_filtered$lat[i])
  if (length(x_index) > 0 && length(y_index) > 0) {
    values_matrix[y_index, x_index] <- df_filtered$value[i]
    upper_matrix[y_index, x_index] <- df_filtered$up[i]
    lower_matrix[y_index, x_index] <- df_filtered$down[i]
  }
}

plot_ly() %>%
  # add_surface(z = values_matrix, x = lon_unique, y = lat_unique,
  #             colorscale = custom_colorscale,
  #             cmin = 22, cmax = 32,
  #             name = 'Values', showlegend = TRUE, showscale = FALSE, opacity = 0.5) %>%
  add_surface(z = upper_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Upper Bound', 
              showlegend = FALSE, showscale = TRUE) %>%
  add_surface(z = lower_matrix, x = lon_unique, y = lat_unique, 
              colorscale = custom_colorscale,
              cmin = 22, cmax = 32,  # Set z limits for scaling
              opacity = 0.7, name = 'Lower Bound', 
              showlegend = FALSE, showscale = FALSE) %>%
  
  add_surface(z = matrix(27, nrow = nrow(values_matrix), ncol = ncol(values_matrix)), 
              x = lon_unique, y = lat_unique, 
              colorscale = list(c(0, 1), c("grey", "grey")), name = 'Temperature = 27', 
              opacity = 0.5, showlegend = FALSE, showscale = FALSE) %>%
  
  layout(scene = list(xaxis = list(title = 'Longitude'),
                      yaxis = list(title = 'Latitude'),
                      zaxis = list(title = 'Temperature', range = c(22, 32)),  
         camera = list(eye = list(x=-2,y=-2,z=1))),
         title = "Temperature - North Tyrrhenian - year 2099 - RCP 8.5")

