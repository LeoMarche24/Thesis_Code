library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)
library(progress)

# Set up functions
source("functions_network.R")
source("functions_simulation.R")
col1 <- "cyan"
col2 <- "darkblue"
third <- rgb(63/255, 128/255, 97/255)
fourth <- rgb(173/255, 38/255, 36/255)
grey <- "grey80"

#### Building the network ####

df <- read.csv('Data/RealData_August_LigurianSea.csv')
df <- df[which(df$years==2006) ,]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x = lon, y = lat, fill=temperature)) +
  geom_raster() + 
  geom_segment(data = df[which(!is.na(df$temperature)) ,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + nord),
               arrow = arrow(length = unit(0.3, "cm")), color = third, alpha=1, 
               linewidth = 2) + 
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = col2, high = col1, name = "Temperature") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text = element_text(size = 15)
  )

# Computing the matrices of proportional influence both upstream and downstream, the lines
# and the distance matrix characterizing the network.
# Note that these functions and this procedure is valid if the data are provided in a
# regular grid and are sort in that sense. Do not use if this is not the case

results <- compute_matrices(df = df)
lines <- results[[1]]
dist <- results[[2]]
PI_in <- results[[3]]
PI_out <- results[[4]]

g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
has_cycles <- !is_acyclic(g)
print("The built network has a cycle: ")
print(has_cycles)

plot_lin_net(lines, df[, 1:2])

# Three objects for the analysis: 
# 1. dist: matrix containing the distances for the edges in the network.
# 2. PI_out: matrix contaning the probability of choosing that edge, markov dynamics.
# 3. PI_in: matrix normalized on the columns
# 3. df: dataframe with the informations on the points of the domain.

#### Weight and distances object ####

# Create the object, it contains all the connections between two points both in the 
# length of the path and in the weight, depending on the model

# Normalized matrix process: object named normalized

inx <- which(!is.na(df$temperature))

dist <- dist[inx, inx]
PI_out <- PI_out[inx, inx]
PI_in <- PI_in[inx, inx]

P <- PI_out
B <- dim(P)[1]
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- P[i,j]/sqrt(sum(PI_out[,j], na.rm = T))
    }
  }
}
normalized <- initialize(dist, P)
updated <- TRUE
order <- 1
max_update <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  normalized <- lapply(normalized, function(p) update(p,dist, P))
  order <- order + 1
}

B <- length(normalized)

# New matrix for the upstream dynamics, object named NewMatrix

P <- PI_in
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
NewMatrix <- initialize(dist, P)
updated <- TRUE
order <- 1
max_update <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  NewMatrix <- lapply(NewMatrix, function(p) update(p,dist, P))
  order <- order + 1
}

# Dynamics ruled by a bistochastic matrix, object named bistochastic

P <- build_bistochastic(PI_out, 1000, 2405, control_bar = T)
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

#### Variogram estimation ####

#### Test the performance of the new estimator

# Normalized-based process - simulate a realization of the process 
# and then evaluate the variogram

sill <- 15
range <- 2e5
cov <- build_covariances(sill, range, normalized, fun = linear_covariance)
simulation <- create_process_covariance(18, cov)
df$Simulation <- NA
df$Simulation[inx] <- simulation
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Simulation)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

variogram <- evaluate_variogram_fcwa(simulation, normalized, 15, 5e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

variogram <- evaluate_variogram_unadjusted(simulation, normalized, 15, 5e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

grid <- seq(0.001, 1, length = 10)
g <- length(grid)
variogram_lambda <- lapply(grid, FUN = function(x) 
  evaluate_variogram_penalization(simulation, normalized, 15, lambda = x))
xx <- seq(0, max(variogram$dist), length = 200)
yy <- linear_kernel(c(sill, range), xx)
plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l', 
     xlab = "Distances", ylab = "Semivariogram")
for (i in 1:g)
  lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff, 
        col = colorRampPalette(c(col2, col1))(g)[i], lwd = 1.5)
lines(variogram$dist, variogram$squared_diff, col = fourth, lwd = 3)

# x11()
# plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l', 
#      xlab = "Distances", ylab = "Semivariogram")
# for (i in 1:g)
#   lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff, 
#         col = colorRampPalette(c(col2, col1))(g)[i], lwd = 1.5)
# lines(variogram$dist, variogram$squared_diff, col = fourth, lwd = 3)
# legend('topleft', legend = c(grid, "unweighted"), 
#        fill = c(colorRampPalette(c(col2, col1))(g), fourth))

# Upstream dynamics ruled by a new matrix process - simulate a realization of the 
# process and then evaluate the variogram

sill <- 15
range <- 1e5
cov <- build_covariances(sill, range, NewMatrix, fun = linear_covariance)
simulation <- create_process_covariance(2405, cov)
df$Simulation <- NA
df$Simulation[inx] <- simulation
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Simulation)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20) 
  )

variogram <- evaluate_variogram(simulation, NewMatrix, 15)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

variogram <- evaluate_variogram_unadjusted(simulation, NewMatrix, 15, 5e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

grid <- seq(0.001, 1, length = 10)
g <- length(grid)
variogram_lambda <- lapply(grid, FUN = function(x) 
  evaluate_variogram_penalization(simulation, NewMatrix, 15, lambda = x))
xx <- seq(0, max(variogram$dist), length = 200)
yy <- exponential_kernel(c(sill, range), xx)
plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l', 
     xlab = "Distances", ylab = "Semivariogram")
for (i in 1:g)
  lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff, 
        col = colorRampPalette(c(col2, col1))(g)[i], lwd = 1.5)
lines(variogram$dist, variogram$squared_diff, col = fourth, lwd = 3)
# 
# x11()
# plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l',
#      xlab = "Distances", ylab = "Semivariogram")
# for (i in 1:g)
#   lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff,
#         col = colorRampPalette(c(col2, col1))(g)[i], lwd = 1.5)
# lines(variogram$dist, variogram$squared_diff, col = fourth, lwd = 3)
# legend('topleft', legend = c(grid, "unweighted"),
#        fill = c(colorRampPalette(c(col2, col1))(g), fourth))

# Bistochastic process - simulate a realization of the 
# process and then evaluate the variogram

sill <- 15
range <- 2e5
cov <- build_covariances(sill, range, bistochastic, fun = spherical_covariance)
simulation <- create_process_covariance(2405, cov)
df$Simulation <- NA
df$Simulation[inx] <- simulation
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Simulation)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20) 
  )

variogram <- evaluate_variogram_unadjusted(simulation, bistochastic, 10, 5e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

grid <- seq(0.01, 1, length = 5)
g <- length(grid)
variogram_lambda <- lapply(grid, FUN = function(x) 
  evaluate_variogram_penalization(simulation, bistochastic, 15, lambda = x))
xx <- seq(0, max(variogram$dist), length = 200)
yy <- spherical_kernel(c(sill, range), xx)
plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l', 
     xlab = "Distances", ylab = "Semivariogram")
for (i in 1:g)
  lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff, 
        col = colorRampPalette(c(col2, col1))(g)[i], lwd = 1.5)
lines(c(0,variogram$dist), c(0,variogram$squared_diff), col = fourth, lwd = 3)

xx <- seq(0, max(unlist(lapply(variogram_lambda, function(x) x$dist))), length = 200)
yy <- spherical_kernel(c(sill, range), xx)

variogram_best <- evaluate_variogram_penalization_best_lambda(simulation, bistochastic, 15)

# x11()
# plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l',
#      xlab = "Distances", ylab = "Semivariogram")
# for (i in 1:g)
#   lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff,
#         col = colorRampPalette(c(grey, third))(g)[i], lwd = 2)
# lines(c(0, variogram$dist), c(0, variogram$squared_diff), col = fourth, lwd = 3)
# legend('topleft', legend = c(grid, "Unweighted", "True semi-variogram"),
#        fill = c(colorRampPalette(c(grey, third))(g), fourth, "black"))

# x11()
# plot(xx, yy, col = 'black', lwd = 2, ylim = c(0, 2*sill), type = 'l',
#      xlab = "Distances (meters)", ylab = "Semivariogram", cex.lab = 2, cex.axis = 2)
# for (i in 1:g)
#   lines(variogram_lambda[[i]]$dist, variogram_lambda[[i]]$squared_diff,
#         col = colorRampPalette(c(grey, third))(g)[i], lwd = 2)
# lines(c(0, variogram$dist), c(0, variogram$squared_diff), col = fourth, lwd = 3)
# lines(c(0, variogram_best$dist), c(0, variogram_best$squared_diff), col = col2, lwd = 3)
# legend("topleft",
#        legend = grid, # The gradient for 'grid' lines
#        fill = colorRampPalette(c(grey, third))(g),
#        bty = "n",     # No border
#        cex = 1.3)     # Increase text size
# legend("top",
#        legend = c("Unweighted", "Best lambda", "True semi-variogram"),
#        fill = c(fourth, col2, "black"),
#        bty = "n",     # No border
#        cex = 1.3,     # Increase text size
#        x.intersp = 0.5,  # Adjusts horizontal spacing
#        inset = c(0.1, 0))  # Slightly shifts to the right

#### LOO - cross validation ####

# Assessing the ability of the model to predict at unsample locations.
# Moreover, a deeper analysis to the variogram estimate.

#### Normalized process ####

B <- length(normalized)
K <- 500
sill <- 15
range <- 2e5
initial_params <- c(sill, range)
cov <- build_covariances(sill, range, normalized, fun = spherical_covariance)
errors_nor_unwe <- matrix(0, nrow = K, ncol = B)
errors_nor_pen <- matrix(0, nrow = K, ncol = B)
errors_new_unwe <- matrix(0, nrow = K, ncol = B)
errors_new_pen <- matrix(0, nrow = K, ncol = B)
errors_eucl <- matrix(0, nrow = K, ncol = B)
variograms_nor <- list('vector', B)
variograms_new <- list('vector', B)
variograms_bist <- list('vector', B)
variograms_unwe <- list('vector', B)
sills <- rep(0, K)
sills_unwe <- rep(0, K)
ranges_nor <- rep(0, K)
ranges_new <- rep(0, K)
ranges_bist <- rep(0, K)
ranges_unwe <- rep(0, K)
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)

for (iter in 1:K)
{
  ma <- create_process_covariance(iter, cov)
  
  variogram <- evaluate_variogram_penalization(ma, normalized, 15, 
                                               return_fuv = T, lambda = lam_nor)
  fuv <- variogram[[1]]
  variogram <- variogram[[2]]
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_sill <- fuv
  estimated_range <- params[2]
  sills[iter] <- estimated_sill
  ranges_nor[iter] <- estimated_range
  covariances_nor <- build_covariances(estimated_sill, estimated_range, normalized, 
                                       fun = spherical_covariance)
  variograms_nor[[iter]] <- variogram
  
  variogram <- evaluate_variogram_penalization(ma, NewMatrix, 15, lambda = lam_new)
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_range <- params[2]
  ranges_new[iter] <- estimated_range
  covariances_new <- build_covariances(estimated_sill, estimated_range, NewMatrix, 
                                       fun = spherical_covariance)
  variograms_new[[iter]] <- variogram
  
  predictions <- loo_predictions(ma, covariances_nor, covariances_new, eucl = T)
  
  aux <- predictions[[1]]
  errors_nor_pen[iter,] <- (aux$preds - ma)^2
  aux <- predictions[[2]]
  errors_new_pen[iter,] <- (aux$preds - ma)^2
  aux <- predictions[[3]]
  errors_eucl[iter,] <- (aux$preds - ma)^2
  
  variogram <- evaluate_variogram_penalization(ma, bistochastic, 15, lambda = lam_bist)
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_range <- params[2]
  ranges_bist[iter] <- estimated_range
  variograms_bist[[iter]] <- variogram
    
  variogram <- evaluate_variogram_unadjusted(ma, normalized, 15, 6e5)
  params <- fit_variogram(variogram, spherical_kernel, initials = initial_params)
  sills_unwe[iter] <- params[1]
  ranges_unwe[iter] <- params[2]
  variograms_unwe[[iter]] <- variogram
  
  covariances_nor <- build_covariances(params[1], params[2], normalized, 
                                       fun = spherical_covariance)
  covariances_new <- build_covariances(params[1], params[2], NewMatrix, 
                                       fun = spherical_covariance)
  predictions <- loo_predictions(ma, covariances_nor, covariances_new, eucl = F)
  errors_nor_unwe[iter, ] <- (predictions[[1]]$preds-ma)^2
  errors_new_unwe[iter, ] <- (predictions[[2]]$preds-ma)^2
  
  pb$tick()
}

# Sill estimate
# Note that the sill estimate for the penalized estimator is equivalent whatever model
# is used

data <- data.frame(
  value = c(sills, sills_unwe),
  group = factor(c(rep("Sill penalized estimator",length(sills)), 
                   rep("Sill unweighted", length(sills_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = sill, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for the sill parameter",
       x = "Value",
       y = "Density") +
  xlim(0, 2*sill) +
  theme(legend.position = "none")

# Range estimate

data <- data.frame(
  value = c(ranges_nor, ranges_new, ranges_bist, ranges_unwe),
  group = factor(c(rep("Range random path process - normalized",length(ranges_nor)), 
                   rep("Range random path process - new matrix", length(ranges_new)),
                   rep("Range random path process - bistochastic", length(ranges_bist)),
                   rep("Ranges unweighted estimator", length(ranges_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = range, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(col1, col2, third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for Three Different Groups",
       x = "Value",
       y = "Density") +
  xlim(0, 5e5) +
  theme(legend.position = "none")

# General look to the shapes of the empirical semi-variograms

xx <- seq(0, 4e5, length = 200)
yy <- spherical_kernel(c(sill, range), xx)
plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_nor, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_new, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_bist, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_unwe, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

# Mean squared error of the empirical semi-variograms

squared_diff <- lapply(variograms_nor, function(x) x$squared_diff)
dists <- lapply(variograms_nor, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_nor <- colMeans(do.call(rbind, dists))
MSE_nor <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                   dists_nor), ncol = dim(variograms)[2], 
                   nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_new, function(x) x$squared_diff)
dists <- lapply(variograms_new, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_new <- colMeans(do.call(rbind, dists))
MSE_new <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
            dists_new), ncol = dim(variograms)[2], 
            nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_bist, function(x) x$squared_diff)
dists <- lapply(variograms_bist, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_bist <- colMeans(do.call(rbind, dists))
MSE_bist <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                             dists_bist), ncol = dim(variograms)[2], 
                                         nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_unwe, function(x) x$squared_diff)
dists <- lapply(variograms_unwe, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_unwe <- colMeans(do.call(rbind, dists))
MSE_unwe <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
           dists_unwe), ncol = dim(variograms)[2], 
           nrow = dim(variograms)[1], byrow = T))^2)

data <- data.frame(
  dist = c(dists_nor, dists_new, dists_bist, dists_unwe),
  MSE = c(MSE_nor, MSE_new, MSE_bist, MSE_unwe),
  type = factor(c(rep("Random process penalized - normalized", length(dists_nor)),
                  rep("Random process penalized - new matrix", length(dists_new)),
                  rep("Random process penalized - bistochastic", length(dists_bist)),
                  rep("Unweighted estimator", length(dists_unwe))))
)
ggplot(data, aes(x = dist, y = MSE, color = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col1, third, col2, fourth)) +
  labs(x = "Distance", y = "MSE") +
  ylim(0, max(data$MSE)) +
  theme_minimal()

# Plot of the errors evaluated via leave one out cross validation

errors_data <- data.frame(
  error = c(rowMeans(errors_new_pen), rowMeans(errors_nor_pen), 
            rowMeans(errors_new_unwe), rowMeans(errors_nor_unwe),
            rowMeans(errors_eucl)),
  process = factor(rep(c("New matrix", "Normalized", 
                         "New matrix unweighted", "Normalized unweighted", "Euclidean"), 
                       each = K))
)

ggplot(errors_data, aes(x = process, y = error, fill = process)) +
  geom_boxplot(outlier.stroke = 1) +
  labs(title = "Boxplots of Errors from Different Processes",
       x = "Process",
       y = "Error") +
  scale_fill_manual(values = c(fourth, third, col1, third, col1)) +
  theme_minimal() +
  theme(legend.position = "none")

# Highlighting the distribution of the errors in the domain:
# A brighter colour stands for a better performance of the normalized process
# which by the way is the process generating the data.
df$errors <- NA

df$errors[inx] <- colMeans(errors_nor_pen-errors_new_pen)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

df$errors[inx] <- colMeans(errors_nor_unwe-errors_new_unwe)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

#### New matrix process ####

B <- length(NewMatrix)
K <- 500
sill <- 15
range <- 2e5
initial_params <- c(sill, range)
cov <- build_covariances(sill, range, NewMatrix, fun = spherical_covariance)
errors_nor_unwe <- matrix(0, nrow = K, ncol = B)
errors_nor_pen <- matrix(0, nrow = K, ncol = B)
errors_new_unwe <- matrix(0, nrow = K, ncol = B)
errors_new_pen <- matrix(0, nrow = K, ncol = B)
variograms_nor <- list('vector', B)
variograms_new <- list('vector', B)
variograms_bist <- list('vector', B)
variograms_unwe <- list('vector', B)
sills <- rep(0, K)
sills_unwe <- rep(0, K)
ranges_nor <- rep(0, K)
ranges_new <- rep(0, K)
ranges_bist <- rep(0, K)
ranges_unwe <- rep(0, K)
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)

for (iter in 1:K)
{
  ma <- create_process_covariance(iter, cov)
  
  variogram <- evaluate_variogram_penalization(ma, normalized, 15, 
                                               return_fuv = T, lambda = lam_nor)
  fuv <- variogram[[1]]
  variogram <- variogram[[2]]
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_sill <- fuv
  estimated_range <- params[2]
  sills[iter] <- estimated_sill
  ranges_nor[iter] <- estimated_range
  covariances_nor <- build_covariances(estimated_sill, estimated_range, normalized, 
                                       fun = spherical_covariance)
  variograms_nor[[iter]] <- variogram
  
  variogram <- evaluate_variogram_penalization(ma, NewMatrix, 15, lambda = lam_new)
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_range <- params[2]
  ranges_new[iter] <- estimated_range
  covariances_new <- build_covariances(estimated_sill, estimated_range, NewMatrix, 
                                       fun = spherical_covariance)
  variograms_new[[iter]] <- variogram
  
  predictions <- loo_predictions(ma, covariances_nor, covariances_new, eucl = F)
  
  aux <- predictions[[1]]
  errors_nor_pen[iter,] <- (aux$preds - ma)^2
  aux <- predictions[[2]]
  errors_new_pen[iter,] <- (aux$preds - ma)^2
  
  variogram <- evaluate_variogram_penalization(ma, bistochastic, 15, lambda = lam_bist)
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_range <- params[2]
  ranges_bist[iter] <- estimated_range
  variograms_bist[[iter]] <- variogram
  
  variogram <- evaluate_variogram_unadjusted(ma, normalized, 15, 6e5)
  params <- fit_variogram(variogram, spherical_kernel, initials = initial_params)
  sills_unwe[iter] <- params[1]
  ranges_unwe[iter] <- params[2]
  variograms_unwe[[iter]] <- variogram
  
  covariances_nor <- build_covariances(params[1], params[2], normalized, 
                                       fun = spherical_covariance)
  covariances_new <- build_covariances(params[1], params[2], NewMatrix, 
                                       fun = spherical_covariance)
  predictions <- loo_predictions(ma, covariances_nor, covariances_new, eucl = F)
  errors_nor_unwe[iter, ] <- (predictions[[1]]$preds-ma)^2
  errors_new_unwe[iter, ] <- (predictions[[2]]$preds-ma)^2
  
  pb$tick()
}

# Sill estimate
# Note that the sill estimate for the penalized estimator is equivalent whatever model
# is used

data <- data.frame(
  value = c(sills, sills_unwe),
  group = factor(c(rep("Sill penalized estimator",length(sills)), 
                   rep("Sill unweighted", length(sills_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = sill, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for the sill parameter",
       x = "Value",
       y = "Density") +
  xlim(0, 2*sill) +
  theme(legend.position = "none")

# Range estimate

data <- data.frame(
  value = c(ranges_nor, ranges_new, ranges_bist, ranges_unwe),
  group = factor(c(rep("Range random path process - normalized",length(ranges_nor)), 
                   rep("Range random path process - new matrix", length(ranges_new)),
                   rep("Range random path process - bistochastic", length(ranges_bist)),
                   rep("Ranges unweighted estimator", length(ranges_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = range, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(col1, col2, third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for Three Different Groups",
       x = "Value",
       y = "Density") +
  xlim(0, 5e5) +
  theme(legend.position = "none")

# General look to the shapes of the empirical semi-variograms

xx <- seq(0, 4e5, length = 200)
yy <- spherical_kernel(c(sill, range), xx)
plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_nor, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_new, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_bist, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_unwe, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

# Mean squared error of the empirical semi-variograms

squared_diff <- lapply(variograms_nor, function(x) x$squared_diff)
dists <- lapply(variograms_nor, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_nor <- colMeans(do.call(rbind, dists))
MSE_nor <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                                                          dists_nor), ncol = dim(variograms)[2], 
                                         nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_new, function(x) x$squared_diff)
dists <- lapply(variograms_new, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_new <- colMeans(do.call(rbind, dists))
MSE_new <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                                                          dists_new), ncol = dim(variograms)[2], 
                                         nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_bist, function(x) x$squared_diff)
dists <- lapply(variograms_bist, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_bist <- colMeans(do.call(rbind, dists))
MSE_bist <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                                                           dists_bist), ncol = dim(variograms)[2], 
                                          nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_unwe, function(x) x$squared_diff)
dists <- lapply(variograms_unwe, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_unwe <- colMeans(do.call(rbind, dists))
MSE_unwe <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                                                           dists_unwe), ncol = dim(variograms)[2], 
                                          nrow = dim(variograms)[1], byrow = T))^2)

data <- data.frame(
  dist = c(dists_nor, dists_new, dists_bist, dists_unwe),
  MSE = c(MSE_nor, MSE_new, MSE_bist, MSE_unwe),
  type = factor(c(rep("Random process penalized - normalized", length(dists_nor)),
                  rep("Random process penalized - new matrix", length(dists_new)),
                  rep("Random process penalized - bistochastic", length(dists_bist)),
                  rep("Unweighted estimator", length(dists_unwe))))
)
ggplot(data, aes(x = dist, y = MSE, color = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col1, third, col2, fourth)) +
  labs(x = "Distance", y = "MSE") +
  ylim(0, max(data$MSE)) +
  theme_minimal()

# Highlighting the distribution of the errors in the domain:
# A brighter colour stands for a better performance of the new matrix process
# which by the way is the process generating the data.
df$errors <- NA

df$errors[inx] <- colMeans(errors_new_pen - errors_nor_pen)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

df$errors[inx] <- colMeans(errors_new_unwe - errors_nor_unwe)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

#### Bistochastic process ####

B <- length(bistochastic)
K <- 500
sill <- 15
range <- 2e5
initial_params <- c(sill, range)
cov <- build_covariances(sill, range, bistochastic, fun = spherical_covariance)
errors_bist_unwe <- matrix(0, nrow = K, ncol = B)
errors_bist_pen <- matrix(0, nrow = K, ncol = B)
errors_new_unwe <- matrix(0, nrow = K, ncol = B)
errors_new_pen <- matrix(0, nrow = K, ncol = B)
variograms_nor <- list('vector', B)
variograms_new <- list('vector', B)
variograms_bist <- list('vector', B)
variograms_unwe <- list('vector', B)
sills <- rep(0, K)
sills_unwe <- rep(0, K)
ranges_nor <- rep(0, K)
ranges_new <- rep(0, K)
ranges_bist <- rep(0, K)
ranges_unwe <- rep(0, K)
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K,
  clear = FALSE,
  width = 60
)
pb$tick(0)

for (iter in 1:K)
{
  ma <- create_process_covariance(iter, cov)
  
  variogram <- evaluate_variogram_penalization_best_lambda(ma, bistochastic, 10, 
                                               return_fuv = T)
  fuv <- variogram[[1]]
  variogram <- variogram[[2]]
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_sill <- fuv
  estimated_range <- params[2]
  sills[iter] <- estimated_sill
  ranges_bist[iter] <- estimated_range
  covariances_bist <- build_covariances(estimated_sill, estimated_range, bistochastic, 
                                       fun = spherical_covariance)
  variograms_bist[[iter]] <- variogram
  
  variogram <- evaluate_variogram_penalization_best_lambda(ma, NewMatrix, 15)
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_range <- params[2]
  ranges_new[iter] <- estimated_range
  covariances_new <- build_covariances(estimated_sill, estimated_range, NewMatrix, 
                                       fun = spherical_covariance)
  variograms_new[[iter]] <- variogram
  
  predictions <- loo_predictions(ma, covariances_bist, covariances_new, eucl = F)
  
  aux <- predictions[[1]]
  errors_bist_pen[iter,] <- (aux$preds - ma)^2
  aux <- predictions[[2]]
  errors_new_pen[iter,] <- (aux$preds - ma)^2
  
  variogram <- evaluate_variogram_penalization_best_lambda(ma, normalized, 15)
  params <- fit_variogram(variogram, spherical_kernel, initial_params)
  estimated_range <- params[2]
  ranges_nor[iter] <- estimated_range
  variograms_nor[[iter]] <- variogram
  
  variogram <- evaluate_variogram_unadjusted(ma, bistochastic, 15, 4e5)
  params <- fit_variogram(variogram, spherical_kernel, initials = initial_params)
  sills_unwe[iter] <- params[1]
  ranges_unwe[iter] <- params[2]
  variograms_unwe[[iter]] <- variogram
  
  covariances_bist <- build_covariances(params[1], params[2], bistochastic, 
                                       fun = spherical_covariance)
  covariances_new <- build_covariances(params[1], params[2], NewMatrix, 
                                       fun = spherical_covariance)
  predictions <- loo_predictions(ma, covariances_bist, covariances_new, eucl = F)
  errors_bist_unwe[iter, ] <- (predictions[[1]]$preds-ma)^2
  errors_new_unwe[iter, ] <- (predictions[[2]]$preds-ma)^2
  
  pb$tick()
}

# Sill estimate
# Note that the sill estimate for the penalized estimator is equivalent whatever model
# is used

data <- data.frame(
  value = c(sills, sills_unwe),
  group = factor(c(rep("Sill penalized estimator",length(sills)), 
                   rep("Sill unweighted", length(sills_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = sill, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for the sill parameter",
       x = "Value",
       y = "Density") +
  xlim(0, 2*sill) +
  theme(legend.position = "none") +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  )

# Range estimate

data <- data.frame(
  value = c(ranges_nor, ranges_new, ranges_bist, ranges_unwe),
  group = factor(c(rep("Range random path process - normalized",length(ranges_nor)), 
                   rep("Range random path process - new matrix", length(ranges_new)),
                   rep("Range random path process - bistochastic", length(ranges_bist)),
                   rep("Ranges unweighted estimator", length(ranges_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = range, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(col1, col2, third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for Three Different Groups",
       x = "Value",
       y = "Density") +
  xlim(0, 5e5) +
  theme(legend.position = "none")

data <- data.frame(
  value = c(ranges_bist, ranges_unwe),
  group = factor(c(rep("Range penalized estimator", length(ranges_bist)),
                   rep("Ranges unweighted estimator", length(ranges_unwe))))
)
# Create the combined histogram
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = range, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(third, fourth)) +
  theme_minimal() +
  labs(title = "Density Estimates for the range parameter",
       x = "Value",
       y = "Density") +
  xlim(0, 5e5) +
  theme(legend.position = "none") +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20) 
  )

# General look to the shapes of the empirical semi-variograms

xx <- seq(0, 4e5, length = 200)
yy <- spherical_kernel(c(sill, range), xx)
plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_nor, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l')
lapply(variograms_new, function(x) lines(x$dist, x$squared_diff, col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l', xlab="Distance (meters)", ylab="Semi-variogram",
     cex.lab=2,cex.axis=2, xlim=c(0,
         max(unlist(lapply(variograms_bist,function(x) max(x$dist))))))
lapply(variograms_bist, function(x) lines(c(0,x$dist), c(0,x$squared_diff), col = grey))
lines(xx, yy, lwd = 5, col = third)

plot(xx, yy, ylim = c(0,2*sill), type = 'l', xlab="Distance (meters)", ylab="Semi-variogram",
     cex.lab=2,cex.axis=2,xlim=c(0,
     max(unlist(lapply(variograms_bist,function(x) max(x$dist))))))
lapply(variograms_unwe, function(x) lines(c(0,x$dist), c(0,x$squared_diff), col = grey))
lines(xx, yy, lwd = 5, col = third)

# Mean squared error of the empirical semi-variograms

squared_diff <- lapply(variograms_nor, function(x) x$squared_diff)
dists <- lapply(variograms_nor, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_nor <- colMeans(do.call(rbind, dists))
MSE_nor <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                                                          dists_nor), ncol = dim(variograms)[2], 
                                         nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_new, function(x) x$squared_diff)
dists <- lapply(variograms_new, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_new <- colMeans(do.call(rbind, dists))
MSE_new <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
                                                          dists_new), ncol = dim(variograms)[2], 
                                         nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_bist, function(x) x$squared_diff)
dists <- lapply(variograms_bist, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_bist <- colMeans(do.call(rbind, dists))
MSE_bist <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
             dists_bist), ncol = dim(variograms)[2], 
          nrow = dim(variograms)[1], byrow = T))^2)

squared_diff <- lapply(variograms_unwe, function(x) x$squared_diff)
dists <- lapply(variograms_unwe, function(x) x$dist)
variograms <- do.call(rbind, squared_diff)
dists_unwe <- colMeans(do.call(rbind, dists))
MSE_unwe <- colMeans((variograms - matrix(spherical_kernel(c(sill, range), 
               dists_unwe), ncol = dim(variograms)[2], 
                                          nrow = dim(variograms)[1], byrow = T))^2)

data <- data.frame(
  dist = c(dists_nor, dists_new, dists_bist, dists_unwe),
  MSE = c(MSE_nor, MSE_new, MSE_bist, MSE_unwe),
  type = factor(c(rep("Random process penalized - normalized", length(dists_nor)),
                  rep("Random process penalized - new matrix", length(dists_new)),
                  rep("Random process penalized - bistochastic", length(dists_bist)),
                  rep("Unweighted estimator", length(dists_unwe))))
)
ggplot(data, aes(x = dist, y = MSE, color = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col1, third, col2, fourth)) +
  labs(x = "Distance", y = "MSE") +
  ylim(0, max(data$MSE)) +
  theme_minimal()

data <- data.frame(
  dist = c(dists_bist, dists_unwe),
  MSE = c(MSE_bist, MSE_unwe),
  Estimator = factor(c(rep("Penalized estimator", length(dists_bist)),
                  rep("Unweighted estimator", length(dists_unwe))))
)
ggplot(data, aes(x = dist, y = MSE, color = Estimator)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col1, col2)) +
  labs(x = "Distance (meters)", y = "MSE") +
  ylim(0, max(data$MSE)) +
  xlim(0, 3e5) +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text = element_text(size = 15)
  )

# Plot of the errors evaluated via leave one out cross validation

errors_data <- data.frame(
  error = c(rowMeans(errors_new_pen), rowMeans(errors_bist_pen), 
            rowMeans(errors_new_unwe), rowMeans(errors_bist_unwe)),
  process = factor(rep(c("New matrix", "Bistochastic", 
                         "New matrix unweighted", "Bistochastic unweighted"), 
                       each = K))
)

ggplot(errors_data, aes(x = process, y = error, fill = process)) +
  geom_boxplot(outlier.stroke = 1) +
  labs(title = "Boxplots of Errors from Different Processes",
       x = "Process",
       y = "Error") +
  scale_fill_manual(values = c(third, col1, third, col1)) +
  theme_minimal() +
  theme(legend.position = "none")

errors_data <- data.frame(
  error = c(rowMeans(errors_bist_pen), rowMeans(errors_bist_unwe)),
  process = factor(rep(c("Bistochastic", "Bistochastic unweighted"), 
                       each = K))
)

ggplot(errors_data, aes(x = process, y = error, fill = process)) +
  geom_boxplot(outlier.stroke = 1) +
  labs(x = "Process",
       y = "Error") +
  scale_fill_manual(values = c(third, third)) +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

# Highlighting the distribution of the errors in the domain:
# A brighter colour stands for a better performance of the bistochastic process
# which by the way is the process generating the data.
df$Errors <- NA

df$Errors[inx] <- colMeans(errors_bist_pen- errors_bist_unwe)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

df$Errors[inx] <- colMeans(errors_bist_unwe-errors_new_unwe)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

df$Errors[inx] <- colMeans(errors_bist_pen)
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Errors)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1) +
  xlab("Longitude")+
  ylab("Latitude") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  )

#### Further Simulations ####

#### Testing the new estimator varying the range in the simulation ####

# Normalized process

B <- length(normalized)
grid <- seq(1e5, 3e5, length = 10)
lam_nor <- seq(0.001, 0.5, length = 10)
g <- length(grid)
l <- length(lam_nor)
ranges_nor <- matrix(0, nrow = l, ncol = g)
up_ranges_nor <- matrix(0, nrow = l, ncol = g)
low_ranges_nor <- matrix(0, nrow = l, ncol = g)
K <- 500
sill <- 1
library(progress)
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K*g,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (r in 1:g)
{
  range <- grid[r]
  initial_params <- c(sill, range)
  cov <- build_covariances(sill, range, normalized, fun = spherical_covariance)
  r_nor <- matrix(0, nrow = K, ncol = l)
  
  for (iter in 1:K)
  {
    ma <- create_process_covariance(iter, cov)
    
    variogram_list <- lapply(lam_nor, FUN = function(x) 
      evaluate_variogram_penalization(ma, normalized, 10, lambda = x))
    
    r_nor[iter,] <- unlist(lapply(variogram_list, function(x) fit_variogram(x, 
                                                                            spherical_kernel, initials = initial_params)[2]))
    
    pb$tick()
  }
  ranges_nor[,r] <- colMeans(r_nor)
  up_ranges_nor[,r] <- apply(r_nor, MARGIN = 2, function(x) quantile(x, probs  = c(.95)))
  low_ranges_nor[,r] <- apply(r_nor, MARGIN = 2, function(x) quantile(x, probs  = c(.05)))
}

df_values <- data.frame(x = grid, y = grid)
df_matrix1 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(ranges_nor), 
                         group = rep(lam_nor, times = g),
                         matrix = 'Point')
df_matrix2 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(up_ranges_nor), 
                         group = rep(lam_nor, times = g),
                         matrix = 'Up')
df_matrix3 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(low_ranges_nor), 
                         group = rep(lam_nor, times = g),
                         matrix = 'Low')
df_matrix <- rbind(df_matrix1, df_matrix2, df_matrix3)
df_matrix$line_type <- df_matrix$matrix

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), linetype = line_type), 
            linewidth = 1) +
  scale_color_manual(values = rep(colorRampPalette(c(col1, col2))(l), 3)) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted")) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - Normalized process",
       x = "Index",
       y = "Value",
       color = "Lambdas") +
  theme_minimal()

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix1, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group)), 
            linewidth = 1) +
  scale_color_manual(values = rep(colorRampPalette(c(col1, col2))(l))) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - Normalized process",
       x = "Index",
       y = "Value",
       color = "Lambdas") +
  theme_minimal()

# The new matrix process

B <- length(NewMatrix)
grid <- seq(1e5, 3e5, length = 10)
lam_new <- seq(0.01, 1, length = 10)
g <- length(grid)
l <- length(lam_new)
ranges_new <- matrix(0, nrow = l, ncol = g)
up_ranges_new <- matrix(0, nrow = l, ncol = g)
low_ranges_new <- matrix(0, nrow = l, ncol = g)
var_new <- matrix(0, nrow = l, ncol = g)
K <- 500
sill <- 1
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K*g,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (r in 1:g)
{
  range <- grid[r]
  initial_params <- c(sill, range)
  cov <- build_covariances(sill, range, NewMatrix, fun = linear_covariance)
  r_new <- matrix(0, nrow = K, ncol = l)
  
  for (iter in 1:K)
  {
    ma <- create_process_covariance(iter, cov)
    
    variogram_list <- lapply(lam_new, FUN = function(x) 
      evaluate_variogram_penalization(ma, NewMatrix, 10, lambda = x))
    
    r_new[iter,] <- unlist(lapply(variogram_list, function(x) fit_variogram(x, 
                                                                            spherical_kernel, initials = initial_params)[2]))
    
    pb$tick()
  }
  ranges_new[,r] <- colMeans(r_new)
  up_ranges_new[,r] <- apply(r_new, MARGIN = 2, function(x) quantile(x, probs  = c(.95)))
  low_ranges_new[,r] <- apply(r_new, MARGIN = 2, function(x) quantile(x, probs  = c(.05)))
  var_new[,r] <- apply(r_new, MARGIN = 2, function(x) var(x))
}

df_values <- data.frame(x = 1:g, y = grid)
df_matrix1 <- data.frame(x = rep(1:g, each = l), 
                         y = as.vector(ranges_new), 
                         group = rep(lam_new, times = g),
                         matrix = 'Point')
df_matrix2 <- data.frame(x = rep(1:g, each = l), 
                         y = as.vector(up_ranges_new), 
                         group = rep(lam_new, times = g),
                         matrix = 'Up')
df_matrix3 <- data.frame(x = rep(1:g, each = l), 
                         y = as.vector(low_ranges_new), 
                         group = rep(lam_new, times = g),
                         matrix = 'Low')
df_matrix <- rbind(df_matrix1, df_matrix2, df_matrix3)
df_matrix$line_type <- df_matrix$matrix

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), linetype = line_type), 
            linewidth = 1) +
  scale_color_manual(values = rep(colorRampPalette(c(col1, col2))(l), 3)) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted")) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - New matrix process",
       x = "Index",
       y = "Value",
       color = "Lambdas",
       linetype = "Processes") +
  theme_minimal()

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix1, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group)), 
            linewidth = 1) +
  scale_color_manual(values = colorRampPalette(c(col1, col2))(l)) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - New matrix process",
       x = "Index",
       y = "Value",
       color = "Lambdas",
       linetype = "Processes") +
  theme_minimal()

df_matrix1 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(ranges_new), 
                         group = rep(lam_new, times = g),
                         vars = sqrt(as.vector(var_new)),
                         matrix = 'Point')
df_matrix1$vars[which(df_matrix1$y>max(grid))] <- 
  min(df_matrix1$vars[which(df_matrix1$y<max(grid))])
df_values <- data.frame(x = grid, y = grid)

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix1, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), size = vars)) +
  scale_color_manual(values = colorRampPalette(c(col1, col2))(l)) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - Bistochastic process",
       x = "Index",
       y = "Value",
       color = "Lambdas",
       size = "Standard deviation") +
  xlab("Range") +
  ylab("Estimate") +
  theme_minimal()

# The bistochastic process

B <- length(bistochastic)
grid <- seq(1e5, 3e5, length = 10)
lam_bist <- c(2.5,2.7,3,4,5)*10^(-2)
g <- length(grid)
l <- length(lam_bist)
ranges_bist <- matrix(0, nrow = l, ncol = g)
up_ranges_bist <- matrix(0, nrow = l, ncol = g)
low_ranges_bist <- matrix(0, nrow = l, ncol = g)
var_bist <- matrix(0, nrow = l, ncol = g)
K <- 500
sill <- 1
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K*g,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (r in 1:g)
{
  range <- grid[r]
  initial_params <- c(sill, range)
  cov <- build_covariances(sill, range, bistochastic, fun = spherical_covariance)
  r_bist <- matrix(0, nrow = K, ncol = l)
  
  eigV <- eigen(cov)
  S <- eigV$vectors %*% diag(sqrt(eigV$values)) %*% t(eigV$vectors)
  
  for (iter in 1:K)
  {
    set.seed(iter)
    ma <- S%*%rnorm(B,0,1)
    
    variogram_list <- lapply(lam_bist, FUN = function(x) 
      evaluate_variogram_penalization(ma, bistochastic, 15, lambda = x))
    
    r_bist[iter,] <- unlist(lapply(variogram_list, function(x) fit_variogram(x, 
                                                                             spherical_kernel, initials = initial_params)[2]))
    
    pb$tick()
  }
  ranges_bist[,r] <- colMeans(r_bist)
  up_ranges_bist[,r] <- apply(r_bist, MARGIN = 2, function(x) quantile(x, probs  = c(.95)))
  low_ranges_bist[,r] <- apply(r_bist, MARGIN = 2, function(x) quantile(x, probs  = c(.05)))
  var_bist[,r] <- apply(r_bist, MARGIN = 2, function(x) var(x))
}

df_values <- data.frame(x = grid, y = grid)
df_matrix1 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(ranges_bist), 
                         group = rep(lam_bist, times = g),
                         matrix = 'Point')
df_matrix2 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(up_ranges_bist), 
                         group = rep(lam_bist, times = g),
                         matrix = 'Up')
df_matrix3 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(low_ranges_bist), 
                         group = rep(lam_bist, times = g),
                         matrix = 'Low')
df_matrix <- rbind(df_matrix1, df_matrix2, df_matrix3)
df_matrix$line_type <- df_matrix$matrix

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), linetype = line_type), 
            linewidth = 1) +
  scale_color_manual(values = rep(colorRampPalette(c(col1, col2))(l), 3)) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted")) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - Bistochastic process",
       x = "Index",
       y = "Value",
       color = "Lambdas",
       linetype = "Processes") +
  theme_minimal()

df_matrix1 <- data.frame(x = rep(grid, each = l), 
                         y = as.vector(ranges_bist), 
                         group = rep(lam_bist, times = g),
                         vars = sqrt(as.vector(var_bist)),
                         matrix = 'Point')
df_matrix1$vars[which(df_matrix1$y>max(grid))] <- 
  min(df_matrix1$vars[which(df_matrix1$y<max(grid))])

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix1, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), size = vars)) +
  scale_color_manual(values = colorRampPalette(c(col1, col2))(l)) +
  ylim(c(0, 3e5)) + 
  labs(title = "Plot of Ranges estimates - Bistochastic process",
       x = "Index",
       y = "Value",
       color = "Lambdas",
       size = "Standard deviation") +
  xlab("Range") +
  ylab("Estimate") +
  theme_minimal()

# x11()
# ggplot() +
#   geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
#   geom_line(data = df_matrix1, 
#             aes(x = x, y = y, group = interaction(matrix, group), 
#                 color = factor(group), size = vars)) +
#   scale_color_manual(values = colorRampPalette(c(col1, col2))(l)) +
#   ylim(c(0, 3e5)) + 
#   labs(title = "Plot of Ranges estimates - Bistochastic process",
#        x = "Index",
#        y = "Value",
#        color = "Lambdas",
#        size = "Variance") +
#   xlab("Range") +
#   ylab("Estimate") +
#   theme_minimal()

#### Comparison of the estimators varying the range ####

# Normalized process

B <- length(normalized)
lam_nor <- lam_nor[which.min(apply(ranges_nor, MARGIN = 1, function(x) sum((x - grid)^2)))]
# lam_nor <- 0.143
grid <- seq(1e5, 3e5, length = 10)
g <- length(grid)
ranges_nor <- rep(0, g)
up_ranges_nor <- rep(0, g)
low_ranges_nor <- rep(0, g)
ranges_unwe <- rep(0, g)
up_ranges_unwe <- rep(0, g)
low_ranges_unwe <- rep(0, g)
ranges_fcwa <- rep(0, g)
up_ranges_fcwa <- rep(0, g)
low_ranges_fcwa <- rep(0, g)
K <- 500
sill <- 1
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K*g,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (r in 1:g)
{
  range <- grid[r]
  initial_params <- c(sill, range)
  cov <- build_covariances(sill, range, normalized, fun = spherical_covariance)
  r_nor <- rep(0, K)
  r_unwe <- rep(0, K)
  r_fcwa <- rep(0, K)
  
  for (iter in 1:K)
  {
    ma <- create_process_covariance(iter, cov)
    
    variogram <- evaluate_variogram_penalization(ma, normalized, 10, lambda = lam_nor)
    r_nor[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    variogram <- evaluate_variogram_fcwa(ma, normalized, 10)
    r_fcwa[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    variogram <- evaluate_variogram_unadjusted(ma, normalized, 10, 5e5)
    r_unwe[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    
    pb$tick()
  }
  ranges_nor[r] <- mean(r_nor)
  up_ranges_nor[r] <- quantile(r_nor, probs  = c(.95))
  low_ranges_nor[r] <- quantile(r_nor, probs  = c(.05))
  ranges_unwe[r] <- mean(r_unwe)
  up_ranges_unwe[r] <- quantile(r_unwe, probs  = c(.95))
  low_ranges_unwe[r] <- quantile(r_unwe, probs  = c(.05))
  ranges_fcwa[r] <- mean(r_fcwa)
  up_ranges_fcwa[r] <- quantile(r_fcwa, probs  = c(.95))
  low_ranges_fcwa[r] <- quantile(r_fcwa, probs  = c(.05))
}

df_values <- data.frame(x = 1:g, y = grid)
df_matrix1 <- data.frame(x = rep(1:g, time = 3), 
                         y = c(ranges_nor, ranges_unwe, ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Point')
df_matrix2 <- data.frame(x = rep(1:g, time = 3), 
                         y = c(up_ranges_nor, up_ranges_unwe, up_ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Up')
df_matrix3 <- data.frame(x = rep(1:g, time = 3), 
                         y = c(low_ranges_nor, low_ranges_unwe, low_ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Low')
df_matrix <- rbind(df_matrix1, df_matrix2, df_matrix3)
df_matrix$line_type <- df_matrix$matrix

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", size = 1) +
  geom_line(data = df_matrix, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), linetype = line_type), 
            linewidth = 1) +
  scale_color_manual(values = rep(c(third, fourth, col1, col2), 3)) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted")) +
  ylim(c(0, 3e5)) + 
  labs(title = "Ranges estimates - Random path process - Normalized matrix",
       x = "Index",
       y = "Value",
       color = "Lambdas") +
  theme_minimal()

# New matrix process

B <- length(NewMatrix)
lam_new <- lam_new[which.min(apply(ranges_new, MARGIN = 1, function(x) sum((x - grid)^2)))]
# lam_new <- 0.15
grid <- seq(1e5, 3e5, length = 10)
ranges_new <- rep(0, g)
up_ranges_new <- rep(0, g)
low_ranges_new <- rep(0, g)
ranges_unwe <- rep(0, g)
up_ranges_unwe <- rep(0, g)
low_ranges_unwe <- rep(0, g)
ranges_fcwa <- rep(0, g)
up_ranges_fcwa <- rep(0, g)
low_ranges_fcwa <- rep(0, g)
K <- 500
sill <- 1
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K*g,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (r in 1:g)
{
  range <- grid[r]
  initial_params <- c(sill, range)
  cov <- build_covariances(sill, range, NewMatrix, fun = spherical_covariance)
  r_new <- rep(0, K)
  r_unwe <- rep(0, K)
  r_fcwa <- rep(0, K)
  
  for (iter in 1:K)
  {
    ma <- create_process_covariance(iter, cov)
    
    variogram <- evaluate_variogram_penalization(ma, NewMatrix, 10, lambda = lam_new)
    r_new[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    variogram <- evaluate_variogram_fcwa(ma, NewMatrix, 10)
    r_fcwa[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    variogram <- evaluate_variogram_unadjusted(ma, NewMatrix, 10, 5e5)
    r_unwe[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    
    pb$tick()
  }
  ranges_new[r] <- mean(r_new)
  up_ranges_new[r] <- quantile(r_new, probs  = c(.95))
  low_ranges_new[r] <- quantile(r_new, probs  = c(.05))
  ranges_unwe[r] <- mean(r_unwe)
  up_ranges_unwe[r] <- quantile(r_unwe, probs  = c(.95))
  low_ranges_unwe[r] <- quantile(r_unwe, probs  = c(.05))
  ranges_fcwa[r] <- mean(r_fcwa)
  up_ranges_fcwa[r] <- quantile(r_fcwa, probs  = c(.95))
  low_ranges_fcwa[r] <- quantile(r_fcwa, probs  = c(.05))
}

df_values <- data.frame(x = 1:g, y = grid)
df_matrix1 <- data.frame(x = rep(1:g, time = 3), 
                         y = c(ranges_new, ranges_unwe, ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Point')
df_matrix2 <- data.frame(x = rep(1:g, time = 3), 
                         y = c(up_ranges_new, up_ranges_unwe, up_ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Up')
df_matrix3 <- data.frame(x = rep(1:g, time = 3), 
                         y = c(low_ranges_new, low_ranges_unwe, low_ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Low')
df_matrix <- rbind(df_matrix1, df_matrix2, df_matrix3)
df_matrix$line_type <- df_matrix$matrix

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", linewidth = 1) +
  geom_line(data = df_matrix, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), linetype = line_type), 
            linewidth = 1) +
  scale_color_manual(values = rep(c(third, fourth, col1), 3)) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted")) +
  ylim(c(0, 3e5)) + 
  labs(title = "Ranges estimates - Random path process - New matrix",
       x = "Index",
       y = "Value",
       color = "Lambdas") +
  theme_minimal()

# Bistochastic process

B <- length(bistochastic)
grid <- seq(5e4, 2e5, length = 10)
g <- length(grid)
ranges_bist <- rep(0, g)
up_ranges_bist <- rep(0, g)
low_ranges_bist <- rep(0, g)
var_bist <- rep(0, g)
ranges_unwe <- rep(0, g)
up_ranges_unwe <- rep(0, g)
low_ranges_unwe <- rep(0, g)
var_unwe <- rep(0, g)
ranges_fcwa <- rep(0, g)
up_ranges_fcwa <- rep(0, g)
low_ranges_fcwa <- rep(0, g)
var_fcwa <- rep(0, g)
K <- 500
sill <- 1
pb <- progress_bar$new(
  format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
  total = K*g,
  clear = FALSE,
  width = 60
)
pb$tick(0)
for (r in 1:g)
{
  range <- grid[r]
  initial_params <- c(sill, range)
  cov <- build_covariances(sill, range, bistochastic, fun = spherical_covariance)
  r_bist <- rep(0, K)
  r_unwe <- rep(0, K)
  r_fcwa <- rep(0, K)
  
  for (iter in 1:K)
  {
    ma <- create_process_covariance(iter, cov)
    
    variogram <- evaluate_variogram_penalization_best_lambda(ma, bistochastic, 15)
    r_bist[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    variogram <- evaluate_variogram_fcwa(ma, bistochastic, 15)
    r_fcwa[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    variogram <- evaluate_variogram_unadjusted(ma, bistochastic, 15, 5e5)
    r_unwe[iter] <- fit_variogram(variogram, spherical_kernel, initials = initial_params)[2]
    
    pb$tick()
  }
  ranges_bist[r] <- median(r_bist)
  up_ranges_bist[r] <- quantile(r_bist, probs  = c(.75))
  low_ranges_bist[r] <- quantile(r_bist, probs  = c(.25))
  var_bist[r] <- var(r_bist)
  ranges_unwe[r] <- median(r_unwe)
  up_ranges_unwe[r] <- quantile(r_unwe, probs  = c(.75))
  low_ranges_unwe[r] <- quantile(r_unwe, probs  = c(.25))
  var_unwe[r] <- var(r_unwe)
  ranges_fcwa[r] <- median(r_fcwa)
  up_ranges_fcwa[r] <- quantile(r_fcwa, probs  = c(.75))
  low_ranges_fcwa[r] <- quantile(r_fcwa, probs  = c(.25))
  var_fcwa[r] <- var(r_fcwa)
}

df_values <- data.frame(x = grid, y = grid)
df_matrix1 <- data.frame(x = rep(grid, time = 3), 
                         y = c(ranges_bist, ranges_unwe, ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Point')
df_matrix2 <- data.frame(x = rep(grid, time = 3), 
                         y = c(up_ranges_bist, up_ranges_unwe, up_ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Up')
df_matrix3 <- data.frame(x = rep(grid, time = 3), 
                         y = c(low_ranges_bist, low_ranges_unwe, low_ranges_fcwa), 
                         group = rep(c("Penalized", "Unweighted", "FCWA"), each = g),
                         matrix = 'Low')
df_matrix <- rbind(df_matrix1, df_matrix2, df_matrix3)
df_matrix$line_type <- df_matrix$matrix

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", linewidth = 1) +
  geom_line(data = df_matrix, 
            aes(x = x, y = y, group = interaction(matrix, group), 
                color = factor(group), linetype = line_type), 
            linewidth = 1) +
  scale_color_manual(values = rep(c(third, fourth, col1), 3)) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted")) +
  ylim(c(0, 3e5)) + 
  labs(title = "Ranges estimates - Random path process - Bistochastic matrix",
       x = "Index",
       y = "Value",
       color = "Lambdas") +
  theme_minimal()

df_matrix <- data.frame(x = rep(grid, time = 2), 
                        y = c(ranges_bist, ranges_unwe), 
                        group = rep(c("Penalized", "Unweighted"), each = g),
                        IQR = c(up_ranges_bist - low_ranges_bist, 
                                up_ranges_unwe - low_ranges_unwe),
                        matrix = 'Point')

ggplot() +
  geom_line(data = df_values, aes(x = x, y = y), color = "black", linewidth = 1) +
  geom_point(data = df_matrix, 
             aes(x = x, y = y, group = interaction(matrix, group), 
                 color = factor(group), size = IQR)) +
  scale_color_manual(values = rep(c(third, fourth, col1), 3)) +
  ylim(c(0, 3e5)) + 
  scale_size_continuous(range = c(1, 5)) + 
  labs(title = "Ranges estimates - Bistochastic matrix",
       x = "Index",
       y = "Value",
       color = "Estimator",
       size = "Interquantile range") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )
