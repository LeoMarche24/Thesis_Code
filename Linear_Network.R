library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)

source("functions_network.R")
source("functions_simulation.R")

col1 <- "cyan"
col2 <- "darkblue"
third <- rgb(63/255, 128/255, 97/255)
fourth <- "darkred"
grey <- "grey80"

df <- read.csv('Data/LigurianSea_Simulation.csv')
B <- dim(df)[1]

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], aes(x = lon, y = lat, fill=temperature)) +
  geom_raster() + 
  geom_segment(data = df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], aes(x = lon, y = lat, xend = lon + east, yend = lat + nord),
               arrow = arrow(length = unit(0.2, "cm")), color = third, alpha=0.8) + 
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = col2, high = col1, name = "Temperature")

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

p <- rep(1, B+1)/B
P <- PI_out
P <- rbind(P, ifelse(colSums(is.na(P)) == nrow(P), 1, NA))
P <- cbind(P, ifelse(rowSums(is.na(P))==ncol(P), 1, NA))
P[which(is.na(P))] <- 0
P[(B+1),] <- P[(B+1),]/sum(P[(B+1),])
for(i in 1:1000)
  p <- p%*%P
df$p <- p[-(B+1)]
df$p <- df$p + p[B+1]/B
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=p)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

plot_lin_net(lines, df[, 1:2])

# Cretion of the distances object. For any pair of points it evaluates, for any possible
# path, the distance of that path and the corresponding probability, according to the
# matrix prob provided in the function.

# Ver Hoef like, i decide in advance which are both the functions of the input and 
# i assume that all the water flows down stream.

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
order <- 1
while (updated)
{
  updated <- FALSE
  print(order)
  distances <- lapply(distances, function(p) update(p,dist, P))
  order <- order + 1
}

# Create a bistochastic matrix and then normalize the column

P <- build_bistochastic(PI_out, 1000, 2405)
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

## Regularity of the distances

dists <- lengths_distribution(distances)
vars <- unlist(lapply(dists, function(x) sd(x)))
means <- unlist(lapply(dists, function(x) mean(x)))

xx <- which(is.na(vars))
plot(means[-xx], vars[-xx])

dists <- lengths_distribution(bistochastic)
vars <- unlist(lapply(dists, function(x) sd(x)))
means <- unlist(lapply(dists, function(x) mean(x)))

xx <- which(is.na(vars))
plot(means[-xx], vars[-xx])

## Try some simulation

sill <- 15
range <- 2e5

cov <- build_covariances(sill, range, distances, fun = spherical_covariance)
simulation <- create_process_covariance(18, cov)

df$Simulation <- simulation
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Simulation)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

res <- evaluate_variogram_penalization(simulation, 
     distances, l = 15, return_fuv = T, lambda = 0.4)

variogram <- res[[2]]
fuv <- res[[1]]

ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

variogram <- evaluate_variogram_unadjusted(simulation, distances, l=15, 5e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

# Try to build it with the bistochastic process

sill <- 15
range <- 2e5

cov <- build_covariances(sill, range, bistochastic, fun = spherical_covariance)
simulation <- create_process_covariance(24, cov)

df$Simulation <- simulation
ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=Simulation)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1) + 
  theme_minimal()

variogram <- evaluate_variogram_penalization(simulation, bistochastic, l = 15, lambda = 0.03)

ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

# Compare the values on the real dataset

#### Temperature ####

inx <- intersect(which(!is.na(df$temperature)), which(!is.na(df$east)))
inx <- intersect(inx, which(!is.na(df$nord)))
set.seed(1807)
test <- sample(inx, 0.1*length(inx))
train <- setdiff(inx, test)

ggplot(df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], aes(x = lon, y = lat, fill=temperature)) +
  geom_raster() + 
  geom_segment(data = df[which(!is.na(df$east) & !is.na(df$nord) & !is.na(df$temperature)) ,], aes(x = lon, y = lat, xend = lon + east, yend = lat + nord),
               arrow = arrow(length = unit(0.2, "cm")), color = third, alpha=0.8) +
  geom_point(data = df[test,], aes(x=lon, y = lat), color = "black", size = 4) +
  scale_fill_gradient(low = col2, high = col1, name = "Temperature")

############
#### OK ####
############

dati_krig <- df[train ,]
coordinates(dati_krig) <- c("lon", "lat")
proj4string(dati_krig) <- CRS("+proj=longlat +datum=WGS84")
dati_krig <- spTransform(dati_krig, CRS("+proj=utm +zone=33 +datum=WGS84"))
svgm <- variogram(temperature ~ 1, dati_krig)
plot(svgm, main = 'Sample Variogram',pch=19)
v.fit <- fit.variogram(svgm, vgm(.6, "Gau", 7e4, 0))
plot(svgm, v.fit, main = 'Fitted variogram')
new <- df[test ,]
coordinates(new) <- c("lon", "lat")
proj4string(new) <- CRS("+proj=longlat +datum=WGS84")
new <- spTransform(new, CRS("+proj=utm +zone=33 +datum=WGS84"))
g.tr <- gstat(formula = temperature ~ 1, data = dati_krig, model=v.fit)
predict_eucl <- predict(g.tr, new)@data
err <- mean((new$temperature-predict_eucl$var1.pred)^2)
errors <- (new$temperature-predict_eucl$var1.pred)^2

#### MY NEW PROCESS ####

variogram <- evaluate_variogram_unadjusted_train(df$temperature, train, distances, 15, 3e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

variogram <- evaluate_variogram_penalization_train(df$temperature,
                                                   distances, train, l = 15, lambda = 0.15)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

initial_params <- c(.15, 1e5)
params <- fit_variogram(variogram, spherical_kernel, initial_params)
sill <- params[1]
range <- params[2]
plot(variogram$dist, variogram$squared_diff, col = col2, pch= 16,
     xlab = "Distances", ylab = "Squared differences")
points(seq(0,5e5, length=200), 
       spherical_kernel(c(sill, range), 
                        seq(0,5e5, length=200)), type='l', col=col1, lwd=4)

cov_new <- build_covariances(sill, range, distances, fun = spherical_covariance)
preds_new <- evaluate_test(df$temperature, cov_new, train, test)
err_new <- mean((df$temperature[test] - preds_new$preds)^2)

variogram <- evaluate_variogram_unadjusted_train(df$temperature, train, bistochastic, 15, 3e5)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

variogram <- evaluate_variogram_penalization_train(df$temperature,
                             bistochastic, train, l = 15, lambda = 0.1)
ggplot() +
  geom_point(data = variogram, aes(x=dist, y=squared_diff, size = np), col = col2) + 
  labs(x = "Distance", y = "Squared differences") +
  theme_minimal() +
  labs(title = "Variogram")

initial_params <- c(.15, 1e5)
params <- fit_variogram(variogram, spherical_kernel, initial_params)
sill <- params[1]
range <- params[2]
plot(variogram$dist, variogram$squared_diff, col = col2, pch= 16,
     xlab = "Distances", ylab = "Squared differences")
points(seq(0,5e5, length=200), 
       spherical_kernel(c(sill, range), 
                        seq(0,5e5, length=200)), type='l', col=col1, lwd=4)

cov_bist <- build_covariances(sill, range, bistochastic, fun = spherical_covariance)
preds_bist <- evaluate_test(df$temperature, cov_bist, train, test)
err_bist <- mean((df$temperature[test] - preds_bist$preds)^2)
