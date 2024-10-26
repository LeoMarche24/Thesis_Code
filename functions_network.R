extract_data <- function(name, lon1, lon2, lat1, lat2)
{
  dataset <- nc_open(name)
  latitude <- ncvar_get(dataset, "latitude")
  longitude <- ncvar_get(dataset, "longitude")
  temperature <- ncvar_get(dataset, "to")
  sal <- ncvar_get(dataset, "so")
  time <- ncvar_get(dataset, "time")
  uo_clean <- ncvar_get(dataset, "ugo")
  vo_clean <- ncvar_get(dataset, "vgo")
  nc_close(dataset)
  uo_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                      lat = rep(latitude, times = length(longitude)),
                      value = c(t(uo_clean)))
  vo_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                      lat = rep(latitude, times = length(longitude)),
                      value = c(t(vo_clean)))
  tem_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                       lat = rep(latitude, times = length(longitude)),
                       value = c(t(temperature)))
  sal_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                       lat = rep(latitude, times = length(longitude)),
                       value = c(t(sal)))
  
  lon_sub <- longitude[which(longitude>lon1 & longitude<lon2)]
  lat_sub <- latitude[which(latitude>lat1 & latitude<lat2)]
  
  coords <- expand.grid(lon_sub, lat_sub)
  B <- dim(coords)[1]
  
  lon_idx <- match(coords[, 1], longitude)
  lat_idx <- match(coords[, 2], latitude)
  rows <- match(paste(longitude[lon_idx], latitude[lat_idx]), paste(uo_df$lon, uo_df$lat))
  temp <- tem_df$value[rows]
  sal_temp <- sal_df$value[rows]
  east_values <- uo_df$value[rows]
  nord_values <- vo_df$value[rows]
  df <- data.frame(
    lon = coords[, 1],
    lat = coords[, 2],
    east = east_values,
    nord = nord_values,
    temperature = temp,
    salinity = sal_temp
  )
  return(df)
}

##################################
##### CREATION OF THE NETWORK ####
##################################

normalize_column <- function(col)
{
  non_na_col <- col[!is.na(col)]
  normalized_col <- non_na_col / sum(non_na_col)
  col[!is.na(col)] <- normalized_col
  return(col)
}

generate_line_id <- function(coords)
{
  paste(coords[, 1], coords[, 2])
}

compute_matrices <- function(df) #This function is ok if the grid is regular - COPERNICUS
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
    if (df$east[i] > abs(df$nord[i]) || is.na(df$nord[i]))
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
    if (negative_nord[i])
    {
      dist_mat[i, (i-len_lon+1)] <- dist_eucl[i, (i-len_lon+1)]
      PI[i, (i-len_lon+1)] <- (sqrt(2)/2)*(df$east[i]+abs(df$nord[i]))
      line2 <- sp::Line(rbind(coord_sub[i,], coord_sub[(i-len_lon+1) ,]))
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
    if (df$nord[i] > abs(df$east[i]) || is.na(df$east[i]))
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
  trues <- which(negative_east)
  for (i in trues)
  {
    temp <- list()
    line <- NULL
    line2 <- NULL
    if (abs(df$east[i]) > abs(df$nord[i]) || is.na(df$nord[i]))
    {
      dist_mat[i, i-1] <- dist_eucl[i, i-1]
      PI[i, i-1] <- abs(df$east[i])
      line <- sp::Line(rbind(coord_sub[i,], coord_sub[i-1 ,]))
    }
    if (positive_nord[i])
    {
      dist_mat[i, (i+len_lon-1)] <- dist_eucl[i, (i+len_lon-1)]
      PI[i, (i+len_lon-1)] <- (sqrt(2)/2)*(abs(df$east[i])+df$nord[i])
      line2 <- sp::Line(rbind(coord_sub[i,], coord_sub[(i+len_lon-1) ,]))
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
    if (abs(df$nord[i]) > abs(df$east[i]) || is.na(df$east[i]))
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
  
  for (i in 1:B)
  {
    for(j in 1:B)
    {
      if (!is.na(PI[i,j]) & !is.na(PI[j,i]))
      {
        if(PI[i,j] > PI[j,i])
        {
          PI[j,i] <- NA
          dist_mat[j,i] <- NA
        }
        else
        {
          PI[i,j] <- NA
          dist_mat[i,j] <- NA
        }
      }
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

build_bistochastic <- function(P,K,seed, control_bar = F) 
  # In input the element with no connections should be Na
{
  diag(P) <- NA
  B <- dim(P)[1]
  start <- which((colSums(is.na(P)) == nrow(P)) & (rowSums(is.na(P)) != nrow(P)))
  end <- which((rowSums(is.na(P)) == nrow(P)) & (colSums(is.na(P)) != nrow(P)))
  
  P_fin <- matrix(0, B,B)
  
  if(control_bar)
  {
    pb <- progress_bar$new(
      format = "[:bar] :percent Elapsed: :elapsedfull Time to finish: :eta",
      total = K,
      clear = FALSE,
      width = 60
    )
    pb$tick(0)
  }
  
  for (k in 1:K)
  {
    P_transition <- matrix(0, B, B)
    set.seed(seed*k)
    s <- sample(start)
    visited <- c()
    while(length(s) > 0)
    {
      P_aux <- matrix(0, B, B)
      i <- s[1]
      s <- s[-1]
      first <- i
      pre <- i
      visited <- c(visited, pre)
      while(!(pre %in% end) & (pre != -1))
      {
        where <- which(!is.na(P[i,]))
        if(length(where)>1)
        {
          i <- sample(where, 1, prob = P[pre,where])
        }
        else
        {
          i <- where
        }
        if (!(i %in% visited))
        {
          P_aux[pre,i] <- 1
          pre <- i
          visited <- c(visited, pre)
        }
        else
        {
          pre <- -1
        }
      }
      
      if (pre != -1)
      {
        P_aux[pre, first] <- 1
        P_transition <- P_transition + P_aux
      }
    }
    
    diag(P_transition)[which(rowSums(P_transition) == 0 & colSums(P_transition) == 0)] <- 1
    
    if (sum(colSums(P_transition)!=rep(1,B)) || sum(rowSums(P_transition)!=rep(1,B)))
      print(c(k, "th iteration is not fine."))
    
    P_fin <- P_fin + P_transition
    if(control_bar)
      pb$tick()
  }
  
  P_fin <- P_fin/K
  
  P_fin[P_fin == 0 | is.na(P_fin) | is.nan(P_fin)] <- NA

  return(P_fin)
}

#############################
#### PLOT OF THE NETWORK ####
#############################

plot_lin_net <- function(lines, coord_sub)
{
  lin_net <- SpatialLines(LinesList = lines)
  lines_coords <- NULL
  k <- 1
  for (i in 1:length(lin_net@lines))
  {
    for (j in 1:length(lin_net@lines[[i]]@Lines))
    {
      new_cord <- as.vector(lin_net@lines[[i]]@Lines[[j]]@coords)
      if (length(new_cord) == 4)
      {
        lines_coords <- rbind(lines_coords, new_cord)
      }
      else
      {
        print(i,j)
        errorCondition("Problem")
      }
    }
  }
  lines_coords <- data.frame(lines_coords)
  names(lines_coords) <- c("x0", "x1", "y0", "y1")
  midpoints_x <- (lines_coords[, "x0"] + lines_coords[, "x1"]) / 2
  midpoints_y <- (lines_coords[, "y0"] + lines_coords[, "y1"]) / 2
  arrow_data <- data.frame(
    x0 = lines_coords$x0,
    y0 = lines_coords$y0,
    x1 = midpoints_x,
    y1 = midpoints_y
  )
  
  # Plot using ggplot2
  ggplot() +
    geom_sf(data = st_as_sf(lin_net), color = grey, linewidth = 2) +
    geom_point(data=data.frame(x=coord_sub[, 1], y=coord_sub[, 2]), aes(x=x, y=y), col='black') + 
    geom_segment(
      data = arrow_data,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
      color = col2, linewidth = .1
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    labs(title = "Linear network")
}

plot_one_point <- function(coord_sub, data)
{
  x1 <- (data$lon[1] + data$lon[2])/2
  y1 <- (data$lat[1] + data$lat[2])/2
  x2 <- (data$lon[1] + data$lon[3])/2
  y2 <- (data$lat[1] + data$lat[3])/2
  arrow_data <- data.frame(
    x0 = rep(data$lon[1], 2),
    y0 = rep(data$lat[1], 2),
    x1 = c(x1, x2),
    y1 = c(y1,y2)
  )
  x1 <- data$lon[2]
  y1 <- data$lat[2]
  x2 <- data$lon[3]
  y2 <- data$lat[3]
  segment_data <- data.frame(
    x0 = rep(data$lon[1], 2),
    y0 = rep(data$lat[1], 2),
    x1 = c(x1, x2),
    y1 = c(y1,y2)
  )
  
  # Plot using ggplot2
  ggplot() +
    geom_point(data=data.frame(x=coord_sub[, 1], y=coord_sub[, 2], 
                               label = c("A", "B", "C", "D", "E", "F", "G", "H", "I")), aes(x=x, y=y), col='black') + 
    geom_text(data=data.frame(x=coord_sub[, 1], y=coord_sub[, 2], 
                               label = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                              hjust = c(1,-1,1,1,-1,1,-1,1,1),
                              vjust = c(1,1,1,1,1,1,-1,-1,-1)), 
              aes(x=x, y=y, label = label, hjust=hjust, vjust=vjust), size = 10) + 
    geom_segment(
      data = segment_data,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      color = grey, linewidth = 2) + 
    geom_segment(
      data = arrow_data,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
      color = col2, linewidth = .5) +
    geom_segment(
      data = data[1 ,],
      aes(x = lon, y = lat, xend = lon + east*5, yend = lat + nord*5),
      arrow = arrow(length = unit(.3, "inches"), type = "closed"),
      color = col1, linewidth = 1) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    labs(title = "Linear network")
}

##########################################
#### CREATION OF THE DISTANCES OBJECT ####
##########################################

initialize <- function(dist, prob)
{
  B <- dim(dist)[1]
  distances <- vector("list", length = B)
  
  for (i in 1:B)
  {
    temp <- vector("list", length = B)
    distances[[i]] <- list(lengths = temp, update = NULL, Id = i)
  }
  
  for (i in 1:B)
  {
    for (j in 1:B)
    {
      distances[[i]]$visited[[j]] <- j
      if ((!is.na(dist[j, i])) & (!is.na(prob[j,i])))
      {
        distances[[i]]$lengths[[j]] <- 
          rbind(distances[[i]]$lengths[[j]], c(dist[j, i], prob[j, i]))
        distances[[i]]$update <- rbind(distances[[i]]$update, 
                                       cbind(j, nrow(distances[[i]]$lengths[[j]])))
      }
    }
  }
  return(distances)
}

update <- function(p, dist, prob) 
{
  list_update <- NULL
  if (length(nrow(p$update)))
  {
    for (row in 1:nrow(p$update))
    {
      i <- p$update[row, 1]
      j <- p$update[row, 2]
      inx <- which(!is.na(dist[, i]))
      if (length(inx))
      {
        for (k in 1:length(inx))
        {
          new_dist <- dist[inx[k],i]+p$lengths[[i]][j, 1]
          new_PI <- p$lengths[[i]][j, 2]*prob[inx[k],i]
          if (new_PI > 1e-7)
          {
            if (is.null(p$lengths[[inx[k]]]))
            {
              p$min[inx[k]] <- new_dist
            }
            p$lengths[[inx[k]]] <- rbind(p$lengths[[inx[k]]], 
                                         c(new_dist, new_PI))
            list_update <- rbind(list_update, cbind(inx[k], nrow(p$lengths[[inx[k]]])))
          }
        }
        updated <<- TRUE
      }
    }
    p$update <- list_update
  }
  return(p)
}
