library(maps)
library(googleway)
library(sp)
library(dplyr)

dat = read.csv("./annual_all_2015.csv")

fix.poly = function(z){
    z = cbind(z$x, z$y)

    v = t(cbind(c(1, which(is.na(z[,1]))+1),
        c(which(is.na(z[,1]))-1, nrow(z))))

    island = NULL
    out = z[v[1,1]:v[2,1],]
    prev.ind = 2
    tmp = as.matrix(v[,-1])
    while (NCOL(tmp) > 0){
        ind = which(sapply(c(tmp), function(x) all(z[v[prev.ind],] == z[x,])))
        if (length(ind) == 0){
            island[[length(island)+1]] = unique(out)
            ind = c(2, 1)
            out = z[tmp[ind[1], ind[2]]:tmp[3-ind[1], ind[2]],]
            prev.ind = which(v == tmp[3-ind[1], ind[2]])
            tmp = as.matrix(tmp[,-ind[2]])
        } else {
            ind = which(tmp == tmp[ind], arr.ind = TRUE)
            out = rbind(out, z[tmp[ind[1], ind[2]]:tmp[3-ind[1], ind[2]],])
            prev.ind = which(v == tmp[3-ind[1], ind[2]])
            tmp = as.matrix(tmp[,-ind[2]])
            }
        }
    if (length(island) > 0){
        island[[length(island)+1]] = unique(out)
        out = island
    } else {
        out = unique(out)
        }
        
    return (out)
    }

### Filter the data
states = c("Maine", "Vermont", "New Hampshire", "Connecticut",
    "Massachusetts", "New Jersey", "Delaware", "Rhode Island",
    "Maryland", "District Of Columbia")
events = c("No Events", "Concurred Events Excluded", "Events Excluded")



### Ozone
ozone.ind = (dat$Parameter.Name == "Ozone") &
    (dat$Sample.Duration == "8-HR RUN AVG BEGIN HOUR") &
    (dat$Pollutant.Standard == "Ozone 8-Hour 2008") &
    (dat$State.Name %in% states) &
    (dat$Event.Type %in% events)
# dat[ozone.ind,]$Arithmetic.Mean * 1000

lon = dat[ozone.ind,]$Longitude
lat = dat[ozone.ind,]$Latitude
lonlat = unique(cbind(lon, lat))
y = double(NROW(lonlat))

# Average over the values where duplicat lon/lat show up
tmp = cbind(lon, lat)
tmp3 = dat[ozone.ind,]$Arithmetic.Mean * 1000
for (j in 1:NROW(lonlat)){
    tmp2 = NULL
    for (i in 1:NROW(tmp)){
        if (all(lonlat[j,] == tmp[i,])){
            tmp2 = c(tmp2, tmp3[i])
            }
        }
    y[j] = mean(tmp2)
    }

ozone.lonlat = lonlat
ozone.y = y




### NO2
no2.ind = (dat$Parameter.Name == "Nitrogen dioxide (NO2)") &
    (dat$Sample.Duration == "1 HOUR") &
    (dat$Pollutant.Standard == "NO2 1-hour") &
    (dat$State.Name %in% states) &
    (dat$Event.Type %in% events)

lon = dat[no2.ind,]$Longitude
lat = dat[no2.ind,]$Latitude
lonlat = unique(cbind(lon, lat))
y = double(NROW(lonlat))

# Average over the values where duplicat lon/lat show up
tmp = cbind(lon, lat)
tmp3 = dat[no2.ind,]$Arithmetic.Mean * 1
for (j in 1:NROW(lonlat)){
    tmp2 = NULL
    for (i in 1:NROW(tmp)){
        if (all(lonlat[j,] == tmp[i,])){
            tmp2 = c(tmp2, tmp3[i])
            }
        }
    y[j] = mean(tmp2)
    }

no2.lonlat = lonlat
no2.y = y




### Particulate Matter
pm25.ind = (dat$Parameter.Name == "PM2.5 - Local Conditions") &
    (dat$Sample.Duration == "24 HOUR") &
    (dat$Pollutant.Standard == "PM25 24-hour 2012") &
    (dat$State.Name %in% states) &
    (dat$Event.Type %in% events)

lon = dat[pm25.ind,]$Longitude
lat = dat[pm25.ind,]$Latitude
lonlat = unique(cbind(lon, lat))
y = double(NROW(lonlat))

# Average over the values where duplicat lon/lat show up
tmp = cbind(lon, lat)
tmp3 = dat[pm25.ind,]$Arithmetic.Mean * 1
for (j in 1:NROW(lonlat)){
    tmp2 = NULL
    for (i in 1:NROW(tmp)){
        if (all(lonlat[j,] == tmp[i,])){
            tmp2 = c(tmp2, tmp3[i])
            }
        }
    y[j] = mean(tmp2)
    }

pm25.lonlat = lonlat
pm25.y = y





### Individual
map("state", xlim = c(-79, -67), ylim = c(38, 48))
points(ozone.lonlat+0.00, pch = 16, col = 'blue')
points(no2.lonlat+0.05, pch = 16, col = 'green')
points(pm25.lonlat-0.05, pch = 16, col = 'red')

### Union
map("state", xlim = c(-79, -67), ylim = c(38, 48))
points(unique(rbind(ozone.lonlat, no2.lonlat, pm25.lonlat)), col = 'blue')

### Intersection
tmp = rbind(ozone.lonlat, no2.lonlat)
tmp = unique(tmp[duplicated(tmp),])
tmp = rbind(tmp, pm25.lonlat)
tmp = unique(tmp[duplicated(tmp),])
points(tmp, col = 'red')


### Knots
loc = unique(rbind(ozone.lonlat, no2.lonlat, pm25.lonlat))
map("state", xlim = c(-79, -67), ylim = c(38, 48))
grid = as.matrix(expand.grid(seq(min(loc[,1]), max(loc[,1]), length = 14)+0.4,
    seq(min(loc[,2]), max(loc[,2]), length = 15)-0.1))


knots = NULL
for (i in 1:length(states)){
    tmp = map("state", states[i])
    tmp = fix.poly(tmp)
    if (is.list(tmp)){
        tmp2 = unlist(lapply(tmp, function(x) pointsInPoly(x, grid)))
        tmp2 = tmp2[!is.na(tmp2)]
    } else {
        tmp2 = pointsInPoly(tmp, grid)
        }
    if (!all(is.na(tmp2)))
        knots = rbind(knots, grid[tmp2,])
    }

### TODO: add another boundary around the knots
knots = expand.grid(seq(0, 1, length = 14), seq(0, 1, length = 19))
plot(knots)
rm.ind = sample(nrow(knots), floor(nrow(knots) * 0.65))
knots = knots[-rm.ind,]
plot(knots, pch = 16)

ukx = sort(unique(knots[,1]))
dx = min(diff(ukx))

uky = sort(unique(knots[,2]))
dy = min(diff(uky))

Ax = t(matrix(c(-1, 0, 1) * dx, 3, 3))
Ay = matrix(c(1, 0, -1) * dy, 3, 3)

border = NULL
for (i in 1:nrow(knots)){
    for (j in 1:9)
        border = rbind(border, c(knots[i,1] + Ax[j], knots[i,2] + Ay[j]))
    }
border = unique(border)
dd = as.matrix(dist(border, diag = TRUE, upper = TRUE))
ind = which(dd < 1e-6 & dd > 0, arr.ind = TRUE)
for (i in 1:nrow(ind))
    border[ind[i, 1],] = border[ind[i, 2],]
border = unique(border)

plot(border, pch = 15, col = 'darkblue')
points(knots, pch = 16, cex = 2)


map("state", xlim = c(-79, -67), ylim = c(38, 48))
points(grid, col = 'red')
points(knots, col = 'green')
nrow(knots)


### All locations with knots
map("state", xlim = c(-79, -67), ylim = c(38, 48))
points(unique(rbind(ozone.lonlat, no2.lonlat, pm25.lonlat)), col = 'blue', pch = 15)
points(knots, col = 'green', pch = 16)


### Finer grid for kriging
loc = unique(rbind(ozone.lonlat, no2.lonlat, pm25.lonlat))
map("state", xlim = c(-79, -67), ylim = c(38, 48))
grid = as.matrix(expand.grid(seq(min(loc[,1])-2, max(loc[,1])+2, length = 120),
    seq(min(loc[,2])-3, max(loc[,2]+1), length = 120)))

krig = NULL
for (i in 1:length(states)){
    tmp = map("state", states[i])
    tmp = fix.poly(tmp)
    if (is.list(tmp)){
        tmp2 = unlist(lapply(tmp, function(x) pointsInPoly(x, grid)))
        tmp2 = tmp2[!is.na(tmp2)]
    } else {
        tmp2 = pointsInPoly(tmp, grid)
        }
    if (!all(is.na(tmp2)))
        krig = rbind(krig, grid[tmp2,])
    }
map("state", xlim = c(-81, -65), ylim = c(36, 50))
points(grid, col = 'red')
points(krig, col = 'green')
nrow(krig)
axis(1); axis(2)



df.lonlat = data.frame(unique(rbind(ozone.lonlat, no2.lonlat, pm25.lonlat)))
df.lonlat = df.lonlat[,2:1]
df.knots = data.frame("lat"=knots[,2], "lon"=knots[,1])
df.krig = data.frame("lat"=krig[,2], "lon"=krig[,1])


### Get elevation
api_key=      # Need to get a Google api key


ele.loc = google_elevation(df.lonlat, key = api_key)$results[,1]
ele.knots = google_elevation(df.knots, key = api_key)$results[,1]

ele.krig = double(nrow(krig))
for (i in 1:ceiling(nrow(krig) / 100)){
    end = min(i * 100, nrow(krig))
    ele.krig[(100*(i-1)+1):end] =
        google_elevation(df.krig[(100*(i-1)+1):end,], key = api_key)$results[,1]
    }

library(rgl)
plot3d(df.lonlat[,2], df.lonlat[,1], ele.loc)
points3d(df.knots[,2], df.knots[,1], ele.knots, col = 'red')
points3d(df.krig[,2], df.krig[,1], ele.krig, col = 'blue')


loc.obs = cbind("lon"=df.lonlat[,2], "lat"=df.lonlat[,1], "alt"=ele.loc / 1000)
loc.knots = cbind("lon"=df.knots[,2], "lat"=df.knots[,1], "alt"=ele.knots / 1000)
loc.krig = cbind("lon"=df.krig[,2], "lat"=df.krig[,1], "alt"=ele.krig / 1000)


### Northing and easting
utm.dat = tbl_df(loc.obs[,1:2])
coordinates(utm.dat) = c("lon", "lat")
proj4string(utm.dat) = CRS("+init=epsg:4326")
utm.dat = as.data.frame(spTransform(utm.dat, CRS(paste("+proj=utm +zone=", 19," ellps=WGS84",sep=''))))
utm.dat = utm.dat / 1000
colnames(utm.dat) = c("kmEast", "kmNorth")
loc.obs = cbind(loc.obs, utm.dat)

utm.dat = tbl_df(loc.knots[,1:2])
coordinates(utm.dat) = c("lon", "lat")
proj4string(utm.dat) = CRS("+init=epsg:4326")
utm.dat = as.data.frame(spTransform(utm.dat, CRS(paste("+proj=utm +zone=", 19," ellps=WGS84",sep=''))))
utm.dat = utm.dat / 1000
colnames(utm.dat) = c("kmEast", "kmNorth")
loc.knots = cbind(loc.knots, utm.dat)

utm.dat = tbl_df(loc.krig[,1:2])
coordinates(utm.dat) = c("lon", "lat")
proj4string(utm.dat) = CRS("+init=epsg:4326")
utm.dat = as.data.frame(spTransform(utm.dat, CRS(paste("+proj=utm +zone=", 19," ellps=WGS84",sep=''))))
utm.dat = utm.dat / 1000
colnames(utm.dat) = c("kmEast", "kmNorth")
loc.krig = cbind(loc.krig, utm.dat)

### Package it all up
tmp = matrix(NA, nrow(loc.obs), 3)
for (i in 1:nrow(loc.obs)){
    for (j in 1:nrow(ozone.lonlat))
        if (all(loc.obs[i,1:2] == ozone.lonlat[j,])) tmp[i, 1] = ozone.y[j]
    for (j in 1:nrow(no2.lonlat))
        if (all(loc.obs[i,1:2] == no2.lonlat[j,])) tmp[i, 2] = no2.y[j]
    for (j in 1:nrow(pm25.lonlat))
        if (all(loc.obs[i,1:2] == pm25.lonlat[j,])) tmp[i, 3] = pm25.y[j]
    }
loc.obs = cbind(loc.obs, tmp)
colnames(loc.obs)[6:8] = c("ozone", "no2", "pm25")



### Write data set for faster load times
write.table(x = loc.obs, file = "./obs.txt", quote = FALSE, row.names = FALSE)
write.table(x = loc.knots, file = "./knots.txt", quote = FALSE, row.names = FALSE)
write.table(x = loc.krig, file = "./krig.txt", quote = FALSE, row.names = FALSE)
