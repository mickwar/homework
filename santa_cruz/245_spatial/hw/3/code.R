library(maps)
dat = read.csv("./annual_all_2015.csv")


### Filter the data
states = c("Maine", "Vermont", "New Hampshire", "Connecticut",
    "Massachusetts", "New Jersey", "Delaware", "Rhode Island",
    "Maryland", "District Of Columbia")
events = c("No Events", "Concurred Events Excluded", "Events Excluded")

param.ind = (dat$Parameter.Name == "Ozone") &
    (dat$Sample.Duration == "8-HR RUN AVG BEGIN HOUR") &
    (dat$Pollutant.Standard == "Ozone 8-Hour 2008")
state.ind = dat$State.Name %in% states
event.ind = dat$Event.Type %in% events

ind = (param.ind & state.ind & event.ind)

lon = dat[ind,]$Longitude
lat = dat[ind,]$Latitude

### Locations and values
lonlat = unique(cbind(lon, lat))
y = double(NROW(lonlat))

### Average over the values where duplicat lon/lat show up
tmp = cbind(lon, lat)
tmp3 = dat[ind,]$Arithmetic.Mean
for (j in 1:NROW(lonlat)){
    tmp2 = NULL
    for (i in 1:NROW(tmp)){
        if (all(lonlat[j,] == tmp[i,])){
            tmp2 = c(tmp2, tmp3[i])
            }
        }
    y[j] = mean(tmp2)
    }


### Plot data points on map
col.vals = (y - min(y)) / diff(range(y))
map("state", xlim = range(lonlat[,1]))
points(lonlat, col = rgb(col.vals, 0, 1-col.vals), pch = 16)

# Possibly a linear trend in the southwest (or northeast) direction

plot(lonlat[,1], y)
plot(lonlat[,2], y)
plot(lm(y ~ lonlat[,1]))
plot(lm(y ~ lonlat[,2]))
