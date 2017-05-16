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

lon = lonlat[,1]
lat = lonlat[,2]

### Get elevation
dat2 = read.csv("./aqs_sites.csv")
ele = double(length(y)) -Inf

l2 = cbind(dat2$Longitude, dat2$Latitude)
for (i in 1:length(y)){
    tmp_dist = Inf
    for (j in 1:nrow(l2)){
        if (any(is.na(l2[j,])))
            next
        comp = sum(abs(lonlat[i,] - l2[j,]))
        if (comp < tmp_dist){
            ele[i] = dat2$Elevation[j]
            tmp_dist = comp
            }
        }
    }

### Write data set for faster load times
write.table(x = data.frame("Ozone"=y, "Longitude"=lon, "Latitude"=lat,
    "Altitude"=ele), file = "./northeast_data.txt", quote = FALSE,
    row.names = FALSE)
