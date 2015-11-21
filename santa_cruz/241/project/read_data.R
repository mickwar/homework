read_data = function(path = "~/downloads/Ta/", truncate = 0.023){
    names = list.files(path)
    out = NULL
    for (i in 1:length(names)){
        # Get contents
        temp.x = readLines(paste0(path, names[i]))

        # Get temperature
        temp.t = temp.x[grep("temperature", temp.x)]
        temp.start = gregexpr("<", temp.t)[[1]]
        temp.end = gregexpr(">", temp.t)[[1]]
        temperature = as.numeric(substr(temp.t, temp.end[1]+1, temp.start[2]-1))

        # Get strain_rate
        temp.t = temp.x[grep("strain_rate", temp.x)]
        temp.start = gregexpr("<", temp.t)[[1]]
        temp.end = gregexpr(">", temp.t)[[1]]
        strain_rate = as.numeric(substr(temp.t, temp.end[1]+1, temp.start[2]-1))

        # Get strain/stress (x/y)
        temp.t = temp.x[grep("xy", temp.x)]
        xy = t(sapply(temp.t, function(x){
            temp.start = gregexpr("<", x)[[1]]
            temp.end = gregexpr(">", x)[[1]]
            x = substr(x, temp.end[1]+1, temp.start[2]-1)
            as.numeric(strsplit(x, ",")[[1]])
            }))
        colnames(xy) = c("x", "y")
        rownames(xy) = NULL

        # Truncate
        xy = xy[which(xy[,1] >= truncate),]

        out[[i]] = list("strain_rate"=strain_rate, "temperature"=temperature,
            "xy"=xy)
        }
    return (out)
    }
