# Monthly values of the The Southern Oscillation Index (SOI) during 
# 1950-1995. This series consists of 540 observations on the SOI, 
# computed as the ``difference of the departure from the long-term 
# monthly mean sea level pressures'' at Tahiti in the South Pacific 
# and Darwin in Northern Australia. The index is one measure of the 
# so-called "El Nino-Southern Oscillation" -- an event of critical 
# importance and interest in climatological studies in recent decades. 
# The fact that most of the observations in the last part of the series 
# take negative values is related to a recent warming in the tropical 
# Pacific, and this pattern has been of key interest in recent studies. 
# A key question of interest is to determine just how unusual this event 
# is, and if it can reasonably be explained by standard "stationary" 
# time series models, or requires models that include drifts/trends 
# that may be related to global climatic change. 

dat = read.table("~/files/repos/data/soi.txt")[,1]

plot(dat, type='l')
