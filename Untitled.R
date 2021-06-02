
tmp <- installed.packages()
installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpkgs, file="installed_old.rda")

library(vctrs)
install.packages('devtools') #assuming it is not already installed
library(devtools)

install_github('andreacirilloac/updateR')

library(updateR)

updateR(admin_password = 'Sedona16')
