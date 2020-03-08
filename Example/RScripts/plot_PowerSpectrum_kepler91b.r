# make the asteroseismology power spectra plots for the paper

library(LSspec)
library(lomb)

# load the data
RGB = read.table(file = "../Data/koi2133.lp.ts")
colnames(RGB) = c("Time", "Flux")

# function to caculate normalization used in lomb package
norm = function(x) var( x - mean(x) )

# calculate the LS periodogram using LSspec
lsRGB = LSspec(x = (RGB$Flux - mean(RGB$Flux)), t = RGB$Time )

# calculate the mtLS periodogram using LSspecMT for three different sets of tapers
K = c(7, 19, 39)
NW = c(4, 10, 20)

RG7 <- LSspecMT(t = RGB$Time, x = RGB$Flux, subtract.mean = TRUE, w=NW[1]/nrow(RGB), k=K[1])

RG19 <- LSspecMT(t = RGB$Time, x = RGB$Flux, subtract.mean = TRUE, w=NW[2]/nrow(RGB), k=K[2])

RG39 <- LSspecMT(t = RGB$Time, x = RGB$Flux, subtract.mean = TRUE, w=NW[3]/nrow(RGB), k=K[3])



# make multipanel plot
yrange = c(1e-13, 1e-6)
xrange= c(0,20)

# adjust line width
linew = 0.2

png("../transit91b_mtLS.png", height=9, width=4, units="in", res=300)
par(oma = c(4, 3, 0, 0))
par(mfrow=c(4,1), mar=c(2,4,1,2))

plot(lsRGB$freq, lsRGB$P, log="y", type="l", xlab="", ylab = "", ylim=yrange, xlim=xrange, lwd=linew)
grid()
text(x = 17.5, y = 5e-7, labels = "LS periodogram")


plot(RG7$freq, RG7$P, log="y", type="l", xlab="", ylab="", ylim=yrange, xlim=xrange, lwd=linew)
grid()
text(x = 17.5, y = 2.5e-7, labels = "MTLS estimate \n NW=4, K=7")


plot(RG19$freq, RG19$P, log="y", type="l", xlab="", ylab="", ylim=yrange, xlim=xrange, lwd=linew*1.5)
grid()
text(x = 17.5, y = 2.5e-7, labels = "MTLS estimate \n NW=10, K=19")

plot(RG39$freq, RG39$P, log="y", type="l", xlab="", ylab="", ylim=yrange, xlim=xrange, lwd=linew*2)
grid()
text(x = 17.5, y = 2.5e-7, labels = "MTLS estimate \n NW=20, K=39")


mtext(text = "Estimated Power Spectral Density", side = 2, outer = TRUE)
mtext(text = expression("Frequency "~(day^-1)), side = 1, outer = TRUE, line=2)

dev.off()

# make zoom-in plot to compare LS and mtLS around solar oscillations

xrange = c(6,13)
yrange = c(2e-10, 2e-7)

png(file = "../transit91b_mtLS_zoom-in.png", width = 6, height = 8, units="in", res=300)  
par(oma = c(4, 3, 0, 0))
par(mfrow=c(2,1), mar=c(2,4,1,2), cex.lab=2)

plot(lsRGB$freq, lsRGB$P, log="y", type="l", xlab="", ylab = "", ylim=yrange, xlim=xrange, lwd=linew)
grid()
legend("topright", legend = "LS periodogram",pch=NA,bty="n")

plot(RG19$freq, RG19$P, log="y", type="l", xlab="", ylab="", ylim=yrange, xlim=xrange, lwd=linew*1.5)
grid()
legend("topright", legend = "MTLS estimate \n NW=10, K=19",pch=NA,bty="n")

mtext(text = "Estimated Power Spectral Density", side = 2, outer = TRUE , cex=1.5)
mtext(text = expression("Frequency "~(day^-1)), side = 1, outer = TRUE, line=2 , cex=1.5)

dev.off()

