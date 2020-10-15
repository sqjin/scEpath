#
# R file to fit a principal curve of data in low-dimensional space
if(!require("princurve")){
  install.packages("princurve")
}
# library(princurve)
require("princurve")

args <- commandArgs()
baseName <- args[6]

#infile1 <- paste(getwd(), baseName, "ydata.txt", sep="/") 
infile1 <- paste(baseName, "ydata.txt", sep="/") 
infile2 <- paste(baseName, "ydataOutCell.txt", sep="/") 
infile3 <- paste(baseName,  "pathLength.txt", sep="/") 

ydata <- read.table(file=infile1, sep="\t",header=FALSE)
ydataOutCell <- read.table(file=infile2, sep="\t",header=FALSE)
pathLength <- read.table(file=infile3, sep="\t",header=FALSE)

y <- data.matrix(ydata)
yOutCell <- data.matrix(ydataOutCell)
pathLength <-data.matrix(pathLength)

# fitpc <- principal.curve(y, smoother="smooth.spline", df=pathLength+1,plot = TRUE, maxit = 50) # This function is deprecated
fitpc <- principal_curve(y, smoother="smooth.spline", df=pathLength+1,plot = TRUE, maxit = 50)

# projectionpc <- get.lam(yOutCell, fitpc$s, fitpc$tag, stretch = 2)# This function is deprecated
projectionpc <- project_to_curve(yOutCell, fitpc$s, stretch = 2)


# pcvout1 <- paste(getwd(), baseName,  "PcurveProjectionValueMainCell.txt", sep="/")
pcvout1 <- paste(baseName,  "PcurveProjectionValueMainCell.txt", sep="/")
pcvout2 <- paste(baseName, "PcurveProjectionValueOutCell.txt", sep="/")
lambdaout1 <- paste(baseName, "PcurveLambdaMainCell.txt", sep="/")
lambdaout2 <- paste(baseName,  "PcurveLambdaOutCell.txt", sep="/")

write.table(file=pcvout1, fitpc$s, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(file=pcvout2, projectionpc$s, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(file=lambdaout1, fitpc$lambda, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(file=lambdaout2, projectionpc$lambda, sep="\t", row.names = FALSE, col.names = FALSE)

# lambda is projection index for each cell,quantify the distance from the beginning of the curve
# s is a matrix corresponding to x, giving their projected value (fitted values) onto the curve.




