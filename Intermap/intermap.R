library(MetaboMate)
library(hastaLaVista)
library(car)

insget.idx=function(range=c(1,5), ppm){
  range=sort(range, decreasing=T);
  which(ppm<=range[1] & ppm>=range[2])}

raw1 <- read.csv2("Fvisit_RawData", sep=",", dec=".", header=TRUE)
country <- factor(raw1[,2])
raw1 <- raw1[,-c(1,2,3),]       #
raw <- raw1
ppm <- read.csv2("ppm.csv", sep=",", dec=".", header=FALSE)
remove(raw1)
# raw <- read.csv2("/home/jul/git/intermap/Xpqm_normF", sep=",", dec=".", header=TRUE)
# raw <- raw[,-1,]
# ppm <- read.csv2("/home/jul/git/intermap/ppmNorm.csv", sep=",", dec=".", header=FALSE)

IDs <- read.csv2("intermapID.csv", sep=",", dec=".", header=TRUE)
ID <- IDs[,1]

# # create json spectra


path = "C:/Users/RLL/Documents/R/win-library/3.6/hastaLaVista/visu/data/json"
for (i in seq_along(ID)) {
   filename = paste0(ID[i],".json")
   hastaLaVista::saveJSON(list("y" = as.numeric(t(raw[i,]))), path, filename)
 }


#ID <- ID
NN <- 100                     #4618 for whole data set 
group <- country[1:NN]
metadata <- data.frame("country"=country[1:NN])
metadata['JcampUrl'] <- paste0('http://127.0.0.1:5474/data/json/', ID[1:NN], '.json')

#x <- matrix(raw[1:50,], dim(raw[1:50,])[1], dim(raw[1:50,])[2])
x_axis <- as.numeric(ppm)
levels(group) <- c('Japan', 'China', 'UK', 'USA')
color = sapply(group, function(x) getColor2(as.character(x)))

d = list()
c <- data.frame(ID = seq_along(ID[1:NN]),
                group = group,
                color = color,
                "_highlight" = seq_along(group) - 1,
                #dataMatrix = I(raw[1:NN,]),#I(matrix( c(rbind(repRow(x_axis, nrow(x)), x)), nrow(x), ncol(x)*2)),
                metadata = I(metadata),
                check.names = FALSE
)
d <- appendData(data = d, variableName = "data", variable = c, type = "table")
remove(c)
d <- appendData(data = d, variableName = "xAxis", variable = x_axis, type = "table")

# v <- new("visualization")
# v@data <- "intermap.data.json"
# v@view <- "dataExplorer_1_1_0.view.json"
# push(v, type="data", d)
# visualize(v)

model <- list()
rangeList <- list(c(0.5, 1.1), 
                  c(1.1, 1.5), 
                  c(1.5, 2), 
                  c(2, 2.475), 
                  c(2.475, 3), 
                  c(3, 3.4), 
                  c(3.4, 4.04), 
                  c(4.04, 4.5595), 
                  c(5.1, 5.5),
                  c(6.5, 6.94),
                  c(6.94, 7.53),
                  c(7.53, 8.16),
                  c(8.16, 8.66),
                  c(8.66, 9),
                  c(9, 9.5))
it <- 1

for (ran in rangeList) {

idx <- get.idx(range(ran), as.numeric(ppm))
pcaName <- paste0("model",rangeList[it])

mod <- MetaboMate::pca(raw[1:NN,idx])
plotscores(mod, an=list(Class=country[1:NN]), title = 'PCA')

ellipse <- dataEllipse(mod@t[,1], mod@t[,2], levels=0.80)

ellipseChart <- data.frame("x" = ellipse[,1],
                           "y" = ellipse[,2],
                           "color" = rep('black', length(ellipse[,1])))

c <- list(name = pcaName,
    scores = mod@t,
    loadings = cov(mod@t, raw[1:NN, idx]),
    loadingsXaxis = x_axis[idx],
    loadingsColor = abs(cor(mod@t, raw[1:NN, idx])),
    ellipses = ellipseChart
)

model[[it]] <- c
remove(c)
it <- it + 1
} ## end of loop

## we compute the pca for the whole spectra

pcaName <- paste0("model"," full spectra")

mod <- MetaboMate::pca(raw[1:NN,], pc = 3)
plotscores(mod, pc = c(1, 3), an=list(Class=country[1:NN]), title = 'PCA')

ellipse <- dataEllipse(mod@t[,1], mod@t[,2], levels=0.80)

ellipseChart <- data.frame("x" = ellipse[,1],
                           "y" = ellipse[,2],
                           "color" = rep('black', length(ellipse[,1])))

c <- list(name = pcaName,
          scores = mod@t,
          loadings = cov(mod@t, raw[1:NN, ]),
          loadingsXaxis = x_axis,
          loadingsColor = abs(cor(mod@t, raw[1:NN, ])),
          ellipses = ellipseChart
)

model[[it]] <- c
remove(c)
d[['model']] <- model

v2 <- new("visualization")
v2@data <- "intermap.data.json"
v2@view <- "modelExplorer_1_0.view.json"
push(v2, type="data", d)                     #this writes data to the disk
visualize(v2)
# rhia link should pop up in GOOGLE CHROME http://127.0.0.1:5474/?viewURL=http://127.0.0.1:5474/view/modelExplorer_1_0.view.json&dataURL=http://127.0.0.1:5474/data/intermap.data.json

## =======================
# if you re use this again, just load the 
library(hastaLaVista)
v2 <- new("visualization")
v2@data <- "intermap.data.json"
v2@view <- "modelExplorer_1_0.view.json"
 
visualize(v2)
