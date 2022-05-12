##############################################################################-#
## The following codes give an illustration of detecting clusters of           #
## exposure-response relationship.                                             #
## These codes can completely replicate our results in the article             #
## Author: Wei Wang                                                            #
##############################################################################-#

###############################################################################
########### load library,function and data ####################################
###############################################################################
library(mvmeta)
library(splines)
library(dlnm)
library(rgdal)
library(tmap)
path <- "your path"
setwd(path)
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%` # pipe operator
source("MMSatScan-function.R")
load("data\\stage1data in motivating example.Rdata")

###############################################################################
########### detecting cluster #################################################
###############################################################################
id <- CINF$citycd
geo <- CINF[,c("POINT_X","POINT_Y")]
n <- nrow(yall)
# maxsize = 0.5 and mcmc = 999 in the article
system.time(
  res <- MMSatScan(id = id,geo = geo,value = yall,slist = Sall,maxsize = 0.1,mcmc = 99)
)
# save("DetectedClustersInMotivatingExample.Rdata")
# load("data\\DetectedClustersInMotivatingExample.Rdata") # The dectected cluster for maxsize = 0.5 and mcmc = 999
res$clusters
data_cluster <- merge(CINF,res$clusters,by.x = "citycd",by.y = "id",all.x = T)
data_cluster$i[is.na(data_cluster$pvalue) | data_cluster$pvalue > 0.05] <- 0
data_cluster <- within(data_cluster,{
  clustertype <- "NA"
  clustertype[i == 1] <- 2 # i = 1 refers to cluster 2 in the article 
  clustertype[i == 2] <- 1 # i = 2 refers to cluster 1 in the article
  clustertype[i == 0] <- 0 # i = 2 refers to noncluster in the article
})
###############################################################################
########## plot map for detected clusters #####################################
###############################################################################
load("data\\mapdata_China_province_project.Rdata")
citypoints <- CINF
citypoints <- merge(citypoints,data_cluster[,c("citycd","clustertype")],by = "citycd")
sp::coordinates(citypoints) <- c("POINT_X","POINT_Y")
citypoints@proj4string <- mapdata@proj4string
mapdata <- sf::st_as_sf(mapdata) # transform sp to sf class
citypoints <- sf::st_as_sf(citypoints)

pdf(path%+%"The distribution of cities.pdf")
tmap_mode("plot")
tm_shape(mapdata) + 
  tm_polygons(col = "grey100") + 
  tm_shape(citypoints) + 
  tm_bubbles(size = 0.2,col = "grey50") +
  tm_add_legend(type = "symbol",col = "grey50",
                labels = "The center of a city",size = 0.8) 
dev.off()

citypoints$col <- "grey"
citypoints$col[citypoints$clustertype == 1] <- "red" 
citypoints$col[citypoints$clustertype == 2] <- "blue" 
n0 <- sum(citypoints$clustertype == 0)
n1 <- sum(citypoints$clustertype == 1)
n2 <- sum(citypoints$clustertype == 2)
pdf("The distribution of clusters.pdf")
tmap_mode("plot")
tm_shape(mapdata) + 
  tm_polygons(col = "grey100") + 
  tm_shape(citypoints) + 
  tm_bubbles(size = 0.15,col = "col") +
  tm_add_legend(type = "symbol",col = c("grey","red","blue"),size = 0.8,
                labels = c("Not cluster ("%+%n0%+%")","Cluster 1 ("%+%n1%+%")","Cluster 2 ("%+%n2%+%")")) 
dev.off()  
###############################################################################
########## plot the exposure-response curve for clusters and noncluster area ##
###############################################################################
xvar <- seq(0,100,by=5)
bvar <- do.call("onebasis",c(list(x=xvar),attr(M.CB,"argvar")))
name <- "The average ERRs for clusters and non-cluster region "

#### for all cities
fit <- mvmeta(yall~1,S = Sall,method="fixed")
cpall <- crosspred(bvar, coef=coef(fit), vcov=vcov(fit), model.link="log",by=5)
x<-seq(0,100,by=5)
tiff(filename=name%+%".tiff", width=15.9, height=8, units="cm",pointsize = 8,res=300)
par(cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(4,4,3.6,0.8))
plot(cpall,ylim=c(0.6,1.8),type="n", ci="n",xlab="Temperature on relative scale (%)",
     ylab="Relative Risk (RR)")
redtrans <- rgb(255, 100, 0, 60, maxColorValue=255) 
points(xvar,cpall$allRRfit, type="l", lwd=1, col="orange")
polygon(c(rev(x),x),c(rev(cpall$allRRhigh),cpall$allRRlow),col=redtrans, border = NA)


#### for cluster
fixed_cluster <- data_cluster$clustertype[match(CINF$citycd,data_cluster$citycd)] %>% as.factor() 
fit <- mvmeta(yall~fixed_cluster,S = Sall,method="fixed")

# for cityies in noncluster area
beta <- fit$coefficients[1,]
covindex <- seq(1,15,3)
covbeta <- vcov(fit)[covindex,covindex]
cpall <- crosspred(bvar, coef=beta, vcov=covbeta, model.link="log",by=5)
lines(cpall,col="grey",lty=4,lwd=1,ci="area",ci.arg=list(density=20,col="grey"))

# for cluster 1
beta <- fit$coefficients[1,] + fit$coefficients[2,]
covindex <- setdiff(1:15,seq(1,15,3)+1) 
a1 <- matrix(0,ncol = 10,nrow = 5)
a1[1,1:2] <- a1[2,3:4] <- a1[3,5:6] <- a1[4,7:8] <- a1[5,9:10] <- 1
covbeta <- vcov(fit)[covindex,covindex]
covbeta <- a1 %*% covbeta %*% t(a1)
cpall <- crosspred(bvar, coef=beta, vcov=covbeta, model.link="log",by=5)
lines(cpall,col="red",lty=4,lwd=1,ci="area",ci.arg=list(density=20,col="red",angle = -45))

# for cityies in cluster 2
beta <- fit$coefficients[1,] + fit$coefficients[3,]
covindex <- setdiff(1:15,seq(1,15,3)+2) 
a1 <- matrix(0,ncol = 10,nrow = 5)
a1[1,1:2] <- a1[2,3:4] <- a1[3,5:6] <- a1[4,7:8] <- a1[5,9:10] <- 1
covbeta <- vcov(fit)[covindex,covindex]
covbeta <- a1 %*% covbeta %*% t(a1)
cpall <- crosspred(bvar, coef=beta, vcov=covbeta, model.link="log",by=5)
redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
points(xvar,cpall$allRRfit, type="l", lwd=1, col="blue")
polygon(c(rev(x),x),c(rev(cpall$allRRhigh),cpall$allRRlow),col=redtrans, border = NA)
legend(x="top",inset =0, legend=c("All cities", "Non-cluster","Cluster 1","Cluster 2"),
       lwd=1.5, lty=1, col=c("orange", "grey","red","blue"), bty="n",ncol=2, cex=1.6)
dev.off()


################################################################################
############ AIC comparison#####################################################
################################################################################
####### meta analysis without incorporating clusters
fit0 <- mvmeta(yall~1,S = Sall,method="ml")
AIC(fit0)
####### meta analysis with incorporating clusters as fixed effects
fixed_cluster <- data_cluster$clustertype[match(CINF$citycd,data_cluster$citycd)] %>% as.factor() 
fit1 <- mvmeta(yall~fixed_cluster,S = Sall,method="ml")
AIC(fit1)

################################################################################
# The following codes implement the simulation study, which will cost much     #
# time due to the large number of simulation data sets. A super computer       #
# may be necessary to replicate our results                                    #
################################################################################
path <- "your path"
setwd(path)
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%` # pipe operator
source("MMSatScan-function.R")
load(path%+%"data\\Simdata.Rdata")
load(path%+%"data\\SimulationScenario.Rdata")

Sall <- Simdata$covbeta
CINF <- Simdata$coordinates
scennames <- names(Simdata)[1:12]
maxsize <- 0.5
core <- parallel::detectCores(logical = T)
for (s in scennames) {
  scen <- s
  nsim <- 1000
  scenres <- vector(mode = "list",length = nsim)
  for (i in 1:nsim) {
    simi <- i
    yall <- t(sapply(Simdata[[scen]],function(x) x[simi,]))
    id <- CINF$citycd
    geo <- CINF[,c("POINT_X","POINT_Y")]
    res <- MMSatScan(id = id,geo = geo,value = yall,slist = Sall,
                     maxsize = maxsize,mcmc = 999,ncore = core)
    scenres[[i]] <- res$clusters
    if(i%%50==0) 
      save(scenres,file = scen%+%"_mws"%+%maxsize%+%".Rdata")
    cat(i,"\n")
  }
}

GetPerformance <- function(object,scen,allscenario,alpha = 0.05){
  x <- allscenario[,scen]
  tc <- allscenario$citycd[which(x != "grey")] # true cluster
  tnc <- allscenario$citycd[which(x == "grey")] # true noncluster
  dc <- object[,"id"][object[,"pvalue"]<=alpha] # detected cluster
  dnc <- setdiff(allscenario$citycd,dc) # detected noncluster
  power <- as.integer(length(dc)>0)
  sensitivity <- ifelse(length(tc)==0,NA,length(intersect(tc,dc))/length(tc)) 
  specificity <- length(intersect(tnc,dnc))/length(tnc)
  ppv <- ifelse(length(dc) == 0,NA,length(intersect(tc,dc))/length(dc)) 
  ydi <- sensitivity + specificity - 1
  misclassificaiton <- 1 - (length(intersect(tc,dc)) + length(intersect(tnc,dnc)))/length(x)
  aa <- c(power,sensitivity,specificity,ppv,ydi,misclassificaiton)
  names(aa) <- c("power","sensitivity","specificity","ppv","ydi","misclassificaiton")
  aa
}
# GetPerformance(object = res$clusters,scen = scen,allscenario = allscenario,alpha=0.05)
allperformance <- NULL
xx <- scennames[1:12]
for (i in xx) {
  scen <- i
  load(scen%+%"_mws0.5.Rdata")
  performance <- t(sapply(scenres, GetPerformance, scen = scen, allscenario = allscenario))
  allperformance <- rbind(allperformance,round(apply(performance, 2, mean,na.rm = T),4)) 
}
rownames(allperformance) <- xx
xlsx::write.xlsx(allperformance,file = "simulation-result.xlsx")




