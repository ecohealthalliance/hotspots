# Project HQ
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())      # Clear all variables
graphics.off()       # Close graphics windows

# Load libraries
library(geosphere) # spherical trigonometry functions for geographic applications
library(foreign) # exporting DBF files

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TODO(EHA) 1 degree analysis

dr<-read.table(file="data/eid08_drivers_19OCT11.csv",sep=",",header=TRUE)

# repeat regression analysis and scaling according to Jones et al 2008 to verify methods
# gets about 100% correspondance to the original values and coef are roughly in the right range

# this was the weight used in the original analysis - a better one would be %land area in each grid cell
weight=dr$landarea/mean(dr$landarea)
abslat=abs(dr$lat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Total GLM model (all drivers and all EIDs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Total model
mt<-glm(anytotr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

# Remove the publication parameter for the predictive model & calculate predictions
# Note: cannot use the predict function because of the missing parameter
pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mt)[[1]]+coef(mt)[[2]]*dr$lndensity+coef(mt)[[4]]*dr$high_pop_g+coef(mt)[[5]]*dr$mamdiv+coef(mt)[[6]]*dr$precip+coef(mt)[[7]]*abslat))) 

# Scale predictions between 1 and 0
P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

# Check correlation between predicted output and original SPSS predictions
cor.test(P_totr1_sc,dr$P_totr1_sc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Zoonotic non-wildlife disease (reponse variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mz<-glm(anyzoor1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mz)[[1]]+coef(mz)[[2]]*dr$lndensity+coef(mz)[[4]]*dr$high_pop_g+coef(mz)[[5]]*dr$mamdiv+coef(mz)[[6]]*dr$precip+coef(mz)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

cor.test(P_totr1_sc,dr$P_zoor1_sc)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wildlife zoonotic disease (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mw<-glm(anywildr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mw)[[1]]+coef(mw)[[2]]*dr$lndensity+coef(mw)[[4]]*dr$high_pop_g+coef(mw)[[5]]*dr$mamdiv+coef(mw)[[6]]*dr$precip+coef(mw)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

cor.test(P_totr1_sc,dr$P_wildr1_s)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Foodborne disease (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mf<-glm(anyfoodr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mf)[[1]]+coef(mf)[[2]]*dr$lndensity+coef(mf)[[4]]*dr$high_pop_g+coef(mf)[[5]]*dr$mamdiv+coef(mf)[[6]]*dr$precip+coef(mf)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

cor.test(P_totr1_sc,dr$P_foodr1_s)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vector-borne diseases (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mv<-glm(anyvectr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mv)[[1]]+coef(mv)[[2]]*dr$lndensity+coef(mv)[[4]]*dr$high_pop_g+coef(mv)[[5]]*dr$mamdiv+coef(mv)[[6]]*dr$precip+coef(mv)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

cor.test(P_totr1_sc,dr$P_vectr1_s)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Drug resistant diseases (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

md<-glm(anydrmr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(md)[[1]]+coef(md)[[2]]*dr$lndensity+coef(md)[[4]]*dr$high_pop_g+coef(md)[[5]]*dr$mamdiv+coef(md)[[6]]*dr$precip+coef(md)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

cor.test(P_totr1_sc,dr$P_drmr1_sc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data mining with new drivers...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# land use, pasture, livestock
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calc new drivers
luc <- dr$crop2000-dr$crop1900  # land-use change crop
pc <- dr$pastur2000-dr$pastur1900  # pasture change

# Total model + landuse + pasture + various livestock
m3<-glm(anytotr1~lndensity+high_pop_g+mamdiv+precip+luc+pc+dr$ctlect2000+dr$pigct2000+dr$shpct2000+dr$goatct2000+dr$bfloct2000 +dr$pltyct2000+GDPcap, family=binomial,weights=weight, data=dr)

# Predictions
pred_totr3<-vector("numeric",nrow(dr))
pred_totr3<-1/(1+exp(-(coef(m3)[[1]]+coef(m3)[[2]]*dr$lndensity+coef(m3)[[3]]*dr$high_pop_g+coef(m3)[[4]]*dr$mamdiv+coef(m3)[[5]]*dr$precip+coef(m3)[[6]]*luc+coef(m3)[[7]]*dr$ctlect2000+coef(m3)[[8]]*dr$pigct2000+coef(m3)[[9]]*dr$shpct2000)))

# Scaling
P_totr3<-(1-1/(1+exp(pred_totr3)))
P_totr3_sc<-(P_totr3-min(P_totr3))/(max(P_totr3)-min(P_totr3))
dr$P_totr3_sc<-P_totr3_sc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Response variable (non-wildlife zoonotic disease)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m2z<-glm(anyzoor1~lndensity+high_pop_g+mamdiv+precip+luc+pc+dr$ctlect2000+dr$pigct2000+dr$shpct2000+dr$goatct2000+ dr$bfloct2000 +dr$pltyct2000+GDPcap, family=binomial,weights=weight, data=dr)

pred_zoor3<-vector("numeric",nrow(dr))
pred_zoor3<-1/(1+exp(-(coef(m2z)[[1]]+coef(m2z)[[2]]*dr$lndensity+coef(m2z)[[3]]*dr$high_pop_g+coef(m2z)[[4]]*dr$mamdiv+coef(m2z)[[5]]*dr$precip+coef(m2z)[[6]]*luc+coef(m2z)[[7]]*dr$ctlect2000+coef(m2z)[[8]]*dr$pigct2000+coef(m2z)[[9]]*dr$shpct2000)))

P_zoor3<-(1-1/(1+exp(pred_zoor3)))
P_zoor3_sc<-(P_zoor3-min(P_zoor3))/(max(P_zoor3)-min(P_zoor3))
dr$P_zoor3_sc<-P_zoor3_sc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Response variable (wildlife zoonotic disease)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m2w<-glm(anywildr1~lndensity+high_pop_g+mamdiv+precip+luc+dr$ctlect2000+dr$pigct2000+dr$shpct2000+GDPcap, family=binomial,weights=weight, data=dr)

pred_wildr3<-vector("numeric",nrow(dr))
pred_wildr3<-1/(1+exp(-(coef(m2w)[[1]]+coef(m2w)[[2]]*dr$lndensity+coef(m2w)[[3]]*dr$high_pop_g+coef(m2w)[[4]]*dr$mamdiv+coef(m2w)[[5]]*dr$precip+coef(m2w)[[6]]*luc+coef(m2w)[[7]]*dr$ctlect2000+coef(m2w)[[8]]*dr$pigct2000+coef(m2w)[[9]]*dr$shpct2000)))

P_wildr3<-(1-1/(1+exp(pred_wildr3)))
P_wildr3_sc<-(P_wildr3-min(P_wildr3))/(max(P_wildr3)-min(P_wildr3))
dr$P_wildr3_sc<-P_wildr3_sc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Response variable (vector disease)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m2v<-glm(anyvectr1~lndensity+high_pop_g+mamdiv+precip+luc+dr$ctlect2000+dr$pigct2000+dr$shpct2000+GDPcap, family=binomial,weights=weight, data=dr)

pred_vectr3<-vector("numeric",nrow(dr))
pred_vectr3<-1/(1+exp(-(coef(m2v)[[1]]+coef(m2v)[[2]]*dr$lndensity+coef(m2v)[[3]]*dr$high_pop_g+coef(m2v)[[4]]*dr$mamdiv+coef(m2v)[[5]]*dr$precip+coef(m2v)[[6]]*luc+coef(m2v)[[7]]*dr$ctlect2000+coef(m2v)[[8]]*dr$pigct2000+coef(m2v)[[9]]*dr$shpct2000)))

P_vectr3<-(1-1/(1+exp(pred_vectr3)))
P_vectr3_sc<-(P_vectr3-min(P_vectr3))/(max(P_vectr3)-min(P_vectr3))
dr$P_vectr3_sc<-P_vectr3_sc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Response variable (foodborne disease)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m2f<-glm(anyfoodr1~lndensity+high_pop_g+mamdiv+precip+luc+dr$ctlect2000+dr$pigct2000+dr$shpct2000+GDPcap, family=binomial,weights=weight, data=dr)

pred_foodr3<-foodor("numeric",nrow(dr))
pred_foodr3<-1/(1+exp(-(coef(m2f)[[1]]+coef(m2f)[[2]]*dr$lndensity+coef(m2f)[[3]]*dr$high_pop_g+coef(m2f)[[4]]*dr$mamdiv+coef(m2f)[[5]]*dr$precip+coef(m2f)[[6]]*luc+coef(m2f)[[7]]*dr$ctlect2000+coef(m2f)[[8]]*dr$pigct2000+coef(m2f)[[9]]*dr$shpct2000)))

P_foodr3<-(1-1/(1+exp(pred_foodr3)))
P_foodr3_sc<-(P_foodr3-min(P_foodr3))/(max(P_foodr3)-min(P_foodr3))
dr$P_foodr3_sc<-P_foodr3_sc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Response variable (drug-resistant disease)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m2d<-glm(anydrmr1~lndensity+high_pop_g+mamdiv+precip+luc+dr$ctlect2000+dr$pigct2000+dr$shpct2000+GDPcap, family=binomial,weights=weight, data=dr)

pred_drmr3<-drmor("numeric",nrow(dr))
pred_drmr3<-1/(1+exp(-(coef(m2d)[[1]]+coef(m2d)[[2]]*dr$lndensity+coef(m2d)[[3]]*dr$high_pop_g+coef(m2d)[[4]]*dr$mamdiv+coef(m2d)[[5]]*dr$precip+coef(m2d)[[6]]*luc+coef(m2d)[[7]]*dr$ctlect2000+coef(m2d)[[8]]*dr$pigct2000+coef(m2d)[[9]]*dr$shpct2000)))

P_drmr3<-(1-1/(1+exp(pred_drmr3)))
P_drmr3_sc<-(P_drmr3-min(P_drmr3))/(max(P_drmr3)-min(P_drmr3))
dr$P_drmr3_sc<-P_drmr3_sc

#write table with new models
write.table(dr,"data/dr_new2.txt", sep=",")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Response variable (non-wildlife disease)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# combined index for livestock
lst<-(dr$goatct2000/max(dr$goatct2000)+dr$bfloct2000/max(dr$bfloct2000)+dr$ctlect2000/max(dr$ctlect2000)+dr$shpct2000/max(dr$shpct2000)+dr$pltyct2000/max(dr$pltyct2000)+dr$pigct2000/max(dr$pigct2000))/max(dr$goatct2000/max(dr$goatct2000)+dr$bfloct2000/max(dr$bfloct2000)+dr$ctlect2000/max(dr$ctlect2000)+dr$shpct2000/max(dr$shpct2000)+dr$pltyct2000/max(dr$pltyct2000)+dr$pigct2000/max(dr$pigct2000))

m3<-glm(anyzoor1~lndensity+high_pop_g+mamdiv+precip+lst+luc+pc+lnpubs+GDPcap, family=binomial,weights=weight, data=dr)

pred_zoor3<-vector("numeric",nrow(dr))
pred_zoor3<-1/(1+exp(-(coef(m3)[[1]]+coef(m3)[[2]]*dr$lndensity+coef(m3)[[3]]*dr$high_pop_g+coef(m3)[[4]]*dr$mamdiv+coef(m3)[[5]]*dr$precip+coef(m3)[[6]]*lst)))

P_zoor3<-(1-1/(1+exp(pred_totr3)))
P_zoor3_sc<-(P_totr3-min(P_totr3))/(max(P_totr3)-min(P_totr3))
dr$P_zoor3_sc<-P_totr3_sc
write.table(drzoo,"data/dr_zoo.txt", sep=",")

# todo(EHA) wildr

# todo(EHA) foodr

# todo(EHA) drmr

# todo(EHA) vectr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new grid selection

# read species information
# only those in the nature paper
sp08<-as.data.frame(read.table("/Users/tiffzoo/Dropbox/Tiff Laptop Files/Manuscripts/Drivers of EIDs/Driver_analysis_March2011/sppnat08.TXT",header=TRUE,sep="\t"))
sp08$nat08=1

#all 475
eid<-as.data.frame(read.table("/Users/tiffzoo/Dropbox/Tiff Laptop Files/Manuscripts/Drivers of EIDs/Driver_analysis_March2011/eha_475.TXT",header=TRUE,sep="\t"))
#add nat08 column to eid to check 1/0 if in jones 2008
eid<-merge(x=eid,y=sp08[,c(1,32)],by="EIDID",all.x=TRUE)

#read in grids with EID locations and weight grid cells for EIDs that occur in >1 grid
eid_grid<-as.data.frame(read.spss(file="/Users/tiffzoo/Dropbox/Tiff Laptop Files/EcoHealth Alliance/CIESIN/Non spatial/eid-grid.sav"))
eid_grid<-eid_grid[,c(1,3,4,5)]
eid_grid$prob<-with(eid_grid,1)

#remove duplications (due to being a 0.25 grid originally)
eid_grid<-eid_grid[!duplicated(eid_grid),]

#count the number of grids for each EID event
ngrid<-with(eid_grid,tapply(prob,EIDID,sum))

#create equal distribution of probabilities among grids of the same EID
for(i in 1:length(unique(eid_grid$EIDID))){
	temp<-which(eid_grid$EIDID==unique(eid_grid$EIDID)[i])
	eid_grid[temp,5]<-1/ngrid[[i]]
	}

#take only those EIDs in the original nature paper
eid_grid$nat08=match(eid_grid[,1],sp08[,1],nomatch=0)
eid_grid_nat=eid_grid[-(which(eid_grid$nat08==0)),]

#calculate weighted number of EIDs per grid cell - this is where it would be best to parcel out by time/pathogen type, etc
eidw<-eid_grid[!duplicated(eid_grid[,c(2,3,4)]),c(2,3,4)]
eidw$n_eid<-with(eid_grid,tapply(prob,seqv1,sum))

#nat08 eids only
eidwn<-eid_grid_nat[!duplicated(eid_grid_nat[,c(2,3,4)]),c(2,3,4)]
eidwn$n_eid<-with(eid_grid_nat,tapply(prob,seqv1,sum))



#compute distance between all points using Meeus method, assuming a WGS84 ellipsoid
d<-as.matrix(distm(cbind(dr$lon,dr$lat)))
d.inv<-1/d
diag(d.inv)<-0
Moran.I(pred_totr,d.inv)