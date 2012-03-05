# HSD Project 
# Collaborators (EHA and CIESIN)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())      # Clear all variables
graphics.off()       # Close graphics windows

# Load libraries

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Todo(EHA) turn all models into a single loop with multipage output

dr<-read.table(file="data/eid08_drivers_19OCT11.csv",sep=",",header=TRUE)

# repeated regression analysis and scaling according to Jones et al 2008 to verify methods
# original values and coef are in the right range

# this was the weight used in the original analysis
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Zoonotic non-wildlife disease (reponse variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mz<-glm(anyzoor1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mz)[[1]]+coef(mz)[[2]]*dr$lndensity+coef(mz)[[4]]*dr$high_pop_g+coef(mz)[[5]]*dr$mamdiv+coef(mz)[[6]]*dr$precip+coef(mz)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

# todo (EHA) Add original values
# cor.test(P_totr1_sc,dr$P_zoor1_sc) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wildlife zoonotic disease (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mw<-glm(anywildr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mw)[[1]]+coef(mw)[[2]]*dr$lndensity+coef(mw)[[4]]*dr$high_pop_g+coef(mw)[[5]]*dr$mamdiv+coef(mw)[[6]]*dr$precip+coef(mw)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

# todo (EHA) Add original values
# cor.test(P_totr1_sc,dr$P_wildr1_s)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Foodborne disease (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mf<-glm(anyfoodr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mf)[[1]]+coef(mf)[[2]]*dr$lndensity+coef(mf)[[4]]*dr$high_pop_g+coef(mf)[[5]]*dr$mamdiv+coef(mf)[[6]]*dr$precip+coef(mf)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

# todo (EHA) Add original values
# cor.test(P_totr1_sc,dr$P_foodr1_s)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vector-borne diseases (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mv<-glm(anyvectr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(mv)[[1]]+coef(mv)[[2]]*dr$lndensity+coef(mv)[[4]]*dr$high_pop_g+coef(mv)[[5]]*dr$mamdiv+coef(mv)[[6]]*dr$precip+coef(mv)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

# todo (EHA) Add original values
# cor.test(P_totr1_sc,dr$P_vectr1_s)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Drug resistant diseases (response variable)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

md<-glm(anydrmr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

pred_totr<-vector("numeric",nrow(dr))
pred_totr<-1/(1+exp(-(coef(md)[[1]]+coef(md)[[2]]*dr$lndensity+coef(md)[[4]]*dr$high_pop_g+coef(md)[[5]]*dr$mamdiv+coef(md)[[6]]*dr$precip+coef(md)[[7]]*abslat))) 

P_totr1<-1-1/(1+exp(pred_totr))
P_totr1_sc<-(P_totr1-min(P_totr1))/(max(P_totr1)-min(P_totr1))

# todo (EHA) Add original values
# cor.test(P_totr1_sc,dr$P_drmr1_sc)

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
