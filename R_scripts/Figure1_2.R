library(glmmTMB)
library(maps)
library(raster)
library(plotrix)
library(emmeans)
library(boot)
library(MonoPoly)

#setwd('path_to_your_repository')

dat<-read.csv('data/combined_africa_olfac.csv')
dat<-dat[!is.na(dat$N),]
dat$pref<-(dat$N_Hu-dat$N_Alt)/(dat$N_Hu+dat$N_Alt)
dat$NR<-(dat$N_Hu+dat$N_Alt)
dat$N0<-dat$N-dat$NR
dat$resp<-dat$NR/dat$N
dat$Day<-factor(dat$Day)
dat$Population<-factor(dat$Population,levels=c("ORL",'T51',"NGO","THI","MIN","KED","PKT","BTT","OGD","OHI","KUM","KIN","BOA","AWK","LBV","LPV","LPF","FCV","ENT","ZIK","KAK","SYL","VMB","RAB","KBO","KWA","SHM","GND","ABK"))

gm0<-glmmTMB(cbind(N_Hu,N_Alt)~(1|Day),data=dat,family='betabinomial')
gm<-glmmTMB(cbind(N_Hu,N_Alt)~Population+(1|Colony)+(1|Day),data=dat,family='betabinomial')
#test glm with population differences in behavior vs null model of only day-to-day variation
anova(gm,gm0)

#extract coefficients and confidence intervals
e<-emmeans(gm,'Population')
conf<-confint(e,type='response')
fitted<-2*conf$prob-1
upf<-2*conf$lower.CL-1
dwnf<-2*conf$upper.CL-1
names(fitted)<-conf$Population
conf$pref<-fitted

###############################################################
pdf(width=8,height=4,file='Day_Effect.pdf',useDingbats=F)
par(mar=c(4,3,2,1),tck=-0.01,mgp=c(1.2,.2,0),mfrow=c(2,1))
coords=c(1,2,5,6,7,8,9,10,13,14,17,18,19,22,25,26,27,28,31,32,35,36,37,38,39,40,41,42,43)
cols<-rainbow(29)[as.numeric(dat$Pop)]
pchs<-as.numeric(dat$Day)
plot(jitter(coords[as.numeric(dat$Pop)],1),dat$pref,cex=0.5,axes=F,ylab='Preference index',pch=pchs,col=cols,xlab='',ylim=c(-1,1))
axis(1,las=2,at=coords,labels=c('ORL','T51',popmeta$code),cex.axis=0.7,lwd=0.5)
axis(1,at=c(7.5,13.5,18,22,26.5,31.5,39),line=2,labels=c('Senegal','Burkina Faso','Ghana','Nigeria','Gabon','Uganda','Kenya'),cex.axis=0.7)
axis(2,at=c(-1,0,1),lwd=0.5)

plot(jitter(as.numeric(dat$Day)),dat$pref,cex=0.5,axes=F,ylab='Preference index',xlab='Trial day',col=cols,pch=pchs)
axis(1,at=1:14,las=1,lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
dev.off()
##############################################################

#metadata for populations
popmeta<-read.csv('data/Colonies_kernel.csv',stringsAsFactors=F)
popmeta<-popmeta[match(names(fitted)[3:length(fitted)],popmeta$code),]
popmeta$pref<-fitted[3:length(fitted)]
popmeta$prob<-conf$prob[-c(1,2)]
popmeta$densCex<-log10(popmeta$dens20+10)/1.7
conf$densCex<-c(1,1,popmeta$densCex)
popmeta$bio15col<- c(colorRampPalette(c('darkgreen','yellow','white'))(262)[popmeta$bio15_20])
conf$bio15col<-c('grey','grey',popmeta$bio15col)
popmeta$prefcol<-colorRampPalette(c('cornflowerblue','red'))(100)[round(popmeta$prob*100)]
conf$prefcol<-colorRampPalette(c('cornflowerblue','red'))(100)[conf$prob*100]
write.table(popmeta,file='data/colonies_fitted.txt',sep='\t',quote=F)
write.table(conf,file='data/colonies_conf.txt',sep='\t',quote=F)

#PCA on bioclim variables
bios<-grep('bio.bio',names(popmeta))
pc.out<-prcomp(popmeta[,bios],scale=T)
popmeta$PC1<-pc.out$x[,1]
popmeta$PC2<-pc.out$x[,2]
popmeta$PC3<-pc.out$x[,3]

######################################################################

fits<-matrix(nrow=17,ncol=6)
for(i in 11:27){
	print(colnames(popmeta[i]))
	print(colnames(popmeta)[i+17])
	print(colnames(popmeta)[i+34])
	lm.out_latlon<-lm(logit(prob)~popmeta[,i]+Latitude*Longitude,data=popmeta)
	lm.out_lin<-lm(logit(prob)~popmeta[,i]+popmeta[,i+17],data=popmeta)
	mp2<-monpol(logit(prob)~popmeta[,i+17],data=popmeta,degree=2,a=0,b=261)
	lm.out_mp2<-lm(logit(prob)~popmeta[,i]+offset(fitted(mp2)),data=popmeta)
	mp3<-monpol(logit(prob)~popmeta[,i+17],data=popmeta,degree=3,a=0,b=261)	
	lm.out_mp3<-lm(logit(prob)~popmeta[,i]+offset(fitted(mp3)),data=popmeta)
	mp4<-monpol(logit(prob)~popmeta[,i+17],data=popmeta,degree=4,a=0,b=261)
	lm.out_mp4<-lm(logit(prob)~popmeta[,i]+offset(fitted(mp4)),data=popmeta)
	lm.out_full<-lm(logit(prob)~popmeta[,i]+offset(fitted(mp3))+popmeta[,i+34],data=popmeta)

	fits[i-10,1]=AIC(lm.out_lin)
	fits[i-10,2]=10-2*logLik(lm.out_mp2)
	fits[i-10,3]=12-2*logLik(lm.out_mp3)
	fits[i-10,4]=14-2*logLik(lm.out_mp4)
	fits[i-10,5]=AIC(lm.out_latlon)
	fits[i-10,6]=16-2*logLik(lm.out_full)
}

lm.out_latlon0<-lm(logit(prob)~Latitude*Longitude,data=popmeta)
bufs<-c(seq(5,50,5),seq(60,100,10),200,300)
pdf(file='Model_selection.pdf',width=4,height=4,useDingbats=F)
cols=rainbow(4,v=0.8)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(3,3,1,1))
plot(bufs,fits[,1],type='l',ylim=c(min(fits),max(fits)),lwd=2,bty='l',axes=F,xlab='Radius of adaptation (km)',ylab='AIC',col=cols[1])
axis(1,at=bufs,las=2,labels=F,tck=-0.01,lwd=0.5)
axis(1,at=c(10,50,100,200,300),tck=0,lwd=0)
axis(2,lwd=0.5)
lines(bufs,fits[,3],type='l',col=cols[3],lwd=2)
lines(bufs,fits[,5],type='l',col='darkgrey',lwd=2)
lines(bufs,fits[,6],type='l',col=cols[4],lwd=2)
abline(v=20,lty=2)
abline(h=AIC(lm.out_latlon0),col='lightgrey',lwd=2)
legend('bottomright',legend=c('Lat x Long','Lat x Long + pop. dens.','Linear bio15 + pop. dens.','3rd degree bio15 + pop. dens.','3rd degree bio15 + linear bio18 + pop. dens.'),fill=c('lightgrey','darkgrey',cols[c(1,3,4)]),bty='n',cex=.7)
dev.off()

bufs[which.min(fits[,5])]
lm.out_latlon<-lm(logit(prob)~popmeta[,'dens45']+Latitude*Longitude,data=popmeta)
summary(lm.out_latlon)
anova(lm(logit(prob)~popmeta[,'dens20']+Latitude*Longitude,data=popmeta),lm.out_latlon0)
anova(lm(logit(prob)~popmeta[,'dens25']+Latitude*Longitude,data=popmeta),lm.out_latlon0)
anova(lm(logit(prob)~popmeta[,'dens30']+Latitude*Longitude,data=popmeta),lm.out_latlon0)
anova(lm(logit(prob)~popmeta[,'dens35']+Latitude*Longitude,data=popmeta),lm.out_latlon0)
anova(lm(logit(prob)~popmeta[,'dens40']+Latitude*Longitude,data=popmeta),lm.out_latlon0)
anova(lm(logit(prob)~popmeta[,'dens45']+Latitude*Longitude,data=popmeta),lm.out_latlon0)
anova(lm(logit(prob)~popmeta[,'dens50']+Latitude*Longitude,data=popmeta),lm.out_latlon0)

########################################################################

fits2<-matrix(nrow=22,ncol=4)
curr=1
for(i in c(bios,grep('PC',colnames(popmeta)))){
	currbio=popmeta[,i]
	print(colnames(popmeta)[i])
	lm.out_lin<-lm(logit(prob)~popmeta[,i]+dens20,data=popmeta)
	mp2<-monpol(logit(prob)~currbio,data=popmeta,degree=2,a=min(currbio),b=max(currbio))
	mp3<-monpol(logit(prob)~currbio,data=popmeta,degree=3,a=min(currbio),b=max(currbio))
	mp4<-monpol(logit(prob)~currbio,data=popmeta,degree=4,a=min(currbio),b=max(currbio))
	lm.out_mp2<-lm(logit(prob)~dens20+offset(fitted(mp2)),data=popmeta)
	lm.out_mp3<-lm(logit(prob)~dens20+offset(fitted(mp3)),data=popmeta)
	lm.out_mp4<-lm(logit(prob)~dens20+offset(fitted(mp4)),data=popmeta)
	fits2[curr,1]=AIC(lm.out_lin)
	fits2[curr,2]=10-2*logLik(lm.out_mp2)
	fits2[curr,3]=12-2*logLik(lm.out_mp3)
	fits2[curr,4]=14-2*logLik(lm.out_mp4)
	curr=curr+1
}

pdf(file='Bioclim_selection.pdf',width=4,height=3,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(4,3,1,1))
plot(fits2[,1],col=cols[1],axes=F,xlab='',ylab='AIC',xlim=c(0,23),ylim=c(20,85))
points(fits2[,2],col=cols[2])
points(fits2[,3],col=cols[3])
points(fits2[,4],col=cols[4])
points(23,fits[4,5],col=cols[1])
ax<-c(paste('bio',1:19,sep=''),'PC1','PC2','PC3','LatLong')
axis(1,at=1:23,labels=ax,las=2,cex.axis=0.8,lwd=0.5)
axis(2,lwd=0.5)
mtext('Bioclimate variable',side=1,line=3)
legend('bottomleft',legend=c('Lat x Long + pop. dens.','Linear + pop. dens.','2nd degree + pop. dens.','3rd degree + pop. dens.','4th degree + pop. dens.'),fill=c('grey',cols),bty='n',cex=.7)
dev.off()

###################################################################################

fits3<-matrix(nrow=18,ncol=4)
curr=1
mp15<-monpol(logit(prob)~bio.bio15,data=popmeta,degree=3,a=0,b=261)
resids<-resid(lm(logit(prob)~dens20+offset(fitted(mp15)),data=popmeta))
for(i in paste('bio.bio',c(1:14,16:19),sep='')){
	currbio=resid(lm(popmeta[,i]~offset(fitted(mp15))+popmeta$dens20))
	lm.out_lin<-lm(logit(prob)~currbio+offset(fitted(mp15))+dens20,data=popmeta)
	mp2<-monpol(resids~currbio,data=popmeta,degree=2,a=min(currbio),b=max(currbio))
	mp3<-monpol(resids~currbio,data=popmeta,degree=3,a=min(currbio),b=max(currbio))
	mp4<-monpol(resids~currbio,data=popmeta,degree=4,a=min(currbio),b=max(currbio))
	lm.out_mp2<-lm(logit(prob)~dens20+offset(fitted(mp15))+offset(fitted(mp2)),data=popmeta)
	lm.out_mp3<-lm(logit(prob)~dens20+offset(fitted(mp15))+offset(fitted(mp3)),data=popmeta)
	lm.out_mp4<-lm(logit(prob)~dens20+offset(fitted(mp15))+offset(fitted(mp4)),data=popmeta)
	fits3[curr,1]=14-2*logLik(lm.out_lin)
	fits3[curr,2]=16-2*logLik(lm.out_mp2)
	fits3[curr,3]=18-2*logLik(lm.out_mp3)
	fits3[curr,4]=20-2*logLik(lm.out_mp4)
	curr=curr+1
}
inds<-c(1:14,16:19)


pdf(file='Bioclim18_selection.pdf',width=4,height=3,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(4,3,1,1))
plot(inds,fits3[,1],ylim=c(30,50),col=cols[1],axes=F,xlab='',ylab='AIC')
points(inds,fits3[,2],col=cols[2])
points(inds,fits3[,3],col=cols[3])
points(inds,fits3[,4],col=cols[4])
abline(h=min(fits2),lty=2)
ax<-c(paste('bio',inds,sep=''))
axis(1,at=inds,labels=ax,las=2,cex.axis=0.8,lwd=0.5)
axis(2,lwd=0.5)
dev.off()

##############################################################################
#ANOVA analysis of models
full18<-(lm(logit(prob)~fitted(mp15)+bio.bio18+dens20,data=popmeta))
full0<-(lm(logit(prob)~Latitude*Longitude+dens20,data=popmeta))
full00<-(lm(logit(prob)~Latitude*Longitude,data=popmeta))
stot<-sum(anova(full18)[,'Sum Sq'])
ss1<-anova(full18)[1,'Sum Sq']
ss3<-anova(full18)[2,'Sum Sq']
ss2<-anova(full18)[3,'Sum Sq']

drop15<-(lm(logit(prob)~dens20+bio.bio18,data=popmeta))
drop20<-(lm(logit(prob)~(fitted(mp15))+bio.bio18,data=popmeta))
drop18<-(lm(logit(prob)~(fitted(mp15))+dens20,data=popmeta))
dropclim<-(lm(logit(prob)~dens20,data=popmeta))

sink(file='additional_bioclim.log')
print('pop density 20km LRT w Lat Long')
anova(lm(logit(prob)~Latitude*Longitude+dens20,data=popmeta),lm(logit(prob)~Latitude*Longitude,data=popmeta))
print('pop density 45km LRT w Lat Long')
anova(lm(logit(prob)~Latitude*Longitude+dens45,data=popmeta),lm(logit(prob)~Latitude*Longitude,data=popmeta))


lm.out<-step(lm(logit(prob)~dens20+(fitted(mp15))+bio.bio18+bio.bio9+bio.bio12,data=popmeta))
summary(lm.out)

drop1(full18)
pchisq(2*logLik(full18)-2*logLik(drop15),df=3,lower=F)
pchisq(2*logLik(full18)-2*logLik(drop18),df=1,lower=F)
pchisq(2*logLik(full18)-2*logLik(drop20),df=1,lower=F)
print('LRT climate')
pchisq(2*logLik(full18)-2*logLik(dropclim),df=4,lower=F)

print('eta sq bio15')
ss1/stot
print('eta sq dens20')
ss2/stot
print('eta sq bio18')
ss3/stot

coefs<-c(coef(full18)[1]+coef(mp15)[1],coef(full18)[2]*coef(mp15)[2:4],coef(full18)[3:4])
names(coefs)<-c('Intercept','bio15','bio15_2','bio15_3','bio18','popdensity_20km')
print("Fitted coefficients for bio15, bio18 and density model")
coefs

popmeta$dryseason<-(coefs[1]+coefs[2]*popmeta$bio.bio15+coefs[3]*popmeta$bio.bio15^2+coefs[4]*popmeta$bio.bio15^3+coefs[5]*popmeta$bio.bio18)
aov.out<-anova(lm(logit(prob)~dryseason+dens20,data=popmeta))

print('eta sq dry season')
aov.out[1,'Sum Sq']/sum(aov.out[,'Sum Sq'])
print('eta sq dens20')
aov.out[2,'Sum Sq']/sum(aov.out[,'Sum Sq'])

sink()


##############################################################

#######################################################################################

climpred<-rowSums(predict(full18,type='terms',terms=c('fitted(mp15)','bio.bio18')))
denspred<-predict(full18,type='terms',terms=c('dens20'))
dens20resid<-resid(lm(logit(prob)~denspred,data=popmeta))+logit(mean(popmeta$prob))
bioresid<-resid(lm(logit(prob)~climpred,data=popmeta))+logit(mean(popmeta$prob))
lm.dens<-lm(bioresid~dens20,data=popmeta)
lm.bio<-lm(dens20resid~dryseason,data=popmeta)
lm.bio18<-lm(dens20resid~bio.bio18,data=popmeta)
lm.bio15.2<-lm(dens20resid[!popmeta$code%in%c('NGO','THI')]~bio.bio15,data=popmeta[!popmeta$code%in%c('NGO','THI'),])
lm.bio18.2<-lm(dens20resid[!popmeta$code%in%c('NGO','THI')]~bio.bio18,data=popmeta[!popmeta$code%in%c('NGO','THI'),])



pdf(file='Climate_pop_behavior.pdf',width=6,height=3,useDingbats=F)
par(mfrow=c(1,2),mgp=c(1.5,.2,0),tck=-0.02,bty='l',mar=c(3,3,1,1))
plot((popmeta$dryseason),2*inv.logit(dens20resid)-1,ylim=c(-1,1),pch=21,bg=popmeta$bio15col,xlab='Dry season intensity index',ylab='Preference index (adj. for pop.)',cex=popmeta$densCex,axes=F,xlim=c(-2.3,1),lwd=0.25)
lines(seq(-3,1,0.1),2*inv.logit(predict(lm.bio,data.frame(dryseason=seq(-3,1,0.1))))-1,lty=2)
# abline(0,1)
summary(lm(logit(popmeta$prob)~fitted(lm.bio)))
text(-1.6,.8,labels=expression(italic(eta^2) == 0.65))
axis(1,at=c(-2,-1,0,1),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
plot(popmeta$dens20,2*inv.logit(bioresid)-1,ylim=c(-1,1),xlim=c(0,2500),pch=21,bg=popmeta$bio15col,xlab= expression("Population density km"^-2),ylab='Preference index (adj. for climate)',cex=popmeta$densCex,axes=F,lwd=0.25)
lines(seq(0,3000,30),2*inv.logit(predict(lm.dens,data.frame(dens20=seq(0,3000,30))))-1,lty=2)
summary(lm.dens)
text(500,.8,labels=expression(italic(eta^2) == 0.18))
axis(1,at=c(0,1000,2000),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
dev.off()

pdf(file='bioclim_pop_behavior.pdf',width=2,height=3,useDingbats=F)
par(mfrow=c(2,1),mgp=c(1.5,.2,0),tck=-0.02,bty='l',mar=c(3,3,1,1))
plot(popmeta$bio.bio15,2*inv.logit(dens20resid)-1,pch=21,bg=popmeta$bio15col,xlab='Precip. seasonality',ylab='',cex=popmeta$densCex/1.5,axes=F,lwd=0.25,ylim=c(-1,1))
axis(1,at=c(50,100,150),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
lines(seq(0,200),2*inv.logit(predict(mp15,data.frame(bio.bio15=seq(0,200))))-1,lty=1)
lines(seq(0,200),2*inv.logit(predict(lm.bio15.2,data.frame(bio.bio15=seq(0,200))))-1,lty=1,col='grey')

plot(popmeta$bio.bio18,2*inv.logit(dens20resid)-1,pch=21,bg=popmeta$bio15col,xlab='Precip. of warmest Q',ylab='',cex=popmeta$densCex/1.5,axes=F,lwd=0.25,ylim=c(-1,1))
axis(1,at=c(0,1000,2000),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
lines(seq(0,5000),2*inv.logit(predict(lm.bio18,data.frame(bio.bio18=seq(0,5000))))-1,lty=1)
lines(seq(0,5000),2*inv.logit(predict(lm.bio18.2,data.frame(bio.bio18=seq(0,5000))))-1,lty=1,col='grey')

dev.off()



pdf(file='Climate_bioclim_pop_behavior_combined.pdf',width=8,height=2.5,useDingbats=F)
par(mgp=c(1.5,.2,0),tck=-0.02,bty='l',mar=c(3,4,1,1))
layout(rbind(c(1,2,3),c(1,2,4)),widths=c(3,3,2))


plot(popmeta$dens20,2*inv.logit(bioresid)-1,ylim=c(-1,1),xlim=c(0,2500),pch=21,bg=popmeta$bio15col,xlab= expression("Population density km"^-2),ylab='Preference index (adj. for climate)',cex=popmeta$densCex,axes=F,lwd=0.25)
lines(seq(0,3000,30),2*inv.logit(predict(lm.dens,data.frame(dens20=seq(0,3000,30))))-1,lty=2)
summary(lm.dens)
text(500,.8,labels=expression(italic(eta^2) == 0.18))
axis(1,at=c(0,500,1500,2500),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
box(bty='l',lwd=0.5)


plot((popmeta$dryseason),2*inv.logit(dens20resid)-1,ylim=c(-1,1),pch=21,bg=popmeta$bio15col,xlab='Dry season intensity index',ylab='Preference index (adj. for pop.)',cex=popmeta$densCex,axes=F,xlim=c(-2.3,1),lwd=0.25)
lines(seq(-3,1,0.1),2*inv.logit(predict(lm.bio,data.frame(dryseason=seq(-3,1,0.1))))-1,lty=2)
# abline(0,1)
summary(lm(logit(popmeta$prob)~fitted(lm.bio)))
text(-1.6,.8,labels=expression(italic(eta^2) == 0.65))
axis(1,at=c(-2,-1,0,1),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
box(bty='l',lwd=0.5)

par(mgp=c(1.5,.2,0),tck=-0.02,bty='l',mar=c(3,5,1,1))

plot(popmeta$bio.bio15,2*inv.logit(dens20resid)-1,pch=21,bg=popmeta$bio15col,xlab='Precip. seasonality',ylab='',cex=popmeta$densCex/1.5,axes=F,lwd=0.25,ylim=c(-1,1))
axis(1,at=c(50,100,150),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
lines(seq(0,200),2*inv.logit(predict(mp15,data.frame(bio.bio15=seq(0,200))))-1,lty=1)
lines(seq(0,200),2*inv.logit(predict(lm.bio15.2,data.frame(bio.bio15=seq(0,200))))-1,lty=1,col='grey')
box(bty='l',lwd=0.5)

plot(popmeta$bio.bio18,2*inv.logit(dens20resid)-1,pch=21,bg=popmeta$bio15col,xlab='Precip. of warmest Q',ylab='',cex=popmeta$densCex/1.5,axes=F,lwd=0.25,ylim=c(-1,1),xlim=c(0,900))
axis(1,at=c(0,400,800),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
lines(seq(0,5000),2*inv.logit(predict(lm.bio18,data.frame(bio.bio18=seq(0,5000))))-1,lty=1)
lines(seq(0,5000),2*inv.logit(predict(lm.bio18.2,data.frame(bio.bio18=seq(0,5000))))-1,lty=1,col='grey')
box(bty='l',lwd=0.5)

dev.off()


#######################################################################################
#alternative GAM-based analysis recovers same model

library(mgcv)
gm.out<-(gam(logit(prob)~dens20+s(bio.bio15,k=5)+s(bio.bio18,k=5)+s(Longitude,Latitude,k=10),data=popmeta,method='REML',select=T))
summary(gm.out)
# gam.check(gm.out)
# par(mfrow=c(1,2))
# plot(gm.out)


######################################################################
##########MAP

r<-raster('wc2.0_bio_10m_15.tif')
pdf(width=6,height=6,file='Map.pdf',useDingbats=F)
par(mar=c(1,1,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(r,xlim=c(-25,55),ylim=c(-40,40),col=colorRampPalette(c('darkgreen','yellow','white'))(200),xlab='',ylab='',bty='n',box=F,axes=F,legend=F,zlim=c(0,262))
segments(popmeta$plotLong,popmeta$plotLat,popmeta$Longitude,popmeta$Latitude,lwd=0.5)
for(i in 1:nrow(popmeta)){
	draw.circle((popmeta$plotLong)[i],(popmeta$plotLat)[i],radius=log10(popmeta$dens20[i]+10)/1.3,col=popmeta$prefcol[i],lwd=0.25)
}
draw.circle(-10,-8,radius=log10(10+10)/1.3,col='black')
draw.circle(-5,-8,radius=log10(10+100)/1.3,col='black')
draw.circle(0,-8,radius=log10(10+1000)/1.3,col='black')
gradient.rect(-15,-17,5,-15,col=colorRampPalette(c('cornflowerblue','red'))(200))
text(c(-15,-5,5),rep(-13,5),c(-1,0,1))
gradient.rect(-15,-28,5,-26,col=colorRampPalette(c('darkgreen','yellow','white'))(200))
text(-5,-30,'Precipitation seasonality')
text(-5,-19,'Preference index')
text(-5,-4,expression("Pop. density km"^-2))
text(-15+(20/261)*seq(0,200,100),rep(-24,5),seq(0,200,100))
text(-10,-8,'10',col='white',cex=0.5)
text(-5,-8,'100',col='white',cex=0.5)
text(0,-8,'1000',col='white',cex=0.5)
dev.off()

######################


pdf(width=6,height=6,file='Map2.pdf',useDingbats=F)
par(mar=c(1,1,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(r,xlim=c(-18,0),ylim=c(0,18),col=colorRampPalette(c('darkgreen','yellow','white'))(200),xlab='',ylab='',bty='n',box=F,axes=F,legend=F,zlim=c(100,200))
plot(r,xlim=c(-18,0),ylim=c(0,18),col=colorRampPalette(c('darkgreen'))(200),xlab='',ylab='',bty='n',box=F,axes=F,legend=F,zlim=c(0,100),add=T)
map('world',add=T)
segments(popmeta$plotLong,popmeta$plotLat,popmeta$Longitude,popmeta$Latitude,lwd=0.5)
gradient.rect(-15,0,-10,1,col=colorRampPalette(c('darkgreen','yellow','white'))(200))
text(-12.5,2,'Precipitation seasonality')
text(-5,-19,'Preference index')
text(-5,-4,expression("Pop. density km"^-2))
text(seq(-15,-10,2.5),rep(1.5,3),seq(100,200,50))
text(-10,-8,'10',col='white',cex=0.5)
text(-5,-8,'100',col='white',cex=0.5)
text(0,-8,'1000',col='white',cex=0.5)
dev.off()



################
psub<-popmeta[match(c('RAB','KBO','GND','ABK','KWA','SHM','KIN','BOA','LPV','LPF','KED','PKT','BTT'),popmeta$code),]
psub$peri<-as.numeric(factor(c(rep(c(2,1),6),1)))
psub$pairing=factor(c('RAB','RAB','GND','GND','KWA','KWA','KIN','KIN','LP','LP','KED','KED','KED'))

pdf(width=2,height=4,file='Habitat.pdf',useDingbats=F)
par(mar=c(4,3,1,1),tck=-0.01,mgp=c(1.2,.2,0))
pp=jitter(psub$peri,0.5)
plot(pp,psub$pref,cex=2,pch=21,bg=psub$prefcol,axes=F,xlab='',ylab='Preference index',ylim=c(-1,1),xlim=c(0.5,2.5),type='n')
points(pp[psub$pairing=='RAB'],psub$pref[psub$pairing=='RAB'], pch=21,bg=psub$prefcol[psub$pairing=='RAB'],cex=psub$densCex[psub$pairing=='RAB']*2.5,lwd=0.25)
segments(pp[c(1)],psub$pref[c(1)],pp[c(2)],psub$pref[c(2)])
points(pp[psub$pairing=='GND'],psub$pref[psub$pairing=='GND'], pch=21,bg=psub$prefcol[psub$pairing=='GND'],cex=psub$densCex[psub$pairing=='GND']*2.5,lwd=0.25)
segments(pp[c(3)],psub$pref[c(3)],pp[c(4)],psub$pref[c(4)])
points(pp[psub$pairing=='KWA'],psub$pref[psub$pairing=='KWA'], pch=21,bg=psub$prefcol[psub$pairing=='KWA'],cex=psub$densCex[psub$pairing=='KWA']*2.5,lwd=0.25)
segments(pp[c(5)],psub$pref[c(5)],pp[c(6)],psub$pref[c(6)])
points(pp[psub$pairing=='KIN'],psub$pref[psub$pairing=='KIN'], pch=21,bg=psub$prefcol[psub$pairing=='KIN'],cex=psub$densCex[psub$pairing=='KIN']*2.5,lwd=0.25)
segments(pp[c(7)],psub$pref[c(7)],pp[c(8)],psub$pref[c(8)])
points(pp[psub$pairing=='LP'],psub$pref[psub$pairing=='LP'], pch=21,bg=psub$prefcol[psub$pairing=='LP'],cex=psub$densCex[psub$pairing=='LP']*2.5,lwd=0.25)
segments(pp[c(9)],psub$pref[c(9)],pp[c(10)],psub$pref[c(10)])
points(pp[psub$pairing=='KED'],psub$pref[psub$pairing=='KED'], pch=21,bg=psub$prefcol[psub$pairing=='KED'],cex=psub$densCex[psub$pairing=='KED']*2.5,lwd=0.25)
segments(pp[c(11,11)],psub$pref[c(11,11)],pp[c(12,13)],psub$pref[c(12,13)])
axis(1,at=c(1,2),labels=c('Forest','Town'),lwd=0,cex.axis=0.7,lwd=0.5)
axis(2,at=c(-1,-.5,0,.5,1),lwd=0.5)
dev.off()

summary(lm(logit(prob)~pairing+peri,data=psub))
############PREFs

pdf(width=8,height=4,file='Prefs.pdf',useDingbats=F)
par(mar=c(4,3,2,1),tck=0.01,mgp=c(1.2,.2,0))
coords=c(1,2,5,6,7,8,9,10,13,14,17,18,19,22,25,26,27,28,31,32,35,36,37,38,39,40,41,42,43)
plot(coords,fitted,pch=19,col=colorRampPalette(c('cornflowerblue','red'))(100)[round(conf$prob*100)],bty='l',ylim=c(-1,1.1),axes=F,xlab='',ylab='Preference index',type='n')
cols= colorRampPalette(c('cornflowerblue','red'))(100)[round((fitted+1)*50)]
for(i in 1:length(fitted)){
	draw.circle(coords[i],fitted[i],radius=conf$densCex[i]/2.5,col=cols[i],border=NA)
}
arrows(coords,dwnf,coords,upf,angle=90,code=3,length=0.02,col=conf$bio15col,lwd=4)
arrows(coords,dwnf,coords,upf,angle=90,code=3,length=0.02)
axis(1,las=2,at=coords,labels=c('ORL','T51',popmeta$code),lwd=0.5)
axis(1,at=c(7.5,13.5,18,22,26.5,31.5,39),line=2,labels=c('Senegal','Burkina Faso','Ghana','Nigeria','Gabon','Uganda','Kenya'),lwd=0.5)
axis(2,lwd=0.5)
dev.off()



#########################################################################################

#Colony-level model for repeatability analysis
gm3<-glmmTMB(cbind(N_Hu,N_Alt)~Colony+(1|Day),data=dat,family='betabinomial')

e3<-emmeans(gm3,'Colony')
conf3<-confint(e3,type='response')

uq<-dat[match(unique(dat$Colony),dat$Colony),]
pairs<-names(table(uq$Population))[table(uq$Population)==2]
uq<-uq[uq$Population%in%pairs,]
uq<-uq[order(uq$Population),]
first<-as.character(uq$Colony[match(unique(uq$Population),uq$Population)])
second<-as.character(uq[!uq$Colony%in%first,'Colony'])
uq<-uq[match(first,uq$Colony),]
fp<-2*(conf3$prob-0.5)
names(fp)<-conf3$Colony

summary(lm(fp[first]~fp[second]))

pdf(width=3,height=3,file='Repeatability.pdf',useDingbats=F)
par(mar=c(4,4,2,1),mgp=c(1.5,.2,0),tck=-0.01)
plot(fp[first],fp[second],xlab='Colony 1 preference index',ylab='Colony 2 preference index',bty='l',pch=21,bg=popmeta$bio15col[match(uq$Population,c(as.character(popmeta$code)))],axes=F,xlim=c(-1,1),ylim=c(-1,1),cex=popmeta$densCex[match(uq$Population,c(as.character(popmeta$code)))],lwd=0.25)
text(-.6,.8,labels=expression(italic(R^2) == 0.60))
axis(1,at=c(-1,0,1),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
abline(0,1,lty=2)
dev.off()


gm4<-glmmTMB(cbind(NR,N0)~Population+(1|Colony)+(1|Day),data=dat,family='betabinomial')
e4<-emmeans(gm4,'Population')
conf4<-confint(e4,type='response')
fittedr<-conf4$prob
lm.out<-(lm(fittedr~fitted))
summary(lm.out)

pdf(width=3,height=3,file='Responserate.pdf',useDingbats=F)
par(mar=c(4,4,2,1),mgp=c(1.5,.2,0),tck=-0.01)
plot(fitted,fittedr,bg=popmeta$bio15col,pch=21,xlab='Preference index',ylab='Response rate',bty='l',axes=F,ylim=c(0,1),xlim=c(-1,1),cex=popmeta$densCex,lwd=0.25)
text(-.6,.8,labels=expression(italic(R^2) == 0.73))
abline(lm.out,lty=2)
axis(1,at=c(-1,0,1),lwd=0.5)
axis(2,at=c(0,0.5,1),lwd=0.5)
dev.off()

############RESPs

pdf(width=12,height=4,file='Responses.pdf',useDingbats=F)
par(mar=c(4,3,2,1),tck=0.01,mgp=c(1.2,.2,0))
coords=c(1,2,5,6,7,8,9,10,13,14,17,18,19,22,25,26,27,28,31,32,35,36,37,38,39,40,41,42,43)
plot(coords,fitted,pch=19,col=colorRampPalette(c('cornflowerblue','red'))(100)[round(conf$prob*100)],bty='l',ylim=c(0,1.1),axes=F,xlab='',ylab='Preference index',type='n',lwd=0.25)
cols= colorRampPalette(c('cornflowerblue','red'))(100)[round((fitted+1)*50)]
for(i in 1:length(fitted)){
	draw.circle(coords[i],fittedr[i],radius=conf$densCex[i]/3,col=cols[i],border=NA)
}
arrows(coords,conf4$lower,coords, conf4$upper,angle=90,code=3,length=0.02,col=conf$bio15col,lwd=4)
arrows(coords, conf4$lower,coords, conf4$upper,angle=90,code=3,length=0.02)
axis(1,las=2,at=coords,labels=c('ORL','T51',popmeta$code),lwd=0.5)
axis(1,lwd=0,at=c(7.5,13.5,18,22,26.5,31.5,39),line=2,labels=c('Senegal','Burkina Faso','Ghana','Nigeria','Gabon','Uganda','Kenya'))
axis(2,lwd=0.5)
dev.off()


##################################################################################

dsub<-dat[dat$Population%in%c('NGO','OGD','AWK','FCV'),1:17]
dsub$Matchup<-paste('HU1vs',dsub$Host,sep='')
alt<-read.csv('data/alternate_host_olfactometer_trials.csv')[,c(1:17,20)]
alt<-rbind(alt,dsub)
alt$Population<-dat$Population[match(alt$Colony,dat$Colony)]
alt$Population[alt$Country=='Nigeria']<-'AWK'
alt$Population<-factor(alt$Population)
alt$pref<-(alt$N_Hu-alt$N_Alt)/(alt$N_Hu+alt$N_Alt)
alt$Matchup<-factor(alt$Matchup,levels=c('HU1vsMoja','HU1vsMbili','HU2vsMoja','HU1vsAfi'))
alt$int<-interaction(alt$Population,alt$Matchup)

agg<-aggregate(cbind(N_Hu,N_Alt)~Population+Matchup,data=alt,FUN=sum)
agg$pref<-(agg$N_Hu-agg$N_Alt)/(agg$N_Hu+agg$N_Alt)

gm5<-glmmTMB(cbind(N_Hu,N_Alt)~Matchup+Population+(1|Colony)+(1|Day),data=alt,family='betabinomial')
e5<-emmeans(gm5,c('Population','Matchup'))
conf5<-confint(e5,type='response')

pdf(width=4,height=4,file='FigAltHost.pdf',useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.2,0),tck=-0.01)
#pos=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19)
pos=c(1:4-.15,1:4-.05,1:4+.05,1:4+.15)
plot(pos,agg$pref,ylim=c(-1.5,1.5),xlab='',ylab='Preference index',axes=F,pch=19,cex=1,type='n')
lines(1:4-.15,agg$pref[1:4],type='b',cex=1.5)
lines(1:4-.05,agg$pref[5:8],type='b',cex=1.5)
lines(1:4+0.05,agg$pref[9:12],type='b',cex=1.5)
lines(1:4+.15,agg$pref[13:16],type='b',cex=1.5)
# arrows(pos,2*conf5$lower-1,pos,2*conf5$upper-1, angle=90,code=3,length=0,lwd=1)
# points(jitter(pos[as.numeric(alt$int)],.0),alt$pref,col='grey',pch=19,cex=0.5)
axis(1,at=1:4,labels=(c('NGO','OGD','AWK','FCV')),las=1,tck=-0.01,lwd=0)
axis(2,at=c(-1,0,1),lwd=0.5)
legend('topright',legend=c('HU1','HU2','GP1','GP2','QU'),bty='n',fill=c('orange','red','steelblue','blue','magenta'))
dev.off()

#################################################################################
