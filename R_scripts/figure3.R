setwd('~/Dropbox/2019 Rose et al AaegVariationAfrica/Figures/Figure3 stuff/')
dat<-read.csv('proj.csv',stringsAsFactors=F)
pm<-read.delim('~/Dropbox/2019 Rose et al AaegVariationAfrica/Figures/Figure1 stuff/colonies_fitted.txt',stringsAsFactors=F)
pm<-pm[match(dat$code,pm$code),]

Intercept =  -3.204739e+00
bio15 = 7.319490e-02
bio15_2 = -9.742814e-04
bio15_3 =  4.322816e-06
bio18 = -9.715850e-04
dens20km  = 5.873513e-04
require("boot")
dat$out2050climOnly <- Intercept + bio15*(dat$bio152050) + bio15_2*(dat$bio152050)^2 +
  bio15_3*(dat$bio152050)^3 + bio18*(dat$bio182050) + dens20km*dat$dens20
dat$out2050climOnly <- 2*inv.logit(dat$out2050climOnly)-1
dat$out2050PopOnly <- Intercept + bio15*(dat$bio.bio15 ) + bio15_2*(dat$bio.bio15)^2 +
  bio15_3*(dat$bio.bio15)^3 + bio18*(dat$bio.bio18) + dens20km*dat$p2050
dat$out2050PopOnly <- 2*inv.logit(dat$out2050PopOnly)-1
dat$out2050<- Intercept + bio15*(dat$bio152050) + bio15_2*(dat$bio152050)^2 +
  bio15_3*(dat$bio152050)^3 + bio18*(dat$bio182050) + dens20km*dat$p2050
dat$out2050 <- 2*inv.logit(dat$out2050)-1
dat$out2000<- Intercept + bio15*(dat$bio.bio15 ) + bio15_2*(dat$bio.bio15 )^2 +
  bio15_3*(dat$bio.bio15 )^3 + bio18*(dat$bio.bio18) + dens20km*dat$dens20
dat$out2000 <- 2*inv.logit(dat$out2000)-1

dat$out2020dryseason <- Intercept + bio15*(dat$bio.bio15) + bio15_2*(dat$bio.bio15)^2 +
  bio15_3*(dat$bio.bio15)^3 + bio18*(dat$bio.bio18)
dat$out2050dryseason <- Intercept + bio15*(dat$bio152050) + bio15_2*(dat$bio152050)^2 +
  bio15_3*(dat$bio152050)^3 + bio18*(dat$bio182050)
dat$prefcol2020<-colorRampPalette(c('cornflowerblue','red'))(100)[round(((dat$out2000+1)/2)*100)]
dat$prefcol2050<-colorRampPalette(c('cornflowerblue','red'))(100)[round(((dat$out2050+1)/2)*100)]
dat$densCex2020<-log10(dat$dens20+10)/1.7
dat$densCex2050<-log10(dat$p2050+10)/1.7

pdf(file='Alt_plot.pdf',width=4,height=4,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(3,3,1,1))
plot(dat$out2020dryseason,dat$dens20,xlim=c(-2.5,1),ylim=c(0,10000),xlab='Dry season intensity',ylab='Population density',type='n')
segments(dat$out2020dryseason, dat$dens20, dat$out2050dryseason, dat$p2050)
points(dat$out2020dryseason, dat$dens20,xlim=c(0,261),ylim=c(0,10000), pch=21,bg=dat$prefcol2020,cex=dat$densCex2020)
points(dat$out2050dryseason, dat$p2050,xlim=c(0,261),ylim=c(0,10000),col=dat$prefcol2050,cex=dat$densCex2050)
dev.off()

clims<-seq(-2.5,1.5,length.out=100)
denss<-seq(0,10000,length.out=100)
mat<-matrix(nrow=100,ncol=100)
for(i in 1:100){
	for(j in 1:100){
		mat[i,j]<-round(100*inv.logit(clims[i]+dens20km*denss[j]))
	}
}

pdf(file='Alt_plot_gradient.pdf',width=4,height=4,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(3,3,1,1))
image(clims,denss,mat,col=colorRampPalette(c('cornflowerblue','red'))(100),xlab='Dry season intensity',ylab='Population density')
arrows(dat$out2020dryseason, dat$dens20, dat$out2050dryseason, dat$p2050,length=0.08)
points(dat$out2020dryseason, dat$dens20,xlim=c(0,261),ylim=c(0,10000), pch=21,bg=pm$prefcol,cex=pm$densCex*1.5)
dev.off()

pdf(file='Alt_plot_gradient_labels.pdf',width=4,height=4,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(3,3,1,1))
image(clims,denss,mat,col=colorRampPalette(c('cornflowerblue','red'))(100),xlab='Dry season intensity',ylab='Population density')
arrows(dat$out2020dryseason, dat$dens20, dat$out2050dryseason, dat$p2050,length=0.08)
points(dat$out2020dryseason, dat$dens20,xlim=c(0,261),ylim=c(0,10000), pch=21,bg=pm$prefcol,cex=pm$densCex*1.5)
text(dat$out2020dryseason, dat$dens20,labels=dat$code,cex=0.5)
dev.off()

pdf(file='Alt_plot_gradient_bgonly.pdf',width=4,height=4,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(3,3,1,1))
image(clims,denss,mat,col=colorRampPalette(c('cornflowerblue','red'))(100),xlab='Dry season intensity',ylab='Population density')
dev.off()

pdf(file='Alt_plot_gradient_circles.pdf',width=4,height=4,useDingbats=F)
par(mgp=c(1.2,.2,0),tck=-0.01,bty='l',mar=c(3,3,1,1))
image(clims,denss,mat,col=rep(NULL,100),xlab='Dry season intensity',ylab='Population density')
arrows(dat$out2020dryseason, dat$dens20, dat$out2050dryseason, dat$p2050,length=0.08)
points(dat$out2020dryseason, dat$dens20,xlim=c(0,261),ylim=c(0,10000), pch=21,bg=pm$prefcol,cex=pm$densCex*1.5)
dev.off()
