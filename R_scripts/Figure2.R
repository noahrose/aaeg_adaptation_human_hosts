library(lfmm)
library(boot)
library(zoo)
library(permute)
library(data.table)
library(emmeans)

#setwd('path_to_your_repository')

sm<-read.delim('data/simple_meta.txt',header=F,stringsAsFactors=F)
colnames(sm)<-c('Sample','Pop')
mm<-sm
mm$Pop<-factor(sm$Pop,levels=c("BKK","SAN","NGO","THI","MIN","KED","PKT","BTT","OGD","OHI","KUM","KIN","BOA","AWK","LBV","LPV","FCV","ENT","KAK","VMB","SHM","KWA","ABK","GND","KBO","RAB","RABd","MASC"))
sm$Pop<-factor(sm$Pop,levels=c("BKK","SAN","NGO","THI","MIN","KED","PKT","BTT","OGD","OHI","KUM","KIN","BOA","AWK","LBV","LPV","FCV","ENT","KAK","VMB","SHM","KWA","ABK","GND","KBO","RAB","RABd"))
m<-na.omit(sm[order(sm$Pop),])
m$Country<-rep('Kenya',nrow(m))
m$Country[m$Pop%in%c('BKK','SAN')]<-'outside Africa'
m$Country[m$Pop%in%c('NGO','THI','MIN','KED','PKT','BTT')]<-'Senegal'
m$Country[m$Pop%in%c('OGD','OHI')]<-'Burkina Faso'
m$Country[m$Pop%in%c('KUM','KIN','BOA')]<-'Ghana'
m$Country[m$Pop%in%c('AWK')]<-'Nigeria'
m$Country[m$Pop%in%c('LBV','LPV','FCV')]<-'Gabon'
m$Country[m$Pop%in%c('ENT')]<-'Uganda'

########################################
aaa<-read.delim('data/verily_unlinked_norepeat.3.Q',header=F,sep=' ')
sm$aaa<-NA
sm$aaa[1:nrow(aaa)]<-aaa[,3]
sm$PopNoOutl<-as.character(sm$Pop)
sm$PopNoOutl[sm$Pop%in%c('RABd')&sm$aaa<0.99]<-NA
sm$PopNoOutl[sm$Pop%in%c('NGO')&sm$aaa<0.2]<-NA
sm$fid<-0
clst<-na.omit(sm[,c('fid','Sample','PopNoOutl')])

write.table(clst,quote=F,sep='\t',row.names=F,col.names=F,file='data/unrelated_nooutl.clst')

#######################################

q2<-read.delim('admixture_analyses/verily_unlinked_norepeat.2.Q',header=F,sep=' ')[order(sm$Pop),][1:nrow(m),]
q3<-read.delim('admixture_analyses/verily_unlinked_norepeat.3.Q',header=F,sep=' ')[order(sm$Pop),][1:nrow(m),]
q4<-read.delim('admixture_analyses/verily_unlinked_norepeat.4.Q',header=F,sep=' ')[order(sm$Pop),][1:nrow(m),]
q5<-read.delim('admixture_analyses/verily_unlinked_norepeat.5.Q',header=F,sep=' ')[order(sm$Pop),][1:nrow(m),]
q6<-read.delim('admixture_analyses/verily_unlinked_norepeat.6.Q',header=F,sep=' ')[order(sm$Pop),][1:nrow(m),]
ord<-order(m$Pop,round(q3[,2]+q3[,1],2),round(q3[,1],2))
rownames(q3)<-sm$Sample[order(sm$Pop)][1:nrow(q3)]

tbl<-table(sm$Pop)
spaces<-c()
for(i in 1:length(tbl)){
	spaces<-c(spaces,1,rep(0,(tbl[i]-1)))
}


pdf('admixture.pdf',width=8,height=3,useDingbats=F)
par(mfrow=c(2,1),mar=c(0,0,0.5,0),oma=c(3,1,0,0),mgp=c(1.5,.5,0),tck=-0.01,bty='l')
bp<-barplot(t(q3[ord,]),names.arg=rep('',nrow(m)),las=2,col=c('#bdd7e7','#08519c','red'),border=NA,space=spaces,axes=F,ylab='',cex.lab=.8)
mtext('K=3',side=2,line=-1,outer=F)
bp<-barplot(t(q6[ord,c(1,5,4,6,2,3)]),names.arg=rep('',nrow(m)),las=2,col=c('#bdd7e7','#6baed6','#3182bd','#08519c','orange','red'),border=NA,space=spaces,axes=F,cex.lab=.8)
mtext('K=6',side=2,line=-1,outer=F)
agg<-aggregate(bp~m$Pop,FUN=mean)
axis(1,at=agg[,2],labels=agg[,1],lwd=0,cex.axis=0.7,las=2)
agg<-aggregate(bp~m$Country,FUN=mean)
axis(1,at=agg[,2],labels=agg[,1],lwd=0,cex.axis=0.7,las=1,line=1.5)
mtext('Ancestry proportion',side=2,line=0,outer=T)
dev.off()

pdf('admixture_k.pdf',width=8,height=8,useDingbats=F)
par(mar=c(2,3,1,1),mgp=c(1.5,0.5,0),tck=-0.01,bty='l',mfrow=c(5,1),oma=c(3,0,2,0))
bp<-barplot(t(q2[ord,]),names.arg=rep('',nrow(m)),las=2,col=c('#08519c','red'),border=NA,space=spaces,axes=F,ylab='Ancestry prop.',cex.lab=.8,main='K=2')
axis(2,at=c(0,1))
bp<-barplot(t(q3[ord,]),names.arg=rep('',nrow(m)),las=2,col=c('#bdd7e7','#08519c','red'),border=NA,space=spaces,axes=F,ylab='Ancestry prop.',cex.lab=.8,main='K=3')
axis(2,at=c(0,1))
bp<-barplot(t(q4[ord,c(2,4,1,3)]),names.arg=rep('',nrow(m)),las=2,col=c('#bdd7e7','#6baed6','#08519c','red'),border=NA,space=spaces,axes=F,ylab='Ancestry prop.',cex.lab=.8,main='K=4')
axis(2,at=c(0,1))
bp<-barplot(t(q5[ord,c(4,5,3,2,1)]),names.arg=rep('',nrow(m)),las=2,col=c('#bdd7e7','#3182bd','#6baed6','#08519c','red'),border=NA,space=spaces,axes=F,ylab='Ancestry prop.',cex.lab=.8,main='K=5')
axis(2,at=c(0,1))
bp<-barplot(t(q6[ord,c(1,4,5,6,2,3)]),names.arg=rep('',nrow(m)),las=2,col=c('#bdd7e7','#3182bd','#6baed6','#08519c','orange','red'),border=NA,space=spaces,axes=F,ylab='Ancestry prop.',cex.lab=.8,main='K=6')
axis(2,at=c(0,1))
agg<-aggregate(bp~m$Pop,FUN=mean)
axis(1,at=agg[,2],labels=agg[,1],lwd=0,cex.axis=0.7,las=2)
agg<-aggregate(bp~m$Country,FUN=mean)
axis(1,at=agg[,2],labels=agg[,1],lwd=0,cex.axis=0.7,las=1,line=1.5)
axis(2,at=c(0,1))
dev.off()


pm<-read.delim('data/colonies_conf.txt',stringsAsFactors=F)
agg<-aggregate(cbind(V1,V2,V3)~m$Pop,FUN=mean,data=q3)
pm$q3.west<-agg[match(pm$Population,agg[,1]),2]
pm$q3.east<-agg[match(pm$Population,agg[,1]),3]
pm$q3.aaa<-agg[match(pm$Population,agg[,1]),4]

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
cols2<-c('#bdd7e7','#3182bd','#6baed6','#08519c','orange','red')[apply(q6[,c(1,4,5,6,2,3)],1,which.max)]
agg<-aggregate(as.numeric(factor(cols2,levels=c('#bdd7e7','#3182bd','#6baed6','#08519c','orange','red')))~as.character(m$Pop),FUN=getmode)
agg[agg[,1]=='SAN',1]<-'ORL'
agg[agg[,1]=='BKK',1]<-'T51'
pm$k6col<-c('#bdd7e7','#3182bd','#6baed6','#08519c','orange','red')[agg[,2]][match(pm$Population,agg[,1])]


pdf('crossValidation.pdf',width=4,height=4,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01,bty='l')
err<-read.delim('admixture_analyses/cv.out',sep=' ',header=F)[c(2:10,1),]
plot(1:10,err[,4],type='b',xlab='K value',ylab='Cross-validation error',axes=F)
axis(1,at=c(1:10),labels=NA,tck=-0.01)
axis(1,at=c(1,3,5,7,9),lwd=0)
axis(2)
dev.off()


lm.out<-lm(logit(prob)~q3.aaa,data=pm)
summary(lm.out)

pdf('admixture_behavior_k6.pdf',width=3,height=3,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01,bty='l')
plot(pm$q3.aaa,pm$pref,cex=2,bg=pm$k6col,pch=21,ylim=c(-1,1),xlim=c(0,.42),xlab='Proportion Aaa',ylab='Preference index',axes=F,lwd=0.25)
axis(1,lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
text(.1,.8,labels=expression(italic(R^2) == 0.76))
lines(seq(0.01,0.99,0.01),2*inv.logit(predict(lm.out,data.frame(q3.aaa=seq(.01,0.99,0.01))))-1,lty=2)
box(bty='l')
dev.off()

e<-read.delim('data/verily_unlinked_norepeat.eigenvec',header=F,sep=' ')[order(sm$Pop),][1:nrow(m),]
agg<-aggregate(cbind(V3,V4,V5)~m$Pop,FUN=mean,data=e)
pm$genoPC1<-agg[match(pm$Population,agg[,1]),2]
pm$genoPC2<-agg[match(pm$Population,agg[,1]),3]
pm$genoPC3<-agg[match(pm$Population,agg[,1]),4]
lm.out2<-step(lm(logit(prob)~genoPC1+genoPC2*genoPC3,data=pm))
cols<-pm$prefcol[match(m$Pop,pm$Population)]
cols[is.na(cols)]<-'grey'
summary(lm.out2)

pdf('PCA2.pdf',width=4,height=4,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01,bty='l')
plot(e[,3],-e[,4],bg=cols2,pch=21,xlab='PC1',ylab='PC2',axes=F,lwd=0.25,cex=1.2)
axis(1,lwd=0.5)
axis(2,lwd=0.5)
dev.off()

##########################

thetas<-list()
for(i in list.files('thetas/',pattern='.*pestPG')){
	print(i)
	curr<-fread(paste('thetas/',i,sep=''))
	curr$plotpos<-curr$WinCenter
	curr$plotpos[curr$Chr=='NC_035108.1']=curr$plotpos[curr$Chr=='NC_035108.1']+380000000
	curr$plotpos[curr$Chr=='NC_035109.1']=curr$plotpos[curr$Chr=='NC_035109.1']+920000000
	thetas[[gsub('.pestPG','',i)]]<-curr
}
names(thetas)[names(thetas)=='Outgroup']<-'MASC'
pis<-sapply(thetas,function(t) sum(t$tP)/sum(t$nSites))
allpis<-(sapply(thetas,function(t) rollsum(t$tP,k=50,na.pad=T)/rollsum(t$nSites,k=50,na.pad=T)))
meanpis<-rowMeans(allpis)
nSites<-rollmean(rowMeans(sapply(thetas,'[[',14)),k=50,na.pad=T)

pdf(file='pi_genome.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(thetas[[1]]$plotpos,meanpis,cex=0.1,axes=F,xlab='Position (Mb)',ylab=expression(pi),ylim=c(-0.002,0.06))
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,0.06,0.03))
points(152e6,.045,lwd=2,lty=2,pch=6)
text(152e6,.052,lwd=2,lty=2,pch=6,labels='M locus',font=3,cex=0.8)
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
dev.off()

######################################################

fstread3<-function(pair,matchpos=NULL,k=50){
	curr<-fread(paste('FST/',pair,'.fst.components.txt',sep=''),header=F)
	curr$plotpos<-curr[,2]+50000
	curr$plotpos[curr[,1]=='NC_035108.1']=curr$plotpos[curr[,1]=='NC_035108.1']+380000000
	curr$plotpos[curr[,1]=='NC_035109.1']=curr$plotpos[curr[,1]=='NC_035109.1']+920000000
	colnames(curr)<-c('CHROM','START','END','num','den','nSites','plotpos')
	curr$FST<-rollmean(curr$num,k=k,na.pad=T)/rollmean(curr$den,k=k,na.pad=T)
	if(!is.null(matchpos)){
		curr<-curr[match(matchpos,curr$plotpos),]		
	}
	return(curr)
}

matchpos=thetas[[1]]$plotpos
BTT.BKK<-fstread3('BTT.BKK',matchpos)
BTT.NGO<-fstread3('BTT.NGO',matchpos)
BTT.THI<-fstread3('BTT.THI',matchpos)
BTT.OGD<-fstread3('BTT.OGD',matchpos)
BTT.OHI<-fstread3('BTT.OHI',matchpos)
BTT.MIN<-fstread3('BTT.MIN',matchpos)
BTT.KED<-fstread3('BTT.KED',matchpos)
MIN.KED<-fstread3('MIN.KED',matchpos)
MIN.OHI<-fstread3('MIN.OHI',matchpos)
MIN.OGD<-fstread3('MIN.OGD',matchpos)
MIN.THI<-fstread3('MIN.THI',matchpos)
MIN.NGO<-fstread3('MIN.NGO',matchpos)
BTT.FCV<-fstread3('BTT.FCV',matchpos)
FCV.BKK<-fstread3('FCV.BKK',matchpos)
FCV.NGO<-fstread3('FCV.NGO',matchpos)
FCV.THI<-fstread3('FCV.THI',matchpos)
FCV.OGD<-fstread3('FCV.OGD',matchpos)
SHM.NGO<-fstread3('NGO.SHM',matchpos)
SHM.THI<-fstread3('THI.SHM',matchpos)
SHM.OGD<-fstread3('OGD.SHM',matchpos)
BTT.SHM<-fstread3('BTT.SHM',matchpos)
FCV.SHM<-fstread3('FCV.SHM',matchpos)
FCV.KED<-fstread3('FCV.KED',matchpos)
FCV.OHI<-fstread3('FCV.OHI',matchpos)


BTT.BKK.pi<-((thetas[['BTT']]$tP/thetas[['BTT']]$nSites)+(thetas[['BKK']]$tP/thetas[['BKK']]$nSites))/2
BTT.NGO.pi<-((thetas[['BTT']]$tP/thetas[['BTT']]$nSites)+(thetas[['NGO']]$tP/thetas[['NGO']]$nSites))/2
BTT.THI.pi<-((thetas[['BTT']]$tP/thetas[['BTT']]$nSites)+(thetas[['THI']]$tP/thetas[['THI']]$nSites))/2
BTT.OGD.pi<-((thetas[['BTT']]$tP/thetas[['BTT']]$nSites)+(thetas[['OGD']]$tP/thetas[['OGD']]$nSites))/2

getPerm<-function(df1,df2,df3,pi,reps=10,k=10){
	perms<-replicate(reps,{
		shuf<-shuffle(length(pi),control=how(blocks=cut(pi,quantile(pi,seq(0,1,0.1),na.rm=T))))
		getPBS(df1[shuf,],df2[shuf,],df3[shuf,],k=k)
	})
	return(perms)
}

getPBS<-function(focal.ref,focal.out,ref.out,k=10){
	T1=-log(1-rollmean(focal.ref$num,k=k,na.pad=T)/rollmean(focal.ref$den,k=k,na.pad=T))
	T2=-log(1-rollmean(focal.out$num,k=k,na.pad=T)/rollmean(focal.out$den,k=k,na.pad=T))
	T3=-log(1-rollmean(ref.out$num,k=k,na.pad=T)/rollmean(ref.out$den,k=k,na.pad=T))
	PBS=(T1+T2-T3)/2
	return(PBS)
}

k=50
ngoPBS<-getPBS(BTT.NGO,FCV.NGO,BTT.FCV,k=k)
thiPBS<-getPBS(BTT.THI,FCV.THI,BTT.FCV,k=k)
ogdPBS<-getPBS(BTT.OGD,FCV.OGD,BTT.FCV,k=k)
bkkPBS<-getPBS(BTT.BKK,FCV.BKK,BTT.FCV,k=k)
bttPBS<-getPBS(BTT.OGD, BTT.FCV,FCV.OGD,k=k)
kedPBS<-getPBS(BTT.KED,FCV.KED,BTT.FCV,k=k)
ohiPBS<-getPBS(BTT.OHI,FCV.OHI,BTT.FCV,k=k)

ngoPBS2<-getPBS(BTT.NGO,SHM.NGO,BTT.SHM,k=k)
thiPBS2<-getPBS(BTT.THI,SHM.THI,BTT.SHM,k=k)
ogdPBS2<-getPBS(BTT.OGD,SHM.OGD,BTT.SHM,k=k)
ngoPBS3<-getPBS(FCV.NGO,SHM.NGO,FCV.SHM,k=k)
thiPBS3<-getPBS(FCV.THI,SHM.THI,FCV.SHM,k=k)
ogdPBS3<-getPBS(FCV.OGD,SHM.OGD,FCV.SHM,k=k)

# permNGO<-getPerm(BTT.NGO,FCV.NGO,BTT.FCV,BTT.NGO.pi,k=k)
# permTHI<-getPerm(BTT.THI,FCV.THI,BTT.FCV,BTT.THI.pi,k=k)
# permOGD<-getPerm(BTT.OGD,FCV.OGD,BTT.FCV,BTT.THI.pi,k=k)
# save(permNGO,permTHI,permOGD,file='perms.Rdata')
load('perms.Rdata')
eNGO<-ecdf(permNGO)
eTHI<-ecdf(permTHI)
eOGD<-ecdf(permOGD)

BTT.NGO$p<-2*pmin((1-eNGO(ngoPBS)),eNGO(ngoPBS))
BTT.NGO$padj<-p.adjust(BTT.NGO$p,method='fdr')
BTT.NGO$outl<-BTT.NGO$padj<0.05&ngoPBS>quantile(eNGO,0.5)
BTT.NGO$outld<-BTT.NGO$padj<0.05&ngoPBS<quantile(eNGO,0.5)

BTT.THI$p<-2*pmin((1-eTHI(thiPBS)),eTHI(thiPBS))
BTT.THI$padj<-p.adjust(BTT.THI$p,method='fdr')
BTT.THI$outl<-BTT.THI$padj<0.05&thiPBS>quantile(eTHI,0.5)
BTT.THI$outld<-BTT.THI$padj<0.05&thiPBS<quantile(eTHI,0.5)

BTT.OGD$p<-2*pmin((1-eOGD(ogdPBS)),eOGD(ogdPBS))
BTT.OGD$padj<-p.adjust(BTT.OGD$p,method='fdr')
BTT.OGD$outl<-BTT.OGD$padj<0.05&ogdPBS>quantile(eOGD,0.5)
BTT.OGD$outld<-BTT.OGD$padj<0.05&ogdPBS<quantile(eOGD,0.5)

r<-rle(BTT.NGO$outl&BTT.THI$outl)
rects2<-(cbind(BTT.NGO$plotpos[(cumsum(r$lengths)-r$lengths)[which(r$values)]]-2500000,rep(0,length(r$length[which(r$values)])), BTT.NGO$plotpos[cumsum(r$lengths)[which(r$values)]]+2500000,rep(0.75,length(r$length[which(r$values)]))))
r<-rle(BTT.NGO$outl&BTT.THI$outl&BTT.OGD$outl)
rects3<-(cbind(BTT.NGO$plotpos[(cumsum(r$lengths)-r$lengths)[which(r$values)]]-2500000,rep(0,length(r$length[which(r$values)])), BTT.NGO$plotpos[cumsum(r$lengths)[which(r$values)]]+2500000,rep(0.75,length(r$length[which(r$values)]))))

chr1<-which(BTT.BKK$CHROM=='NC_035107.1')
chr2<-which(BTT.BKK$CHROM=='NC_035108.1')
chr3<-which(BTT.BKK$CHROM=='NC_035109.1')

pdf(file='FST_pairwise.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,BTT.BKK$FST,cex=0.1,ylim=c(0,0.8),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(F[ST]~vs.~BTT),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos,BTT.BKK$FST,cex=0.1,col='darkred')
points(BTT.BKK$plotpos,BTT.NGO$FST,cex=0.1,col='darkorange')
points(BTT.BKK$plotpos,BTT.THI$FST,cex=0.1,col='darkcyan')
points(BTT.BKK$plotpos,BTT.OGD$FST,cex=0.1,col='darkmagenta')
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,0.6,0.3))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
legend('topright',fill=c('darkblue','darkorange','darkcyan','darkmagenta'),legend=c('BKK','NGO','THI','OGD'),ncol=4,cex=0.7,bty='n')
dev.off()

pdf(file='PBS_pairwise.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,ngoPBS,cex=0.1,ylim=c(0,.75),col='darkorange',axes=F, xlab='Position (Mb)',ylab='PBS',type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],bkkPBS[chr1],cex=0.1,col='darkgrey',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],ngoPBS[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],thiPBS[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],ogdPBS[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],bkkPBS[chr2],cex=0.1,col='darkgrey',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ngoPBS[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],thiPBS[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ogdPBS[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],bkkPBS[chr3],cex=0.1,col='darkgrey',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ngoPBS[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],thiPBS[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ogdPBS[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,.8,0.2))
legend('topright',fill=c('darkgrey','darkorange','darkcyan','darkmagenta'),legend=c('BKK','NGO','THI','OGD'),ncol=1,cex=0.65,bty='n')
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
dev.off()

pdf(file='PBS_pairwise_KED_OHI.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,ngoPBS,cex=0.1,ylim=c(0,.75),col='darkorange',axes=F, xlab='Position (Mb)',ylab='PBS',type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],kedPBS[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],ohiPBS[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],kedPBS[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ohiPBS[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],kedPBS[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ohiPBS[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,.8,0.2))
legend('topright',fill=c('darkcyan','darkorange'),legend=c('KED','OHI'),ncol=1,cex=0.65,bty='n')
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
dev.off()


pdf(file='PBS_pairwise_nobkk.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,ngoPBS,cex=0.1,ylim=c(0,.45),col='darkorange',axes=F, xlab='Position (Mb)',ylab='PBS',type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],ngoPBS[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],thiPBS[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],ogdPBS[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ngoPBS[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],thiPBS[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ogdPBS[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ngoPBS[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],thiPBS[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ogdPBS[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,.8,0.2))
legend('topright',fill=c('darkorange','darkcyan','darkmagenta'),legend=c('NGO','THI','OGD'),ncol=1,cex=0.7,bty='n')
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
dev.off()

pdf(file='PBS_pairwise_nobkk_BTT.SHM.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,ngoPBS,cex=0.1,ylim=c(0,.45),col='darkorange',axes=F, xlab='Position (Mb)',ylab='PBS',type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],ngoPBS2[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],thiPBS2[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],ogdPBS2[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ngoPBS2[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],thiPBS2[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ogdPBS2[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ngoPBS2[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],thiPBS2[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ogdPBS2[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,.8,0.2))
legend('topright',fill=c('darkorange','darkcyan','darkmagenta'),legend=c('NGO','THI','OGD'),ncol=1,cex=0.7,bty='n')
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
dev.off()

pdf(file='PBS_pairwise_nobkk_FCV.SHM.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,ngoPBS,cex=0.1,ylim=c(0,.45),col='darkorange',axes=F, xlab='Position (Mb)',ylab='PBS',type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],ngoPBS3[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],thiPBS3[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],ogdPBS3[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ngoPBS3[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],thiPBS3[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],ogdPBS3[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ngoPBS3[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],thiPBS3[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],ogdPBS3[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
axis(2,at=seq(0,.8,0.2))
legend('topright',fill=c('darkorange','darkcyan','darkmagenta'),legend=c('NGO','THI','OGD'),ncol=1,cex=0.7,bty='n')
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
dev.off()


pdf(file='FST_pairwise2.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,BTT.KED$FST,cex=0.1,ylim=c(0,0.8),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(F[ST]~vs.~BTT),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos,BTT.KED$FST,cex=0.1,col='darkblue')
points(BTT.BKK$plotpos,BTT.OHI$FST,cex=0.1,col='darkorange')
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.6,0.3))
legend('topright',fill=c('darkblue','darkorange'),legend=c('KED','OHI'),ncol=2,cex=0.7,bty='n')
dev.off()

pdf(file='FST_pairwise2_MIN.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,MIN.KED$FST,cex=0.1,ylim=c(0,0.75),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(F[ST]~vs.~MIN),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],BTT.MIN$FST[chr1],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],MIN.KED$FST[chr1],cex=0.1,col='darkgreen',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],MIN.OHI$FST[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],BTT.MIN$FST[chr2],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],MIN.KED$FST[chr2],cex=0.1,col='darkgreen',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],MIN.OHI$FST[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],BTT.MIN$FST[chr3],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],MIN.KED$FST[chr3],cex=0.1,col='darkgreen',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],MIN.OHI$FST[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.6,0.3))
legend('topright',fill=c('darkblue','darkorange','darkgreen'),legend=c('KED','OHI','BTT'),ncol=3,cex=0.7,bty='n')
dev.off()

pdf(file='FST_pairwise_MIN.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,MIN.KED$FST,cex=0.1,ylim=c(0,0.8),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(F[ST]~vs.~MIN),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1],MIN.NGO$FST[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],MIN.THI$FST[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1],MIN.OGD$FST[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],MIN.NGO$FST[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],MIN.THI$FST[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2],MIN.OGD$FST[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],MIN.NGO$FST[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],MIN.THI$FST[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3],MIN.OGD$FST[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.6,0.3))
legend('topright',fill=c('darkorange','darkcyan','darkmagenta'),legend=c('NGO','THI','OGD'),ncol=3,cex=0.7,bty='n')
dev.off()


##############################

dxyread<-function(file,fstmatch,k=50){
	curr<-fread(file)
	curr$plotpos<-curr[,2]+50000
	curr$plotpos[curr[,1]=='NC_035108.1']=curr$plotpos[curr[,1]=='NC_035108.1']+380000000
	curr$plotpos[curr[,1]=='NC_035109.1']=curr$plotpos[curr[,1]=='NC_035109.1']+920000000
	curr$dxy<-rollsum(curr[,4],k=k,na.pad=T)/rollsum(fstmatch$nSites[match(curr$plotpos,fstmatch$plotpos)],k=k,na.pad=T)
	curr<-curr[match(fstmatch$plotpos,curr$plotpos),]
	return(curr)
}

dxys<-list()
for(i in list.files('FST/',pattern='.*\\.dxy$')){
	dxys[[i]]<-dxyread(paste('FST/',i,sep=''),BTT.BKK)
}
extradxys<-list()
for(i in list.files('FST/extradxy/',pattern='.*\\.dxy$')){
	extradxys[[i]]<-dxyread(paste('FST/extradxy/',i,sep=''),BTT.BKK)
}

meandxy<-rowMeans(sapply(dxys,'[[',6))

pdf(file='pi_pops.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,allpis[,'BKK'],cex=0.1,ylim=c(0,0.06),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(pi),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1], allpis[chr1,'BKK'],cex=0.1,col='darkred',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'NGO'],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'THI'],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'OGD'],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'BTT'],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'BKK'],cex=0.1,col='darkred',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'NGO'],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'THI'],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'OGD'],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'BTT'],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'BKK'],cex=0.1,col='darkred',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'NGO'],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'THI'],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'OGD'],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'BTT'],cex=0.1,col='darkblue',type='l',lwd=2)

axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.05,0.025))
legend('topright',fill=c('darkred','darkorange','darkcyan','darkmagenta','darkblue'),legend=c('BKK','NGO','THI','OGD','BTT'),ncol=5,cex=0.7,bty='n')
dev.off()

pdf(file='dxy_pairwise.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,dxys[['BTT.BKK.dxy']]$dxy,cex=0.1,ylim=c(0,0.05),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(d[xy]~vs.~BTT),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1], dxys[['BTT.BKK.dxy']]$dxy[chr1],cex=0.1,col='darkred',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], dxys[['BTT.NGO.dxy']]$dxy[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], dxys[['BTT.THI.dxy']]$dxy[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], dxys[['BTT.OGD.dxy']]$dxy[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], extradxys[['BTT.OHI.dxy']]$dxy[chr1],cex=0.1,col='green',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], extradxys[['BTT.KED.dxy']]$dxy[chr1],cex=0.1,col='darkgreen',type='l',lwd=2)

points(BTT.BKK$plotpos[chr2], dxys[['BTT.BKK.dxy']]$dxy[chr2],cex=0.1,col='darkred',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], dxys[['BTT.NGO.dxy']]$dxy[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], dxys[['BTT.THI.dxy']]$dxy[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], dxys[['BTT.OGD.dxy']]$dxy[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], extradxys[['BTT.OHI.dxy']]$dxy[chr2],cex=0.1,col='green',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], extradxys[['BTT.KED.dxy']]$dxy[chr2],cex=0.1,col='darkgreen',type='l',lwd=2)

points(BTT.BKK$plotpos[chr3], dxys[['BTT.BKK.dxy']]$dxy[chr3],cex=0.1,col='darkred',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], dxys[['BTT.NGO.dxy']]$dxy[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], dxys[['BTT.THI.dxy']]$dxy[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], dxys[['BTT.OGD.dxy']]$dxy[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], extradxys[['BTT.OHI.dxy']]$dxy[chr3],cex=0.1,col='green',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], extradxys[['BTT.KED.dxy']]$dxy[chr3],cex=0.1,col='darkgreen',type='l',lwd=2)

axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.05,0.025))
legend('topright',fill=c('darkred','darkorange','darkcyan','darkmagenta','green','darkgreen'),legend=c('BKK','NGO','THI','OGD','OHI','KED'),ncol=6,cex=0.7,bty='n')
dev.off()

pdf(file='Supplementary Figures/dxy_pairwise_MIN.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,dxys[['BTT.BKK.dxy']]$dxy,cex=0.1,ylim=c(0,0.05),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(d[xy]~vs.~MIN),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1], extradxys[['BTT.MIN.dxy']]$dxy[chr1],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], extradxys[['MIN.KED.dxy']]$dxy[chr1],cex=0.1,col='blue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], extradxys[['MIN.OHI.dxy']]$dxy[chr1],cex=0.1,col='lightblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], extradxys[['BTT.MIN.dxy']]$dxy[chr2],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], extradxys[['MIN.KED.dxy']]$dxy[chr2],cex=0.1,col='blue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], extradxys[['MIN.OHI.dxy']]$dxy[chr2],cex=0.1,col='lightblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], extradxys[['BTT.MIN.dxy']]$dxy[chr3],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], extradxys[['MIN.KED.dxy']]$dxy[chr3],cex=0.1,col='blue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], extradxys[['MIN.OHI.dxy']]$dxy[chr3],cex=0.1,col='lightblue',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.05,0.025))
legend('topright',fill=c('darkblue','blue','lightblue'),legend=c('BTT','KED','OHI'),ncol=3,cex=0.7,bty='n')
dev.off()

pdf(file='Supplementary Figures/pi_pops_MIN.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(BTT.BKK$plotpos,allpis[,'BKK'],cex=0.1,ylim=c(0,0.06),col='darkblue',axes=F, xlab='Position (Mb)',ylab=expression(pi),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(BTT.BKK$plotpos[chr1], allpis[chr1,'MIN'],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'KED'],cex=0.1,col='blue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'OHI'],cex=0.1,col='lightblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr1], allpis[chr1,'BTT'],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'MIN'],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'KED'],cex=0.1,col='blue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'OHI'],cex=0.1,col='lightblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr2], allpis[chr2,'BTT'],cex=0.1,col='darkblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'MIN'],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'KED'],cex=0.1,col='blue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'OHI'],cex=0.1,col='lightblue',type='l',lwd=2)
points(BTT.BKK$plotpos[chr3], allpis[chr3,'BTT'],cex=0.1,col='darkblue',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='red',cex=0.5)
axis(2,at=seq(0,0.05,0.025))
legend('topright',fill=c('darkmagenta','darkblue','blue','lightblue'),legend=c('MIN','BTT','KED','OHI'),ncol=5,cex=0.7,bty='n')
dev.off()


#######

pdf(file='Supplementary Figures/dxy_pi_norm.pdf',width=8,height=3,useDingbats=F)
par(mfrow=c(3,4),mar=c(2,5,1,2),mgp=c(1,.1,0),tck=-0.01,bty='l',font.main=1)
cols=rep('grey',nrow(dxys[['BTT.BKK.dxy']]))
outls<-BTT.NGO$outl&BTT.THI$outl&BTT.OGD$outl
outls2<-BTT.NGO$outl&BTT.THI$outl
cols[outls2&thetas[[1]][,Chr]=='NC_035107.1']<-'darkred'
cols[outls2&thetas[[1]][,Chr]=='NC_035108.1']<-'darkblue'
cols[outls2&thetas[[1]][,Chr]=='NC_035109.1']<-'darkorange'
cols[outls&thetas[[1]][,Chr]=='NC_035107.1']<-'red'
cols[outls&thetas[[1]][,Chr]=='NC_035108.1']<-'blue'
cols[outls&thetas[[1]][,Chr]=='NC_035109.1']<-'orange'



##########

plot(allpis[,'BKK']/meanpis,dxys[['BTT.BKK.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab=expression(~Normalized~d[xy]),main='BKK vs BTT',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'BKK']/meanpis)[outls],(dxys[['BTT.BKK.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'BKK']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BTT.BKK.dxy']]$dxy/meandxy,na.rm=T),lty=2)
legend('bottomright',fill=c('red','blue','orange'),legend=c('Chr. 1','Chr. 2','Chr. 3'),bty='n',cex=0.8)

plot(allpis[,'NGO']/meanpis,dxys[['BTT.NGO.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab='',main='NGO vs BTT',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'NGO']/meanpis)[outls],(dxys[['BTT.NGO.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'NGO']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BTT.NGO.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'THI']/meanpis,dxys[['BTT.THI.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab='',main='THI vs BTT',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'THI']/meanpis)[outls],(dxys[['BTT.THI.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'THI']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BTT.THI.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'OGD']/meanpis,dxys[['BTT.OGD.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab='',main='OGD vs BTT',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'OGD']/meanpis)[outls],(dxys[['BTT.OGD.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'OGD']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BTT.OGD.dxy']]$dxy/meandxy,na.rm=T),lty=2)

###############

plot(allpis[,'BKK']/meanpis,dxys[['FCV.BKK.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab=expression(~Normalized~d[xy]),main='BKK vs FCV',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'BKK']/meanpis)[outls],(dxys[['FCV.BKK.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'BKK']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['FCV.BKK.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'NGO']/meanpis,dxys[['FCV.NGO.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab='',main='NGO vs FCV',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'NGO']/meanpis)[outls],(dxys[['FCV.NGO.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'NGO']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['FCV.NGO.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'THI']/meanpis,dxys[['FCV.THI.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab='',main='THI vs FCV',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'THI']/meanpis)[outls],(dxys[['FCV.THI.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'THI']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['FCV.THI.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'OGD']/meanpis,dxys[['FCV.OGD.dxy']]$dxy/meandxy,cex=0.1,xlab='',ylab='',main='OGD vs FCV',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'OGD']/meanpis)[outls],(dxys[['FCV.OGD.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'OGD']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['FCV.OGD.dxy']]$dxy/meandxy,na.rm=T),lty=2)

##########

plot(allpis[,'BTT']/meanpis,dxys[['BTT.FCV.dxy']]$dxy/meandxy,cex=0.1,xlab=expression(Normalized~pi),ylab=expression(~Normalized~d[xy]),main='FCV vs BTT',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'BTT']/meanpis)[outls],(dxys[['BTT.FCV.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'BTT']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BTT.FCV.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'NGO']/meanpis,dxys[['BKK.NGO.dxy']]$dxy/meandxy,cex=0.1,xlab=expression(Normalized~pi),ylab='',main='NGO vs BKK',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'NGO']/meanpis)[outls],(dxys[['BKK.NGO.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'NGO']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BKK.NGO.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'THI']/meanpis,dxys[['BKK.THI.dxy']]$dxy/meandxy,cex=0.1,xlab=expression(Normalized~pi),ylab='',main='THI vs BKK',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'THI']/meanpis)[outls],(dxys[['BKK.THI.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'THI']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BKK.THI.dxy']]$dxy/meandxy,na.rm=T),lty=2)

plot(allpis[,'OGD']/meanpis,dxys[['BKK.OGD.dxy']]$dxy/meandxy,cex=0.1,xlab=expression(Normalized~pi),ylab='',main='OGD vs BKK',col=cols,ylim=c(0.5,1.3),xlim=c(0.3,1.3))
points((allpis[,'OGD']/meanpis)[outls],(dxys[['BKK.OGD.dxy']]$dxy/meandxy)[outls],col=cols[outls],cex=0.3,pch=19)
abline(v=mean(allpis[,'OGD']/meanpis,na.rm=T),lty=2)
abline(h=mean(dxys[['BKK.OGD.dxy']]$dxy/meandxy,na.rm=T),lty=2)

dev.off()

#################################################################

am<-read.delim('data/africa_only_meta.txt',header=F,stringsAsFactors=F)
colnames(am)<-c('Sample','Pop')
am$Pop<-factor(am$Pop,levels=c("NGO","THI","MIN","KED","PKT","BTT","OGD","OHI","KUM","KIN","BOA","AWK","LBV","LPV","FCV","ENT","KAK","VMB","SHM","KWA","ABK","GND","KBO","RAB"))

load('pcadapt.Rdata')
sites<-fread('data/verily_maf05_africa.sites.gz')
sites$plotpos<-sites[,2]+c(0,380e6,920e6)[unlist(sites[,1])]
pdf(file='-- FINAL subplots/pcadapt.pdf',width=8,height=2,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(sites$plotpos[pvalues[,2]<0.01],-log10(pvalues[pvalues[,2]<0.01,2]),cex=0.1,axes=F,xlab='Position (Mb)',ylab=expression(-log[10]~p),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(sites$plotpos[pvalues[,2]<0.01],-log10(pvalues[pvalues[,2]<0.01,2]),cex=0.1)
padj<-p.adjust(pvalues[,2],method='bonf')
abline(h=-log10(0.05/length(pvalues[,2])),col='blue')
points(sites$plotpos[padj<0.05],-log10(pvalues[padj<0.05,2]),cex=0.1,col='red')
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(152e6,0,lwd=2,lty=2,pch=16,col='blue',cex=0.5)
points(380e6+230e6,0,lwd=2,lty=2,pch=16,col='blue',cex=0.5)
points(920e6+198e6,0,lwd=2,lty=2,pch=16,col='blue',cex=0.5)
axis(2)
dev.off()

lagg<-aggregate(scores~Pop,FUN=mean,data=am)
lagg$pref<-logit(pm$prob[match(lagg$Pop,pm$Population)])
pmm<-pm[match(lagg$Pop,pm$Population),]
summary(lm(pref~V2,data=lagg))
lm.out<-(lm(pref~V2,data=lagg))
pdf('-- FINAL subplots/pca_behavior.pdf',width=4,height=4,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(-lagg[,3],2*inv.logit(lagg$pref)-1,cex=1.5*pmm$densCex,bg=pmm$bio15col,pch=21,ylim=c(-1,1),xlab='-PC2',ylab='Preference index',axes=F,lwd=0.25)
axis(1,lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
text(.1,.8,labels=expression(italic(R^2) == 0.86))
lines(seq(-.1,.2,0.01),2*inv.logit(predict(lm.out,data.frame(V2=seq(.1,-.2,-0.01))))-1,lty=2)
dev.off()

cs<-as.data.frame(sites[padj<0.05,c(1,2,2)])
cs[,2]<-cs[,2]-1
cs[cs[,1]==1,1]<-'NC_035107.1'
cs[cs[,1]==2,1]<-'NC_035108.1'
cs[cs[,1]==3,1]<-'NC_035109.1'
write.table(cs,file='data/pcadapt_cand.bed',sep='\t',quote=F,row.names=F,col.names=F)

############################################


fdread<-function(trio,fstmatch){
	curr<-read.delim(paste('FST/',trio,'.Outgroup.fd',sep=''),header=F)
	curr$plotpos<-curr[,2]+50000
	curr$plotpos[curr[,1]=='NC_035108.1']=curr$plotpos[curr[,1]=='NC_035108.1']+380000000
	curr$plotpos[curr[,1]=='NC_035109.1']=curr$plotpos[curr[,1]=='NC_035109.1']+920000000
	curr$fd<-rollsum(curr[,4],k=k,na.pad=T)/rollsum(curr[,6],k=k,na.pad=T)
	curr$fd[curr$fd<0]<-0
	curr$fdwin<-curr[,4]/curr[,6]
	curr$fdwin[curr$fdwin<0]<-0
	curr<-curr[match(fstmatch$plotpos,curr$plotpos),]
	return(curr)
}
NGO.fd<-fdread('BTT.NGO.BKK',BTT.BKK)
THI.fd<-fdread('BTT.THI.BKK',BTT.BKK)
OGD.fd<-fdread('BTT.OGD.BKK',BTT.BKK)
MIN.fd<-fdread('BTT.MIN.BKK',BTT.BKK)


pdf(file='-- FINAL subplots/Fd.pdf',width=8,height=2.5,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(NGO.fd$plotpos,NGO.fd$fd,cex=0.1,ylim=c(-0.02,1),col='darkorange',axes=F, xlab='Position (Mb)',ylab=expression(italic(f[d])),type='n')
rect(rects2[,1],rects2[,2],rects2[,3],rects2[,4],border=F,col=grey(0.9))
rect(rects3[,1],rects3[,2],rects3[,3],rects3[,4],border=F,col=grey(0.75))
points(NGO.fd$plotpos[chr1],NGO.fd$fd[chr1],cex=0.1,col='darkorange',type='l',lwd=2)
points(THI.fd$plotpos[chr1],THI.fd$fd[chr1],cex=0.1,col='darkcyan',type='l',lwd=2)
points(OGD.fd$plotpos[chr1],OGD.fd$fd[chr1],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(NGO.fd$plotpos[chr2],NGO.fd$fd[chr2],cex=0.1,col='darkorange',type='l',lwd=2)
points(THI.fd$plotpos[chr2],THI.fd$fd[chr2],cex=0.1,col='darkcyan',type='l',lwd=2)
points(OGD.fd$plotpos[chr2],OGD.fd$fd[chr2],cex=0.1,col='darkmagenta',type='l',lwd=2)
points(NGO.fd$plotpos[chr3],NGO.fd$fd[chr3],cex=0.1,col='darkorange',type='l',lwd=2)
points(THI.fd$plotpos[chr3],THI.fd$fd[chr3],cex=0.1,col='darkcyan',type='l',lwd=2)
points(OGD.fd$plotpos[chr3],OGD.fd$fd[chr3],cex=0.1,col='darkmagenta',type='l',lwd=2)
axis(1,at=seq(0,300000000,100000000),labels=seq(0,300,100))
axis(1,at=seq(380000000,780000000,100000000),labels=seq(0,400,100))
axis(1,at=seq(920000000,1320000000,100000000),labels=seq(0,400,100))
points(281795001,.75,lwd=2,lty=2,pch=6)
text(281795001,.85,lwd=2,lty=2,pch=6,labels='Or4',font=3)
points(152e6,-0.02,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(380e6+230e6,-0.02,lwd=2,lty=2,pch=16,col='red',cex=0.5)
points(920e6+198e6,-0.02,lwd=2,lty=2,pch=16,col='red',cex=0.5)

axis(2,at=seq(0,0.6,0.3),cex.axis=0.8)
legend('topright',fill=c('darkorange','darkcyan','darkmagenta'),legend=c('NGO','THI','OGD'),ncol=3,cex=0.7,bty='n')

dev.off()



############################
library(plotrix)

sc<-read.csv('data/128_Threshold_Characters_pop.csv',stringsAsFactors=F)
lm.sc<-lm(logit(X.Area/100)~Pop,data=sc)
e<-emmeans(lm.sc,'Pop')
conf<-confint(e,type='response')
conf<-conf[match(pmg$Pop,conf$Pop),]
pmg<-pm[pm$Population%in%c('T51','ORL',levels(m$Pop)),]
pmg$percArea<-conf$response[match(pmg$Population,conf$Pop)]

lm.out<-(lm(logit(prob)~percArea,data=pmg))
summary(lm.out)

pdf('scaling_behavior_k6.pdf',width=3,height=3,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(pmg$percArea,pmg$pref,bg=pmg$k6col,pch=21,ylim=c(-1,1),xlab='Percent white scales',ylab='Preference index',axes=F,lwd=0.25,cex=2)
axis(1,at=c(.2,.4,.6),labels=c(20,40,60),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
text(20,.8,labels=expression(italic(R^2) == 0.81))
lines(seq(0,1,0.01),2*inv.logit(predict(lm.out,data.frame(percArea=seq(0,1,.01))))-1,lty=2)
box(bty='l')
dev.off()

pmg2<-pmg[pmg$Population%in%c(levels(m$Pop)),]

lm.out2<-(lm(logit(prob)~percArea,data=pmg2))
summary(lm.out2)

pmg3<-pmg[pmg$Population%in%m$Pop[m$Country!='Senegal'],]
lm.out3<-(lm(logit(prob)~percArea,data=pmg3))
summary(lm.out3)



pdf('scaling_behavior_k6_noaaa.pdf',width=3,height=3,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(pmg2$percArea,pmg2$pref,bg=pmg2$k6col,pch=21,ylim=c(-1,1),xlab='Percent white scales',ylab='Preference index',axes=F,lwd=0.25,cex=2)
axis(1,at=c(.2,.4,.6),labels=c(20,40,60),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
text(20,.8,labels=expression(italic(R^2) == 0.78))
lines(seq(0,1,0.01),2*inv.logit(predict(lm.out2,data.frame(percArea=seq(0,1,.01))))-1,lty=2)
box(bty='l')
dev.off()

pmg$k6col[pmg$Population%in%c('T51','ORL')]<-'grey'
pdf('scaling_behavior_k6_aaa_heldout.pdf',width=3,height=3,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),tck=-0.01)
plot(pmg$percArea,pmg$pref,bg=pmg$k6col,pch=21,ylim=c(-1,1),xlab='Percent white scales',ylab='Preference index',axes=F,lwd=0.25,cex=2)
axis(1,at=c(.2,.4,.6),labels=c(20,40,60),lwd=0.5)
axis(2,at=c(-1,0,1),lwd=0.5)
text(20,.8,labels=expression(italic(R^2) == 0.78))
lines(seq(0,1,0.01),2*inv.logit(predict(lm.out2,data.frame(percArea=seq(0,1,.01))))-1,lty=2)
box(bty='l')
dev.off()


pdf(width=8,height=1.5,file='Scaling.pdf',useDingbats=F)
par(mar=c(2,3,1,1),tck=-0.01,mgp=c(1.2,.2,0))
plot(conf$response,pch=21,bg=pmg$prefcol,bty='l',ylim=c(0,0.85),axes=F,xlab='',ylab='% white scales',lwd=0.25,cex=1.5)
arrows(1:26,(conf$lower),1:26,(conf$upper),angle=90,code=3,length=0)
axis(1,las=2,at=1:26,labels=c(pmg$Population),lwd=0.5,cex.axis=0.7)
# axis(1,lwd=0,at=c(7.5,13.5,18,22,26.5,31.5,39),line=2,labels=c('Senegal','Burkina Faso','Ghana','Nigeria','Gabon','Uganda','Kenya'))
axis(2,lwd=0.5,at=c(0,0.8),labels=c(0,80))
axis(2,at=seq(0,0.75,0.05),labels=NA)
axis(2,at=c(0,0.06,0.29,0.65,0.82),labels=NA,tck=0.02)
dev.off()

###########################

indseq<-seq(25,nrow(NGO.fd),50)
fisher.test(NGO.fd$fd[indseq]>quantile(NGO.fd$fd[indseq],probs=0.9,na.rm=T),outls2[indseq])
fisher.test(THI.fd$fd[indseq]>quantile(THI.fd$fd[indseq],probs=0.9,na.rm=T),outls2[indseq])
fisher.test(OGD.fd$fd[indseq]>quantile(OGD.fd$fd[indseq],probs=0.9,na.rm=T),outls2[indseq])