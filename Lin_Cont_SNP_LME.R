## The below function is to run linear mixed effects models with a vector of numeric covariates as the dependent variable.  
## The covs here are SNPs coded as 0,1,2 for an additive model.
## This is only for a random intercept model, most likely in our situation the family ID number shared among siblings. 

total=read.csv("Z:/Jackie/skin tone/skin/skin tone/skin_total.csv", header=T)

library (nlme)
Lin_Cont_SNP_LME=function(data,outcome,snp,random){
	
	data=na.omit(data[,c(outcome, snp, random)])
	colnames(data)<-c("outcome", "snp", "random")
	
	results=rep(NA,13)
	g=lme(outcome~snp,random= ~1|random, data=data, method="ML", na.action=na.omit)
	g.red=lme(outcome~1, random=~1|random, data=data, method="ML",  na.action=na.omit)
	
	results[1]=n=dim(na.omit(data[,c("outcome","snp")]))[1]
	results[2]=LR=2*(summary(g)$logLik - summary(g.red)$logLik)
	results[3]=LR_P=1-pchisq((2*(summary(g)$logLik - summary(g.red)$logLik)),1,ncp=0)
	results[4:5]=betas=summary(g)$tTable[,1]
	results[6:7]=betase=summary(g)$tTable[,2]
	results[8:9]=betat=summary(g)$tTable[,4]
	results[10:11]=betap=summary(g)$tTable[,5]
	results[12]=AIC=summary(g)$AIC
	results[13]=BIC=summary(g)$BIC
	
	return(results)

}


## Call function
path="./skin tone/"
data=total
pooled_sex="M_F"
project="time test"
dataset="skin_total"
date="05082008"
outcome="sysbp"
random="random"    
snp=c("rs8179526","rs10177833","rs673","rs672")

temp.results=matrix(NA,length(snp),13)

for(i in 1:length(snp)){
	temp.results[i,]=Lin_Cont_SNP_LME(data,outcome,snp[i],random)

}
temp.results=cbind(snp,temp.results)
colnames(temp.results)=c("SNP","N","LR","LR_P","Int_E","B_SNP_E","Int_SE","B_SNP_SE","Int_t","B_SNP_t","Int_P","B_SNP_P","AIC","BIC")
output = paste(path, outcome," SNP LME Stat ",pooled_sex," ",project," ",dataset," ",date,".csv",sep="")
write.table(temp.results,output,row.names=F,sep=",")
