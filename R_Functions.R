# Name: FreqHW
# Function definition: FreqHW statistics about allelic frequency, genotype frequency, and Hardy-Weinberg equilibrium
# Arguments: data set, SNP names, output file name 
  
# Results: 	Column label    		Description
	     	#SNP				SNP name	
		#Signal				Percent of subjects with genotype
		#Heterozygosity			2pq	
		#i				Major allele
		#j				Minor allele
		#Ni				Number of major alleles
		#Nj				Number of minor alleles
		#RelFreqi			Frequency of major allele
		#RelFreqj			Frequency of minor allele
		#Nii				Number of subjects who are homozygous for major allele
		#Nij				Number of heterozygous subjects
		#Njj				Number of subjects who are homozygous for minor allele
		#RelFreqii			Frequency of major homozygotes
		#RelFreqij			Frequency of heterozygotes
		#RelFreqjj			Frequency of minor homozygotes
		#HWEChiP			P-value of Hardy-Weinberg equilibrium

# Must have genetics library installed to run program (To install select Packages -> Install package(s) from CRAN -> genetics)

library(genetics)
#######################################
#definition of function FreqHW. 
#Do not change
######################################
FreqHW=function(test.data,snps,output){

  len=length(snps)
  results=matrix(NA,len,17)
 
 for (i in 1:len)
 {
 
	temp=genotype(test.data[,snps[i]],sep="")
	results[i,1]=signal=summary(temp)$n.typed/summary(temp)$n.total
	results[i,2]=heterozygosity=summary(temp)$Hu
	results[i,3]=major.allele.name=summary(temp)$allele.names[1]
	results[i,4]=major.allele.name=summary(temp)$allele.names[2]
	
	if(summary(temp)$nallele==2){
	 
		results[i,5]=major.allele.count=HWE.exact(temp)$par[1]
		results[i,6]=minor.allele.count=HWE.exact(temp)$par[2]
		results[i,7]=major.allele.freq=major.allele.count/(major.allele.count+minor.allele.count)
		results[i,8]=minor.allele.freq=minor.allele.count/(major.allele.count+minor.allele.count)

		results[i,9]=major.homozygotes.count=HWE.exact(temp)$stat[1]
		results[i,10]=heterozygotes.count=HWE.exact(temp)$stat[2]
		results[i,11]=minor.homozygotes.count=HWE.exact(temp)$stat[3]
		results[i,12]=major.homozygotes.freq=major.homozygotes.count/sum(HWE.exact(temp)$stat)
		results[i,13]=heterozygotes.freq=heterozygotes.count/sum(HWE.exact(temp)$stat)
		results[i,14]=minor.homozygotes.freq=minor.homozygotes.count/sum(HWE.exact(temp)$stat)

		 if(min(HWE.exact(temp)$stat)<=5){
			results[i,16]=HWE.exact(temp)$p.value
			results[i,17]=HWE.exact(temp)$p.value

		 }
		else{
		 	results[i,15]=HWE.chisq(temp)$p.value
		 	results[i,17]=results[i,15]
		 }
	} #end of "if(summary(temp)$nallele==2)"
	
	else{
		results[i,5]=major.allele.count=summary(temp)$allele.freq[1,1]
		results[i,6]=minor.allele.count=0
		results[i,7]=major.allele.freq=1
		results[i,8]=minor.allele.freq=0
		
		results[i,9]=major.homozygotes.count=summary(temp)$genotype.freq[1,1]
		results[i,10]=heterozygotes.count=0
		results[i,11]=minor.homozygotes.count=0
		results[i,12]=major.homozygotes.freq=1
		results[i,13]=heterozygotes.freq=0
		results[i,14]=minor.homozygotes.freq=0
	}	
	
}
results=cbind(snps, results)

colnames(results)=c("SNP","Signal","Heterozygosity","i","j", "Ni","Nj","RelFreqi","RelFreqj","Nii","Nij","Njj","RelFreqii","RelFreqij","RelFreqjj","HWEChiP","FisherP","TestP")

write.table(results,output,row.names=F,sep=",")
}

FreqHW(ds1,c("rs4646116","rs740387"),"c:/temp/try.csv")



library(gee)
#########################################################################
# Function: GEE_Cont_SNP_Univ(gee(outcome~snp[i]))
# outcome: continous outcome varible
# snp[i]: Genotypes of snp[i], has 3 or 2 levels.
###06-13-2007 Add error handler to the function, For some snps the gee got error: , Can not find the #########reason#### 
###
library(gee)
GEE_Cont_SNP_Univ=function(data,outcome,snps, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(snps)
	temp.results=matrix(NA,len,22)
	for (i in 1:len){
		results=rep(NA,22)
		if (sum(table(data[,snps[i]])>0)==2 & min(table(data[,snps[i]])[table(data[,snps[i]])>0])>1 ) {

			a=try(gee(data[,outcome]~data[,snps[i]],id=pedid, data=data, family="quasi", corstr="exchangeable"),TRUE)
			if (attributes(a)[1]!="try-error"){
				g=gee(data[,outcome]~data[,snps[i]],id=pedid, data=data, family="quasi", corstr="exchangeable")
				results[1:length(table(data[,snps[i]]))]= dimnames(table(data[,snps[i]]))[[1]]
				results[4]=dim(na.omit(data[,c(outcome,snps[i])]))[1]
				results[5:6]=betas=summary(g)$coef[,1]
				results[8:9]=betase=summary(g)$coef[,2]
				results[11:12]=naivez=summary(g)$coef[,3]
				results[14:15]=robustse=summary(g)$coef[,4]
				results[17:18]=robustz=summary(g)$coef[,5]
				results[20:21]=2*(1-pnorm(abs(as.numeric(results[17:18]))))
			}
		}
		else if (sum(table(data[,snps[i]])>0)==3) {


			a=try(gee(data[,outcome]~data[,snps[i]],id=pedid, data=data, family="quasi", corstr="exchangeable"),TRUE)
			if (attributes(a)[1]!="try-error"){
				g=gee(data[,outcome]~data[,snps[i]],id=pedid, data=data, family="quasi", corstr="exchangeable")
				results[1:3]= dimnames(table(data[,snps[i]]))[[1]]
				results[4]=dim(na.omit(data[,c(outcome,snps[i])]))[1]
				results[5:7]=betas=summary(g)$coef[,1]
				results[8:10]=betase=summary(g)$coef[,2]
				results[11:13]=naivez=summary(g)$coef[,3]
				results[14:16]=robustse=summary(g)$coef[,4]
				results[17:19]=robustz=summary(g)$coef[,5]
				results[20:22]=2*(1-pnorm(abs(as.numeric(results[17:19]))))
			}
		}
				
		temp.results[i,]=results
		print(i)
	}
	
	temp.results=cbind(snps,temp.results)
	colnames(temp.results)=c("SNP","SNPii","SNPij","SNPjj","N","Int_E","B_SNPij_E","B_SNPjj_E","Naive_Int_SE","Naive_B_SNPij_SE","Naive_B_SNPjj_SE","Naive_Int_Z","Naive_B_SNPij_Z","Naive_B_SNPjj_Z","Robust_Int_SE","Robust_B_SNPij_SE","Robust_B_SNPjj_SE","Robust_Int_Z","Robust_B_SNPij_Z","Robust_B_SNPjj_Z","Robust_Int_P","Robust_B_SNPij_P","Robust_B_SNPjj_P")
	write.table(temp.results,output,row.names=F,sep=",")

}

#######################################################################
# Function: GEE_Cat_SNP_Univ(gee(outcome~snp[i]))
# outcome: catigorical outcome varible
# snp[i]: Genotypes of snp[i], has 3 or 2 levels.

GEE_Cat_SNP_Univ=function(data,outcome,snps,output){
	data=data[is.na(data[,outcome])==F,]
	len=length(snps)
	temp.results=matrix(NA,len,22)
	for (i in 1:len){
		results=rep(NA,22)
		if (sum(table(data[,snps[i]])>0)==2 & min(table(data[,outcome],data[,snps[i]]))>0 ) {
			if(sum(is.na(glm(data[,outcome]~data[,snps[i]],family=binomial,data=data,na.action=na.omit)$coef))==0){
				g=gee(data[,outcome]~data[,snps[i]],id=netid, data=data, family=binomial, corstr="exchangeable")
				results[1:length(table(data[,snps[i]]))]= dimnames(table(data[,snps[i]]))[[1]]
				results[4]=dim(na.omit(data[,c(outcome,snps[i])]))[1]
				results[5:6]=betas=summary(g)$coef[,1]
				results[8:9]=betase=summary(g)$coef[,2]
				results[11:12]=naivez=summary(g)$coef[,3]
				results[14:15]=robustse=summary(g)$coef[,4]
				results[17:18]=robustz=summary(g)$coef[,5]
				results[20:21]=2*(1-pnorm(abs(as.numeric(results[17:18]))))
			}
		}
		else if (sum(table(data[,snps[i]])>0)==3) {
			if(sum(is.na(glm(data[,outcome]~data[,snps[i]],family=binomial, data=data,na.action=na.omit)$coef))==0){

				g=gee(data[,outcome]~data[,snps[i]],id=netid, data=data, family=binomial, corstr="exchangeable")
				results[1:3]= dimnames(table(data[,snps[i]]))[[1]]
				results[4]=dim(na.omit(data[,c(outcome,snps[i])]))[1]
				results[5:7]=betas=summary(g)$coef[,1]
				results[8:10]=betase=summary(g)$coef[,2]
				results[11:13]=naivez=summary(g)$coef[,3]
				results[14:16]=robustse=summary(g)$coef[,4]
				results[17:19]=robustz=summary(g)$coef[,5]
				results[20:22]=2*(1-pnorm(abs(as.numeric(results[17:19]))))
			}
		}
				
		temp.results[i,]=results
	}
	temp.results=cbind(snps,temp.results)
	colnames(temp.results)=c("SNP","SNPii","SNPij","SNPjj","N","Int_E","B_SNPij_E","B_SNPjj_E","Naive_Int_SE","Naive_B_SNPij_SE","Naive_B_SNPjj_SE","Naive_Int_Z","Naive_B_SNPij_Z","Naive_B_SNPjj_Z","Robust_Int_SE","Robust_B_SNPij_SE","Robust_B_SNPjj_SE","Robust_Int_Z","Robust_B_SNPij_Z","Robust_B_SNPjj_Z","Robust_Int_P","Robust_B_SNPij_P","Robust_B_SNPjj_P")
	write.table(temp.results,output,row.names=F,sep=",")

}
#########################################################################
# Function: GEE_Cont_Cov_Univ(gee(outcome~covs[i]))
# outcome: continous outcome varible
# covs[i]: covariate

GEE_Cont_Cov_Univ=function(data,outcome,covs,output){
	data=data[is.na(data[,outcome])==F,]
	len=length(covs)
	temp.results=matrix(NA,len,13)
	for (i in 1:len){
		results=rep(NA,13)
		results=rep(NA,13)
		if (covs[i]!=outcome){
			g=gee(data[,outcome]~data[,covs[i]],id=netid, data=data, family=quasi, corstr="exchangeable")
			results[1]=dim(na.omit(data[,c(outcome,covs[i])]))[1]
			results[2:3]=betas=summary(g)$coef[,1]
			results[4:5]=betase=summary(g)$coef[,2]
			results[6:7]=naivez=summary(g)$coef[,3]
			results[8:9]=robustse=summary(g)$coef[,4]
			results[10:11]=robustz=summary(g)$coef[,5]
			results[12:13]=2*(1-pnorm(abs(as.numeric(summary(g)$coef[,5]))))
		}
		temp.results[i,]=results
	}
	temp=rep(outcome,len)
	temp.results=cbind(temp,covs,temp.results)
	colnames(temp.results)=c("Outcome","Cov","N","Int_E","B_Cov_E","Naive_Int_SE","Naive_B_Cov_SE","Naive_Int_Z","Naive_B_Cov_Z","Robust_Int_SE","Robust_B_Cov_SE","Robust_Int_Z","Robust_B_Cov_Z","Robust_Int_P","Robust_B_Cov_P")
	write.table(temp.results,output,row.names=F,sep=",")

}

#########################################################################
# Function: GEE_Cat_Cov_Univ(gee(outcome~covs[i]))
# outcome: catigorical outcome varible
# covs[i]: covariate

GEE_Cat_Cov_Univ=function(data,outcome,covs, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(covs)
	temp.results=matrix(NA,len,13)
	for (i in 1:len){
		results=rep(NA,13)
		if (covs[i]!=outcome ){
			g=gee(data[,outcome]~data[,covs[i]],id=netid, data=data, family=binomial, corstr="exchangeable")
			results[1]=dim(na.omit(data[,c(outcome,covs[i])]))[1]
			results[2:3]=betas=summary(g)$coef[,1]
			results[4:5]=betase=summary(g)$coef[,2]
			results[6:7]=naivez=summary(g)$coef[,3]
			results[8:9]=robustse=summary(g)$coef[,4]
			results[10:11]=robustz=summary(g)$coef[,5]
			results[12:13]=2*(1-pnorm(abs(as.numeric(summary(g)$coef[,5]))))
		}
		temp.results[i,]=results
	}
	temp=rep(outcome,len)
	temp.results=cbind(temp,covs,temp.results)
	colnames(temp.results)=c("Outcome","Cov","N","Int_E","B_Cov_E","Naive_Int_SE","Naive_B_Cov_SE","Naive_Int_Z","Naive_B_Cov_Z","Robust_Int_SE","Robust_B_Cov_SE","Robust_Int_Z","Robust_B_Cov_Z","Robust_Int_P","Robust_B_Cov_P")
	write.table(temp.results,output,row.names=F,sep=",")
}

#########################################################################
# Function: GEE_SNP_Cov_Int(gee(outcome~covs[j]*snps[i]))
# outcome: continous outcome varible
# covs[j]: covariate 
# snps[i]: Genotypes of snp[i], has 3 or 2 levels.

GEE_SNP_Cov_Int=function(data,outcome,snps,covs, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(snps)
	len2=length(covs)
	colname=matrix(NA,1,38)
	colname[1,]=c("SNP","Cov","SNPii","SNPij","SNPjj","Nii","Nij","Njj","IntE","BCovE","BSNPijE","BSNPjjE","BCovSNPijE","BCovSNPjjE","NaiveIntSE","NaiveBCovSE","NaiveBSNPijSE","NaiveBSNPjjSE","NaiveBCovSNPijSE","NaiveBCovSNPijSE","NaiveIntZ","NaiveBCovZ","NaiveBSNPijZ","NaiveBSNPjjZ","NaiveBCovSNPijZ","NaiveBCovSNPjjZ","RobustIntSE","RobustBCovSE","RobustBSNPijSE","RobustBSNPjjSE","RobustBCovSNPijSE","RobustBCovSNPjjSE","RobustIntZ","RobustBCovZ","RobustBSNPijz","RobustBSNPjjZ","RobustBCovSNPijZ","RobustBCovSNPjjZ")
	write.table(colname,output,col.names=F,row.names=F,sep=",")

	for (i in 1:len){
		temp.results=matrix(NA,len2,38)
		for (j in 1:len2){
			temp.data=na.omit(data[,c(outcome, covs[j],snps[i])])
			results=rep(NA,38)
			results[1]=snps[i]
			results[2]=covs[j]
			if (sum(table(data[,snps[i]])>0)==2) {
				if(sum(is.na(lm(data[,outcome]~data[,covs[j]]*data[,snps[i]],data=data,na.action=na.omit,singular.ok=T)$coef))>0){}
				else{
					g=gee(data[,outcome]~data[,covs[j]]+data[,snps[i]]+data[,covs[j]]*data[,snps[i]],id=netid, data=data, family="quasi", corstr="exchangeable")
					results[6:(5+length(table(temp.data[,snps[i]])))]=summary(temp.data[,snps[i]])
					results[3:(2+length(table(temp.data[,snps[i]])))]= dimnames(table(temp.data[,snps[i]]))[[1]]
					results[9:11]=summary(g)$coef[1:3,1]
					results[13]=summary(g)$coef[4,1]
					results[15:17]=summary(g)$coef[1:3,2]
					results[19]=summary(g)$coef[4,2]
					results[21:23]=summary(g)$coef[1:3,3]
					results[25]=summary(g)$coef[4,3]
					results[27:29]=summary(g)$coef[1:3,4]
					results[31]=summary(g)$coef[4,4]
					results[32:35]=summary(g)$coef[1:3,5]
					results[37]=summary(g)$coef[4,5]

			
				}
			}
			else if (sum(table(data[,snps[i]])>0)==3) {
				if(sum(is.na(lm(data[,outcome]~data[,covs[j]]*data[,snps[i]],data=data,na.action=na.omit,singular.ok=T)$coef))>0){}
				else{
					g=gee(data[,outcome]~data[,covs[j]]+data[,snps[i]]+data[,covs[j]]*data[,snps[i]],id=netid, data=data, family="quasi", corstr="exchangeable")
					results[6:8]=summary(temp.data[,snps[i]])
					results[3:5]= dimnames(table(temp.data[,snps[i]]))[[1]]
					results[9:14]=betas=summary(g)$coef[,1]
					results[15:20]=betase=summary(g)$coef[,2]
					results[21:26]=naivez=summary(g)$coef[,3]
					results[27:32]=robustse=summary(g)$coef[,4]
					results[33:38]=robustz=summary(g)$coef[,5]

				}
			}

			temp.results[j,]=results
		}
		write.table(temp.results,output,row.names=F,col.names=F,sep=",",append=T)
		

	}
	
	
}

#########################################################################
# Function: GEE_Cont_Cov_poly
# Models: gee(outcome~covs[i]); gee(outcome~covs[i]+cos[i]^2); gee(outcome~covs[i]+cos[i]^2)covs[i]^3)
# outcome: continous outcome varible
# covs[i]: covariate 

GEE_Cont_Cov_poly=function(data,outcome,covs, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(covs)
	temp.results=matrix(NA,len*3,26)
	for (i in 1:len){
		results=rep(NA,26)
		g=gee(data[,outcome]~data[,covs[i]]+I(data[,covs[i]]^2)+I(data[,covs[i]]^3),id=netid, data=data, family="quasi", corstr="exchangeable")
		results[1]=paste(outcome,"~",covs[i],"+",covs[i],"^2+",covs[i],"^3",sep='')
		results[2]=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[3:6]=betas=summary(g)$coef[,1]
		results[7:10]=betase=summary(g)$coef[,2]
		results[11:14]=naivez=summary(g)$coef[,3]
		results[15:18]=robustse=summary(g)$coef[,4]
		results[19:22]=robustz=summary(g)$coef[,5]
		results[23:26]=2*(1-pnorm(abs(as.numeric(results[19:22]))))
		temp.results[(i*3),]=results
		
		results=rep(NA,26)
		g=gee(data[,outcome]~data[,covs[i]]+I(data[,covs[i]]^2),id=netid, data=data, family="quasi", corstr="exchangeable")
		results[1]=paste(outcome,"~",covs[i],"+",covs[i],"^2",sep='')
		results[2]=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[3:5]=betas=summary(g)$coef[,1]
		results[7:9]=betase=summary(g)$coef[,2]
		results[11:13]=naivez=summary(g)$coef[,3]
		results[15:17]=robustse=summary(g)$coef[,4]
		results[19:21]=robustz=summary(g)$coef[,5]
		results[23:25]=2*(1-pnorm(abs(as.numeric(results[19:21]))))
		temp.results[(i*3-1),]=results
		
		results=rep(NA,26)
		g=gee(data[,outcome]~data[,covs[i]],id=netid, data=data, family="quasi", corstr="exchangeable")
		results[1]=paste(outcome,"~",covs[i],sep='')
		results[2]=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[3:4]=betas=summary(g)$coef[,1]
		results[7:8]=betase=summary(g)$coef[,2]
		results[11:12]=naivez=summary(g)$coef[,3]
		results[15:16]=robustse=summary(g)$coef[,4]
		results[19:20]=robustz=summary(g)$coef[,5]
		results[23:24]=2*(1-pnorm(abs(as.numeric(results[19:20]))))
		temp.results[(i*3-2),]=results
	}
	
	colnames(temp.results)=c("Model","N","Int_E","B_Lin_E","B_Quad_E","B_Cubic_E","Naive_Int_SE","Naive_B_Lin_SE","Naive_B_Quad_SE","Naive_B_Cubic_SE","Naive_Int_Z","Naive_B_Lin_Z","Naive_B_Quad_Z","Naive_B_Cubic_Z","Robust_Int_SE","Robust_B_Lin_SE","Robust_B_Quad_SE","Robust_B_Cubic_SE","Robust_Int_Z","Robust_B_Lin_Z","Robust_B_Quad_Z","Robust_B_Cubic_Z","Robust_Int_P","Robust_B_Lin_P","Robust_B_Quad_P","Robust_B_Cubic_P")
	write.table(temp.results,output,row.names=F,sep=",")

}


#########################################################################
# Function: GEE_Cov_Pair(gee(outcome~covs[j]+covs[i]))
# outcome: continous outcome varible
# covs[j],[covs[j]: covariate 

GEE_Cov_Pair=function(data,outcome,covs,output){
	len=length(covs)
	temp.results=matrix(NA,(len*(len-1)/2),27)
	index=0
	for (i in 1:(len-1)){
		for (j in (i+1):len){
		results=rep(NA,27)
		results[1]=covs[i]
		results[2]=covs[j]
			g=gee(data[,outcome]~data[,covs[i]]*data[,covs[j]],id=netid, data=data, family=quasi, corstr="exchangeable")
			results[3]=dim(na.omit(data[,c(outcome,covs[i])]))[1]
			results[4:7]=betas=summary(g)$coef[,1]
			results[8:11]=betase=summary(g)$coef[,2]
			results[12:15]=naivez=summary(g)$coef[,3]
			results[16:19]=robustse=summary(g)$coef[,4]
			results[20:23]=robustz=summary(g)$coef[,5]
			results[24:27]=2*(1-pnorm(abs(as.numeric(summary(g)$coef[,5]))))
			index=index+1
			temp.results[index,]=results
		}		
	}
	
	colnames(temp.results)=c("Cov1","Cov2","N","Int_E","B_Cov1_E","B_Cov2_E","B_Cov1_2_E","Naive_Int_SE","Naive_B_Cov1_SE","Naive_B_Cov2_SE","Naive_B_Cov1_2_SE","Naive_Int_Z","Naive_B_Cov1_Z","Naive_B_Cov2_Z","Naive_B_Cov1_2_Z","Robust_Int_SE","Robust_B_Cov1_SE","Robust_B_Cov2_SE","Robust_B_Cov1_2_SE","Robust_Int_Z","Robust_B_Cov1_Z","Robust_B_Cov2_Z","Robust_B_Cov1_2_Z","Robust_Int_P","Robust_B_Cov1_P","Robust_B_Cov2_P","Robust_B_Cov1_2_P")

	write.table(temp.results,output,row.names=F,sep=",")
}

################################################################

Lin_Cont_SNP_Univ=function(data, outcome,snp, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(snp)
	temp.results=matrix(NA,len,19)
	for (i in 1:len){
		results=rep(NA,19)
		if (sum(table(data[,snp[i]])>0)==2 & min(table(data[,snp[i]])[table(data[,snp[i]])>0])>1) {
			g=lm(data[,outcome]~data[,snp[i]],data=data,na.action=na.omit, singular.ok=T)
			results[1]=substr(rownames(summary(g)$coef)[-1], nchar(rownames(summary(g)$coef)[-1])-1,nchar(rownames(summary(g)$coef)[-1]))
			results[3:4]=betas=summary(g)$coef[,1]
			results[6:7]=betase=summary(g)$coef[,2]
			results[9:10]=betap=summary(g)$coef[,4]
			results[12]=rsq=summary(g)$r.sq
			results[13]=f=summary(g)$fstat[1]
			results[14]=df1=summary(g)$fstat[2]
			results[15]=df2=summary(g)$fstat[3]
			results[16]=pmodel=1-pf(summary(g)$fstatistic[1], summary(g)$fstatistic[2], summary(g)$fstatistic[3], ncp=0)
			results[17]=sse=sum(g$res^2)
			results[18]=sst=sum((na.omit(data[,c(outcome,snp[i])])[,1]-mean(na.omit(data[,c(outcome,snp[i])])[,1]))^2)
			results[19]=n=dim(na.omit(data[,c(outcome,snp[i])]))[1]
		}
		else if (sum(table(data[,snp[i]])>0)==3) {

			g=lm(data[,outcome]~data[,snp[i]],na.action=na.omit, singular.ok=T)
			results[1:2]=substr(rownames(summary(g)$coef)[-1], nchar(rownames(summary(g)$coef)[2])-1,nchar(rownames(summary(g)$coef)[2]))
			results[3:5]=betas=summary(g)$coef[,1]
			results[6:8]=betase=summary(g)$coef[,2]
			results[9:11]=betap=summary(g)$coef[,4]
			results[12]=rsq=summary(g)$r.sq
			results[13]=f=summary(g)$fstat[1]
			results[14]=df1=summary(g)$fstat[2]
			results[15]=df2=summary(g)$fstat[3]
			results[16]=pmodel=1-pf(summary(g)$fstatistic[1], summary(g)$fstatistic[2], summary(g)$fstatistic[3], ncp=0)
			results[17]=sse=sum(g$res^2)
			results[18]=sst=sum((na.omit(data[,c(outcome,snp[i])])[,1]-mean(na.omit(data[,c(outcome,snp[i])])[,1]))^2)
			results[19]=n=dim(na.omit(data[,c(outcome,snp[i])]))[1]
		}

		temp.results[i,]=results
	}
	temp.results=cbind(snp, temp.results)
	colnames(temp.results)=c("SNP","SNPij","SNPjj","Int_E","B_SNPij_E","B_SNPjj_E","Int_SE","B_SNPij_SE","B_SNPjj_SE","Int_P","B_SNPij_P","B_SNPjj_P","Rsq","Fstat","Df1","Df2", "ModelP","SSE","SST","N")
	write.table(temp.results,output,row.names=F,sep=",")
}

Logit_SNP_Univ(ds1,"BPControlled",c("rs1801278","rs4762","rs2010501"),"c:/temp/try2.csv")

Logit_SNP_Univ=function(data, outcome,snp, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(snp)
	temp.results=matrix(NA,len,21)
	for (i in 1:len){
		results=rep(NA,21)
		if (sum(table(data[,snp[i]])>0)==2 & min(table(data[,snp[i]])[table(data[,snp[i]])>0])>1) {
			g=glm(data[,outcome]~data[,snp[i]]+exf_age+GENDER,data=data,na.action=na.omit, family=binomial)

			results[1]=substr(rownames(summary(g)$coef)[-1], nchar(rownames(summary(g)$coef)[-1])-1,nchar(rownames(summary(g)$coef)[-1]))
			results[3:6]=betas=summary(g)$coef[,1]
			results[8:11]=betase=summary(g)$coef[,2]
			results[13:16]=betap=summary(g)$coef[,4]
			results[18]=summary(g)$aic
			results[19]=logLik(g)
			results[20]=n=dim(na.omit(data[,c(outcome,snp[i])]))[1]
			results[21]=anova(g,test="Chisq")$P[2]
		}
		else if (sum(table(data[,snp[i]])>0)==3) {

			g=glm(data[,outcome]~data[,snp[i]]+exf_age+GENDER,data=data,na.action=na.omit, family=binomial)

			results[1:2]=substr(rownames(summary(g)$coef)[-1], nchar(rownames(summary(g)$coef)[2])-1,nchar(rownames(summary(g)$coef)[2]))
			results[3:7]=betas=summary(g)$coef[,1]
			results[8:12]=betase=summary(g)$coef[,2]
			results[13:17]=betap=summary(g)$coef[,4]
			results[18]=summary(g)$aic
			results[19]=logLik(g)
			results[20]=n=dim(na.omit(data[,c(outcome,snp[i])]))[1]
			results[21]=anova(g,test="Chisq")$P[2]
		}

		temp.results[i,]=results
	}
	temp.results=cbind(snp, temp.results)
	colnames(temp.results)=c("SNP","SNPij","SNPjj","Int_E","B_SNPij_E","B_SNPjj_E","B_exf_age_E","B_Gender_E","Int_SE","B_SNPij_SE","B_SNPjj_SE","B_exf_age_SE","B_Gender_SE","Int_P","B_SNPij_P","B_SNPjj_P","B_exf_age_P","B_Gender_P","AIC","NegLogLike","N","ModelP")
	write.table(temp.results,output,row.names=F,sep=",")
}

Logit_SNP_Univ(ds1,"BPControlled",c("rs1801278","rs635311","rs2010501"),"c:/temp/try8.csv")

Lin_Cont_Cov_Univ=function(data,outcome,covs, output){
	data=data[is.na(data[,outcome])==F,]
	len=length(covs)
	temp.results=matrix(NA,len,14)
	for (i in 1:len){
		results=rep(NA,14)
		g=lm(data[,outcome]~data[,covs[i]],na.action=na.omit, singular.ok=T)
		results[1]=n=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[2]=pmodel=1-pf(summary(g)$fstatistic[1], summary(g)$fstatistic[2], summary(g)$fstatistic[3], ncp=0)
		results[3:4]=betas=summary(g)$coef[,1]
		results[5:6]=betase=summary(g)$coef[,2]
		results[7:8]=betap=summary(g)$coef[,4]
		results[9]=rsq=summary(g)$r.sq
		results[10]=f=summary(g)$fstat[1]
		results[11]=df1=summary(g)$fstat[2]
		results[12]=df2=summary(g)$fstat[3]
		results[13]=sse=sum(g$res^2)
		results[14]=sst=sum((na.omit(data[,c(outcome,covs[i])])[,1]-mean(na.omit(data[,c(outcome,covs[i])])[,1]))^2)
				
	
		temp.results[i,]=results
	}
	
	temp.results=cbind(covs,temp.results)
	colnames(temp.results)=c("Cov","N","ModelP","Int_E","B_Cov_E","Int_SE","B_Cov_SE","Int_P","B_Cov_P","Rsq","Fstat","Df1","Df2", "SSE","SST")
	write.table(temp.results,output,row.names=F,sep=",")

}

Lin_Cont_Cov_Poly=function(data,outcome,covs, output){
	len=length(covs)
	temp.results=matrix(NA,len*3,21)
	for (i in 1:len){
		results=rep(NA,21)
		g=lm(data[,outcome]~data[,covs[i]]+I(data[,covs[i]]^2)+I(data[,covs[i]]^3),na.action=na.omit, singular.ok=T)
		results[1]=paste(outcome,"~",covs[i],"+",covs[i],"^2+",covs[i],"^3",sep='')
		results[21]=n=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[2]=pmodel=1-pf(summary(g)$fstatistic[1], summary(g)$fstatistic[2], summary(g)$fstatistic[3], ncp=0)
		results[3:6]=betas=summary(g)$coef[,1]
		results[7:10]=betase=summary(g)$coef[,2]
		results[11:14]=betap=summary(g)$coef[,4]
		results[15]=rsq=summary(g)$r.sq
		results[16]=f=summary(g)$fstat[1]
		results[17]=df1=summary(g)$fstat[2]
		results[18]=df2=summary(g)$fstat[3]
		results[19]=sse=sum(g$res^2)
		results[20]=sst=sum((na.omit(data[,c(outcome,covs[i])])[,1]-mean(na.omit(data[,c(outcome,covs[i])])[,1]))^2)
		temp.results[(i*3),]=results
		
		results=rep(NA,21)
		g=lm(data[,outcome]~data[,covs[i]]+I(data[,covs[i]]^2),na.action=na.omit, singular.ok=T)
		results[1]=paste(outcome,"~",covs[i],"+",covs[i],"^2",sep='')
		results[21]=n=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[2]=pmodel=1-pf(summary(g)$fstatistic[1], summary(g)$fstatistic[2], summary(g)$fstatistic[3], ncp=0)
		results[3:5]=betas=summary(g)$coef[,1]
		results[7:9]=betase=summary(g)$coef[,2]
		results[11:13]=betap=summary(g)$coef[,4]
		results[15]=rsq=summary(g)$r.sq
		results[16]=f=summary(g)$fstat[1]
		results[17]=df1=summary(g)$fstat[2]
		results[18]=df2=summary(g)$fstat[3]
		results[19]=sse=sum(g$res^2)
		results[20]=sst=sum((na.omit(data[,c(outcome,covs[i])])[,1]-mean(na.omit(data[,c(outcome,covs[i])])[,1]))^2)
		temp.results[(i*3-1),]=results
		
		results=rep(NA,21)
		g=lm(data[,outcome]~data[,covs[i]],na.action=na.omit, singular.ok=T)
		results[1]=paste(outcome,"~",covs[i],sep='')
		results[21]=n=dim(na.omit(data[,c(outcome,covs[i])]))[1]
		results[2]=pmodel=1-pf(summary(g)$fstatistic[1], summary(g)$fstatistic[2], summary(g)$fstatistic[3], ncp=0)
		results[3:4]=betas=summary(g)$coef[,1]
		results[7:8]=betase=summary(g)$coef[,2]
		results[11:12]=betap=summary(g)$coef[,4]
		results[15]=rsq=summary(g)$r.sq
		results[16]=f=summary(g)$fstat[1]
		results[17]=df1=summary(g)$fstat[2]
		results[18]=df2=summary(g)$fstat[3]
		results[19]=sse=sum(g$res^2)
		results[20]=sst=sum((na.omit(data[,c(outcome,covs[i])])[,1]-mean(na.omit(data[,c(outcome,covs[i])])[,1]))^2)
		temp.results[(i*3-2),]=results
	}
	
	colnames(temp.results)=c("Model","ModelP","Int_E","B_Lin_E","B_Quad_E","B_Cubic_E","Int_SE","B_Lin_SE","B_Quad_SE","B_Cubic_SE","Int_P","B_Lin_P","B_Quad_P","B_Cubic_P","Rsq","Fstat","Df1","Df2", "SSE","SST","N")
	write.table(temp.results,output,row.names=F,sep=",")

}

