## This is an archived version of the software originally used. Not that previously the sdY was not normalised so it needed to be estimated.
## current script is in this dir  computePPI_nosdY.R

library(data.table)
library(GenomicRanges)

W=0.15
pi_i=10^-4
block_file='0.1cM.regions.new.bed'
in.file<-'mpv.tsv'
out.file<-'mpv.pmi'

sdY.est <- function(vbeta, maf, n) {
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover)
  return(sqrt(coef(m)[['oneover']]))
}

logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}



## exploratory data analysis - Will said that plt was the best to pilot
## the software on.
tmp<-fread(in.file,stringsAsFactors=FALSE)

## estimate the # of samples

N<-ceiling(max((with(head(tmp,n=1000),all_AA+all_AB+all_BB))))

## prior variance or should this be 0.15 for quant trait ?

## we can compute Z score by beta/standard error and
## variance of beta is standard error of beta squared

tmp<-tmp[,c('Z','V'):=list(frequentist_add_beta_1/frequentist_add_se_1,
	frequentist_add_se_1^2),by=alternate_ids][]
## compute the shrinkage ratio r (ratio of prior variance to total variance)

##r = W/W+V -- we need to estimate W the prior variance. Where we estimate 
## V using allele frequencies and sample sizes we scale V so that sdY=1
## thus we don't have to bother with it. However in this case we do have 
## to bother with it.

## if have standard error then we cannot assume the sdY is scaled so as to be 1
## we have to estimate it.

sdY<-with(tmp,sdY.est(V,all_maf,N))

## we can for ease assume that n is fixed for a trait
## but we should check to make sure.

## I think that we wish to estimate this over the whole dataset

sd.prior<-W * sdY
tmp$r<-with(tmp,sd.prior^2 / (sd.prior^2 + V))
tmp$lABF<-with(tmp,0.5 * (log(1-r) + (r * Z^2)))

snp.gr<-with(tmp,GRanges(seqnames=Rle(chromosome),ranges=IRanges(start=position,end=position),id=alternate_ids))

blocks<-fread(block_file,sep="\t",header=FALSE,stringsAsFactors=FALSE)
setnames(blocks,c('chr','start','end','id'))
blocks.gr<-with(blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end)))

ol<-as.matrix(findOverlaps(snp.gr,blocks.gr))
block<-integer(length=length(snp.gr))
block[ol[,1]]<-ol[,2]
tmp$block<-block

tmp<-tmp[,sBF:=logsum(lABF + log(pi_i)),by=block][]
tmp<-tmp[,ppi:=exp(lABF + log(pi_i))/(exp(sBF) + 1),by=block][]

## as a  test this should be positive i.e. smaller p gives larger ppi
## cor(-log(tmp$frequentist_add_pvalue),tmp$ppi)

## spit out in a format that the downstream software can recognise.

out<-with(tmp,data.table(chr=chromosome,start=position-1,end=position,id=rsid,r2=0,imp.snp=0,MAF=all_maf,pval=frequentist_add_pvalue,ppi=ppi,norm.minus.lp=0)))

out<-out[order(out$chr,out$start),]

write.table(out,file=out.file,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

## using blockshifter 
