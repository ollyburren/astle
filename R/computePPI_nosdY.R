library(data.table)
library(GenomicRanges)

pi_i<-10^-4

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")
DATA_DIR<-file.path(GRPATH,'CHIGP/DATA')
script.dir <- file.path(GRPATH,'CHIGP/R')
## read in library functions
source(file.path(script.dir,'common.R'))

logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

args<-list(
  W=0.15,
  block.file='/Volumes/stats/STATS/PROJECTS/FGWAS_IC/support/0.1cM.regions.new.bed',
  in.file='/Volumes/stats/WILL_ASTLE/version_2_raw/plt_common.tsv',
  out.file='/Volumes/stats/WILL_ASTLE/PPI_TAB/plt_common.ppi'
)

numerics=c('W')

if(!interactive())
  args<-getArgs(verbose=TRUE,numeric=numerics)

## exploratory data analysis - Will said that plt was the best to pilot
## the software on.
tmp<-fread(args[['in.file']],stringsAsFactors=FALSE)

## prior variance or should this be 0.15 for quant trait ?

## we can compute Z score by beta/standard error and
## variance of beta is standard error of beta squared

tmp<-tmp[,c('Z','V'):=list(Effect/StdErr,
	StdErr^2),by=SNP][]

## Will has told us that the sd of Y is 1 for all SNPs therefore we no longer need to estimate it.
sdY<-1
sd.prior<-args[['W']] * sdY
tmp$r<-with(tmp,sd.prior^2 / (sd.prior^2 + V))
tmp$lABF<-with(tmp,0.5 * (log(1-r) + (r * Z^2)))

snp.gr<-with(tmp,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,end=BP),id=SNP))

blocks<-fread(args[['block.file']],sep="\t",header=FALSE,stringsAsFactors=FALSE)
setnames(blocks,c('chr','start','end','id'))
blocks.gr<-with(blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end)))

ol<-as.matrix(findOverlaps(snp.gr,blocks.gr))
block<-integer(length=length(snp.gr))
block[ol[,1]]<-ol[,2]
tmp$block<-block

tmp<-tmp[,sBF:=logsum(lABF + log(pi_i)),by=block][]
tmp<-tmp[,ppi:=exp(lABF + log(pi_i))/(exp(sBF) + 1),by=block][]

## spit out in a format that the downstream software can recognise.

out<-with(tmp,data.table(chr=CHR,start=BP-1,end=BP,rsid=ID,maf=0,pval=signif(10^(minus_log10P*-1),digits=3),ppi=ppi,imp.snp.pos=0,imp.r2=0))
out<-out[order(out$chr,out$start),]
write.table(out,file=args[['out.file']],sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

