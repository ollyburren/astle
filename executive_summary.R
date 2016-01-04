library(data.table)

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")
DATA_DIR<-file.path(GRPATH,'astle/DATA')
script.dir <- file.path(GRPATH,'CHIGP/R')
## read in library functions
source(file.path(script.dir,'common.R'))

score.dir=file.path(DATA_DIR,'out/geneScore/')
score.thresh<-0.8

t.f<-function(x){
  top.hyp<-unlist(sub("_gene_score","",names(which.max(x))))
  s.x<-sort(x,decreasing = TRUE)
  top.score<-signif(s.x[1],digits=4)
  support<-signif(s.x[1]/s.x[2],digits=4)
  #data.frame(hypothesis=top.hyp,top.score=top.score,support=support)
  data.frame(top.tissue.hypothesis=top.hyp,top.score=top.score,tissue.support=support)
}


exec.summary<-function(f){
  message(paste("Processing",f))
  scores<-fread(f,stringsAsFactors = FALSE)
  scores<-subset(scores,all_gene_score>=score.thresh)
  ## for each left which is the top scoring type and how does this compare to the next best scoring ?

  h.classes<-names(scores)[7:length(names(scores))]
  h.classes<-h.classes[h.classes!='all_gene_score']

  ft<-scores[,h.classes,with=FALSE]

  ## do the data table way when working

  ft<-as.matrix(ft)
  rownames(ft)<-scores$ensg

  res<-do.call("rbind",apply(ft,1,t.f))
  res$ensg<-scores$ensg
  res<-data.table(res[order(res$top.score,res$tissue.support,decreasing = TRUE),])
  setkey(res,ensg)

  out<-scores[,c(1:6,24:26),with=FALSE]
  setkey(out,ensg)

  out<-out[res]
  out<-out[order(out$all_gene_score,out$top.score,decreasing=TRUE),]
  out
}

gs.f<-list.files(path=score.dir,pattern="*.tab",full.names = TRUE)

gs<-lapply(gs.f,exec.summary)

out.dir<-file.path(DATA_DIR,'out','gs_summary')

dev.null<-lapply(gs,function(o){
  fname<-paste(unique(o$disease),'tab',sep=".")
  write.table(o,file=file.path(out.dir,fname),row.names = FALSE,quote=FALSE,sep="\t")
})

