#' Utils: Experimental sketch
#'
#' @description Experimental sketching function
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @name sketch
#' @param trat Vector with factor A levels
#' @param trat1 Vector with levels of factor B (Set to NULL if not factorial or psub)
#' @param trat2 Vector with levels of factor C (Set to NULL if not factorial)
#' @param r Number of repetitions
#' @param design Experimental design (see note)
#' @param pos Repeat position (line or column),
#' @param color.sep Color box
#' @param ID plot Add only identification in sketch
#' @param label.x text in x
#' @param label.y text in y
#' @param labelsize Label size
#' @param legendsize Title legend size
#' @param axissize Axis size
#' @param add.streets.x Adds streets by separating treatments in row or column. The user must supply a numeric vector grouping the rows or columns that must be together. See the example.
#' @param add.streets.y Adds streets by separating treatments in row or column. The user must supply a numeric vector grouping the rows or columns that must be together. See the example.
#' @param export.csv Save table template based on sketch in csv
#' @param comment.caption Add comment in caption
#' @importFrom utils write.csv
#' @keywords croqui
#' @keywords experimental
#' @return Returns an experimental sketch according to the specified design.
#' @note The sketches have only a rectangular shape, and the blocks (in the case of randomized blocks) can be in line or in a column.
#' @note For the design argument, you can choose from the following options:
#'   \describe{
#'   \item{\code{design="DIC"}}{Completely randomized design}
#'   \item{\code{design="DBC"}}{Randomized block design}
#'   \item{\code{design="DQL"}}{Latin square design}
#'   \item{\code{design="FAT2DIC"}}{DIC experiments in double factorial}
#'   \item{\code{design="FAT2DBC"}}{DBC experiments in double factorial}
#'   \item{\code{design="FAT3DIC"}}{DIC experiments in triple factorial}
#'   \item{\code{design="FAT3DBC"}}{DBC experiments in triple factorial}
#'   \item{\code{design="PSUBDIC"}}{DIC experiments in split-plot}
#'   \item{\code{design="PSUBDBC"}}{DBC experiments in split-plot}
#'   \item{\code{design="PSUBSUBDBC"}}{DBC experiments in split-split-plot}
#'   \item{\code{design="STRIP-PLOT"}}{Strip-plot DBC experiments}
#'   }
#' @note For the color.sep argument, you can choose from the following options:
#'   \describe{
#'   \item{\code{design="DIC"}}{use "all" or "none"}
#'   \item{\code{design="DBC"}}{use "all","bloco" or "none"}
#'   \item{\code{design="DQL"}}{use "all", "column", "line" or "none"}
#'   \item{\code{design="FAT2DIC"}}{use "all", "f1", "f2" or "none"}
#'   \item{\code{design="FAT2DBC"}}{use "all", "f1", "f2", "block" or "none"}
#'   \item{\code{design="FAT3DIC"}}{use "all", "f1", "f2", "f3" or "none"}
#'   \item{\code{design="FAT3DBC"}}{use "all", "f1", "f2", "f3", "block" or "none"}
#'   \item{\code{design="PSUBDIC"}}{use "all", "f1", "f2" or "none"}
#'   \item{\code{design="PSUBDBC"}}{use "all", "f1", "f2", "block" or "none"}
#'   \item{\code{design="PSUBSUBDBC"}}{use "all", "f1", "f2", "f3", "block" or "none"}
#'   }
#'
#' @references
#' Mendiburu, F., & de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#' @export
#' @examples
#' Trat=paste("Tr",1:6)
#'
#' #=============================
#' # Completely randomized design
#' #=============================
#' sketch(Trat,r=3)
#' sketch(Trat,r=3,pos="column")
#' sketch(Trat,r=3,color.sep="none")
#' sketch(Trat,r=3,color.sep="none",ID=TRUE)
#' sketch(Trat,r=3,pos="column",add.streets.x=c(1,1,2,2,3,3))
#'
#' #=============================
#' # Randomized block design
#' #=============================
#' sketch(Trat, r=3, design="DBC")
#' sketch(Trat, r=3, design="DBC",pos="column")
#' sketch(Trat, r=3, design="DBC",pos="column",add.streets.x=c(1,1,2))
#' #=============================
#' # Completely randomized experiments in double factorial
#' #=============================
#' sketch(trat=c("A","B"),
#'        trat1=c("A","B","C"),
#'        design = "FAT2DIC",
#'        r=3)
#'
#' sketch(trat=c("A","B"),
#'        trat1=c("A","B","C"),
#'        design = "FAT2DIC",
#'        r=3,
#'        pos="column")

sketch=function(trat,
                trat1=NULL,
                trat2=NULL,
                r,
                design="DIC",
                pos="line",
                color.sep="all",
                ID=FALSE,
                add.streets.y=NA,
                add.streets.x=NA,
                label.x="",
                label.y="",
                axissize=12,
                legendsize=12,
                labelsize=4,
                export.csv=FALSE,
                comment.caption=NULL){
  requireNamespace("ggplot2")
  design.crd <-function(trt,r,serie=2,seed=0,kinds="Super-Duper",randomization=TRUE)
  {
    number<-0
    if(serie>0) number<-10^serie
    junto<-data.frame(trt,r)
    junto<-junto[order(junto[,1]),]
    TR<-as.character(junto[,1])
    r<-as.numeric(junto[,2])
    y <- rep(TR[1], r[1])
    tr <- length(TR)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    parameters<-list(design="crd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
    for (i in 2:tr) y <- c(y, rep(TR[i], r[i]))
    trat<-y
    if(randomization)trat <- sample(y, length(y), replace = FALSE)
    plots <- number+1:length(trat)
    dca<-data.frame(plots, trat)
    dca[,1]<-as.numeric(dca[,1])
    xx<-dca[order(dca[,2],dca[,1]),]
    r1<-seq(1,r[1])
    for (i in 2:length(r)) {
      r1<-c(r1,seq(1,r[i]))
    }
    yy<-data.frame(plots=xx[,1],r=r1,xx[,2])
    book<-yy[order(yy[,1]),]
    rownames(book)<-rownames(yy)
    names(book)[3]<-c(paste(deparse(substitute(trt))))
    outdesign<-list(parameters=parameters,book=book)
    return(outdesign)
  }

  design.rcbd <-function (trt, r,serie=2,seed=0,
                          kinds="Super-Duper",first=TRUE,
                          continue=FALSE,randomization=TRUE){
    number<-10
    if(serie>0) number<-10^serie
    ntr <- length(trt)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    parameters<-list(design="rcbd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
    mtr <-trt
    if(randomization)mtr <- sample(trt, ntr, replace = FALSE)
    block <- c(rep(1, ntr))
    for (y in 2:r) {
      block <- c(block, rep(y, ntr))
      if(randomization)mtr <- c(mtr, sample(trt, ntr, replace = FALSE))
    }
    if(randomization){
      if(!first) mtr[1:ntr]<-trt
    }
    plots <- block*number+(1:ntr)
    book <- data.frame(plots, block = as.factor(block), trt = as.factor(mtr))
    names(book)[3] <- c(paste(deparse(substitute(trt))))
    names(book)[3]<-c(paste(deparse(substitute(trt))))
    if(continue){
      start0<-10^serie
      if(serie==0) start0<-0
      book$plots<-start0+1:nrow(book)
    }
    outdesign<-list(parameters=parameters,sketch=matrix(book[,3], byrow = TRUE, ncol = ntr),book=book)
    return(outdesign)
  }

  design.lsd <-function (trt,serie=2,seed=0,
                         kinds="Super-Duper",first=TRUE,
                         randomization=TRUE){
    number<-10
    if(serie>0) number<-10^serie
    r <- length(trt)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    parameters<-list(design="lsd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
    a <- 1:(r * r)
    dim(a) <- c(r, r)
    for (i in 1:r) {
      for (j in 1:r) {
        k <- i + j - 1
        if (k > r)
          k <- i + j - r - 1
        a[i, j] <- k
      }
    }
    m<-2:r
    if(randomization)m<-sample(2:r,r-1)
    a<-a[,c(1,m)]
    if(randomization){
      if (first) {
        m<-sample(1:r,r)
        a<-a[m,]
      }}
    trat<-trt[a]
    columna <- rep(gl(r, 1), r)
    fila <- gl(r, r)
    fila <- as.character(fila)
    fila <- as.numeric(fila)
    plots <- fila*number+(1:r)
    book <- data.frame(plots, row = as.factor(fila), col = as.factor(columna),
                       trat = as.factor(trat))
    names(book)[4] <- c(paste(deparse(substitute(trt))))
    outdesign<-list(parameters=parameters,sketch=matrix(book[,4], byrow = TRUE, ncol = r),book=book)
    return(outdesign)
  }

  design.split <-function (trt1, trt2,r=NULL, design=c("rcbd","crd","lsd"),serie = 2, seed = 0, kinds = "Super-Duper",
                           first=TRUE,randomization=TRUE){
    n1<-length(trt1)
    n2<-length(trt2)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    design <- match.arg(design)
    number<-10^serie +1
    if (design == "crd") {
      plan<-design.crd(trt1,r,serie, seed, kinds,randomization)
      k<-3
    }
    if (design == "rcbd"){
      plan<-design.rcbd(trt1,r,serie, seed, kinds, first,randomization)
      k<-3
    }
    if (design == "lsd") {
      plan<-design.lsd(trt1,serie, seed, kinds, first,randomization)
      r<-n1
      k<-4
    }
    book<-plan$book
    parameters<-plan$parameters
    names(parameters)[2]<-"trt1"
    parameters$applied<-parameters$design
    parameters$design<-"split"
    parameters$trt2<-trt2
    j<-0
    B<-list()
    for(i in c(1,7,2,8,3:6)){
      j<-j+1
      B[[j]]<-parameters[[i]]
      names(B)[j]<-names(parameters)[i]
    }
    nplot<-nrow(book)
    d<-NULL
    if(randomization){
      for(i in 1:nplot)d<-rbind(d,sample(trt2,n2))
    }
    else{
      d<-rbind(d,trt2[1:n2])
    }
    aa<-data.frame(book,trt2=d[,1])
    for(j in 2:n2) aa<-rbind(aa,data.frame(book,trt2=d[,j]))
    aa<-aa[order(aa[,1]),]
    splots<-rep(gl(n2,1),nplot)
    book <- data.frame(plots=aa[,1],splots,aa[,-1])
    rownames(book)<-1:(nrow(book))
    names(book)[k+1] <- c(paste(deparse(substitute(trt1))))
    names(book)[k+2] <- c(paste(deparse(substitute(trt2))))
    outdesign<-list(parameters=B,book=book)
    return(outdesign)
  }

  design.ab <-function(trt, r=NULL,serie=2,design=c("rcbd","crd","lsd"),seed=0,kinds="Super-Duper",
                       first=TRUE,randomization=TRUE ){
    design <- match.arg(design)
    if( design=="rcbd" | design=="crd") posicion <- 3
    else posicion <- 4
    serie<-serie; seed<-seed; kinds<-kinds; first<-first;
    ntr<-length(trt)
    fact<-NULL
    tr0<-1:trt[1]
    k<-0
    a<-trt[1];b<-trt[2]
    for(i in 1:a){
      for(j in 1:b){
        k<-k+1
        fact[k]<-paste(tr0[i],j)
      }
    }

    if(ntr >2) {
      for(m in 3:ntr){
        k<-0
        tr0<-fact
        fact<-NULL
        a<-a*b
        b<-trt[m]
        for(i in 1:a){
          for(j in 1:b){
            k<-k+1
            fact[k]<-paste(tr0[i],j)
          }
        }
      }
    }
    if(design=="rcbd")plan<-design.rcbd(trt=fact, r, serie, seed, kinds, first,randomization )
    if(design=="crd")plan<-design.crd(trt=fact, r, serie, seed, kinds,randomization)
    if(design=="lsd")plan<-design.lsd(trt=fact, serie, seed, kinds, first,randomization )
    parameters<-plan$parameters
    parameters$applied<-parameters$design
    parameters$design<-"factorial"
    plan<-plan$book
    trt<-as.character(plan[,posicion])
    nplan<-nrow(plan)
    A<-rep(" ",nplan*ntr)
    dim(A)<-c(nplan,ntr)
    colnames(A)<-LETTERS[1:ntr]

    for(i in 1:nplan) {
      A[i,]<-unlist(strsplit(trt[i], " "))
    }
    A<-as.data.frame(A)
    book<-data.frame(plan[,1:(posicion-1)],A)
    outdesign<-list(parameters=parameters,book=book)
    return(outdesign)
  }

  #=================
  if(design=="DIC" | design=="dic"){sort=design.crd(trat,r,serie=0)
  data=sort$book
  data$x=rep(1:length(unique(data$trat)),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(data$trat)))
  data$y=factor(data$y,unique(data$y))
  x=data$x
  y=data$y
  if(color.sep=="all"){separate=data$trat}
  if(color.sep=="none"){separate=rep("white",e=length(data$trat))}

  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y))+
    geom_tile(aes(fill=separate),color="black")+
    labs(x=label.x,y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=y,fill=separate))+
    geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}

  if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
    if(pos=="line"){
      data$y=factor(data$y,levels = rev(unique(data$y)))
      graph=ggplot(data,aes(x=x,y=y,fill=separate))+
        geom_tile(color="black")+
        labs(x=label.x,y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}
    if(pos=="column"){
      graph=ggplot(data,aes(y=x,x=y,fill=separate))+
        geom_tile(color="black")+
        labs(x=label.x,y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}}

  if(color.sep=="none"){graph=graph+
    scale_fill_manual(values = "white",label="plots")+
    labs(fill="")}
  if(ID==FALSE){graph=graph+geom_text(aes(label=trat),size=labelsize)}
  if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(trat)),size=labelsize)}
  tabela=data.frame("ID"=data$plots,
                    "trat"=data$trat)}

  #=================
  if(design=="DBC" | design=="dbc"){sort=design.rcbd(trat,r,serie=0)
  data=sort$book
  data$x=rep(1:length(unique(data$trat)),r)
  data$x=factor(data$x,unique(data$x))
  x=data$x
  block=data$block
  if(color.sep=="all"){separate=data$trat}
  if(color.sep=="block"){separate=data$block}
  if(color.sep=="none"){separate=rep("white",e=length(data$trat))}

  if(pos=="line"){graph=ggplot(data,aes(x=x,y=block,fill=separate))+
    geom_tile(color="black")+labs(y="Block",x=label.x,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize))}
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=separate))+
    geom_tile(color="black")+labs(y=label.y,x="Block",fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize))}

  if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
    data$block=factor(data$block,levels = rev(unique(data$block)))
    if(pos=="line"){
      graph=ggplot(data,aes(x=x,y=block,fill=separate))+
        geom_tile(color="black")+
        labs(x=label.x,y="Block",fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.x = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}
    if(pos=="column"){
      data$block=factor(data$block,levels = unique(data$block))
      graph=ggplot(data,aes(y=x,x=block,fill=separate))+
        geom_tile(color="black")+
        labs(x="Block",y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}}

  if(color.sep=="none"){graph=graph+
    scale_fill_manual(values = "white",label="plots")+
    labs(fill="")}
  if(color.sep=="block"){graph=graph+labs(fill="block")}
  if(ID==FALSE){graph=graph+geom_text(aes(label=trat),size=labelsize)}
  if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(trat)),size=labelsize)}
  tabela=data.frame("ID"=data$plots,
                    "block"=data$block,
                    "trat"=data$trat)}

  #=================
  if(design=="DQL" | design=="dql"){sort=design.lsd(trat,r,serie=0)
  data=sort$book
  if(color.sep=="all"){separate=data$trat}
  if(color.sep=="line"){separate=data$row}
  if(color.sep=="column"){separate=data$col}
  if(color.sep=="none"){separate=rep("white",e=length(data$trat))}

  graph=ggplot(data,aes(x=row,y=col,fill=separate))+
    geom_tile(color="black")+labs(x="Row",y="Column",fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

  if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
    data$row=factor(data$row,levels = rev(unique(data$row)))
    if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
      ruas1=rep(add.streets.x,e=(length(data$row)/length(add.streets.x)))
      data$ruas1=ruas1
      graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
    if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
      ruas2=rep(rev(add.streets.y),(length(data$row)/length(add.streets.y)))
      data$ruas2=ruas2
      graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
    if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
      ruas1=rep(add.streets.x,e=(length(data$row)/length(add.streets.x)))
      data$ruas1=ruas1
      ruas2=rep(rev(add.streets.y),(length(data$row)/length(add.streets.y)))
      data$ruas2=ruas2
      graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
    graph=graph+
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            line = element_blank())}

  if(color.sep=="none"){graph=graph+
    scale_fill_manual(values = "white",label="plots")+
    labs(fill="")}
  if(color.sep=="line"){graph=graph+labs(fill="Line")}
  if(color.sep=="column"){graph=graph+labs(fill="Column")}

  if(ID==FALSE){graph=graph+geom_text(aes(label=trat),size=labelsize)}
  if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(trat)),size=labelsize)}
  tabela=data.frame("ID"=data$plots,
                    "line"=data$row,
                    "column"=data$col,
                    "trat"=data$trat)}

  #=================
  if(design=="PSUBDIC"  | design=="psubdic"){sort=design.split(trat,trat1,r,design = "crd",serie=0)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(paste(data$trat,data$trat1))))
  data$y=factor(data$y,unique(data$y))
  x=data$x
  y=data$y
  if(color.sep=="all"){separate=paste(data$trat,data$trat1)}
  if(color.sep=="f1"){separate=data$trat}
  if(color.sep=="f2"){separate=data$trat1}
  if(color.sep=="none"){separate=rep("white",e=length(data$trat))}

  if(pos=="column"){graph=ggplot(data,aes(x=y,y=x,fill=separate))+
    geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=separate))+
    geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}

  if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
    if(pos=="line"){
      data$y=factor(data$y,levels = rev(unique(data$y)))
      graph=ggplot(data,aes(x=x,y=y,fill=separate))+
        geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.x = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}
    if(pos=="column"){
      data$y=factor(data$y,levels = unique(data$y))
      graph=ggplot(data,aes(y=x,x=y,fill=separate))+
        geom_tile(color="black")+
        labs(x=label.x,y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}}
  if(color.sep=="none"){graph=graph+
    scale_fill_manual(values = "white",label="plots")+
    labs(fill="")}
  if(ID==FALSE){graph=graph+geom_text(aes(label=paste(trat,trat1)),size=labelsize)}
  if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(paste(trat,trat1))),size=labelsize)}
  tabela=data.frame("ID"=1:length(data$plots),
                    "plot"=data$trat,
                    "split_plot"=data$trat1,
                    "Repetition"=data$r)}

  #================
  if(design=="PSUBDBC" | design=="psubdbc"){sort=design.split(trat,trat1,r,design = "rcbd",serie=0)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  x=data$x
  block=data$block
  if(color.sep=="all"){separate=paste(data$trat,data$trat1)}
  if(color.sep=="block"){separate=data$block}
  if(color.sep=="f1"){separate=data$trat}
  if(color.sep=="f2"){separate=data$trat1}
  if(color.sep=="none"){separate=rep("white",e=length(data$trat))}

  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=separate))+
    geom_tile(color="black")+labs(x="Block",y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
  if(pos=="line"){graph=ggplot(data,aes(y=block,x=x,fill=separate))+
    geom_tile(color="black")+labs(x=label.x,y="Block",fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
  if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
    if(pos=="line"){
      data$y=factor(data$block,levels = rev(unique(data$block)))
      graph=ggplot(data,aes(x=x,y=y,fill=separate))+
        geom_tile(color="black")+labs(x=label.x,y="Block",fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.x = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}
    if(pos=="column"){
      data$y=factor(data$block,levels = unique(data$block))
      graph=ggplot(data,aes(y=x,x=y,fill=separate))+
        geom_tile(color="black")+
        labs(x="Block",y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}}

  if(color.sep=="none"){graph=graph+
    scale_fill_manual(values = "white",label="plots")+
    labs(fill="")}
  if(color.sep=="block"){graph=graph+labs(fill="block")}
  if(ID==FALSE){graph=graph+geom_text(aes(label=paste(trat,trat1)),size=labelsize)}
  if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(paste(trat,trat1))),size=labelsize)}
  tabela=data.frame("ID"=1:length(data$plots),
                    "plot"=data$trat,
                    "split_plot"=data$trat1,
                    "Block"=data$block)}

  faixas=function(t1,t2,r){
    faixas.design=function (trt1, trt2, r, serie = 2, seed = 0, kinds = "Super-Duper",
                            randomization = TRUE){
      number <- 10
      if (serie > 0)
        number <- 10^serie
      n1 <- length(trt1)
      n2 <- length(trt2)
      if (seed == 0) {
        genera <- runif(1)
        seed <- .Random.seed[3]}
      set.seed(seed, kinds)
      a <- trt1[1:n1]
      b <- trt2[1:n2]
      if (randomization) {
        a <- sample(trt1, n1)
        b <- sample(trt2, n2)}
      fila <- rep(b, n1)
      columna <- a[gl(n1, n2)]
      block <- rep(1, n1 * n2)
      if (r > 1) {
        for (i in 2:r) {
          a <- trt1[1:n1]
          b <- trt2[1:n2]
          if (randomization){
            a <- sample(trt1, n1)
            b <- sample(trt2, n2)}
          fila <- c(fila, rep(b, n1))
          columna <- c(columna, a[gl(n1, n2)])
          block <- c(block, rep(i, n1 * n2))
        }
      }
      parameters <- list(design = "strip", trt1 = trt1, trt2 = trt2, r = r, serie = serie, seed = seed, kinds = kinds)
      plots <- block * number + 1:(n1 * n2)
      book <- data.frame(plots, block = as.factor(block), column = as.factor(columna),
                         row = as.factor(fila))
      names(book)[3] <- c(paste(deparse(substitute(trt1))))
      names(book)[4] <- c(paste(deparse(substitute(trt2))))
      outdesign <- list(parameters = parameters, book = book)
      return(outdesign)
    }
    outdesign <-faixas.design(t1,t2,r, serie=2,seed=45,kinds ="Super-Duper") # seed = 45
    book <-outdesign$book # field book
    book$block=factor(book$block,levels = unique(book$block))
    graphs=as.list(1:length(levels(book$block)))
    for(i in 1:length(levels(book$block))){
      d1=book[book$block==levels(book$block)[i],]
      d1$t1=factor(d1$t1,unique(d1$t1))
      d1$t2=factor(d1$t2,unique(d1$t2))
      graphs[[i]]=ggplot(d1,aes(x=t1,y=t2,fill=paste(t1,t2)))+geom_tile(color="black",show.legend = FALSE)+
        facet_wrap(~paste("Block",block))+
        ylab("")+xlab("")+
        geom_text(aes(label=paste(t1,t2)))+
        theme_classic()+theme(axis.line = element_blank(),
                              axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize),
                              strip.text = element_text(size=12))}
    requireNamespace("cowplot")
    graph=do.call("plot_grid", c(graphs, ncol=length(levels(book$block))))
    print(graph)}
  if(design=="STRIP-PLOT"  | design=="stripplot"){(graph=faixas(trat,trat1,r))}
  #=================
  if(design=="FAT2DIC"  | design=="fat2dic"){sort=design.ab(c(length(trat),length(trat1)),r,design = "crd",serie=0)
  sort$book$A=as.factor(sort$book$A)
  sort$book$B=as.factor(sort$book$B)
  levels(sort$book$A)=trat
  levels(sort$book$B)=trat1
  sort$book$trat=paste(sort$book$A,sort$book$B)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(paste(data$trat,data$trat1))))
  data$y=factor(data$y,unique(data$y))
  A=data$A
  B=data$B
  x=data$x
  y=data$y
  if(color.sep=="all"){separate=paste(data$A,data$B)}
  if(color.sep=="f1"){separate=data$A}
  if(color.sep=="f2"){separate=data$B}
  if(color.sep=="none"){separate=rep("white",e=length(data$A))}

  if(pos=="column"){graph=ggplot(data,aes(x=y,y=x,fill=separate))+
    geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=separate))+
    geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
    theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}

  if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
    if(pos=="line"){
      data$y=factor(data$y,levels = rev(unique(data$y)))
      graph=ggplot(data,aes(x=x,y=y,fill=separate))+
        geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.x = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}
    if(pos=="column"){
      data$y=factor(data$y,levels = unique(data$y))
      graph=ggplot(data,aes(y=x,x=y,fill=separate))+
        geom_tile(color="black")+
        labs(x=label.x,y=label.y,fill="Treatments")+
        theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
      if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
        ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
      if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
        ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
        data$ruas1=ruas1
        ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
        data$ruas2=ruas2
        graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
      graph=graph+
        theme(axis.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              line = element_blank())}}
  if(color.sep=="none"){graph=graph+
    scale_fill_manual(values = "white",label="plots")+
    labs(fill="")}
  if(ID==FALSE){graph=graph+geom_text(aes(label=paste(A,B)),size=labelsize)}
  if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(paste(A,B))),size=labelsize)}
  tabela=data.frame("ID"=1:length(data$plots),
                    "Factor 1"=data$A,
                    "Factor 2"=data$B)}

  #=================
  if(design=="FAT2DBC"  | design=="fat2dbc"){
    sort=design.ab(c(length(trat),length(trat1)),r,design = "rcbd",serie=0)
    sort$book$A=as.factor(sort$book$A)
    sort$book$B=as.factor(sort$book$B)
    levels(sort$book$A)=trat
    levels(sort$book$B)=trat1
    sort$book$trat=paste(sort$book$A,sort$book$B)
    data=sort$book
    data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
    data$x=factor(data$x,unique(data$x))
    A=data$A
    B=data$B
    x=data$x
    block=data$block
    if(color.sep=="all"){separate=paste(data$A,data$B)}
    if(color.sep=="block"){separate=data$block}
    if(color.sep=="f1"){separate=data$A}
    if(color.sep=="f2"){separate=data$B}
    if(color.sep=="none"){separate=rep("white",e=length(data$A))}

    if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=separate))+
      geom_tile(color="black")+labs(x="Block",y=label.y,fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
    if(pos=="line"){graph=ggplot(data,aes(y=block,x=x,fill=separate))+
      geom_tile(color="black")+labs(x=label.x,y="Block",fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}

    if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
      if(pos=="line"){
        data$y=factor(data$block,levels = rev(unique(data$block)))
        graph=ggplot(data,aes(x=x,y=y,fill=separate))+
          geom_tile(color="black")+labs(x=label.x,y="Block",fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}
      if(pos=="column"){
        data$y=factor(data$block,levels = unique(data$block))
        graph=ggplot(data,aes(y=x,x=y,fill=separate))+
          geom_tile(color="black")+
          labs(x="Block",y=label.y,fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.y = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}}
    if(color.sep=="none"){graph=graph+
      scale_fill_manual(values = "white",label="plots")+
      labs(fill="")}
    if(color.sep=="block"){graph=graph+labs(fill="block")}
    if(ID==FALSE){graph=graph+geom_text(aes(label=paste(A,B)),size=labelsize)}
    if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(paste(A,B))),size=labelsize)}
    tabela=data.frame("ID"=1:length(data$plots),
                      "Factor 1"=data$A,
                      "Factor 2"=data$B,
                      "Block"=data$block)}

  #######################################################
  if(design=="FAT3DIC"  | design=="fat3dic"){
    trat=expand.grid(trat,trat1,trat2)
    tr=paste(trat$Var1,"@#",trat$Var2,"@#",trat$Var3)
    trats=rep(tr,r)
    sorteio=sample(trats)
    sortd=data.frame(t(matrix(unlist(strsplit(sorteio,"@#")),nrow=3)))
    sorteio=paste(sortd$X1,sortd$X2,sortd$X3)
    x=rep(1:(length(sorteio)/r),r)
    y=rep(1:r,e=(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))
    if(color.sep=="all"){separate=sorteio}
    if(color.sep=="f1"){separate=sortd$X1}
    if(color.sep=="f2"){separate=sortd$X2}
    if(color.sep=="f3"){separate=sortd$X3}
    if(color.sep=="none"){separate=rep("white",e=length(sorteio))}

    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=separate))+
      geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=separate))+
      geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
    if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
      if(pos=="line"){
        data$y=factor(data$y,levels = rev(unique(data$y)))
        graph=ggplot(data,aes(x=x,y=y,fill=separate))+
          geom_tile(color="black")+labs(x=label.x,y=label.y,fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}
      if(pos=="column"){
        data$y=factor(data$y,levels = unique(data$y))
        graph=ggplot(data,aes(y=x,x=y,fill=separate))+
          geom_tile(color="black")+
          labs(x=label.x,y=label.y,fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.y = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}}
    if(color.sep=="none"){graph=graph+
      scale_fill_manual(values = "white",label="plots")+
      labs(fill="")}
    if(ID==FALSE){graph=graph+geom_text(aes(label=sorteio),size=labelsize)}
    if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(sorteio)),size=labelsize)}
    tabela=data.frame("ID"=1:length(data$plots),
                      "Factor 1"=sortd$X1,
                      "Factor 2"=sortd$X2,
                      "Factor 3"=sortd$X3)}

  if(design=="FAT3DBC"  | design=="fat3dbc"){
    trat=expand.grid(trat,trat1,trat2)
    tr=paste(trat$Var1,"@#",trat$Var2,"@#",trat$Var3)

    sorteio=matrix(NA,ncol=length(tr),nrow=r)
    for(i in 1:r){
      sorteio[i,]=sample(tr)
    }
    sorteio=as.vector(sorteio)
    sortd=data.frame(t(matrix(unlist(strsplit(sorteio,"@#")),nrow=3)))
    sorteio=paste(sortd$X1,sortd$X2,sortd$X3)

    x=rep(1:(length(sorteio)/r),e=r)
    y=rep(1:r,(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))
    if(color.sep=="all"){separate=sorteio}
    if(color.sep=="f1"){separate=sortd$X1}
    if(color.sep=="f2"){separate=sortd$X2}
    if(color.sep=="f3"){separate=sortd$X3}
    if(color.sep=="none"){separate=rep("white",e=length(sorteio))}

    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=separate))+
      geom_tile(color="black")+labs(x="block",y=label.y,fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=separate))+
      geom_tile(color="black")+labs(x=label.x,y="block",fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}

    if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
      if(pos=="line"){
        data$y=factor(data$y,levels = rev(unique(data$y)))
        graph=ggplot(data,aes(x=x,y=y,fill=separate))+
          geom_tile(color="black")+labs(x=label.x,y="Block",fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(add.streets.y,(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(add.streets.y,(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}
      if(pos=="column"){
        data$y=factor(data$y,levels = unique(data$y))
        graph=ggplot(data,aes(y=x,x=y,fill=separate))+
          geom_tile(color="black")+
          labs(x="Block",y=label.y,fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(rev(add.streets.y),e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(rev(add.streets.y),e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.y = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}}
    if(color.sep=="none"){graph=graph+
      scale_fill_manual(values = "white",label="plots")+
      labs(fill="")}
    if(color.sep=="block"){graph=graph+labs(fill="block")}
    if(ID==FALSE){graph=graph+geom_text(aes(label=sorteio),size=labelsize)}
    if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(sorteio)),size=labelsize)}
  }

  if(design=="PSUBSUBDBC"  | design=="psubsubdbc"){
    sorteio=list()
    for(i in 1:r){
      sorteio[[i]]=sample(trat)}
    nv1=length(trat)
    nv2=length(trat1)
    nv3=length(trat2)
    sorteiof1=rep(unlist(sorteio),e=nv2*nv3)

    sorteiof2=list()
    for(i in 1:(r*nv1)){
      sorteiof2[[i]]=sample(trat1)}
    sorteiof2=rep(unlist(sorteiof2),e=nv3)

    sorteiof3=list()
    for(i in 1:(r*nv1*nv2)){
      sorteiof3[[i]]=sample(trat2)}
    sorteiof3=unlist(sorteiof3)
    data.frame(sorteiof1,sorteiof2,sorteiof3)
    tr=paste(sorteiof1," x ",sorteiof2," x ",sorteiof3)
    sorteio=as.vector(tr)
    x=rep(1:(length(sorteio)/r),r)
    y=rep(1:r,e=(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))
    if(color.sep=="all"){separate=sorteio}
    if(color.sep=="f1"){separate=sorteiof1}
    if(color.sep=="f2"){separate=sorteiof2}
    if(color.sep=="f3"){separate=sorteiof3}
    if(color.sep=="none"){separate=rep("white",e=length(sorteio))}

    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=separate))+
      geom_tile(color="black")+labs(x="block",y=label.y,fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=separate))+
      geom_tile(color="black")+labs(x=label.x,y="block",fill="Treatments")+
      theme_classic()+theme(axis.line = element_blank())+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))}

    if(is.na(add.streets.x[1])==FALSE | is.na(add.streets.y[1])==FALSE){
      if(pos=="line"){
        data$y=factor(data$y,levels = rev(unique(data$y)))
        graph=ggplot(data,aes(x=x,y=y,fill=separate))+
          geom_tile(color="black")+labs(x=label.x,y="Block",fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(add.streets.y,e=(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}
      if(pos=="column"){
        data$y=factor(data$y,levels = unique(data$y))
        graph=ggplot(data,aes(y=x,x=y,fill=separate))+
          geom_tile(color="black")+
          labs(x="Block",y=label.y,fill="Treatments")+
          theme_classic()+theme(axis.text=element_text(size=axissize),legend.text=element_text(size=axissize),legend.title=element_text(size=legendsize))

        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==TRUE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          graph=graph+facet_grid(~data$ruas1,scales="free",space = "free")}
        if(is.na(add.streets.y[1])==FALSE & is.na(add.streets.x[1])==TRUE){
          ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~.,scales="free",space = "free")}
        if(is.na(add.streets.x[1])==FALSE & is.na(add.streets.y[1])==FALSE){
          ruas1=rep(add.streets.x,e=(length(data$x)/length(add.streets.x)))
          data$ruas1=ruas1
          ruas2=rep(rev(add.streets.y),(length(data$x)/length(add.streets.y)))
          data$ruas2=ruas2
          graph=graph+facet_grid(data$ruas2~data$ruas1,scales="free",space = "free")}
        graph=graph+
          theme(axis.text.y = element_blank(),
                strip.background = element_blank(),
                strip.text = element_blank(),
                line = element_blank())}}
    if(color.sep=="none"){graph=graph+
      scale_fill_manual(values = "white",label="plots")+
      labs(fill="")}
    if(color.sep=="block"){graph=graph+labs(fill="block")}
    if(ID==FALSE){graph=graph+geom_text(aes(label=sorteio),size=labelsize)}
    if(ID==TRUE){graph=graph+geom_text(aes(label=1:length(sorteio)),size=labelsize)}
    tabela=data.frame("ID"=1:length(data$plots),
                      "plot"=sorteiof1,
                      "split_plot"=sorteiof2,
                      "split_split_plot"=sorteiof3)}
  if(isTRUE(ID)==TRUE){print(data)}
  if(export.csv==TRUE){write.csv(tabela,"dataset.csv")}
  #=================
  if(design!="STRIP-PLOT"){print(graph+labs(caption = comment.caption))}
}

