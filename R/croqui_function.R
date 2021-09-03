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
#' @param design Experimental design ("dic", "dbc", "dql","psubdic","psubdbc","fat2dic","fat2dbc")
#' @param pos Repeat position (line or column)
#' @keywords croqui
#' @keywords experimental
#' @return Returns an experimental sketch according to the specified design.
#' @note The sketches have only a rectangular shape, and the blocks (in the case of randomized blocks) can be in line or in a column.
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
#'
#' #=============================
#' # Randomized block design
#' #=============================
#' sketch(Trat, r=3, design="dbc")
#' sketch(Trat, r=3, design="dbc",pos="column")
#'
#' #=============================
#' # Completely randomized experiments in double factorial
#' #=============================
#' sketch(trat=c("A","B"),
#'        trat1=c("A","B","C"),
#'        design = "fat2dic",
#'        r=3)
#'
#' sketch(trat=c("A","B"),
#'        trat1=c("A","B","C"),
#'        design = "fat2dic",
#'        r=3,
#'        pos="column")

sketch=function(trat,
                trat1=NULL,
                trat2=NULL,
                r,
                design="dic",
                pos="line"){
  requireNamespace("ggplot2")
  requireNamespace("gridExtra")
  requireNamespace("grid")

  #=================
  if(design=="dic"){sort=design.crd(trat,r,serie=0)
  data=sort$book
  data$x=rep(1:length(unique(data$trat)),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(data$trat)))
  data$y=factor(data$y,unique(data$y))
  x=data$x
  y=data$y
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=trat))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=y,fill=trat))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}}

  #=================
  if(design=="dbc"){sort=design.rcbd(trat,r,serie=0)
  data=sort$book
  data$x=rep(1:length(unique(data$trat)),r)
  data$x=factor(data$x,unique(data$x))
  x=data$x
  block=data$block
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=block,fill=trat))+
    geom_tile(color="black")+labs(y="Block",x="",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=trat))+
    geom_tile(color="black")+labs(y="",x="Block",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}}

  #=================
  if(design=="dql"){sort=design.lsd(trat,r,serie=0)
  data=sort$book
  graph=ggplot(data,aes(x=row,y=col,fill=trat))+
    geom_tile(color="black")+labs(x="Row",y="Column",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}

  #=================
  if(design=="psubdic"){sort=design.split(trat,trat1,r,design = "crd",serie=0)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(paste(data$trat,data$trat1))))
  data$y=factor(data$y,unique(data$y))
  x=data$x
  y=data$y

  if(pos=="column"){graph=ggplot(data,aes(x=y,y=x,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}}

  #=================
  if(design=="psubdbc"){sort=design.split(trat,trat1,r,design = "rcbd",serie=0)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  x=data$x
  block=data$block
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="Block",y="",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(y=block,x=x,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="",y="Block",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}}

  #=================
  if(design=="fat2dic"){sort=design.ab(c(length(trat),length(trat1)),r,design = "crd",serie=0)
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
  if(pos=="column"){graph=ggplot(data,aes(x=y,y=x,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}}

  #=================
  if(design=="fat2dbc"){sort=design.ab(c(length(trat),length(trat1)),r,design = "rcbd",serie=0)
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
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="Block",y="",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(y=block,x=x,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="",y="Block",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}}

  if(design=="fat3dic"){
    trat=expand.grid(trat,trat1,trat2)
    tr=paste(trat$Var1,trat$Var2,trat$Var3)
    trats=rep(tr,r)
    sorteio=sample(trats)
    x=rep(1:(length(sorteio)/r),r)
    y=rep(1:r,e=(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))

    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=sorteio))+
      geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=sorteio))+
      geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
    print(graph)}

  if(design=="fat3dbc"){
    trat=expand.grid(trat,trat1,trat2)
    tr=paste(trat$Var1,trat$Var2,trat$Var3)

    sorteio=matrix(NA,ncol=length(tr),nrow=r)
    for(i in 1:r){
      sorteio[i,]=sample(tr)
    }
    sorteio=as.vector(sorteio)
    x=rep(1:(length(sorteio)/r),e=r)
    y=rep(1:r,(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))
    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=sorteio))+
      geom_tile(color="black")+labs(x="block",y="",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=sorteio))+
      geom_tile(color="black")+labs(x="",y="block",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
  }
  #=================
  print(graph)
}

