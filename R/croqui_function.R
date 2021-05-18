#' Utils: Experimental sketch
#'
#' @description Experimental sketching function
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @name sketch
#' @param trat Vector with factor A levels
#' @param trat1 Vector with levels of factor B (Set to NULL if not factorial or psub)
#' @param r Number of repetitions
#' @param design Experimental design ("dic", "dbc", "dql","psubdic","psubdbc","fat2dic","fat2dbc")
#' @param pos Repeat position (line or column)
#' @param ncol Default is NA. Warning!!! Use only in a completely randomized design
#' @keywords croqui
#' @keywords experimental
#' @return Returns an experimental sketch according to the specified design.
#' @note The sketches have only a rectangular shape, and the blocks (in the case of randomized blocks) can be in line or in a column.
#' @export
#' @examples
#' Trat=paste("Treatments",1:6)
#' sketch(Trat,r=3)
#' set.seed(1)
#' sketch(Trat,r=3)

sketch=function(trat,
                trat1=NULL,
                r,
                design="dic",
                pos="line",
                ncol=NA){
  if(is.na(ncol)==TRUE){ncol=r}
  requireNamespace("ggplot2")
  requireNamespace("gridExtra")
  requireNamespace("grid")
  if(design=="dic"){sort=design.crd(trat,r,serie=0)}
  if(design=="dbc"){sort=design.rcbd(trat,r,serie=0)}
  if(design=="dql"){sort=design.lsd(trat,r,serie=0)}
  if(design=="psubdic"){sort=design.split(trat,trat1,r,design = "crd",serie=0)}
  if(design=="psubdbc"){sort=design.split(trat,trat1,r,design = "rcbd",serie=0)}
  if(design=="fat2dic"){sort=design.ab(c(length(trat),length(trat1)),r,design = "crd",serie=0)
  sort$book$trat=paste(sort$book$A,sort$book$B)}
  if(design=="fat2dbc"){sort=design.ab(c(length(trat),length(trat1)),r,design = "rcbd",serie=0)
  sort$book$trat=paste(sort$book$A,sort$book$B)}
  if(pos=="column"){sort$book$trat=as.vector(matrix(sort$book$trat,nrow = r,byrow=TRUE))}
  if(pos=="line"){sort$book$trat=as.vector(t(matrix(sort$book$trat,nrow = r,byrow=TRUE)))}

  if(design=="psubdic" & pos=="column"){sort$book$trat=as.vector(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE))}
  if(design=="psubdbc" & pos=="column"){sort$book$trat=as.vector(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE))}
  if(design=="fat2dic" & pos=="column"){sort$book$trat=as.vector(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE))}
  if(design=="fat2dbc" & pos=="column"){sort$book$trat=as.vector(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE))}

  if(design=="psubdic" & pos=="line"){sort$book$trat=as.vector(t(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE)))}
  if(design=="psubdbc" & pos=="line"){sort$book$trat=as.vector(t(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE)))}
  if(design=="fat2dic" & pos=="line"){sort$book$trat=as.vector(t(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE)))}
  if(design=="fat2dbc" & pos=="line"){sort$book$trat=as.vector(t(matrix(paste(sort$book$trat,sort$book$trat1),nrow = r,byrow=TRUE)))}

  sort$book$trat=as.factor(sort$book$trat)

  if(pos=="column"){ncol=ncol}
  if(pos=="line"){ncol=length(levels(sort$book$trat))}
  gs <- lapply(sort$book$trat, function(ii)
    grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5,lwd=3)),
             textGrob(ii)))
  grid.arrange(grobs=gs, ncol=ncol)
}
