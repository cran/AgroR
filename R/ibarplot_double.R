#' Invert letters for two factor chart
#' @description invert uppercase and lowercase letters in graph for factorial scheme the subdivided plot with significant interaction
#'
#' @param analysis FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC object
#' @export
#' @return Return column chart for two factors
#' @examples
#' data(covercrops)
#' attach(covercrops)
#' a=FAT2DBC(A, B, Bloco, Resp, ylab=expression("Yield"~(Kg~"100 m"^2)),
#' legend = "Cover crops",alpha.f = 0.3,family = "serif")
#' ibarplot.double(a)



ibarplot.double=function(analysis){
  lista=analysis[[1]]$plot
  data=analysis[[1]]$plot$graph
  media=data$media
  desvio=data$desvio
  f1=data$f1
  f2=data$f2
  letra=data$letra
  letra1=data$letra1
  requireNamespace("ggplot2")
  ggplot(data,aes(x=f1,y=media,fill=f2))+
    lista$theme+
    ylab(lista$ylab)+
    geom_col(position = position_dodge(width = 0.9),color="black")+
    geom_errorbar(aes(ymin=media-desvio,ymax=media+desvio),
                  position = position_dodge(width = 0.9),width=0.3)+
    geom_text(aes(y=media+desvio+lista$sup,label=paste(round(media,lista$dec),
                                                       toupper(letra),
                                                       tolower(letra1),sep = "")),
              position = position_dodge(width = 0.9),family=lista$family)+
    theme(axis.text = element_text(size = lista$textsize,family = lista$family,color="black"),
          axis.title = element_text(size = lista$textsize,family = lista$family,color="black"),
          legend.text = element_text(size = lista$textsize,family = lista$family,color="black"),
          legend.title = element_text(size = lista$textsize,family = lista$family,color="black"))
}
