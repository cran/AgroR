#' Graph: Bar graph for one factor with facets
#'
#' @description This is a function of the bar graph for one factor with facets
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param facet vector with facets
#' @param theme ggplot2 theme
#' @param fill fill bars
#' @param horiz horizontal bar or point (\emph{default} is FALSE)
#' @param geom graph type (columns or segments)
#' @param width.bar width of the error bars of a regression graph.
#' @param pointsize Point size
#' @param facet.background Color background in facet
#' @export
#' @return Returns a bar chart for one factor
#'
#' @examples
#' library(AgroR)
#' data("laranja")
#' a=with(laranja, DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,sup = 10,family = "serif",
#'      ylab = "Number of fruits/plants"))
#' barfacet(a,c("S1","S1","S1","S1","S1",
#'              "S2","S2","S3","S3"))

barfacet=function(model,
                  facet=NULL,
                  theme=theme_bw(),
                  horiz=FALSE,
                  geom="bar",
                  fill="lightblue",
                  pointsize=4.5,
                  width.bar=0.15,
                  facet.background="gray80"){
  requireNamespace("ggplot2")
  data=model[[1]]$data
  media=data$media
  desvio=data$desvio
  trats=data$trats
  limite=data$limite
  letra=data$letra
  groups=data$groups
  if(is.null(facet[1])==FALSE){
    data$trats=as.character(data$trats)
    fac=factor(facet,unique(facet))
    nomes=unique(facet)
    comp=tapply(fac,fac,length)
    n=length(levels(fac))
    graph=as.list(1:n)
    data$fac=fac
    sup=model[[1]]$plot$sup
    if(geom=="point" & horiz==FALSE){
      graph=ggplot(data,
                   aes(x=trats,
                       y=media))+
        theme+
        geom_errorbar(aes(ymin=media-desvio,
                          ymax=media+desvio),width=model[[1]]$plot$width.bar)+
        geom_point(fill=fill,shape=21,size=pointsize,color="black")+
        geom_text(aes(y=media+desvio+sup,
                      x=trats,
                      label = letra),vjust=0)+
        labs(x=model[[1]]$labels$x,
             y=model[[1]]$labels$y)+
        facet_grid(~fac,scales = "free", space='free')+
        theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
              strip.text = element_text(size=model[[1]]$plot$textsize),
              legend.position = "none",
              # axis.text.x = element_text(angle = model[[1]]$plot$angle),
              strip.background = element_rect(fill=facet.background))+
        ylim(layer_scales(model[[1]])$y$range$range*1.1)}
    if(geom=="bar" & horiz==FALSE){
      graph=ggplot(data,
                   aes(x=trats,
                       y=media))+
        theme+
        geom_col(fill=fill,size=0.3,color="black")+
        geom_errorbar(aes(ymin=media-desvio,
                          ymax=media+desvio),width=model[[1]]$plot$width.bar)+
        geom_text(aes(y=media+desvio+sup,
                      x=trats,
                      label = letra),vjust=0)+
        labs(x=model[[1]]$labels$x,
             y=model[[1]]$labels$y)+
        facet_grid(~fac,scales = "free", space='free')+
        theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
              strip.text = element_text(size=model[[1]]$plot$textsize),
              legend.position = "none",
              # axis.text.x = element_text(angle = model[[1]]$plot$angle),
              strip.background = element_rect(fill=facet.background))+
        ylim(layer_scales(model[[1]])$y$range$range*1.1)}}
  if(geom=="point" & horiz==TRUE){
    graph=ggplot(data,
                 aes(y=trats,
                     x=media))+
      theme+
      geom_errorbar(aes(xmin=media-desvio,
                        xmax=media+desvio),width=model[[1]]$plot$width.bar)+
      geom_point(fill=fill,shape=21,size=pointsize,color="black")+
      geom_text(aes(x=media+desvio+sup,
                    y=trats,
                    label = letra),hjust=0)+
      labs(y=model[[1]]$labels$x,
           x=model[[1]]$labels$y)+
      facet_grid(fac,scales = "free", space='free')+
      theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
            strip.text = element_text(size=model[[1]]$plot$textsize),
            legend.position = "none",strip.background = element_rect(fill=facet.background))+
      xlim(layer_scales(model[[1]])$y$range$range*1.1)}
  if(geom=="bar" & horiz==TRUE){
    graph=ggplot(data,
                 aes(y=trats,
                     x=media))+
      theme+
      geom_col(fill=fill,size=0.3,color="black")+
      geom_errorbar(aes(xmin=media-desvio,
                        xmax=media+desvio),width=model[[1]]$plot$width.bar)+
      geom_text(aes(x=media+desvio+sup,
                    y=trats,
                    label = letra),hjust=0)+
      labs(y=model[[1]]$labels$x,
           x=model[[1]]$labels$y)+
      facet_grid(fac,scales = "free", space='free')+
      theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
            strip.text = element_text(size=model[[1]]$plot$textsize),
            legend.position = "none",strip.background = element_rect(fill=facet.background))+
      xlim(layer_scales(model[[1]])$y$range$range*1.1)}
graph}
