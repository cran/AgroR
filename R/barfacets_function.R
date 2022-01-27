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
#' @param geom graph type (columns or segments)
#' @param width.bar width of the error bars of a regression graph.
#' @param pointsize Point size
#' @export
#' @return Returns a bar chart for one factor
#'
#' @examples
#' library(AgroR)
#' data("laranja")
#' a=with(laranja, DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,
#'      ylab = "Number of fruits/plants"))
#' barfacet(a,c("S1","S1","S1","S1","S1",
#'              "S2","S2","S3","S3"))

barfacet=function(model,
                  facet=NULL,
                  theme=theme_bw(),
                  geom="bar",
                  fill="lightblue",
                  pointsize=4.5,
                  width.bar=0.15){
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

      if(geom=="point"){
      graph=ggplot(data,
                   aes(x=trats,
                       y=media))+
        theme+
        geom_errorbar(aes(ymin=media-desvio,
                          ymax=media+desvio),width=width.bar)+
        geom_point(fill=fill,shape=21,size=pointsize,color="black")+
        geom_text(aes(y=media+desvio+1/15*media,
                      x=trats,
                      label = letra),vjust=0)+
        labs(x=model[[1]]$labels$x,
             y=model[[1]]$labels$y)+
        facet_grid(~fac,scales = "free", space='free')+
        theme(axis.text = element_text(size=12,color="black"),
              strip.text = element_text(size=12),
              legend.position = "none")+
        ylim(layer_scales(model[[1]])$y$range$range)}
      if(geom=="bar"){
        graph=ggplot(data,
                     aes(x=trats,
                         y=media))+
          theme+
          geom_col(fill=fill,size=0.3,color="black")+
          geom_errorbar(aes(ymin=media-desvio,
                            ymax=media+desvio),width=0.2)+
          geom_text(aes(y=media+desvio+1/15*media,
                        x=trats,
                        label = letra),vjust=0)+
          labs(x=model[[1]]$labels$x,
               y=model[[1]]$labels$y)+
          facet_grid(~fac,scales = "free", space='free')+
          theme(axis.text = element_text(size=12,color="black"),
                strip.text = element_text(size=12),
                legend.position = "none")+
          ylim(layer_scales(model[[1]])$y$range$range)}}
  graph}
