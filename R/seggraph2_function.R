#' Graph: Point graph for one factor model 2
#'
#' @description This is a function of the point graph for one factor
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param theme ggplot2 theme
#' @param horiz Horizontal Column (\emph{default} is TRUE)
#' @param pointsize Point size
#' @param pointshape Format point (default is 16)
#' @param vjust vertical adjusted
#' @export
#' @return Returns a point chart for one factor
#' @seealso \link{radargraph}, \link{barplot_positive}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @examples
#' data("laranja")
#' a=with(laranja, DBC(trat, bloco, resp,
#'        mcomp = "sk",angle=45,
#'        ylab = "Number of fruits/plants"))
#' seg_graph2(a,horiz = FALSE)

seg_graph2=function(model,
                    theme=theme_gray(),
                    pointsize=4,
                    pointshape=16,
                    horiz=TRUE,
                    vjust=-0.6){
  requireNamespace("ggplot2")
  data=model[[1]]$data
  media=data$media
  desvio=data$desvio
  trats=data$trats
  limite=data$limite
  letra=data$letra
  groups=data$groups
  sup=model[[1]]$plot$sup
  textsize=model[[1]]$plot$textsize
  if(horiz==TRUE){
  graph=ggplot(data,aes(y=trats,
                          x=media))+theme+
      geom_errorbar(aes(xmin=media-desvio,
                        xmax=media+desvio),width=0,size=0.8)+
      geom_point(size=pointsize,shape=16, fill="black", color="black")+
      geom_text(aes(x=media,
                    y=trats,
                    label = letra),vjust=vjust,family=model[[1]]$plot$family)+
      labs(y=model[[1]]$labels$x,
           x=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=12,color="black"),
            strip.text = element_text(size=12),
            legend.position = "none")+
      scale_y_discrete(limits=trats)+
      xlim(layer_scales(model[[1]])$y$range$range)}
  if(horiz==FALSE){
    graph=ggplot(data,aes(x=trats,
                          y=media))+theme+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),width=0,size=0.8)+
      geom_point(size=pointsize,shape=16, fill="black", color="black")+
      geom_text(aes(y=media,
                     x=trats,
                     label = letra),vjust=vjust,angle=90,family=model[[1]]$plot$family)+
      labs(x=model[[1]]$labels$x,
           y=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=textsize,color="black"),
            strip.text = element_text(size=textsize),
            legend.position = "none")+
      scale_x_discrete(limits=trats)+
      ylim(layer_scales(model[[1]])$y$range$range)}
  graph
}

