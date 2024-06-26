#' Graph: Bar graph for one factor
#'
#' @description This is a function of the bar graph for one factor
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param fill fill bars
#' @param horiz Horizontal Column (\emph{default} is TRUE)
#' @param width.col Width Column
#' @param axis.0 If TRUE causes the columns or bars to start just above the axis line.
#' @export
#' @return Returns a bar chart for one factor
#' @seealso \link{radargraph}, \link{barplot_positive}, \link{plot_TH}, \link{plot_TH1}, \link{corgraph}, \link{spider_graph}, \link{line_plot}, \link{plot_cor}, \link{plot_interaction}, \link{plot_jitter}, \link{seg_graph}, \link{TBARPLOT.reverse}
#' @examples
#' data("laranja")
#'a=with(laranja, DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,
#'      ylab = "Number of fruits/plants"))
#'bar_graph(a,horiz = FALSE)


bar_graph=function(model,
                   fill="lightblue",
                   horiz=TRUE,
                   width.col=0.9,
                   axis.0=FALSE){
  requireNamespace("ggplot2")
  data=model[[1]]$data
  media=data$media
  desvio=data$desvio
  trats=data$trats
  limite=data$limite
  letra=data$letra
  groups=data$groups
  sup=model[[1]]$plot$sup
  if(horiz==TRUE){
    graph=ggplot(data,aes(y=trats,
                          x=media))+
      model[[1]]$theme+
      geom_col(size=0.3,fill=fill,
               color="black",width = width.col)+
      geom_errorbar(aes(xmin=media-desvio,
                        xmax=media+desvio),width=model[[1]]$plot$width.bar)+
      geom_text(aes(x=media+desvio+sup,
                     y=trats,
                     label = letra),hjust=0)+
      labs(y=model[[1]]$labels$x,
           x=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
            strip.text = element_text(size=model[[1]]$plot$textsize),
            legend.position = "none")+
      scale_y_discrete(limits=trats)+
      xlim(layer_scales(model[[1]])$y$range$range*1.1)
    if(axis.0==TRUE){graph=graph+scale_x_continuous(expand = c(0,0))}}
  if(horiz==FALSE){
    graph=ggplot(data,aes(x=trats,
                          y=media))+
      model[[1]]$theme+
      geom_col(fill=fill,size=0.3,color="black",width = width.col)+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),width=model[[1]]$plot$width.bar)+
      geom_text(aes(y=media+desvio+sup,
                     x=trats,
                     label = letra),vjust=0)+
      labs(x=model[[1]]$labels$x,
           y=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
            strip.text = element_text(model[[1]]$plot$textsize),
            legend.position = "none")+
      scale_x_discrete(limits=trats)+
      ylim(layer_scales(model[[1]])$y$range$range*1.1)
    if(axis.0==TRUE){graph=graph+scale_y_continuous(expand = c(0,0))}}
  graph
}
