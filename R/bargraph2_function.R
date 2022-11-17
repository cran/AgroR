#' Graph: Bar graph for one factor model 2
#'
#' @description This is a function of the bar graph for one factor
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param fill Fill bars
#' @param point.color Point color
#' @param point.size Point size
#' @param point.shape Format point
#' @param text.color Text color
#' @param label.color Label color
#' @param bar.color Errorbar color
#' @param title.size Title size
#' @param y.text Y-axis height for x-axis legend
#' @param add.info Add other information
#' @param y.info Y-axis height for other information
#' @param color.info Color text information
#' @export
#' @return Returns a bar chart for one factor
#' @seealso \link{radargraph}, \link{barplot_positive}, \link{plot_TH}, \link{plot_TH1}, \link{corgraph}, \link{spider_graph}, \link{line_plot}, \link{plot_cor}, \link{plot_interaction}, \link{plot_jitter}, \link{seg_graph}, \link{TBARPLOT.reverse}
#' @examples
#' data("laranja")
#'a=with(laranja, DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,sup = 10,
#'      family = "serif",
#'      ylab = "Number of fruits/plants"))
#'bar_graph2(a)
#'bar_graph2(a,fill="darkblue",point.color="orange",text.color='white')

bar_graph2=function(model,
                    point.color="black",
                    point.size=2,
                    point.shape=16,
                    text.color="black",
                    label.color="black",
                    bar.color="black",
                    title.size=14,
                    y.text=0,
                    add.info=NA,
                    y.info=0,
                    color.info="black",
                    fill="lightblue"){
  requireNamespace("ggplot2")
  data=model[[1]]$data
  media=data$media
  desvio=data$desvio
  trats=data$trats
  limite=data$limite
  letra=data$letra
  groups=data$groups
  sup=model[[1]]$plot$sup

  graph=ggplot(data,aes(x=trats,
                        y=media))+
    model[[1]]$theme+
    geom_col(fill=fill,size=0.3,color="black")+
    geom_errorbar(aes(ymin=media-desvio,
                      ymax=media+desvio),color=bar.color,width=0)+
    geom_point(size=point.size,color=point.color,fill=point.color)+
    geom_text(aes(y=media+desvio+sup,
                  x=trats,
                  label = letra),vjust=0,size=model[[1]]$plot$labelsize,color=label.color,family=model[[1]]$plot$family)+
    geom_text(aes(y=y.text,
                  x=trats,
                  label = trats),hjust=0,angle=90,size=model[[1]]$plot$labelsize,color=text.color,family=model[[1]]$plot$family)+
    labs(x=model[[1]]$labels$x,
         y=model[[1]]$labels$y)+
    theme(axis.text = element_text(size=model[[1]]$plot$textsize,color="black"),
          axis.title = element_text(size=title.size,color="black"),
          axis.text.x = element_blank(),
          strip.text = element_text(size=model[[1]]$plot$textsize),
          legend.position = "none")+
    scale_x_discrete(limits=trats)+
    ylim(layer_scales(model[[1]])$y$range$range*1.1)
  if(is.na(add.info[1])==FALSE){
    graph=graph+geom_text(aes(y=y.info,x=trats,label=add.info),hjust=0,
                          size=model[[1]]$plot$labelsize,color=color.info,
                          family=model[[1]]$plot$family)
  }
  graph
}
