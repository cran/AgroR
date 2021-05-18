#' Graph: Interaction plot
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Performs an interaction graph from an output of the FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC commands.
#' @param a FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC object
#' @export
#' @import ggplot2
#' @importFrom crayon green
#' @importFrom crayon bold
#' @importFrom crayon italic
#' @importFrom crayon red
#' @importFrom crayon blue
#' @import stats
#' @return Returns an interaction graph with averages and letters from the multiple comparison test
#' @examples
#' data(cloro)
#' attach(cloro)
#' a=FAT2DIC(f1, f2, resp)
#' plot_interaction(a)

plot_interaction=function(a){
  data=a[[1]]
  requireNamespace("ggplot2")
  graph=ggplot(data$data,
               aes(x=data$data[,1],
                   y=data$data$media,
                   color=data$data[,2],
                   group=data$data[,2]))+
    geom_point(show.legend = TRUE, size=3)+
    geom_errorbar(aes(ymax=data$data$media+data$data$desvio,
                      ymin=data$data$media-data$data$desvio),
                  width=0.05,size=0.8)+
    geom_line(size=0.8)+
    geom_label(aes(label=data$data$numero),show.legend = FALSE)+
    data$theme+
    labs(caption=data$labels$caption,
         color=data$labels$fill,
         x=data$labels$x,
         y=data$labels$y)
  print(graph)
  graphs=list(graph)
}
