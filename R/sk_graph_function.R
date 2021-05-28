#' Graph: Scott-Knott graphics
#'
#' @description This is a function of the bar graph for the Scott-Knott test
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param dec Number of cells
#' @param horiz Horizontal Column (\emph{default} is TRUE)
#' @param transf If the data has been transformed (\emph{default} is FALSE)
#' @export
#' @return Returns a bar chart with columns separated by color according to the Scott-Knott test
#' @seealso \link{radargraph}, \link{barplot_positive}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @examples
#' data("laranja")
#' attach(laranja)
#'a=DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,
#'      ylab = "Number of fruits/plants")
#'sk_graph(a,dec = 0,horiz = FALSE)


sk_graph=function(model,
                  dec=3,
                  horiz=TRUE,
                  transf=FALSE){
  requireNamespace("ggplot2")

  data=model[[1]]$data
  if(transf==FALSE){data=data[,c(5,1,2)]}
  if(transf==TRUE){data=data[,c(6,3,2)]}
  if(horiz==TRUE){
  graph=ggplot(data,aes(y=as.vector(data[,1]),x=as.vector(data[,2])))+
    model[[1]]$theme+
    geom_col(aes(fill=as.vector(data[,3])),
             size=0.3,
             color="black")+
    geom_label(aes(x=as.vector(data[,2])+1/15*as.vector(data[,2]),
                   y=as.vector(data[,1]),
                   label = paste(round(as.vector(data[,2]),dec),as.vector(data[,3]))),
               fill="lightyellow",hjust=0)+
    labs(y=model[[1]]$labels$x,
         x=model[[1]]$labels$y)+
    theme(axis.text = element_text(size=12,color="black"),
          strip.text = element_text(size=12),
          legend.position = "none")+
    scale_y_discrete(limits=data$Tratamentos)+
    xlim(layer_scales(model[[1]])$y$range$range)}
  if(horiz==FALSE){
    graph=ggplot(data,aes(x=as.vector(data[,1]),y=as.vector(data[,2])))+
      model[[1]]$theme+
      geom_col(aes(fill=as.vector(data[,3])),size=0.3,
               color="black")+
      geom_label(aes(y=as.vector(data[,2])+1/15*as.vector(data[,2]),
                     x=as.vector(data[,1]),
                     label = paste(round(as.vector(data[,2]),dec),as.vector(data[,3]))),
                 fill="lightyellow",vjust=0)+
      labs(x=model[[1]]$labels$x,
           y=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=12,color="black"),
            strip.text = element_text(size=12),
            legend.position = "none")+
      scale_x_discrete(limits=data$Tratamentos)+
      ylim(layer_scales(model[[1]])$y$range$range)}
  graph
}
