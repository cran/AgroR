#' utils: group graphs of the output of simple experiments in dic, dbc or dql
#'
#' @description group graphs of the output of simple experiments into dic, dbc or dql. It is possible to group up to 6 graphs in different arrangements (see model argument)
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param output Vector with the outputs of the DIC, DBC or DQL functions
#' @param model Graph arrangement model, see in detail.
#' @export
#' @details
#' - `type1`: Two graphs next to each other
#' - `type2`: Two graphs one below the other
#' - `type3`: Three graphs, two top and one centered below
#' - `type4`: Three graphs one below the other
#' - `type5`: Four graphs, two at the top and two at the bottom
#' - `type6`: Four graphs one below the other
#' - `type7`: Five graphs, two at the top, two in the middle and one centered at the bottom
#' - `type8`: Five graphs, three at the top, two centered at the bottom
#' - `type9`: Six graphs, three at the top, three centered at the bottom
#' - `type10`: Six graphs, two at the top, two in the middle and two at the bottom
#' @return returns grouped graphs
#' @examples
#' data("pomegranate")
#' attach(pomegranate)
#' a=DIC(trat, WL, geom = "point", ylab = "WL")
#' b=DIC(trat, SS, geom = "point", ylab="SS")
#' c=DIC(trat, AT, geom = "point", ylab = "AT")
#' grid.onefactor(c(a,b),model = "type1")
#' grid.onefactor(c(a,b),model = "type2")
#' grid.onefactor(c(a,b,c),model = "type3")
#' grid.onefactor(c(a,b,c),model = "type4")

grid.onefactor=function(output,model="type1"){
  requireNamespace("gridExtra")
  requireNamespace("ggplot2")
  arrangeplot=ggplot()
  fonte=theme(title=element_text(family = output[[1]]$plot$family))
  if(model=="type1"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,ncol=2)}
  if(model=="type2"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,ncol=1)}
  if(model=="type3"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              layout_matrix=rbind(c(1,1,2,2),
                                                                  c(NA,3,3,NA)))}
  if(model=="type4"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              ncol=1)}
  if(model=="type5"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              output[[4]]+labs(title="D)")+fonte,ncol=2)}
  if(model=="type6"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              output[[4]]+labs(title="D)")+fonte,ncol=1)}
  if(model=="type7"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              output[[4]]+labs(title="D)")+fonte,
                                              output[[5]]+labs(title="E)")+fonte,
                                              layout_matrix=rbind(c(1,1,2,2),
                                                                  c(3,3,4,4),
                                                                  c(NA,5,5,NA)))}
  if(model=="type8"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              output[[4]]+labs(title="D)")+fonte,
                                              output[[5]]+labs(title="E)")+fonte,
                                              layout_matrix=rbind(c(1,1,2,2,3,3),
                                                                  c(NA,4,4,5,5,NA)))}
  if(model=="type9"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              output[[4]]+labs(title="D)")+fonte,
                                              output[[5]]+labs(title="E)")+fonte,
                                              output[[6]]+labs(title="F)")+fonte,
                                              layout_matrix=rbind(c(1,1,2,2),
                                                                  c(3,3,4,4),
                                                                  c(5,5,6,6)))}
  if(model=="type10"){arrangeplot=grid.arrange(output[[1]]+labs(title="A)")+fonte,
                                              output[[2]]+labs(title="B)")+fonte,
                                              output[[3]]+labs(title="C)")+fonte,
                                              output[[4]]+labs(title="D)")+fonte,
                                              output[[5]]+labs(title="E)")+fonte,
                                              output[[6]]+labs(title="F)")+fonte,
                                              layout_matrix=rbind(c(1,1,2,2,3,3),
                                                                  c(4,4,5,5,6,6)))}
}

