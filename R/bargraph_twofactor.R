#' Graph: Group FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC functions column charts
#'
#' @description Groups two or more column charts exported from FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC function
#' @param analysis List with DIC, DBC or DQL object
#' @param labels Vector with the name of the facets
#' @param ocult.facet Hide facets
#' @param ocult.box Hide box
#' @param facet.size Font size facets
#' @param ylab Y-axis name
#' @param width.bar Width bar
#' @param sup Number of units above the standard deviation or average bar on the graph
#'
#' @return Returns a column chart grouped by facets
#'
#' @export
#' @examples
#' library(AgroR)
#' data(corn)
#' a=with(corn, FAT2DIC(A, B, Resp, quali=c(TRUE, TRUE),ylab="Heigth (cm)"))
#' b=with(corn, FAT2DIC(A, B, Resp, mcomp="sk", quali=c(TRUE, TRUE),ylab="Heigth (cm)"))
#' bargraph_twofactor(analysis = list(a,b), labels = c("One","Two"),ocult.box = TRUE)

bargraph_twofactor=function(analysis,
                       labels=NULL,
                       ocult.facet=FALSE,
                       ocult.box=FALSE,
                       facet.size=14,
                       ylab=NULL,
                       width.bar=0.3,
                       sup=NULL){
  requireNamespace("ggplot2")
  results=as.list(1:length(analysis))
  for(i in 1:length(analysis)){
    if(is.null(labels)==TRUE){analysis[[i]][[2]]$plot$graph$facet=rep(i,
                                                                      e=nrow(analysis[[i]][[2]]$plot$graph))}else{
                                                                        analysis[[i]][[2]]$plot$graph$facet=rep(labels[i],e=nrow(analysis[[i]][[2]]$plot$graph))}
    results[[i]]=analysis[[i]][[2]]$plot$graph}
  tabela=do.call("rbind",results)
  if(is.null(sup)==TRUE){sup=0.1*mean(tabela$media)}

  media=tabela$media
  desvio=tabela$desvio
  f1=tabela$f1
  f2=tabela$f2
  numero=tabela$numero
  letra=tabela$letra
  graph=ggplot(tabela,aes(y=media,x=f1,fill=f2))+
    geom_col(color="black",position = position_dodge(width = 0.9))+
    facet_grid(~facet,scales = "free", space='free')+
    geom_errorbar(aes(ymin=media-desvio,ymax=media+desvio),width=width.bar,position = position_dodge(width=0.9))+
    geom_text(aes(y=media+desvio+sup,x=f1,label=numero),position = position_dodge(width=0.9))+xlab("")+
    analysis[[1]][[2]]$theme+theme(strip.text = element_text(size=facet.size))
  if(is.null(ylab)==TRUE){graph=graph+ylab(analysis[[1]][[2]]$plot$ylab)}else{graph=graph+ylab(ylab)}
  if(ocult.facet==TRUE){graph=graph+theme(strip.text = element_blank())}
  if(ocult.box==TRUE){graph=graph+theme(strip.background = element_blank())}
  list(graph)[[1]]}
