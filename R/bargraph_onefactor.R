#' Graph: Group DIC, DBC and DQL functions column charts
#'
#' @description Groups two or more column charts exported from DIC, DBC or DQL function
#' @param analysis List with DIC, DBC or DQL object
#' @param labels Vector with the name of the facets
#' @param ocult.facet Hide facets
#' @param ocult.box Hide box
#' @param facet.size Font size facets
#' @param ylab Y-axis name
#' @param width.bar Width error bar
#' @param width.col Width Column
#' @param sup Number of units above the standard deviation or average bar on the graph
#'
#' @return Returns a column chart grouped by facets
#'
#' @export
#' @examples
#' library(AgroR)
#' data("laranja")
#' a=with(laranja, DBC(trat, bloco, resp, ylab = "Number of fruits/plants"))
#' b=with(laranja, DBC(trat, bloco, resp,  ylab = "Number of fruits/plants"))
#' c=with(laranja, DBC(trat, bloco, resp, ylab = "Number of fruits/plants"))
#' bargraph_onefactor(analysis = list(a,b,c), labels = c("One","Two","Three"),ocult.box = TRUE)

bargraph_onefactor=function(analysis,
                            labels=NULL,
                            ocult.facet=FALSE,
                            ocult.box=FALSE,
                            facet.size=14,
                            ylab=NULL,
                            width.bar=0.3,
                            width.col=0.9,
                            sup=NULL){
  requireNamespace("ggplot2")
  results=as.list(1:length(analysis))
  for(i in 1:length(analysis)){
    if(is.null(labels)==TRUE){analysis[[i]][[1]]$plot$dadosm$facet=rep(i,
                                                                       e=nrow(analysis[[i]][[1]]$plot$dadosm))}else{
                                                                         analysis[[i]][[1]]$plot$dadosm$facet=rep(labels[i],e=nrow(analysis[[i]][[1]]$plot$dadosm))}
    results[[i]]=analysis[[i]][[1]]$plot$dadosm}
  tabela=do.call("rbind",results)
  if(is.null(sup)==TRUE){sup=0.1*mean(tabela$media)}
  media=tabela$media
  desvio=tabela$desvio
  trats=tabela$trats
  letra=tabela$letra
  graph=ggplot(tabela,aes(y=media,x=trats))+
    geom_col(color="black",fill="lightblue",width = width.col)+
    facet_grid(~facet,scales = "free", space='free')+
    geom_errorbar(aes(ymin=media-desvio,ymax=media+desvio),width=width.bar)+
    geom_text(aes(y=media+desvio+sup,x=trats,label=letra))+xlab("")+
    analysis[[1]][[1]]$theme+theme(strip.text = element_text(size=facet.size))
  if(is.null(ylab)==TRUE){graph=graph+ylab(analysis[[1]][[1]]$plot$ylab)}else{graph=graph+ylab(ylab)}
  if(ocult.facet==TRUE){graph=graph+theme(strip.text = element_blank())}
  if(ocult.box==TRUE){graph=graph+theme(strip.background = element_blank())}
  list(graph)[[1]]}

