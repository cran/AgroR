#' Analysis: Method to evaluate similarity of experiments based on QMres
#'
#' @description This function presents a method to evaluate similarity of experiments based on a matrix of QMres of all against all. This is used as a measure of similarity and applied in clustering.
#' @param qmres Vector containing mean squares of residuals or output from list DIC or DBC function
#' @param information Option to choose the return type. `matrix`, `bar` or `cluster`
#' @param method.cluster Grouping method
#' @return Returns a residual mean square ratio matrix, bar graph with ratios sorted in ascending order, or cluster analysis.
#' @keywords Joint analysis
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @export
#' @examples
#' qmres=c(0.344429, 0.300542, 0.124833, 0.04531, 0.039571, 0.011812, 0.00519)
#' jointcluster(qmres,information = "cluster")
#' jointcluster(qmres,information = "matrix")
#' jointcluster(qmres,information = "bar")
#'
#' data(mirtilo)
#' m=lapply(unique(mirtilo$exp),function(x){
#'   m=with(mirtilo[mirtilo$exp==x,],DBC(trat,bloco,resp))})
#' jointcluster(m)

jointcluster=function(qmres,
                      information="matrix",
                      method.cluster="ward.D"){
  if(is.list(qmres)==TRUE){
    if(nrow(qmres[[1]][[1]]$plot$anava)==2){
      qmres1=numeric(0)
      for(i in 1:length(qmres)){
        qmres1[i]=qmres[[i]][[1]]$plot$anava$QM[2]}}
    if(nrow(qmres[[1]][[1]]$plot$anava)==3){
      qmres1=numeric(0)
      for(i in 1:length(qmres)){
        qmres1[i]=qmres[[i]][[1]]$plot$anava$QM[3]}}
    qmres=qmres1}
  resp=qmres
  requireNamespace("ggplot2")
  matriztodos=function(resp){
  expe=paste("Exp",1:length(resp))
  dados=expand.grid(expe,expe)
  Var1=dados$Var1
  Var2=dados$Var2
  dados$resp1=rep(resp,length(resp))
  dados$resp2=rep(resp,e=length(resp))
  dados$resp=dados$resp1/dados$resp2
  graph=ggplot(dados,aes(x=Var1,y=Var2,fill=resp))+
    geom_tile(color = "gray50", linewidth = 1) +
    scale_x_discrete(position = "top") +
    scale_fill_distiller(palette = "RdBu", direction = -1) +
    ylab("Numerator") + xlab("Denominator") +
    geom_label(aes(label = format(resp, digits = 2)), fill = "white") + labs(fill = "ratio") +
    theme(axis.text = element_text(size = 12, color = "black"),
          legend.text = element_text(size = 12), axis.ticks = element_blank(),
          panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  graph}
  matrizmaiores=function(resp){
    expe=paste("Exp",1:length(resp))
    dados=expand.grid(expe,expe)
    dados$resp1=rep(resp,length(resp))
    dados$resp2=rep(resp,e=length(resp))
    dados$respAB=dados$resp1/dados$resp2
    dados$respBA=dados$resp2/dados$resp1
    dados$respAB[dados$respAB>dados$respBA]
    dados
    n=c()
    for(i in 2:length(resp)){
      n[[i-1]]=rep(i:length(resp))+length(resp)*(i-2)}
    unlist(n)
    dados=dados[unlist(n),]
    dados$maior=pmax(dados$respAB,dados$respBA)
    dados$combinacao=paste(dados$Var1,dados$Var2)
    combinacao=dados$combinacao
    maior=dados$maior
    datam=dados[order(dados$maior),]
    datam$combinacao=factor(datam$combinacao,datam$combinacao)
    graph=ggplot(datam,
           aes(y=combinacao,x=maior))+
      geom_col(fill="lightblue",color="black")+
      geom_vline(xintercept = 7,lty=1)+
      labs(y="Combination",x="Ratio MSr/MSr")+
      theme(axis.text = element_text(size = 12, color = "black"),
            legend.text = element_text(size = 12), axis.ticks = element_blank(),
            panel.background = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    graph}
  similar=function(resp,k=2,method){
    expe=paste("Exp",1:length(resp))
    dados=expand.grid(expe,expe)
    dados$resp1=rep(resp,length(resp))
    dados$resp2=rep(resp,e=length(resp))
    dados$respAB=dados$resp1/dados$resp2
    dados$respBA=dados$resp2/dados$resp1
    dados$maior=pmax(dados$respAB,dados$respBA)
    matriz=matrix(dados$maior,ncol=length(resp))
    rownames(matriz)=paste("Exp",1:length(resp))
    colnames(matriz)=paste("Exp",1:length(resp))
    matriz=hclust(d = as.dist(matriz),method = method)
    plot(matriz,ylab="Ratio MSr/MSr",main="Cluster experiment",xlab="")}
  if(information=="matrix"){print(matriztodos(qmres))}
  if(information=="bar"){print(matrizmaiores(qmres))}
  if(information=="cluster"){similar(qmres,method = method.cluster)}
  }


