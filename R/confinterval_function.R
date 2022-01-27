#' Utils: Interval of confidence for groups
#'
#' @description Calculates confidence interval for groups
#' @param resp numeric vector with responses
#' @param group vector with groups or list with two factors
#' @param alpha confidence level of the interval
#' @param type lower or upper range
#' @export
#' @return returns a numeric vector with confidence interval grouped by treatment.
#' @examples
#'
#' #===================================
#' # One factor
#' #===================================
#'
#' dados=rnorm(100,10,1)
#' trat=rep(paste("T",1:10),10)
#' confinterval(dados,trat)
#'
#' #===================================
#' # Two factor
#' #===================================
#' f1=rep(c("A","B"),e=50)
#' f2=rep(paste("T",1:5),e=10,2)
#' confinterval(dados,list(f1,f2))

confinterval=function(resp,
                      group,
                      alpha=0.95,
                      type="upper"){
  if(is.list(group)==FALSE){
    lower=c()
    upper=c()
    for(i in 1:length(unique(group))){
    group=factor(group,unique(group))
    ic=t.test(resp[group==levels(group)[i]])
    lower[i]=ic$conf.int[1]
    upper[i]=ic$conf.int[2]}
    names(lower)=levels(group)
    names(upper)=levels(group)}
  if(is.list(group)==TRUE){
    f1=group[[1]]
    f2=group[[2]]
    group=paste(f1,f2)
    f1=factor(f1,unique(f1))
    f2=factor(f2,unique(f2))
    lower=c()
    upper=c()
    for(i in 1:length(unique(group))){
    group=factor(group,unique(group))
    ic=t.test(resp[group==levels(group)[i]])
    lower[i]=ic$conf.int[1]
    upper[i]=ic$conf.int[2]}
  lower=matrix(lower,nrow=length(levels(f2)))
  colnames(lower)=levels(f1)
  rownames(lower)=levels(f2)
  upper=matrix(upper,nrow=length(levels(f2)))
  colnames(upper)=levels(f1)
  rownames(upper)=levels(f2)}
  if(type=="upper"){saida=upper}
  if(type=="lower"){saida=lower}
  saida}
