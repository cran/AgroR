#' Utils: Dunnett's Test Summary
#' @export
#' @description Performs a summary in table form from a list of Dunnett's test outputs
#' @param variable List object Dunnett test
#' @param colnames Names of column
#' @param info Information of table
#' @return A summary table from Dunnett's test is returned
#'
#' @examples
#' library(AgroR)
#' data("pomegranate")
#' a=with(pomegranate,dunnett(trat=trat,resp=WL,control="T1"))
#' b=with(pomegranate,dunnett(trat=trat,resp=SS,control="T1"))
#' c=with(pomegranate,dunnett(trat=trat,resp=AT,control="T1"))
#' d=with(pomegranate,dunnett(trat=trat,resp=ratio,control="T1"))
#' summarise_dunnett(list(a,b,c,d))

summarise_dunnett=function(variable, colnames=NA, info="sig"){
  variaveis=length(variable)
  nomes=rownames(variable[[1]]$data)
  if(is.na(colnames[1])==TRUE){colnames=rep(paste("Var",1:variaveis))}
  datas=data.frame(matrix(rep(NA,variaveis*length(nomes)),ncol=variaveis))
  for(i in 1:variaveis){
    if(info=="sig"){datas[,i]=variable[[i]]$data$sig}
    if(info=="t"){datas[,i]=variable[[i]]$data$t.value}
    if(info=="p"){datas[,i]=variable[[i]]$data$p.value}
    if(info=="estimate"){datas[,i]=variable[[i]]$data$Estimate}
    if(info=="IC.lower"){datas[,i]=variable[[i]]$data$IC.lwr}
    if(info=="IC.upper"){datas[,i]=variable[[i]]$data$IC.upr}
  }
  rownames(datas)=nomes
  colnames(datas)=colnames
  datas
}
