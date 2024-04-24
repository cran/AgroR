#' Utils: Summary of Analysis of Variance and Test of Means for Joint analysis
#' @description Summarizes the output of the analysis of variance and the multiple comparisons test for completely randomized (DIC) and randomized block (DBC) designs for Joint analysis with qualitative factor.
#' @author Gabriel Danilo Shimizu
#' @param analysis List with the analysis outputs of the conjdic and conjdbc functions
#' @param design Type of experimental project (DIC or DBC)
#' @param info Analysis of variance information (can be "p", "f", "QM" or "SQ")
#' @note The column names in the final output are imported from the ylab argument within each function.
#' @note This function is only for declared qualitative factors. In the case of a quantitative factor and the other qualitative in projects with two factors, this function will not work.
#' @import knitr
#' @export
#' @examples
#' library(AgroR)
#' data(mirtilo)
#' set.seed(1); resp1=rnorm(36,10,4)
#' set.seed(4); resp2=rnorm(36,10,3)
#' set.seed(8); resp3=rnorm(36,100,40)
#' type1=with(mirtilo, conjdbc(trat, bloco, exp, resp, ylab = "var1"))
#' type2=with(mirtilo, conjdbc(trat, bloco, exp, resp1, ylab = "var2"))
#' type3=with(mirtilo, conjdbc(trat, bloco, exp, resp2, ylab = "var3"))
#' type4=with(mirtilo, conjdbc(trat, bloco, exp, resp3, ylab = "var4"))
#' summarise_conj(analysis = list(type1,type2,type3,type4))

summarise_conj=function(analysis,design="DBC",info="p"){
  if(design=="DIC"){
    qmres=c()
    nomes=c()
    df1=c();ss1=c();qm1=c();fv1=c();pvalue1=c()
    df2=c();ss2=c();qm2=c();fv2=c();pvalue2=c()
    df3=c();ss3=c();qm3=c();fv3=c();pvalue3=c()
    df4=c();ss4=c();qm4=c();fv4=c();pvalue4=c()
    comp=data.frame(matrix(rep(rep("-",
                                   length(unique(analysis[[1]][[1]]$plot$trat))),length(analysis)),
                           ncol = length(analysis)))
    for(i in 1:length(analysis)){
      nomes[i]=analysis[[i]][[1]]$plot$ylab
      qmres[i]=analysis[[i]][[1]]$plot$qmresmedio
      ano=analysis[[i]][[1]]$plot$datas

      ##########################################
      # teste medias geral
      if(qmres[i]<analysis[[i]][[1]]$plot$homog.value & ano$`Pr(>F)`[3]>analysis[[i]][[1]]$plot$alpha.f){
        compa=analysis[[i]][[1]]$plot$dadosm[,"letra"]}else{
          compa=rep("-",e=length(unique(analysis[[i]][[1]]$plot$trat)))
        }
      comp[,i]=compa

      df1[i]=ano$Df[1];df2[i]=ano$Df[2];df3[i]=ano$Df[3];df4[i]=ano$Df[4]
      ss1[i]=ano$`Sum Sq`[1];ss2[i]=ano$`Sum Sq`[2];ss3[i]=ano$`Sum Sq`[3];ss4[i]=ano$`Sum Sq`[4]
      qm1[i]=ano$`Mean Sq`[1];qm2[i]=ano$`Mean Sq`[2];qm3[i]=ano$`Mean Sq`[3];qm4[i]=ano$`Mean Sq`[4]
      fv1[i]=ano$`F value`[1];fv2[i]=ano$`F value`[2];fv3[i]=ano$`F value`[3]
      pvalue1[i]=ano$`Pr(>F)`[1];pvalue2[i]=ano$`Pr(>F)`[2];pvalue3[i]=ano$`Pr(>F)`[3]}
    #"p", "f", "QM" or "SQ"
    if(info=="p"){data=rbind(pvalue1, pvalue2, pvalue3, round(qmres,2))
    rownames(data)=c("Trat","Exp","Exp:Trat","Ratio QMres")}
    if(info=="f"){data=rbind(fv1,fv2,fv3,round(qmres,2))
    rownames(data)=c("Trat","Exp","Exp:Trat","Ratio QMres")}
    if(info=="QM"){data=rbind(qm1, qm2, qm3, qm4,round(qmres,2))
    rownames(data)=c("Trat","Exp","Exp:Trat","Average residue","Ratio QMres")}
    if(info=="SQ"){data=rbind(ss1, ss2, ss3, ss4,round(qmres,2))
    rownames(data)=c("Trat","Exp","Exp:Trat","Average residue","Ratio QMres")}
    data=data.frame(data)
    colnames(data)=nomes
    rownames(comp)=rownames(analysis[[1]][[1]]$plot$dadosm)
    colnames(comp)=nomes
  }
  if(design=="DBC"){
    qmres=c()
    nomes=c()
    df1=c();ss1=c();qm1=c();fv1=c();pvalue1=c()
    df2=c();ss2=c();qm2=c();fv2=c();pvalue2=c()
    df3=c();ss3=c();qm3=c();fv3=c();pvalue3=c()
    df4=c();ss4=c();qm4=c();fv4=c();pvalue4=c()
    df5=c();ss5=c();qm5=c();fv5=c();pvalue5=c()
    comp=data.frame(matrix(rep(rep("-",
                   length(unique(analysis[[1]][[1]]$plot$trat))),length(analysis)),
                   ncol = length(analysis)))
    for(i in 1:length(analysis)){
    nomes[i]=analysis[[i]][[1]]$plot$ylab
    qmres[i]=analysis[[i]][[1]]$plot$qmresmedio
    ano=analysis[[i]][[1]]$plot$datas

    ##########################################
    # teste medias geral
    if(qmres[i]<analysis[[i]][[1]]$plot$homog.value & ano$`Pr(>F)`[4]>analysis[[i]][[1]]$plot$alpha.f){
      compa=analysis[[i]][[1]]$plot$dadosm[,"letra"]}else{
        compa=rep("-",e=length(unique(analysis[[i]][[1]]$plot$trat)))
      }
    comp[,i]=compa
    ##########################################
    # teste de medias especifico
    # if(qmres[i]>7 | ano$`Pr(>F)`[4]<0.05){
    #   compa1=analysis[[i]][[1]]$plot$tukey
    #   nexp=length(compa1)
    #   compara1=as.list(1:nexp)
    #   for(j in 1:nexp){compara1[j]=compa1[j]}}
    # comp[,i]=compa
    ##########################################

    df1[i]=ano$Df[1];df2[i]=ano$Df[2];df3[i]=ano$Df[3];df4[i]=ano$Df[4];df5[i]=ano$Df[5]
    ss1[i]=ano$`Sum Sq`[1];ss2[i]=ano$`Sum Sq`[2];ss3[i]=ano$`Sum Sq`[3];ss4[i]=ano$`Sum Sq`[4];ss5[i]=ano$`Sum Sq`[5]
    qm1[i]=ano$`Mean Sq`[1];qm2[i]=ano$`Mean Sq`[2];qm3[i]=ano$`Mean Sq`[3];qm4[i]=ano$`Mean Sq`[4];qm5[i]=ano$`Mean Sq`[5]
    fv1[i]=ano$`F value`[1];fv2[i]=ano$`F value`[2];fv3[i]=ano$`F value`[3];fv4[i]=ano$`F value`[4];fv5[i]=ano$`F value`[5]
    pvalue1[i]=ano$`Pr(>F)`[1];pvalue2[i]=ano$`Pr(>F)`[2];pvalue3[i]=ano$`Pr(>F)`[3];pvalue4[i]=ano$`Pr(>F)`[4];pvalue5[i]=ano$`Pr(>F)`[5]}
    #"p", "f", "QM" or "SQ"
    if(info=="p"){data=rbind(pvalue1, pvalue2, pvalue3, pvalue4, round(qmres,2))
                  rownames(data)=c("Trat","Exp","Block/Exp","Exp:Trat","Ratio QMres")}
    if(info=="f"){data=rbind(fv1,fv2,fv3,fv4,round(qmres,2))
    rownames(data)=c("Trat","Exp","Block/Exp","Exp:Trat","Ratio QMres")}
    if(info=="QM"){data=rbind(qm1, qm2, qm3, qm4, qm5,round(qmres,2))
    rownames(data)=c("Trat","Exp","Block/Exp","Exp:Trat","Average residue","Ratio QMres")}
    if(info=="SQ"){data=rbind(ss1, ss2, ss3, ss4, ss5,round(qmres,2))
    rownames(data)=c("Trat","Exp","Block/Exp","Exp:Trat","Average residue","Ratio QMres")}
    data=data.frame(data)
    colnames(data)=nomes
    rownames(comp)=rownames(analysis[[1]][[1]]$plot$dadosm)
    colnames(comp)=nomes
    # list("Anova"=data,
    #      "Post-Hoc"=comp)
  }
  list("Anova"=data,
       "Post-Hoc"=comp)
}

