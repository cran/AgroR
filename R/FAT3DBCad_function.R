#' Analysis: DBC experiments in triple factorial with aditional
#'
#' @description Analysis of an experiment conducted in a randomized block design in a triple factorial scheme with one aditional control using analysis of variance of fixed effects.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param f3 Numeric or complex vector with factor 3 levels
#' @param block Numerical or complex vector with blocks
#' @param response Numerical vector containing the response of the experiment.
#' @param responseAd Numerical vector containing the aditional response
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{qualitative})
#' @param names.fat Allows labeling the factors 1, 2 and 3.
#' @param grau Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with three elements.
#' @param grau12 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 2, in the case of interaction f1 x f2 and qualitative factor 2 and quantitative factor 1.
#' @param grau21 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 1, in the case of interaction f1 x f2 and qualitative factor 1 and quantitative factor 2.
#' @param grau13 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 3, in the case of interaction f1 x f3 and qualitative factor 3 and quantitative factor 1.
#' @param grau31 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 1, in the case of interaction f1 x f3 and qualitative factor 1 and quantitative factor 3.
#' @param grau23 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 3, in the case of interaction f2 x f3 and qualitative factor 3 and quantitative factor 2.
#' @param grau32 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 2, in the case of interaction f2 x f3 and qualitative factor 2 and quantitative factor 3.
#' @param grau123 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 1, in the case of interaction f1 x f2 x f3 and quantitative factor 1.
#' @param grau213 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 2, in the case of interaction f1 x f2 x f3 and quantitative factor 2.
#' @param grau312 Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 3, in the case of interaction f1 x f2 x f3 and quantitative factor 3.
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab.factor Provide a vector with two observations referring to the x-axis name of factors 1, 2 and 3, respectively, when there is an isolated effect of the factors. This argument uses `parse`.
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param transf Applies data transformation (\emph{default} is 1; for log consider 0; `angular` for angular transformation)
#' @param constant Add a constant for transformation (enter value)
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param geom Graph type (columns or segments)
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angulo x-axis scale text rotation
#' @param textsize Font size
#' @param labelsize Label size
#' @param dec Number of cells
#' @param family Font family
#' @param ad.label Aditional label
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param point This function defines whether the point must have all points ("all"), mean ("mean"), standard deviation (\emph{default} - "mean_sd") or mean with standard error ("mean_se") if quali= FALSE. For quali=TRUE, `mean_sd` and `mean_se` change which information will be displayed in the error bar.
#' @param angle.label label angle
#' @note The order of the chart follows the alphabetical pattern. Please use `scale_x_discrete` from package ggplot2, `limits` argument to reorder x-axis. The bars of the column and segment graphs are standard deviation.
#' @return The analysis of variance table, the Shapiro-Wilk error normality test, the Bartlett homogeneity test of variances, the Durbin-Watson error independence test, multiple comparison test (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. The column chart for qualitative treatments is also returned.For significant triple interaction only, no graph is returned.
#' @note The function does not perform multiple regression in the case of two or more quantitative factors. The bars of the column and segment graphs are standard deviation.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Practical Nonparametrics Statistics. W.J. Conover, 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#'
#' Ferreira, E. B., Cavalcanti, P. P., and Nogueira, D. A. (2014). ExpDes: an R package for ANOVA and experimental designs. Applied Mathematics, 5(19), 2952.
#'
#' Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#'
#' @keywords DIC
#' @keywords Factorial
#' @export
#' @examples
#' library(AgroR)
#' data(enxofre)
#' respAd=c(2000,2400,2530,2100)
#' attach(enxofre)
#' with(enxofre, FAT3DBC.ad(f1, f2, f3, bloco, resp, respAd))

FAT3DBC.ad = function(f1,
                      f2,
                      f3,
                      block,
                      response,
                      responseAd,
                      norm = "sw",
                      alpha.f = 0.05,
                      alpha.t = 0.05,
                      quali = c(TRUE, TRUE, TRUE),
                      mcomp = 'tukey',
                      transf = 1,
                      constant = 0,
                      names.fat = c("F1", "F2", "F3"),
                      ylab = "Response",
                      xlab = "",
                      xlab.factor=c("F1","F2","F3"),
                      sup = NA,
                      grau=c(NA,NA,NA), # isolado e interação tripla
                      grau12=NA, # F1/F2
                      grau13=NA, # F1/F3
                      grau23=NA, # F2/F3
                      grau21=NA, # F2/F1
                      grau31=NA, # F3/F1
                      grau32=NA, # F3/F2
                      grau123=NA,
                      grau213=NA,
                      grau312=NA,
                      fill = "lightblue",
                      theme = theme_classic(),
                      ad.label="Additional",
                      angulo = 0,
                      errorbar = TRUE,
                      addmean = TRUE,
                      family = "sans",
                      dec = 3,
                      geom = "bar",
                      textsize = 12,
                      labelsize=4,
                      point="mean_sd",
                      angle.label = 0) {
  if(is.na(sup==TRUE)){sup=0.2*mean(response)}
  if(angle.label==0){hjust=0.5}else{hjust=0}
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")
  if(transf==1){resp=response+constant}else{if(transf!="angular"){resp=((response+constant)^transf-1)/transf}}
  # if(transf==1){resp=response+constant}else{resp=((response+constant)^transf-1)/transf}
  if(transf==0){resp=log(response+constant)}
  if(transf==0.5){resp=sqrt(response+constant)}
  if(transf==-0.5){resp=1/sqrt(response+constant)}
  if(transf==-1){resp=1/(response+constant)}
  if(transf=="angular"){resp=asin(sqrt((response+constant)/100))}

  if(transf==1){respAd=responseAd+constant}else{if(transf!="angular"){respAd=((responseAd+constant)^transf-1)/transf}}
  # if(transf==1){respAd=responseAd+constant}else{respAd=((responseAd+constant)^transf-1)/transf}
  if(transf==0){respAd=log(responseAd+constant)}
  if(transf==0.5){respAd=sqrt(responseAd+constant)}
  if(transf==-0.5){respAd=1/sqrt(responseAd+constant)}
  if(transf==-1){respAd=1/(responseAd+constant)}
  if(transf=="angular"){respAd=asin(sqrt((responseAd+constant)/100))}

  resp1=resp

  ordempadronizado=data.frame(f1,f2,f3,block,response,resp)
  organiz=data.frame(f1,f2,f3,block,response,resp)
  organiz=organiz[order(organiz$block),]
  organiz=organiz[order(organiz$f3),]
  organiz=organiz[order(organiz$f2),]
  organiz=organiz[order(organiz$f1),]
  f1=organiz$f1
  f2=organiz$f2
  f3=organiz$f3
  block=organiz$block
  response=organiz$response
  resp=organiz$resp
  fator1=f1
  fator2=f2
  fator3=f3
  fator1a=fator1
  fator2a=fator2
  fator3a=fator3
  bloco=block
  fac.names=names.fat
  fatores<-data.frame(fator1,fator2,fator3)
  Fator1<-factor(fator1,levels=unique(fator1));
  Fator2<-factor(fator2,levels=unique(fator2));
  Fator3<-factor(fator3,levels=unique(fator3))
  nv1<-length(summary(Fator1)); nv2<-length(summary(Fator2)); nv3<-length(summary(Fator3))
  J<-(length(resp))/(nv1*nv2*nv3)
  lf1<-levels(Fator1); lf2<-levels(Fator2); lf3<-levels(Fator3)
  bloco=as.factor(bloco)

  # =================================
  ## Anova
  # =================================
  anava<-aov(resp~Fator1*Fator2*Fator3+bloco)
  anavaF3<-anova(anava)
  anovaF3=anavaF3
  colnames(anovaF3)=c("GL","SQ","QM","Fcal","p-value")

  anavares<-aov(resp~as.factor(f1)*
                  as.factor(f2)*
                  as.factor(f3)+
                  as.factor(block),data = ordempadronizado)
  respad=anava$residuals/sqrt(anavaF3$`Mean Sq`[9])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}
  Ids=ifelse(respad>3 | respad<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(respad,Ids),aes(y=respad,x=1:length(respad)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(respad),label=1:length(respad),color=Ids,size=labelsize)+
    scale_x_continuous(breaks=1:length(respad))+
    theme_classic()+theme(axis.text.y = element_text(size=textsize),
                          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)
  print(residplot)

  col1<-numeric(0)
  for(i in 1:c(nv1*nv2*nv3)) {
    col1<-c(col1, rep(i,J))
  }
  col1<-c(col1,rep('ad',J))
  col2<-c(bloco,rep(1:J))
  col3<-c(resp,respAd)
  tabF3ad<-data.frame("TRAT"=col1, "BLOCO"=col2, "RESP2"=col3)
  TRAT<-factor(tabF3ad[,1])
  BLOCO<-factor(tabF3ad[,2])
  anava1<-aov(tabF3ad[,3]~ BLOCO + TRAT)
  anavaTr<-anova(anava1)
  SQAd=anavaTr[2,2]-sum(anavaF3[c(1,2,3,5,6,7,8),2])
  DfAd=1
  QMAd=SQAd
  FAd=QMAd/anavaTr[3,3]
  pAd=1-pf(FAd,DfAd,anavaTr[3,3])
  DfE=anavaTr[3,1]
  SQE=anavaTr[3,2]
  QME=anavaTr[3,3]
  Fef=anavaF3$`Mean Sq`[1:8]/anavaTr[3,3]
  pef=1-pf(Fef,anavaF3$Df[1:8],DfE)
  anavaF3=rbind(anavaF3[1:8,],
        "Factorial vs Aditional"=c(DfAd,SQAd,QMAd,FAd,pAd),
        "Residuals"=c(DfE,SQE,QME,NA,NA))
  anavaF3[1:8,4]=Fef
  anavaF3[1:8,5]=pef
  anavaF3[4,]=anavaTr[1,]
  anovaF3=anavaF3

  # =================================
  ## Saída inicial
  # =================================

  #Teste de normalidade
  norm1<-shapiro.test(anava$residuals)
  cat(green(bold("\n------------------------------------------")))
  cat(green(bold("\nNormality of errors")))
  cat(green(bold("\n------------------------------------------")))
  print(norm1)
  message(if(norm1$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors do not follow a normal distribution"})

  # Teste de homogeneidade das variâncias
  homog1=bartlett.test(anava$residuals~paste(Fator1,Fator2,Fator3))
  cat(green(bold("\n------------------------------------------")))
  cat(green(bold("\nHomogeneity of Variances")))
  cat(green(bold("\n------------------------------------------")))
  print(homog1)
  message(if(homog1$p.value[1]>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})

  # Independencia dos erros
  indep=dwtest(anava)
  cat(green(bold("\n------------------------------------------")))
  cat(green(bold("\nIndependence from errors")))
  cat(green(bold("\n------------------------------------------")))
  print(indep)
  message(if(indep$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors are not independent"})
  cat("\n")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV (%) = ",round(sqrt(anavaF3$`Mean Sq`[10])/mean(c(resp,respAd),na.rm=TRUE)*100,2)))
  cat(paste("\nMean Factorial = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian Factorial = ",round(median(response,na.rm=TRUE),4)))
  cat(paste("\nMean Aditional = ",round(mean(responseAd,na.rm=TRUE),4)))
  cat(paste("\nMedian Aditional = ",round(median(responseAd,na.rm=TRUE),4)))
  cat("\nPossible outliers = ", out)
  cat("\n")
  cat(green(bold("\n------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n------------------------------------------\n")))
  anava1=as.matrix(data.frame(anovaF3))
  colnames(anava1)=c("Df","Sum Sq","Mean.Sq","F value","Pr(F)" )
  rownames(anava1)=c(names.fat[1],
                     names.fat[2],
                     names.fat[3],
                     "Block",
                     paste(names.fat[1],"\u00D7",names.fat[2]),
                     paste(names.fat[1],"\u00D7",names.fat[3]),
                     paste(names.fat[2],"\u00D7",names.fat[3]),
                     paste(names.fat[1],"\u00D7",names.fat[2],"\u00D7",names.fat[3]),
                     "Factorial \u00D7 Aditional",
                     "Residuals")
  print(anava1,na.print = "")
  cat("\n")

  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
    message("\n Your analysis is not valid, suggests using a non-parametric test and try to transform the data\n")}else{}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){
    message("\n Your analysis is not valid, suggests using the function FATDIC.art\n")}else{}
  message(if(transf !=1){blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})

  if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[8,5]>alpha.f) {
    graficos=list(1,2,3)
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold('Non-significant interaction: analyzing the simple effects')))
    cat(green(bold("\n------------------------------------------\n")))
    fatores<-data.frame('fator 1'=fator1,'fator 2' = fator2,'fator 3' = fator3)

    for(i in 1:3){
      if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))
        if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],
                                       anavaF3[10,1],anavaF3[10,3],alpha.t)
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="sk"){
          nrep=table(fatores[,i])[1]
          medias=sort(tapply(resp,fatores[i],mean, na.rm=TRUE),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3[10,1],
                        nrep = nrep,
                        QME = anavaF3[10,3],
                        alpha = alpha.t)
          letra1=data.frame(resp=medias,groups=sk)
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="duncan"){
          ad=data.frame(Fator1,Fator2,Fator3)
          letra <- duncan(resp,fatores[,i],
                          anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="lsd"){
          ad=data.frame(Fator1,Fator2,Fator3)
          letra <- LSD(resp,fatores[,i],
                       anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        teste=if(mcomp=="tukey"){"Tukey HSD"}else{
          if(mcomp=="sk"){"Scott-Knott"}else{
            if(mcomp=="lsd"){"LSD-Fischer"}else{
              if(mcomp=="duncan"){"Duncan"}}}}
        cat(green(italic(paste("Multiple Comparison Test:",teste,"\n"))))
        print(letra1)
        cat(green(bold("\n------------------------------------------\n")))
        #=====================================================
        if(point=="mean_sd"){desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)]}
        if(point=="mean_se"){desvio=(tapply(response, c(fatores[i]), sd, na.rm=TRUE)/
                                       sqrt(tapply(response, c(fatores[i]), length)))[rownames(letra1)]}
        dadosm=data.frame(letra1,
                          media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                          desvio=desvio)
        dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
        dadosm$limite=dadosm$media+dadosm$desvio
        dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]
        if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
        if(addmean==FALSE){dadosm$letra=dadosm$groups}
        media=dadosm$media
        desvio=dadosm$desvio
        Tratamentos=dadosm$Tratamentos
        letra=dadosm$letra
        if(geom=="bar"){grafico=ggplot(dadosm,
                                       aes(x=Tratamentos,
                                           y=media))
        if(fill=="trat"){grafico=grafico+
          geom_col(aes(fill=Tratamentos),
                   color=1)}else{grafico=grafico+
                     geom_col(aes(fill=Tratamentos),
                              fill=fill,color=1)}
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=media+sup+
                          if(sup<0){-desvio}else{desvio},
                        label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
        if(errorbar==FALSE){grafico=grafico+
          geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,
                        aes(ymin=media-desvio,
                            ymax=media+desvio,color=1),
                        color="black",width=0.3)}
        grafico=grafico+theme+
          ylab(ylab)+
          xlab(parse(text = xlab.factor[i]))+
          theme(text = element_text(size=textsize,color="black", family = family),
                axis.text = element_text(size=textsize,color="black", family = family),
                axis.title = element_text(size=textsize,color="black", family = family))}

        if(geom=="point"){grafico=ggplot(dadosm,
                                         aes(x=Tratamentos,
                                             y=media))
        if(fill=="trat"){grafico=grafico+
          geom_point(aes(color=Tratamentos))}
        else{grafico=grafico+
          geom_point(aes(color=Tratamentos),color=fill,size=4)}
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                        label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
        if(errorbar==FALSE){grafico=grafico+
          geom_text(aes(y=media+sup,
                        label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,
                        aes(ymin=media-desvio,
                            ymax=media+desvio,color=1),
                        color="black",width=0.3)}
        grafico=grafico+theme+
          ylab(ylab)+
          xlab(parse(text = xlab.factor[i]))+
          theme(text = element_text(size=textsize,color="black", family = family),
                axis.text = element_text(size=textsize,color="black", family = family),
                axis.title = element_text(size=textsize,color="black", family = family))+
          geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=TRUE)),lty=2)+
          scale_color_manual(values = "black")+labs(color="")}
        grafico=grafico+
          geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=TRUE)),lty=2)+
          scale_color_manual(values = "black")+labs(color="")
        print(grafico)
      }

      if(anavaF3[i,5]>alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))
        mean.table<-mean_stat(response,fatores[,i],mean)
        colnames(mean.table)<-c('Niveis','Medias')
        print(mean.table)
        grafico=NA}

      if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))

        dose=as.numeric(as.vector(unlist(fatores[,i])))
        grafico=polynomial(dose,resp,grau = grau[i],
                           DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = xlab,point = point)[[1]]
        cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))
        cat(green(bold("\n------------------------------------------")))}
      graficos[[1]]=residplot
      graficos[[i+1]]=grafico
    }
  }

    if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[2],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' within the combination of levels ', fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))
    des<-aov(resp~Fator2/Fator1+Fator3+Fator2+Fator2:Fator3+Fator1:Fator2:Fator3+bloco)
    l<-vector('list',nv2)
    names(l)<-names(summary(Fator2))
    v<-numeric(0)
    for(j in 1:nv2) {
      for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
      l[[j]]<-v
      v<-numeric(0)}
    des1<-summary(des,split=list('Fator2:Fator1'=l))[[1]]
    des1[nrow(des1),]=c(DfE,SQE,QME,NA,NA)
    des1$`F value`=c(des1$`Mean Sq`[1:nrow(des1)-1]/des1$`Mean Sq`[nrow(des1)],NA)
    des1$`Pr(>F)`=c(1-pf(des1$`F value`[1:nrow(des1)-1],
                       des1$Df[1:nrow(des1)-1],
                       des1$Df[nrow(des1)]),NA)
    des1a=des1[-c(1,2,3,4,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
    #============================
    rn<-numeric(0)
    for (j in 1:nv2) {
      rn <- c(rn, paste(paste(names.fat[1], ":", names.fat[2],
                              sep = ""), lf2[j]))
    }
    rownames(des1a)=rn
    #============================
    print(des1a)

    if(quali[1]==TRUE & quali[2]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
          tukeygrafico[[i]]=tukey$groups[levels(trati),2]
          ordem[[i]]=rownames(tukey$groups[levels(trati),])
        }
        letra=unlist(tukeygrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "duncan"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
          duncangrafico[[i]]=duncan$groups[levels(trati),2]
          ordem[[i]]=rownames(duncan$groups[levels(trati),])
        }
        letra=unlist(duncangrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "lsd"){
        lsdgrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
          lsdgrafico[[i]]=lsd$groups[levels(trati),2]
          ordem[[i]]=rownames(lsd$groups[levels(trati),])
        }
        letra=unlist(lsdgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "sk"){
        skgrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3$Df[10],
                        nrep = nrep,
                        QME = anavaF3$`Mean Sq`[10],
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
          skgrafico[[i]]=sk[levels(trati),2]
          ordem[[i]]=rownames(sk[levels(trati),])
        }
        letra=unlist(skgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}}

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], " inside of the level of ",fac.names[1])
    cat(green(bold("\n------------------------------------------\n")))
    des<-aov(resp~Fator1/Fator2+Fator3+Fator1+Fator1:Fator3+Fator1:Fator2:Fator3+bloco)
    l<-vector('list',nv1)
    names(l)<-names(summary(Fator1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
      l[[j]]<-v
      v<-numeric(0)}
    des1<-summary(des,split=list('Fator1:Fator2'=l))[[1]]
    des1[nrow(des1),]=c(DfE,SQE,QME,NA,NA)
    des1$`F value`=c(des1$`Mean Sq`[1:nrow(des1)-1]/des1$`Mean Sq`[nrow(des1)],NA)
    des1$`Pr(>F)`=c(1-pf(des1$`F value`[1:nrow(des1)-1],
                         des1$Df[1:nrow(des1)-1],
                         des1$Df[nrow(des1)]),NA)
    des1a=des1[-c(1,2,3,4,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
    #============================
    rn<-numeric(0)
    for (j in 1:nv1) {
      rn <- c(rn, paste(paste(names.fat[2], ":", names.fat[1],
                              sep = ""), lf1[j]))
    }
    rownames(des1a)=rn
    #============================
    print(des1a)

    if(quali[1]==TRUE & quali[2]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
          tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
        }
        letra1=unlist(tukeygrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "duncan"){
        duncangrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
          duncangrafico1[[i]]=duncan$groups[levels(trati),2]
        }
        letra1=unlist(duncangrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "lsd"){
        lsdgrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
          lsdgrafico1[[i]]=lsd$groups[levels(trati),2]
        }
        letra1=unlist(lsdgrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "sk"){
        skgrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3$Df[10],
                        nrep = nrep,
                        QME = anavaF3$`Mean Sq`[10],
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
          skgrafico1[[i]]=sk[levels(trati),2]
        }
        letra1=unlist(skgrafico1)
        letra1=toupper(letra1)}}

    if(quali[1]==TRUE & quali[2]==TRUE){
      f1=rep(levels(Fator1),e=length(levels(Fator2)))
      f2=rep(unique(as.character(Fator2)),length(levels(Fator1)))
      f1=factor(f1,levels = unique(f1))
      f2=factor(f2,levels = unique(f2))
      media=tapply(resp,paste(Fator1,Fator2), mean, na.rm=TRUE)[unique(paste(f1,f2))]
      if(point=="mean_sd"){desvio=tapply(response,paste(Fator1,Fator2), sd, na.rm=TRUE)}
      if(point=="mean_se"){desvio=tapply(response,paste(Fator1,Fator2), sd, na.rm=TRUE)/
        sqrt(tapply(response,paste(Fator1,Fator2), length))}
      desvio=desvio[unique(paste(f1,f2))]

      graph=data.frame(f1=f1,
                       f2=f2,
                       media,
                       desvio,
                       letra,letra1,
                       numero=format(media,digits = dec))
      numero=paste(graph$numero,graph$letra,graph$letra1,sep="")
      graph$numero=numero
      colint=ggplot(graph, aes(x=f2,
                               y=media,
                               fill=f1))+
        ylab(ylab)+xlab(xlab)+
        theme+
        labs(fill=fac.names[1])+
        geom_col(position = "dodge",color="black")+
        geom_errorbar(aes(ymin=media-desvio,
                          ymax=media+desvio),
                      width=0.3,position = position_dodge(width=0.9))+
        geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                      label=numero),
                  position = position_dodge(width=0.9),angle=angle.label, hjust=hjust,size=labelsize)+
        theme(text=element_text(size=textsize,family=family),
              axis.text = element_text(size=textsize,color="black",family=family),
              axis.title = element_text(size=textsize,color="black",family=family))+
        geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=T)),lty=2)+
        scale_color_manual(values = "black")+labs(color="")
      colint1=colint
      print(colint)
      letras=paste(graph$letra,graph$letra1,sep="")
      matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
      rownames(matriz)=levels(Fator1)
      colnames(matriz)=levels(Fator2)
      cat(green(bold("\n------------------------------------------\n")))
      cat(green(bold("Final table")))
      cat(green(bold("\n------------------------------------------\n")))
      print(matriz)
      cat("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")")
    }
    if(quali[1]==FALSE | quali[2]==FALSE){
      if(quali[1]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv1) {
            trati=fatores[, 2][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
          for (i in 1:nv1) {
            trati=fatores[, 2][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(duncan$groups)}}
        if (mcomp == "lsd"){
          for (i in 1:nv1) {
            trati=fatores[, 2][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(lsd$groups)}}
        if (mcomp == "sk"){
          for (i in 1:nv1) {
            trati=fatores[, 2][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            nrep=table(trati)[1]
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(sk)}}}
      if(quali[1]==FALSE){
        Fator1a=fator1a
        colint1=polynomial2(Fator1a,
                            response,
                            Fator2,
                            grau = grau12,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2])}
      if(quali[2]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv2) {
            trati=fatores[, 1][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
          for (i in 1:nv2) {
            trati=fatores[, 1][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(duncan$groups)}}
        if (mcomp == "lsd"){
          for (i in 1:nv2) {
            trati=fatores[, 1][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(lsd$groups)}}
        if (mcomp == "sk"){
          for (i in 1:nv2) {
            trati=fatores[, 1][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            nrep=table(trati)[1]
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(sk)}}
      }
      if(quali[2]==FALSE){
        Fator2a=fator2a
        colint1=polynomial2(Fator2a,
                            response,
                            Fator1,
                            grau = grau21,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2])}
      cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial2\" command\n"))}

    #Checar o Fator3
    if(anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f) {


      i<-3
      {
        #Para os fatores QUALITATIVOS, teste de Tukey
        if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
          cat(green(bold("\n------------------------------------------\n")))
          cat(green(italic('Analyzing the simple effects of the factor ',fac.names[i])))
          cat(green(bold("\n------------------------------------------\n")))
          cat(fac.names[i])
          if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="sk"){

            nrep=table(fatores[,i])[1]
            medias=sort(tapply(resp,fatores[i],mean, na.rm=TRUE),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3[10,1],
                          nrep = nrep,
                          QME = anavaF3[10,3],
                          alpha = alpha.t)
            letra1=data.frame(resp=medias,groups=sk)
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="duncan"){
            ad=data.frame(Fator1,Fator2,Fator3)
            letra <- duncan(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="lsd"){
            ad=data.frame(Fator1,Fator2,Fator3)
            letra <- LSD(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          print(letra1)
          cat(green(bold("\n------------------------------------------\n")))
          if(point=="mean_sd"){desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)]}
          if(point=="mean_se"){desvio=(tapply(response, c(fatores[i]), sd, na.rm=TRUE)/
                                         sqrt(tapply(response, c(fatores[i]), length)))[rownames(letra1)]}
          dadosm=data.frame(letra1,
                            media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                            desvio=desvio)
          dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
          dadosm$limite=dadosm$media+dadosm$desvio
          dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]

          if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
          if(addmean==FALSE){dadosm$letra=dadosm$groups}
          media=dadosm$media
          desvio=dadosm$desvio
          Tratamentos=dadosm$Tratamentos
          letra=dadosm$letra

          grafico=ggplot(dadosm,
                         aes(x=Tratamentos,
                             y=media))
          if(fill=="trat"){grafico=grafico+
            geom_col(aes(fill=Tratamentos),color=1)}
          else{grafico=grafico+
            geom_col(aes(fill=Tratamentos),fill=fill,color=1)}
          if(errorbar==TRUE){grafico=grafico+
            geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                          label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
          if(errorbar==FALSE){grafico=grafico+
            geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
          if(errorbar==TRUE){grafico=grafico+
            geom_errorbar(data=dadosm,
                          aes(ymin=media-desvio,
                              ymax=media+desvio,color=1), color="black",width=0.3)}
          grafico1=grafico+theme+
            ylab(ylab)+
            xlab(parse(text = xlab.factor[3]))+
            theme(text = element_text(size=textsize,color="black", family = family),
                  axis.text = element_text(size=textsize,color="black", family = family),
                  axis.title = element_text(size=textsize,color="black", family = family))+
            geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=T)),lty=2)+
            scale_color_manual(values = "black")+labs(color="")
          print(grafico1)}

      }

      #Para os fatores QUANTITATIVOS, regressao
      if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
        cat(green(bold("\n------------------------------------------\n")))
        cat('Analyzing the simple effects of the factor ',fac.names[3])
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        grafico1=polynomial(resp, fatores[,i],grau=grau[i],
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = parse(text = xlab.factor[3]),point = point)[[1]]
        cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))}
    }
  }

  #####################################################################################################
  #Interacao Fator1*Fator3       + fator2
  #####################################################################################################
  if(anavaF3[8,5]>alpha.f && anavaF3[6,5]<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------")))

    #Desdobramento de FATOR 1 dentro do niveis de FATOR 3
    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' within the combination of levels ', fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))
    des<-aov(resp~Fator3/Fator1+Fator2+Fator3+Fator2:Fator3+Fator1:Fator2:Fator3+bloco)
    l<-vector('list',nv3)
    names(l)<-names(summary(Fator3))
    v<-numeric(0)
    for(j in 1:nv3) {
      for(i in 0:(nv1-2)) v<-cbind(v,i*nv3+j)
      l[[j]]<-v
      v<-numeric(0)
    }
    des1<-summary(des,split=list('Fator3:Fator1'=l))[[1]]
    des1[nrow(des1),]=c(DfE,SQE,QME,NA,NA)
    des1$`F value`=c(des1$`Mean Sq`[1:nrow(des1)-1]/des1$`Mean Sq`[nrow(des1)],NA)
    des1$`Pr(>F)`=c(1-pf(des1$`F value`[1:nrow(des1)-1],
                         des1$Df[1:nrow(des1)-1],
                         des1$Df[nrow(des1)]),NA)
    des1a=des1[-c(1,2,3,4,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
    #============================
    rn<-numeric(0)
    for (j in 1:nv3) {
      rn <- c(rn, paste(paste(names.fat[1], ":", names.fat[3],
                              sep = ""), lf3[j]))
    }
    rownames(des1a)=rn
    #============================
    print(des1a)

    # Teste de Tukey
    if(quali[1]==TRUE & quali[3]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 1][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          tukeygrafico[[i]]=tukey$groups[levels(trati),2]
          ordem[[i]]=rownames(tukey$groups[levels(trati),])
        }
        letra=unlist(tukeygrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "duncan"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 1][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          duncangrafico[[i]]=duncan$groups[levels(trati),2]
          ordem[[i]]=rownames(duncan$groups[levels(trati),])
        }
        letra=unlist(duncangrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "lsd"){
        lsdgrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 1][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          lsdgrafico[[i]]=lsd$groups[levels(trati),2]
          ordem[[i]]=rownames(lsd$groups[levels(trati),])
        }
        letra=unlist(lsdgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "sk"){
        skgrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 1][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3$Df[10],
                        nrep = nrep,
                        QME = anavaF3$`Mean Sq`[10],
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          skgrafico[[i]]=sk[levels(trati),2]
          ordem[[i]]=rownames(sk[levels(trati),])
        }
        letra=unlist(skgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}}

    # Desdobramento de F3 dentro de F1

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], " inside of the level of ",fac.names[1])
    cat(green(bold("\n------------------------------------------\n")))
    des<-aov(resp~Fator1/Fator3+Fator1+Fator2+Fator2:Fator1+Fator1:Fator2:Fator3+bloco)
    l<-vector('list',nv1)
    names(l)<-names(summary(Fator1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv3-2)) v<-cbind(v,i*nv1+j)
      l[[j]]<-v
      v<-numeric(0)
    }
    des1<-summary(des,split=list('Fator1:Fator3'=l))[[1]]
    des1[nrow(des1),]=c(DfE,SQE,QME,NA,NA)
    des1$`F value`=c(des1$`Mean Sq`[1:nrow(des1)-1]/des1$`Mean Sq`[nrow(des1)],NA)
    des1$`Pr(>F)`=c(1-pf(des1$`F value`[1:nrow(des1)-1],
                         des1$Df[1:nrow(des1)-1],
                         des1$Df[nrow(des1)]),NA)
    des1a=des1[-c(1,2,3,4,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
    #============================
    rn<-numeric(0)
    for (j in 1:nv1) {
      rn <- c(rn, paste(paste(names.fat[3], ":", names.fat[1],
                              sep = ""), lf1[j]))
    }
    rownames(des1a)=rn
    #============================
    print(des1a)

    #-------------------------------------
    # Teste de Tukey
    #-------------------------------------
    if(quali[1]==TRUE & quali[3]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 3][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
        }
        letra1=unlist(tukeygrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "duncan"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 3][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          duncangrafico[[i]]=duncan$groups[levels(trati),2]
          ordem[[i]]=rownames(duncan$groups[levels(trati),])
        }
        letra=unlist(duncangrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "lsd"){
        lsdgrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 3][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          lsdgrafico[[i]]=lsd$groups[levels(trati),2]
          ordem[[i]]=rownames(lsd$groups[levels(trati),])
        }
        letra=unlist(lsdgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "sk"){
        skgrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 3][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3$Df[10],
                        nrep = nrep,
                        QME = anavaF3$`Mean Sq`[10],
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          skgrafico1[[i]]=sk[levels(trati),2]
        }
        letra1=unlist(skgrafico1)
        letra1=toupper(letra1)}}

    # -----------------------------
    # Gráfico de colunas
    #------------------------------
    if(quali[1]==TRUE & quali[3]==TRUE){
      f1=rep(levels(Fator1),e=length(levels(Fator3)))
      f3=rep(unique(as.character(Fator3)),length(levels(Fator1)))
      f1=factor(f1,levels = unique(f1))
      f3=factor(f3,levels = unique(f3))
      media=tapply(response,paste(Fator1,Fator3), mean, na.rm=TRUE)[unique(paste(f1,f3))]
      if(point=="mean_sd"){desvio=tapply(response,paste(Fator1,Fator3), sd, na.rm=TRUE)}
      if(point=="mean_se"){desvio=tapply(response,paste(Fator1,Fator3), sd, na.rm=TRUE)/
        sqrt(tapply(response,paste(Fator1,Fator3), length))}
      desvio=desvio[unique(paste(f1,f3))]

      graph=data.frame(f1=f1,
                       f3=f3,
                       media,
                       desvio,
                       letra,
                       letra1,
                       numero=format(media,digits = dec))
      numero=paste(graph$numero,graph$letra,graph$letra1,sep="")
      graph$numero=numero
      colint=ggplot(graph,
                    aes(x=f3,
                        y=media,
                        fill=f1))+
        geom_col(position = "dodge",color="black")+
        ylab(ylab)+
        xlab(xlab)+
        theme+
        labs(fill=fac.names[1])+
        geom_errorbar(aes(ymin=media-desvio,
                          ymax=media+desvio),
                      width=0.3,
                      position = position_dodge(width=0.9))+
        geom_text(aes(y=media+desvio+sup,
                      label=numero),
                  position = position_dodge(width=0.9),angle=angle.label, hjust=hjust,size=labelsize)+
        theme(text=element_text(size=textsize,family=family),
              axis.text = element_text(size=textsize,color="black",family=family),
              axis.title = element_text(size=textsize,color="black",family=family))+
        geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=T)),lty=2)+
        scale_color_manual(values = "black")+labs(color="")
      colint2=colint
      print(colint)
      letras=paste(graph$letra,graph$letra1,sep="")
      matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
      rownames(matriz)=levels(Fator1)
      colnames(matriz)=levels(Fator3)
      cat(green(bold("\n------------------------------------------\n")))
      cat(green(bold("Final table")))
      cat(green(bold("\n------------------------------------------\n")))
      print(matriz)
      cat("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")")
    }
    if(quali[1]==FALSE | quali[3]==FALSE){
      if(quali[1]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv1) {
            trati=fatores[, 3][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
          for (i in 1:nv3) {
            trati=fatores[, 3][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(duncan$groups)}}
        if (mcomp == "lsd"){
          for (i in 1:nv3) {
            trati=fatores[, 3][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(lsd$groups)}}
        if (mcomp == "sk"){
          for (i in 1:nv1) {
            trati=fatores[, 3][Fator1 == lf1[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator1 == lf1[i]]
            nrep=table(trati)[1]
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(sk)}}}
      if(quali[1]==FALSE){
        Fator1a=fator1a
        colint2=polynomial2(Fator1a,
                            response,
                            Fator3,
                            grau = grau13,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2])}
      if(quali[3]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv3) {
            trati=fatores[, 1][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Sum Sq`[10],alpha.t)
            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
          for (i in 1:nv3) {
            trati=fatores[, 1][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(duncan$groups)}}
        if (mcomp == "lsd"){
          for (i in 1:nv3) {
            trati=fatores[, 1][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(lsd$groups)}}
        if (mcomp == "sk"){
          for (i in 1:nv3) {
            trati=fatores[, 1][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            nrep=table(trati)[1]
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(sk)}}}
      if(quali[3]==FALSE){
        Fator3a=fator3a
        colint2=polynomial2(Fator3a,
                            response,
                            Fator1,
                            grau = grau31,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2])}
      cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial2\" command\n"))
    }

    if(anavaF3[5,5]>alpha.f && anavaF3[7,5]>alpha.f) {


      i<-2
      {
        #Para os fatores QUALITATIVOS, teste de Tukey
        if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
          cat(green(bold("\n------------------------------------------\n")))
          cat(green(italic('Analyzing the simple effects of the factor ',fac.names[i])))
          cat(green(bold("\n------------------------------------------\n")))
          cat(fac.names[i])
          if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="sk"){
            nrep=table(fatores[,i])[1]
            medias=sort(tapply(resp,fatores[i],mean, na.rm=TRUE),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3[10,1],
                          nrep = nrep,
                          QME = anavaF3[10,3],
                          alpha = alpha.t)
            letra1=data.frame(resp=medias,groups=sk)
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="duncan"){
            ad=data.frame(Fator1,Fator2,Fator3)
            letra <- duncan(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="lsd"){
            ad=data.frame(Fator1,Fator2,Fator3)
            letra <- LSD(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          print(letra1)
          cat(green(bold("\n------------------------------------------")))
          if(point=="mean_sd"){desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)]}
          if(point=="mean_se"){desvio=(tapply(response, c(fatores[i]), sd, na.rm=TRUE)/
                                         sqrt(tapply(response, c(fatores[i]), length)))[rownames(letra1)]}
          dadosm=data.frame(letra1,
                            media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                            desvio=desvio)
          dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
          dadosm$limite=dadosm$media+dadosm$desvio
          dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]
          if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
          if(addmean==FALSE){dadosm$letra=dadosm$groups}
          media=dadosm$media
          desvio=dadosm$desvio
          Tratamentos=dadosm$Tratamentos
          letra=dadosm$letra

          grafico=ggplot(dadosm,
                         aes(x=Tratamentos,
                             y=media))
          if(fill=="trat"){grafico=grafico+
            geom_col(aes(fill=Tratamentos),color=1)}
          else{grafico=grafico+
            geom_col(aes(fill=Tratamentos),fill=fill,color=1)}
          if(errorbar==TRUE){grafico=grafico+
            geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
          if(errorbar==FALSE){grafico=grafico+
            geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
          if(errorbar==TRUE){grafico=grafico+
            geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                          ymax=media+desvio,color=1),
                          color="black",width=0.3)
          grafico2=grafico+theme+
            ylab(ylab)+
            xlab(parse(text = xlab.factor[2]))+
            theme(text = element_text(size=textsize,color="black", family = family),
                  axis.text = element_text(size=textsize,color="black", family = family),
                  axis.title = element_text(size=textsize,color="black", family = family))+
            geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=TRUE)),lty=2)+
            scale_color_manual(values = "black")+labs(color="")
          print(grafico2)}
        }

        #Para os fatores QUANTITATIVOS, regressao
        if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
          cat(green(bold("\n------------------------------------------\n")))
          cat('Analyzing the simple effects of the factor ',fac.names[2])
          cat(green(bold("\n------------------------------------------\n")))
          cat(fac.names[i])
          grafico2=polynomial(resp, fatores[,i],grau=grau[i],
                              DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = parse(text = xlab.factor[2]),point = point)[[1]]
          cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))
        }

        cat('\n')
      }
    }
  }

  if(anavaF3[8,5]>alpha.f && anavaF3[7,5]<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' within the combination of levels ', fac.names[3])
    cat("\n-------------------------------------------------\n")
    des<-aov(resp~Fator3/Fator2+Fator1+Fator3+Fator1:Fator3+Fator1:Fator2:Fator3+bloco)
    l<-vector('list',nv3)
    names(l)<-names(summary(Fator3))
    v<-numeric(0)
    for(j in 1:nv3) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv3+j)
      l[[j]]<-v
      v<-numeric(0)
    }
    des1<-summary(des,split=list('Fator3:Fator2'=l))[[1]]
    des1[nrow(des1),]=c(DfE,SQE,QME,NA,NA)
    des1$`F value`=c(des1$`Mean Sq`[1:nrow(des1)-1]/des1$`Mean Sq`[nrow(des1)],NA)
    des1$`Pr(>F)`=c(1-pf(des1$`F value`[1:nrow(des1)-1],
                         des1$Df[1:nrow(des1)-1],
                         des1$Df[nrow(des1)]),NA)
    des1a=des1[-c(1,2,3,4,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
    #============================
    rn<-numeric(0)
    for (j in 1:nv3) {
      rn <- c(rn, paste(paste(names.fat[2], ":", names.fat[3],
                              sep = ""), lf3[j]))
    }
    rownames(des1a)=rn
    #============================
    print(des1a)

    # Teste de Tukey
    if(quali[2]==TRUE & quali[3]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 2][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          tukeygrafico[[i]]=tukey$groups[levels(trati),2]
          ordem[[i]]=rownames(tukey$groups[levels(trati),])
        }
        letra=unlist(tukeygrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "duncan"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 2][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          duncangrafico[[i]]=duncan$groups[levels(trati),2]
          ordem[[i]]=rownames(duncan$groups[levels(trati),])
        }
        letra=unlist(duncangrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "lsd"){
        lsdgrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 2][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          lsdgrafico[[i]]=lsd$groups[levels(trati),2]
          ordem[[i]]=rownames(lsd$groups[levels(trati),])
        }
        letra=unlist(lsdgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
      if (mcomp == "sk"){
        skgrafico=c()
        ordem=c()
        for (i in 1:nv3) {
          trati=fatores[, 2][Fator3 == lf3[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator3 == lf3[i]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3$Df[10],
                        nrep = nrep,
                        QME = anavaF3$`Mean Sq`[10],
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          skgrafico[[i]]=sk[levels(trati),2]
          ordem[[i]]=rownames(sk[levels(trati),])
        }
        letra=unlist(skgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}}

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], " inside of the level of ",fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))
    cat("\n")
    des<-aov(resp~Fator2/Fator3+Fator1+Fator2+Fator1:Fator2+Fator1:Fator2:Fator3+bloco)
    l<-vector('list',nv2)
    names(l)<-names(summary(Fator2))
    v<-numeric(0)
    for(j in 1:nv2) {
      for(i in 0:(nv3-2)) v<-cbind(v,i*nv2+j)
      l[[j]]<-v
      v<-numeric(0)
    }
    des1<-summary(des,split=list('Fator2:Fator3'=l))[[1]]
    des1[nrow(des1),]=c(DfE,SQE,QME,NA,NA)
    des1$`F value`=c(des1$`Mean Sq`[1:nrow(des1)-1]/des1$`Mean Sq`[nrow(des1)],NA)
    des1$`Pr(>F)`=c(1-pf(des1$`F value`[1:nrow(des1)-1],
                         des1$Df[1:nrow(des1)-1],
                         des1$Df[nrow(des1)]),NA)
    des1a=des1[-c(1,2,3,4,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
    #============================
    rn<-numeric(0)
    for (j in 1:nv2) {
      rn <- c(rn, paste(paste(names.fat[3], ":", names.fat[2],
                              sep = ""), lf2[j]))
    }
    rownames(des1a)=rn
    #============================
    print(des1a)

    if(quali[2]==TRUE & quali[3]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico1=c()
        for (i in 1:nv2) {
          trati=fatores[, 3][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
        }
        letra1=unlist(tukeygrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "duncan"){
        duncangrafico1=c()
        for (i in 1:nv2) {
          trati=fatores[, 3][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          duncangrafico1[[i]]=duncan$groups[levels(trati),2]
        }
        letra1=unlist(duncangrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "lsd"){
        lsdgrafico1=c()
        for (i in 1:nv2) {
          trati=fatores[, 3][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
          lsdgrafico1[[i]]=lsd$groups[levels(trati),2]
        }
        letra1=unlist(lsdgrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "sk"){
        skgrafico1=c()
        for (i in 1:nv2) {
          trati=fatores[, 3][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = anavaF3$Df[10],
                        nrep = nrep,
                        QME = anavaF3$`Mean Sq`[10],
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          skgrafico1[[i]]=sk[levels(trati),2]
        }
        letra1=unlist(skgrafico1)
        letra1=toupper(letra1)}}

    if(quali[2]==TRUE & quali[3]==TRUE){
      f2=rep(levels(Fator2),e=length(levels(Fator3)))
      f3=rep(unique(as.character(Fator3)),length(levels(Fator2)))
      f2=factor(f2,levels = unique(f2))
      f3=factor(f3,levels = unique(f3))
      media=tapply(response,paste(Fator2,Fator3), mean, na.rm=TRUE)[unique(paste(f2,f3))]
      if(point=="mean_sd"){desvio=tapply(response,paste(Fator2,Fator3), sd, na.rm=TRUE)}
      if(point=="mean_se"){desvio=tapply(response,paste(Fator2,Fator3), sd, na.rm=TRUE)/
        sqrt(tapply(response,paste(Fator2,Fator3), length))}
      desvio=desvio[unique(paste(f2,f3))]

      graph=data.frame(f2=f2,
                       f3=f3,
                       media,
                       desvio,
                       letra,letra1,
                       numero=format(media,digits = dec))
      numero=paste(graph$numero,graph$letra,graph$letra1,sep="")
      graph$numero=numero
      colint=ggplot(graph, aes(x=f3,
                               y=media,
                               fill=f2))+
        geom_col(position = "dodge",color="black")+
        ylab(ylab)+xlab(xlab)+
        theme+
        labs(fill=fac.names[2])+
        geom_errorbar(aes(ymin=media-desvio,
                          ymax=media+desvio),
                      width=0.3,position = position_dodge(width=0.9))+
        geom_text(aes(y=media+desvio+sup,
                      label=numero),
                  position = position_dodge(width=0.9),angle=angle.label, hjust=hjust,size=labelsize)+
        theme(text=element_text(size=textsize,family=family),
              axis.text = element_text(size=textsize,color="black",family=family),
              axis.title = element_text(size=textsize,color="black",family=family))+
        geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=T)),lty=2)+
        scale_color_manual(values = "black")+labs(color="")
      colint3=colint
      print(colint)
      letras=paste(graph$letra,graph$letra1,sep="")
      matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),
                                 ncol = length(levels(Fator2)))))
      rownames(matriz)=levels(Fator2)
      colnames(matriz)=levels(Fator3)
      cat(green(bold("\n------------------------------------------\n")))
      cat(green(bold("Final table")))
      cat(green(bold("\n------------------------------------------\n")))
      print(matriz)
      cat("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")")
    }
    if(quali[2]==FALSE | quali[3]==FALSE){
      if(quali[2]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv2) {
            trati=fatores[, 3][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
          for (i in 1:nv2) {
            trati=fatores[, 3][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(duncan$groups)}}
        if (mcomp == "lsd"){
          for (i in 1:nv2) {
            trati=fatores[, 3][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(lsd$groups)}}
        if (mcomp == "sk"){
          for (i in 1:nv2) {
            trati=fatores[, 3][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator2 == lf2[i]]
            nrep=table(trati)[1]
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(sk)}}}
      if(quali[2]==FALSE){
        Fator2a=fator2a
        colint3=polynomial2(Fator2a,
                            response,
                            Fator3,
                            grau = grau23,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2])}
      if(quali[3]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv3) {
            trati=fatores[, 2][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            tukey=TUKEY(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(tukey)}}
        if (mcomp == "duncan"){
          for (i in 1:nv3) {
            trati=fatores[, 2][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            duncan=duncan(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(duncan)}}
        if (mcomp == "lsd"){
          for (i in 1:nv3) {
            trati=fatores[, 2][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            lsd=LSD(respi,trati,anavaF3$Df[10],anavaF3$`Mean Sq`[10],alpha.t)
            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(lsd)}}
        if (mcomp == "sk"){
          for (i in 1:nv3) {
            trati=fatores[, 2][Fator3 == lf3[i]]
            trati=factor(trati,levels = unique(trati))
            respi=resp[Fator3 == lf3[i]]
            nrep=table(trati)[1]
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
            cat("\n----------------------\n")
            print(sk)}}}
      if(quali[3]==FALSE){
        Fator3a=fator3a
        colint3=polynomial2(Fator3a,
                            response,
                            Fator2,
                            grau = grau32,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            DFres= anavaF3[10,1],SSq = anavaF3[10,2])}

      cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial2\" command\n"))
    }

    if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f) {

      i<-1
      {
        #Para os fatores QUALITATIVOS, teste de Tukey
        if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
          cat(green(bold("\n------------------------------------------\n")))
          cat(green(italic('Analyzing the simple effects of the factor ',fac.names[i])))
          cat(green(bold("\n------------------------------------------\n")))
          cat(fac.names[i])
          if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="sk"){
            nrep=table(fatores[,i])[1]
            medias=sort(tapply(resp,fatores[i],mean, na.rm=TRUE),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3[10,1],
                          nrep = nrep,
                          QME = anavaF3[10,3],
                          alpha = alpha.t)
            letra1=data.frame(resp=medias,groups=sk)
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="duncan"){
            ad=data.frame(Fator1,Fator2,Fator3)
            letra <- duncan(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          if(mcomp=="lsd"){
            ad=data.frame(Fator1,Fator2,Fator3)
            letra <- LSD(resp,fatores[,i],anavaF3[10,1],anavaF3[10,3], alpha=alpha.t)
            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
          print(letra1)
          cat(green(bold("\n------------------------------------------")))
          if(point=="mean_sd"){desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)]}
          if(point=="mean_se"){desvio=(tapply(response, c(fatores[i]), sd, na.rm=TRUE)/
                                         sqrt(tapply(response, c(fatores[i]), length)))[rownames(letra1)]}
          dadosm=data.frame(letra1,
                            media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                            desvio=desvio)
          dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
          dadosm$limite=dadosm$media+dadosm$desvio
          dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]
          if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
          if(addmean==FALSE){dadosm$letra=dadosm$groups}
          media=dadosm$media
          desvio=dadosm$desvio
          Tratamentos=dadosm$Tratamentos
          letra=dadosm$letra

          grafico=ggplot(dadosm,aes(x=Tratamentos,
                                    y=media))
          if(fill=="trat"){grafico=grafico+
            geom_col(aes(fill=Tratamentos),color=1)}
          else{grafico=grafico+
            geom_col(aes(fill=Tratamentos),fill=fill,color=1)}
          if(errorbar==TRUE){grafico=grafico+
            geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                          label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
          if(errorbar==FALSE){grafico=grafico+
            geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust,size=labelsize)}
          if(errorbar==TRUE){grafico=grafico+
            geom_errorbar(data=dadosm,
                          aes(ymin=media-desvio,
                              ymax=media+desvio,color=1),
                          color="black",width=0.3)
          grafico3=grafico+theme+
            ylab(ylab)+
            xlab(parse(text = xlab.factor[1]))+
            theme(text = element_text(size=textsize,color="black", family = family),
                  axis.text = element_text(size=textsize,color="black", family = family),
                  axis.title = element_text(size=textsize,color="black", family = family))+
            geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=TRUE)),lty=2)+
            scale_color_manual(values = "black")+labs(color="")
          print(grafico3)}
        }
        if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
          cat(green(bold("\n------------------------------------------\n")))
          cat('\nAnalyzing the simple effects of the factor ',fac.names[1],'\n')
          cat(green(bold("\n------------------------------------------\n")))
          cat(fac.names[i])
          grafico3=polynomial(resp, fatores[,i],grau=grau[i],
                              DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = parse(text = xlab.factor[1]),point = point)[[1]]
          cat(green("\nTo edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))
        }

        cat('\n')
      }
    }
  }

  if(anavaF3[8,5]<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("\nInteraction",paste(fac.names[1],'*',fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction\n")))
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' within the combination of levels ', fac.names[2], 'and',fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))

    m1=aov(resp~(Fator2*Fator3)/Fator1+bloco)
    anova(m1)
    pattern <- c(outer(levels(Fator2), levels(Fator3),
                       function(x,y) paste("Fator2",x,":Fator3",y,":",sep="")))
    des.tab <- sapply(pattern, simplify=FALSE,
                      grep, x=names(coef(m1)[m1$assign==5]))
    des1.tab <- summary(m1, split = list("Fator2:Fator3:Fator1" = des.tab))
    des1.tab[[1]][nrow(des1.tab[[1]]),]=c(DfE,SQE,QME,NA,NA)
    des1.tab[[1]]$`F value`=c(des1.tab[[1]]$`Mean Sq`[1:nrow(des1.tab[[1]])-1]/
                                des1.tab[[1]]$`Mean Sq`[nrow(des1.tab[[1]])],NA)
    des1.tab[[1]]$`Pr(>F)`=c(1-pf(des1.tab[[1]]$`F value`[1:nrow(des1.tab[[1]])-1],
                         des1.tab[[1]]$Df[1:nrow(des1.tab[[1]])-1],
                         des1.tab[[1]]$Df[nrow(des1.tab[[1]])]),NA)
    desd=des1.tab[[1]][-c(1,2,3,4,5),]
    desd=data.frame(desd[-length(rownames(desd)),])
    # rownames(desd)=cbind(paste("Fator2:",rep(levels(Fator2),length(levels(Fator3))),
    #                            "Fator3:",rep(levels(Fator3),e=length(levels(Fator2)))))
    rownames(desd)=cbind(paste(names.fat[2],":",rep(levels(Fator2),length(levels(Fator3))),
                               names.fat[3],":",rep(levels(Fator3),e=length(levels(Fator2)))))
    colnames(desd)=c("Df",  "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    print(desd)

    ii<-0
    for(i in 1:nv2) {
      for(j in 1:nv3) {
        ii<-ii+1
        if(quali[1]==TRUE){
          cat('\n\n',fac.names[1],' inside of each level of ',lf2[i],' of ',fac.names[2],' and ',lf3[j],' of ',fac.names[3],"\n")
          if(mcomp=='tukey'){tukey=TUKEY(y = resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                         trt = fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                                         DFerror = anavaF3[10,1],
                                         MSerror = anavaF3[10,3],
                                         alpha.t)
          tukey=tukey$groups;colnames(tukey)=c("resp","letters")
          if(transf !=1){tukey$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                            fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(tukey)]}
          print(tukey)}
          if(mcomp=='duncan'){duncan=duncan(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                            fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                                            anavaF3[10,1],
                                            anavaF3[10,3],
                                            alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("resp","letters")
          if(transf !=1){duncan$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                             fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(duncan)]}
          print(duncan)}
          if(mcomp=='lsd'){lsd=LSD(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                   fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                                   anavaF3[10,1],
                                   anavaF3[10,3],
                                   alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("resp","letters")
          if(transf !=1){lsd$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                          fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(lsd)]}
          print(lsd)}
          if(mcomp=='sk'){
            fat= fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]]
            fat1=factor(fat,unique(fat))
            levels(fat1)=1:length(levels(fat1))
            resp1=resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]]
            nrep=table(fat1)[1]
            medias=sort(tapply(respi,fat1,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            sk=sk[as.character(unique(fat1)),]
            rownames(sk)=unique(fat)
            if(transf !=1){sk$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                           fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(sk)]}
            print(sk)}
        }
        if(quali[1]==FALSE){
          cat('\n',fac.names[1],' within the combination of levels ',lf2[i],' of  ',fac.names[2],' and ',lf3[j],' of  ',fac.names[3],"\n")
          polynomial(fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                     resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],grau=grau123,
                     DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = xlab,point = point)[[1]]}
      }
    }

    cat('\n\n')

    cat("\n------------------------------------------\n")
    cat("Analyzing ", fac.names[2], ' within the combination of levels ', fac.names[1], 'and',fac.names[3])
    cat("\n------------------------------------------\n")
    m1=aov(resp~(Fator1*Fator3)/Fator2+bloco)
    anova(m1)
    pattern <- c(outer(levels(Fator1), levels(Fator3),
                       function(x,y) paste("Fator1",x,":Fator3",y,":",sep="")))
    des.tab <- sapply(pattern, simplify=FALSE,
                      grep, x=names(coef(m1)[m1$assign==5]))
    des1.tab <- summary(m1, split = list("Fator1:Fator3:Fator2" = des.tab))
    des1.tab[[1]][nrow(des1.tab[[1]]),]=c(DfE,SQE,QME,NA,NA)
    des1.tab[[1]]$`F value`=c(des1.tab[[1]]$`Mean Sq`[1:nrow(des1.tab[[1]])-1]/
                                des1.tab[[1]]$`Mean Sq`[nrow(des1.tab[[1]])],NA)
    des1.tab[[1]]$`Pr(>F)`=c(1-pf(des1.tab[[1]]$`F value`[1:nrow(des1.tab[[1]])-1],
                                  des1.tab[[1]]$Df[1:nrow(des1.tab[[1]])-1],
                                  des1.tab[[1]]$Df[nrow(des1.tab[[1]])]),NA)

    desd=des1.tab[[1]][-c(1,2,3,4,5),]
    desd=data.frame(desd[-length(rownames(desd)),])
    # rownames(desd)=cbind(paste("Fator1:",rep(levels(Fator1),length(levels(Fator3))),
    #                            "Fator3:",rep(levels(Fator3),e=length(levels(Fator1)))))
    rownames(desd)=cbind(paste(names.fat[1],":",rep(levels(Fator1),length(levels(Fator3))),
                               names.fat[3],":",rep(levels(Fator3),e=length(levels(Fator1)))))
    colnames(desd)=c("Df",  "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    print(desd)

    ii<-0
    for(k in 1:nv1) {
      for(j in 1:nv3) {
        ii<-ii+1
        if(quali[2]==TRUE){
          cat('\n\n',fac.names[2],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf3[j],' of ',fac.names[3],'\n')
          if(mcomp=='tukey'){tukey=TUKEY(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                         fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                                         anavaF3[10,1],
                                         anavaF3[10,3],
                                         alpha.t)
          tukey=tukey$groups;colnames(tukey)=c("resp","letters")
          if(transf !=1){tukey$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                            fatores[,2][Fator1==lf1[k]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(tukey)]}
          print(tukey)}
          if(mcomp=='duncan'){duncan=duncan(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                            fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                                            anavaF3[10,1],
                                            anavaF3[10,3],
                                            alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("resp","letters")
          if(transf !=1){duncan$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                             fatores[,2][Fator1==lf1[k]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(duncan)]}

          print(duncan)}
          if(mcomp=='lsd'){lsd=LSD(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                   fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                                   anavaF3[10,1],
                                   anavaF3[10,3],
                                   alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("resp","letters")
          if(transf !=1){lsd$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                          fatores[,2][Fator1==lf1[k]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(lsd)]}
          print(lsd)}
          if(mcomp=='sk'){
            fat=fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]]
            fat1=factor(fat,unique(fat))
            levels(fat1)=1:length(levels(fat1))
            resp1=resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]]
            nrep=table(fat1)[1]
            medias=sort(tapply(respi,fat1,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            sk=sk[as.character(unique(fat1)),]
            rownames(sk)=unique(fat)
            if(transf !=1){sk$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                           fatores[,2][Fator1==lf1[k]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(sk)]}
            print(sk)}

        }
        if(quali[2]==FALSE){
          cat('\n\n',fac.names[2],' within the combination of levels ',lf1[k],' of  ',fac.names[1],' and ',lf3[j],' of  ',fac.names[3],'\n')
          polynomial(fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                     resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],grau=grau213,
                     DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = xlab,point = point)[[1]]}
      }
    }

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' within the combination of levels ', fac.names[1], 'and',fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))

    m1=aov(resp~(Fator1*Fator2)/Fator3+bloco)
    anova(m1)
    pattern <- c(outer(levels(Fator1), levels(Fator2),
                       function(x,y) paste("Fator1",x,":Fator2",y,":",sep="")))
    des.tab <- sapply(pattern, simplify=FALSE,
                      grep, x=names(coef(m1)[m1$assign==5]))
    des1.tab <- summary(m1, split = list("Fator1:Fator2:Fator3" = des.tab))
    des1.tab[[1]][nrow(des1.tab[[1]]),]=c(DfE,SQE,QME,NA,NA)
    des1.tab[[1]]$`F value`=c(des1.tab[[1]]$`Mean Sq`[1:nrow(des1.tab[[1]])-1]/
                                des1.tab[[1]]$`Mean Sq`[nrow(des1.tab[[1]])],NA)
    des1.tab[[1]]$`Pr(>F)`=c(1-pf(des1.tab[[1]]$`F value`[1:nrow(des1.tab[[1]])-1],
                                  des1.tab[[1]]$Df[1:nrow(des1.tab[[1]])-1],
                                  des1.tab[[1]]$Df[nrow(des1.tab[[1]])]),NA)
    desd=des1.tab[[1]][-c(1,2,3,4,5),]
    desd=data.frame(desd[-length(rownames(desd)),])
    # rownames(desd)=cbind(paste("Fator1:",rep(levels(Fator1),length(levels(Fator2))),
    #                            "Fator2:",rep(levels(Fator2),e=length(levels(Fator1)))))
    rownames(desd)=cbind(paste(names.fat[1],":",rep(levels(Fator1),length(levels(Fator2))),
                               names.fat[2],":",rep(levels(Fator2),e=length(levels(Fator1)))))
    colnames(desd)=c("Df",  "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    print(desd)

    ii<-0
    for(k in 1:nv1) {
      for(i in 1:nv2) {
        ii<-ii+1
        if(quali[3]==TRUE){
          cat('\n\n',fac.names[3],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf2[i],' of ',fac.names[2],'\n')
          if(mcomp=='tukey'){tukey=TUKEY(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                         fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                         anavaF3[10,1],
                                         anavaF3[10,3],
                                         alpha.t)
          tukey=tukey$groups;colnames(tukey)=c("resp","letters")
          if(transf !=1){tukey$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                            fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                            mean, na.rm=TRUE)[rownames(tukey)]}
          print(tukey)}
          if(mcomp=='duncan'){duncan=duncan(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                            fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                            anavaF3[10,1],
                                            anavaF3[10,3],
                                            alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("resp","letters")
          if(transf !=1){duncan$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                             fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                             mean, na.rm=TRUE)[rownames(duncan)]}
          print(duncan)}
          if(mcomp=='lsd'){lsd=LSD(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                   fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                   anavaF3[10,1],
                                   anavaF3[10,3],
                                   alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("resp","letters")
          if(transf !=1){lsd$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                          fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                          mean, na.rm=TRUE)[rownames(lsd)]}
          print(lsd)}
          if(mcomp=='sk'){
            fat=fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
            fat1=factor(fat,unique(fat))
            levels(fat1)=1:length(levels(fat1))
            resp1=resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
            nrep=table(fat1)[1]
            medias=sort(tapply(respi,fat1,mean),decreasing = TRUE)
            sk=scottknott(means = medias,
                          df1 = anavaF3$Df[10],
                          nrep = nrep,
                          QME = anavaF3$`Mean Sq`[10],
                          alpha = alpha.t)
            sk=data.frame(respi=medias,groups=sk)
            colnames(sk)=c("resp","letters")
            sk=sk[as.character(unique(fat1)),]
            rownames(sk)=unique(fat)
            if(transf !=1){sk$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                           fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                           mean, na.rm=TRUE)[rownames(sk)]}
            print(sk)}

        }
        if(quali[3]==FALSE){
          cat('\n\n',fac.names[3],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf2[i],' of ',fac.names[2],'\n')
          polynomial(fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                     resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],grau=grau312,
                     DFres= anavaF3[10,1],SSq = anavaF3[10,2],ylab = ylab,xlab = xlab,point = point)[[1]]}
      }
    }

  }

  if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[8,5]>alpha.f){
    if(anavaF3[1,5]<=alpha.f | anavaF3[2,5]<=alpha.f | anavaF3[3,5]<=alpha.f){
      graficos}else{graficos=NA}}
  if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f){
    graficos=list(residplot,colint1)
    if(anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[3,5]<=alpha.f){
      graficos=list(residplot,colint1,grafico1)}
    graficos}
  if(anavaF3[8,5]>alpha.f && anavaF3[6,5]<=alpha.f){
    graficos=list(residplot,colint2)
    if(anavaF3[5,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[2,5]<=alpha.f){
      graficos=list(residplot,colint2,grafico2)}
    graficos}
  if(anavaF3[8,5]>alpha.f && anavaF3[7,5]<=alpha.f){
    graficos=list(residplot,colint3)
    if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f && anavaF3[1,5]<=alpha.f){
      graficos=list(residplot,colint3,grafico3)}}
  if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f && anavaF3[6,5]<=alpha.f){
    graficos=list(residplot,colint1,colint2)
    graficos}
  if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f && anavaF3[7,5]<=alpha.f){
    graficos=list(residplot,colint1,colint3)
    graficos}
  if(anavaF3[8,5]>alpha.f && anavaF3[6,5]<=alpha.f && anavaF3[7,5]<=alpha.f){
    graficos=list(residplot,colint2,colint3)
    graficos}
  if(anavaF3[8,5]<=alpha.f){graficos=list(residplot)}
  graficos=graficos
}

