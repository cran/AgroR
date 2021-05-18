#' Analysis: Joint analysis of experiments in completely randomized design
#'
#' @description Function of the AgroR package for joint analysis of experiments conducted in a completely randomized design with a qualitative or quantitative factor with balanced data.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param repet Numerical or complex vector with repetitions
#' @param local Numeric or complex vector with locations or times
#' @param response Numerical vector containing the response of the experiment.
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param title Graph title
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param dec Number of cells
#' @param color When the columns are different colors (Set fill-in argument as "trat")
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angulo x-axis scale text rotation
#' @param textsize Font size
#' @param family Font family
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @return Returns the assumptions of the analysis of variance, the assumption of the joint analysis by means of a QMres ratio matrix, the analysis of variance, the multiple comparison test or regression.
#' @keywords DIC
#' @keywords Joint Analysis
#' @references
#'
#' FERREIRA, P. V. Estatistica experimental aplicada a agronomia. Edufal, 2018.
#'
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel & Torry & Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Practical Nonparametrics Statistics. W.J. Conover, 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentação em Genética e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#' @export
#' @examples
#' library(AgroR)
#' data(mirtilo)
#' attach(mirtilo)
#' conjdic(trat, bloco, exp, resp)

conjdic=function(trat,
                 repet,
                 local,
                 response,
                 transf=1,
                 norm="sw",
                 homog="bt",
                 mcomp="tukey",
                 quali=TRUE,
                 alpha.f=0.05,
                 alpha.t=0.05,
                 grau=NA,
                 theme=theme_bw(),
                 ylab="response",
                 title="",
                 xlab="",
                 color="rainbow",
                 fill="lightblue",
                 angulo=0,
                 textsize=12,
                 dec=3,
                 family="sans",
                 errorbar=TRUE){
  sup=0.2*mean(response, na.rm=TRUE)
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("lmtest")
  if(transf=="1"){resp=response}else{resp=(response^transf-1)/transf}
  if(transf=="0"){resp=log(response)}
  if(transf=="0.5"){resp=sqrt(response)}
  if(transf=="-0.5"){resp=1/sqrt(response)}
  if(transf=="-1"){resp=1/response}
  tratamento=as.factor(trat)
  bloco=as.factor(repet)
  local=as.factor(local)
  a = anova(aov(resp ~ local + tratamento + local:tratamento))[c(3:4), ]
  b = summary(aov(resp ~ local + local:bloco + tratamento + Error(local/tratamento)))
  c = aov(resp ~ local + local:bloco + tratamento + local:tratamento)
  dados=data.frame(resp,response,tratamento,local,bloco)
  anova=c()
  tukey=c()
  graficos=list()
  qmres=data.frame(QM=NA)
  for(i in 1:length(levels(local))){
    anova[[i]]=anova(aov(resp~tratamento, data=dados[dados$local==levels(dados$local)[i],]))
    qm=anova[[i]]$`Mean Sq`[2]
    qmres[i,]=c(qm)
    qmres=as.vector(qmres)
    names(anova)[i]=levels(local)[i]
    aov1=aov(resp~tratamento, data=dados[dados$local==levels(dados$local)[i],])}
  matriza=matrix(rep(qmres[[1]],e=length(qmres[[1]])),
           ncol=length(qmres[[1]]))
  matrizb=matrix(rep(qmres[[1]],length(qmres[[1]])),
           ncol=length(qmres[[1]]))
  ratio=matriza/matrizb
  rownames(ratio)=levels(local)
  colnames(ratio)=levels(local)
  razao=data.frame(resp=c(ratio),
                   var1=rep(rownames(ratio),e=length(rownames(ratio))),
                   var2=rep(colnames(ratio),length(colnames(ratio))))
  ratioplot=ggplot(razao,
                   aes(x=razao$var2,
                       y=razao$var1,
                       fill=razao$resp))+
    geom_tile(color="gray50",size=1)+
    scale_x_discrete(position = "top")+
    scale_fill_distiller(palette = "RdBu",direction = 1)+
    ylab("Numerator")+xlab("Denominator")+
    geom_label(aes(label=format(razao$resp,digits=2)),fill="white")+
    labs(fill="ratio")+
    theme(axis.text = element_text(size=12,color="black"),
          legend.text = element_text(size=12),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(caption = "The ratio must be less than 7 (Ferreira et al., 2018)",
         title="Matrix of average square of the residue")
  print(ratioplot)

  (QMRES=as.vector(qmres$QM))
  (qmresmedio=max(QMRES)/min(QMRES))

  d=aov(resp~tratamento*local+bloco)
  if(norm=="sw"){norm1 = shapiro.test(d$res)}
  if(norm=="li"){norm1=nortest::lillie.test(d$residuals)}
  if(norm=="ad"){norm1=nortest::ad.test(d$residuals)}
  if(norm=="cvm"){norm1=nortest::cvm.test(d$residuals)}
  if(norm=="pearson"){norm1=nortest::pearson.test(d$residuals)}
  if(norm=="sf"){norm1=nortest::sf.test(d$residuals)}

  if(homog=="bt"){
    homog1 = bartlett.test(d$res ~ trat)
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = leveneTest(d$res~trat)
    statistic=homog1$`F value`[1]
    phomog=homog1$`Pr(>F)`[1]
    method="Levene's Test (center = median)(F)"
    names(homog1)=c("Df", "F value","p.value")}

  indep = dwtest(d)

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Normality of errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  normal=data.frame(Method=paste(norm1$method,"(",names(norm1$statistic),")",sep=""),
                    Statistic=norm1$statistic,
                    "p-value"=norm1$p.value)
  rownames(normal)=""
  print(normal)
  cat("\n")
  message(if(norm1$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors do not follow a normal distribution"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Homogeneity of Variances")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  homoge=data.frame(Method=method,
                    Statistic=statistic,
                    "p-value"=phomog)
  rownames(homoge)=""
  print(homoge)
  cat("\n")
  message(if(homog1$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous\n"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Independence from errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  indepe=data.frame(Method=paste(indep$method,"(",
                                 names(indep$statistic),")",sep=""),
                    Statistic=indep$statistic,
                    "p-value"=indep$p.value)
  rownames(indepe)=""
  print(indepe)
  cat("\n")
  message(if(indep$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors are not independent"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Test Homogeneity of experiments")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(qmresmedio)
  message(blue("\nBased on the analysis of variance and homogeneity of experiments, it can be concluded that: "))
  if(qmresmedio<7 && a$`Pr(>F)`[1]>alpha.f){
    message(black("The experiments can be analyzed together"))}else{
    message("Experiments cannot be analyzed together (Separate by experiment)")}
  cat("\n\n")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Anova location and treatment interaction")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(a)
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of variances isolated")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(b)
  cat(green(bold("\n-----------------------------------------------------------------\n")))

  if(a$`Pr(>F)`[1] < alpha.f | qmresmedio > 7){
    for(i in 1:length(levels(local))){
      if(a$`Pr(>F)`[1] < alpha.f && quali==TRUE | qmresmedio > 7 && quali==TRUE){
        anova[[i]]=anova(aov(resp~tratamento, data=dados[dados$local==levels(dados$local)[i],]))
        qm=anova[[i]]$`Mean Sq`[2]
        qmres[i,]=c(qm)
        qmres=as.vector(qmres)
        names(anova)[i]=levels(local)[i]
        aov1=aov(resp~tratamento, data=dados[dados$local==levels(dados$local)[i],])
        if(quali==TRUE){
          if(mcomp=="tukey"){tukey[[i]]=HSD.test(aov1,"tratamento",alpha = alpha.t)$groups
          comp=HSD.test(aov1,"tratamento")$groups}
          if(mcomp=="duncan"){tukey[[i]]=duncan.test(aov1,"tratamento",alpha = alpha.t)$groups
          comp=duncan.test(aov1,"tratamento")$groups}
          if(mcomp=="lsd"){tukey[[i]]=LSD.test(aov1,"tratamento",alpha = alpha.t)$groups
          comp=LSD.test(aov1,"tratamento")$groups}
          if(mcomp=="sk"){
            anova=anova(aov1)
            data=dados[dados$local==levels(dados$local)[i],]
            tukey[[i]]=sk_triple(resp,tratamento,anova$Df[2],anova$`Sum Sq`[2],alpha = alpha.t)
            comp=sk_triple(resp,tratamento,anova$Df[2],anova$`Sum Sq`[2],alpha = alpha.t)}

          if(transf=="1"){}else{tukey[[i]]$respo=with(dados[dados$local==levels(dados$local)[i],],
                                                      tapply(response, tratamento,
                                                             mean, na.rm=TRUE))[rownames(comp)]}
          names(tukey)[i]=levels(local)[i]

          dadosm=data.frame(comp,
                            media=with(dados[dados$local==levels(dados$local)[i],],
                                       tapply(response, tratamento, mean, na.rm=TRUE))[rownames(comp)],
                            desvio=with(dados[dados$local==levels(dados$local)[i],],
                                        tapply(response, tratamento, sd, na.rm=TRUE))[rownames(comp)])
          dadosm$Tratamentos=rownames(dadosm)
          dadosm$limite=dadosm$media+dadosm$desvio
          dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)
          grafico=ggplot(dadosm,
                         aes(x=dadosm$Tratamentos,
                             y=dadosm$media))+
            geom_col(aes(fill=dadosm$Tratamentos),fill=fill,color=1)+
            theme_bw()+
            ylab(ylab)+
            xlab(xlab)+ylim(0,1.5*max(dadosm$limite))+
            geom_errorbar(aes(ymin=dadosm$media-dadosm$desvio,
                              ymax=dadosm$media+dadosm$desvio),
                          color="black",width=0.3)+
            geom_text(aes(y=dadosm$media+dadosm$desvio+sup,
                          label=dadosm$letra))+
            theme(text = element_text(size=textsize,color="black", family = family),
                  axis.title = element_text(size=textsize,color="black", family = family),
                  axis.text = element_text(size=textsize,color="black", family = family),
                  legend.position = "none")
          print(tukey)}
        }
      if(a$`Pr(>F)`[1] < alpha.f && quali==FALSE | qmresmedio > 7 && quali==FALSE){
        data=dados[dados$local==levels(dados$local)[i],]
        dose1=as.numeric(as.character(data$tratamento))
        resp=data$response
        grafico=polynomial(dose1,
                           resp,grau = grau,
                           textsize=textsize,
                           family=family,
                           ylab=ylab,
                           xlab=xlab,
                           theme=theme,
                           posi="top",
                           se=errorbar)}
      graficos[[i]]=grafico}
    }

  if(a$`Pr(>F)`[1] > alpha.f && qmresmedio < 7){
    if(quali==TRUE){
    if(mcomp=="tukey"){
      tukeyjuntos=HSD.test(resp,tratamento,a$Df[1], a$`Mean Sq`[1], alpha = alpha.t)
      if(transf!="1"){tukeyjuntos$groups$repo=tapply(response, tratamento, mean,
                                                     na.rm=TRUE)[rownames(tukeyjuntos$groups)]}
      tukeyjuntos=tukeyjuntos$groups}
    if(mcomp=="duncan"){
      tukeyjuntos=duncan.test(resp,tratamento,a$Df[1], a$`Mean Sq`[1], alpha = alpha.t)
      if(transf!="1"){tukeyjuntos$groups$repo=tapply(response, tratamento, mean,
                                                     na.rm=TRUE)[rownames(tukeyjuntos$groups)]}
      tukeyjuntos=tukeyjuntos$groups}
    if(mcomp=="lsd"){
      tukeyjuntos=LSD.test(resp,tratamento,a$Df[1], a$`Mean Sq`[1], alpha = alpha.t)
      if(transf!="1"){tukeyjuntos$groups$repo=tapply(response, tratamento, mean,
                                                     na.rm=TRUE)[rownames(tukeyjuntos$groups)]}
      tukeyjuntos=tukeyjuntos$groups}
    if(mcomp=="sk"){
      tukeyjuntos=sk_triple(resp,tratamento,a$Df[1], a$`Sum Sq`[1], alpha = alpha.t)
      colnames(tukeyjuntos)=c("resp","groups")
      if(transf!="1"){tukeyjuntos$respo=tapply(response, tratamento, mean,
                                               na.rm=TRUE)[rownames(tukeyjuntos)]}}

    # ================================
    # data.frame para grafico
    # ================================
    dadosm=data.frame(tukeyjuntos,
                      media=tapply(response, tratamento, mean, na.rm=TRUE)[rownames(tukeyjuntos)],
                      desvio=tapply(response, tratamento, sd, na.rm=TRUE)[rownames(tukeyjuntos)])
    dadosm$Tratamentos=rownames(dadosm)
    dadosm$limite=dadosm$media+dadosm$desvio
    dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)

    grafico1=ggplot(dadosm,
                    aes(x=dadosm$Tratamentos,y=dadosm$media))
    if(fill=="trat"){grafico1=grafico1+
      geom_col(aes(fill=dadosm$Tratamentos),color=1)}
    else{grafico1=grafico1+
      geom_col(aes(fill=dadosm$Tratamentos),fill=fill,color=1)}
    if(errorbar==TRUE){grafico1=grafico1+
      geom_text(aes(y=dadosm$media+sup+if(sup<0){-dadosm$desvio}else{dadosm$desvio},label=dadosm$letra),family=family)}
    if(errorbar==FALSE){grafico1=grafico1+
      geom_text(aes(y=dadosm$media+sup,label=dadosm$letra),family=family)}
    if(errorbar==TRUE){grafico1=grafico1+
      geom_errorbar(data=dadosm,aes(ymin=dadosm$media-dadosm$desvio,
                                    ymax=dadosm$media+dadosm$desvio,color=1),
                    color="black", width=0.3)}
    grafico1=grafico1+theme+
      ylab(ylab)+
      xlab(xlab)+
      theme(text = element_text(size=textsize,color="black",family=family),
            axis.text = element_text(size=textsize,color="black",family=family),
            axis.title = element_text(size=textsize,color="black",family=family),
            legend.position = "none")
    if(angulo !=0){grafico1=grafico1+theme(axis.text.x=element_text(hjust = 1.01,angle = angulo))}
    print(tukeyjuntos)
    print(grafico1)
    graficos=list(grafico1)
    }
    if(quali==FALSE){grafico1=polynomial(as.numeric(as.character(tratamento)),
                           response,grau = grau,
                           textsize=textsize,
                           family=family,
                           ylab=ylab,
                           xlab=xlab,
                           theme=theme,
                           posi="top",
                           se=errorbar)
    graficos=list(grafico1)
    }
    }
  cat(if(transf=="1"){}else{blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})

  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 && homog1$p.value<0.05){cat(red("\n \nWarning!!! Your analysis is not valid, suggests using a non-parametric test and try to transform the data"))}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){cat(red("\n \nWarning!!! Your analysis is not valid, suggests using a non-parametric test"))}
  graph=as.list(graficos)
}
