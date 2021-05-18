#' Analysis: Completely randomized design
#'
#' @description Statistical analysis of experiments conducted in a completely randomized and balanced design with a factor considering the fixed model.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param response Numerical vector containing the response of the experiment.
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param transf Applies data transformation (\emph{default} is 1; for log consider 0)
#' @param test "parametric" - Parametric test or "noparametric" - non-parametric test
#' @param p.adj Method for adjusting p values for Kruskal-Wallis ("none","holm","hommel", "hochberg", "bonferroni", "BH", "BY", "fdr")
#' @param geom Graph type (columns, boxes or segments)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param CV Plotting the coefficient of variation and p-value of Anova (\emph{default} is TRUE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param textsize Font size
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angle x-axis scale text rotation
#' @param family Font family
#' @param dec Number of cells
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param posi Legend position
#' @import ggplot2
#' @import stats
#' @importFrom ScottKnott SK
#' @importFrom ScottKnott SK.nest
#' @importFrom crayon green
#' @importFrom crayon bold
#' @importFrom crayon italic
#' @importFrom crayon red
#' @importFrom crayon blue
#' @importFrom nortest lillie.test
#' @importFrom nortest ad.test
#' @importFrom nortest cvm.test
#' @importFrom nortest pearson.test
#' @importFrom nortest sf.test
#' @importFrom lmtest dwtest
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom graphics abline
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggrepel geom_text_repel
#' @importFrom emmeans emmeans
#' @importFrom multcomp cld
#' @importFrom ARTool artlm
#' @importFrom ARTool art
#' @importFrom lme4 lmer
#' @importFrom grid grobTree
#' @importFrom grid rectGrob
#' @importFrom grid textGrob
#' @importFrom grid gpar
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics par
#' @importFrom utils read.table
#' @importFrom stringr str_trim
#' @importFrom reshape2 melt
#' @importFrom Hmisc rcorr
#' @importFrom cowplot plot_grid
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
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
#'
#' Mendiburu, F., & de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#'
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett or Levene), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. Non-parametric analysis can be used by the Kruskal-Wallis test. The column, segment or box chart for qualitative treatments is also returned. The function also returns a standardized residual plot.
#' @keywords DIC
#' @keywords Experimental
#' @seealso \link{DBC} \link{DQL}
#' @export
#' @examples
#' library(AgroR)
#' data(pomegranate)
#' attach(pomegranate)
#' DIC(trat, WL) # tukey
#' DIC(trat, WL, mcomp = "sk")
#' DIC(trat, WL, mcomp = "duncan")
#' DIC(trat, WL, test = "noparametric")
#' DIC(trat, WL, transf = 0)
#' DIC(trat, WL, geom="point")
#' DIC(trat, WL, ylab = "Weight loss (%)", xlab="Treatments")
#' data("phao")
#' attach(phao)
#' DIC(dose,comp,quali=FALSE,grau=2)

######################################################################################
## Analise de variancia para experimentos em DIC
######################################################################################
DIC <- function(trat,
                response,
                norm="sw",
                homog="bt",
                mcomp="tukey",
                quali=TRUE,
                alpha.f=0.05,
                alpha.t=0.05,
                grau=1,
                transf=1,
                test="parametric",
                p.adj="holm",
                geom="bar",
                theme=theme_bw(),
                sup=NA,
                CV=TRUE,
                ylab="Response",
                xlab="",
                fill="lightblue",
                angle=0,
                family="sans",
                textsize=12,
                dec=3,
                addmean=TRUE,
                errorbar=TRUE,
                posi="top"){
  if(is.na(sup==TRUE)){sup=0.1*mean(response)}
  requireNamespace("gridExtra")
  requireNamespace("nortest")
  requireNamespace("lmtest")
  requireNamespace("ScottKnott")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  if(test=="parametric"){
  if(transf=="1"){resp=response}else{resp=(response^transf-1)/transf}
  if(transf=="0"){resp=log(response)}
  if(transf=="0.5"){resp=sqrt(response)}
  if(transf=="-0.5"){resp=1/sqrt(response)}
  if(transf=="-1"){resp=1/response}
  trat=as.factor(trat)
  a = anova(aov(resp ~ trat))
  aa = summary(aov(resp ~ trat))
  b = aov(resp ~ trat)
  anava=a
  colnames(anava)=c("GL","SQ","QM","Fcal","p-value")
  if(norm=="sw"){norm1 = shapiro.test(b$res)}
  if(norm=="li"){norm1=nortest::lillie.test(b$residuals)}
  if(norm=="ad"){norm1=nortest::ad.test(b$residuals)}
  if(norm=="cvm"){norm1=nortest::cvm.test(b$residuals)}
  if(norm=="pearson"){norm1=nortest::pearson.test(b$residuals)}
  if(norm=="sf"){norm1=nortest::sf.test(b$residuals)}
  if(homog=="bt"){
    homog1 = bartlett.test(b$res ~ trat)
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = leveneTest(b$res~trat)
    statistic=homog1$`F value`[1]
    phomog=homog1$`Pr(>F)`[1]
    method="Levene's Test (center = median)(F)"
    names(homog1)=c("Df", "F value","p.value")}
  indep = dwtest(b)
  plot(b$residuals/sqrt(a$`Mean Sq`[2]),
       ylab = "Standardized Residuals",
       las = 1,pch = 16,col = "blue")
  abline(h=c(0,3,-3),lty=c(1,2,2),col="red")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Normality of errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  normal=data.frame(Method=paste(norm1$method,"(",names(norm1$statistic),")",sep=""),
                    Statistic=norm1$statistic,
                    "p-value"=norm1$p.value)
  rownames(normal)=""
  print(normal)
  cat("\n")
  message(if(norm1$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors do not follow a normal distribution"})
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Homogeneity of Variances")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  homoge=data.frame(Method=method,
                    Statistic=statistic,
                    "p-value"=phomog)
  rownames(homoge)=""
  print(homoge)
  cat("\n")
  message(if(homog1$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level,hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected.Therefore, the variances are not homogeneous"})
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Independence from errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  indepe=data.frame(Method=paste(indep$method,"(",
                                 names(indep$statistic),")",sep=""),
                    Statistic=indep$statistic,
                    "p-value"=indep$p.value)
  rownames(indepe)=""
  print(indepe)
  cat("\n")
  message(if(indep$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected.Therefore, errors are not independent"})
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(anava)
  cat("\n\n")
  message(if (a$`Pr(>F)`[1]<alpha.f){
    black("As the calculated p-value, it is less than the 5% significance level.The hypothesis H0 of equality of means is rejected. Therefore, at least two treatments differ")}
      else {"As the calculated p-value is greater than the 5% significance level, H0 is not rejected"})
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  if(quali==TRUE){cat(green(bold("Multiple Comparison Test")))}else{cat(green(bold("Regression")))}
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  if(quali==TRUE){
  if(mcomp=="tukey"){
    letra <- HSD.test(b, "trat", alpha=alpha.t)
    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
  if(mcomp=="sk"){
    letra=SK(b,"trat",sig.level=alpha.t)
    letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
    #letra1$resp=as.numeric(as.character(letra1$resp))
    }
  if(mcomp=="duncan"){
    letra <- duncan.test(b, "trat", alpha=alpha.t)
    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
  if(mcomp=="lsd"){
      letra <- LSD.test(b, "trat", alpha=alpha.t)
      letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
  media = tapply(response, trat, mean, na.rm=TRUE)
  if(transf=="1"){letra1}else{letra1$respO=media[rownames(letra1)]}
  print(if(a$`Pr(>F)`[1]<0.05){letra1}else{"H0 is not rejected"})
  cat("\n")
  message(if(transf=="1"){}else{blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})
  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
    message("\nYour analysis is not valid, suggests using a non-parametric test and try to transform the data")
    }
  else{}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){cat(red("\nWarning!!! Your analysis is not valid, suggests using a non-parametric test"))}else{}
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV (%) = ",round(sqrt(a$`Mean Sq`[2])/mean(resp)*100,2)))
  cat(paste("\nR-squared = ",round(a$`Mean Sq`[1]/(a$`Mean Sq`[2]+a$`Mean Sq`[1]),2)))
  cat(paste("\nMean = ",round(mean(resp),4)))
  cat(paste("\nMedian = ",round(median(resp),4)))
  dadosm=data.frame(letra1,
                    media=tapply(response, trat, mean, na.rm=TRUE)[rownames(letra1)],
                    desvio=tapply(response, trat, sd, na.rm=TRUE)[rownames(letra1)])
  dadosm$trats=factor(rownames(dadosm),levels = unique(trat))
  dadosm$limite=dadosm$media+dadosm$desvio
  dadosm=dadosm[unique(as.character(trat)),]
  if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
  if(addmean==FALSE){dadosm$letra=dadosm$groups}
  trats=dadosm$trats
  limite=dadosm$limite
  media=dadosm$media
  desvio=dadosm$desvio
  letra=dadosm$letra
  if(geom=="bar"){grafico=ggplot(dadosm,aes(x=trats,y=media))
    if(fill=="trat"){grafico=grafico+
      geom_col(aes(fill=trats),color=1)}
    else{grafico=grafico+
      geom_col(aes(fill=trats),fill=fill,color=1)}
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                  label=letra),family=family)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,label=letra),family=family)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                  ymax=media+desvio,color=1),
                  color="black",width=0.3)}}
  if(geom=="point"){grafico=ggplot(dadosm,aes(x=trats,
                                              y=media))
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                  label=letra),family=family)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,
                  label=letra),family=family)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,
                  aes(ymin=media-desvio,
                      ymax=media+desvio,color=1),
                  color="black",width=0.3)}
  if(fill=="trat"){grafico=grafico+
    geom_point(aes(color=trats),size=5)}
  else{grafico=grafico+
    geom_point(aes(color=trats),
               color="black",
               fill=fill,shape=21,size=5)}}
  if(geom=="box"){
  datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),
                    response)
  dadosm2=data.frame(letra1,
                     superior=tapply(response, trat, mean, na.rm=TRUE)[rownames(letra1)])
  dadosm2$trats=rownames(dadosm2)
  dadosm2=dadosm2[unique(as.character(trat)),]
  dadosm2$limite=dadosm$media+dadosm$desvio
  dadosm2$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)
  trats=dadosm2$trats
  limite=dadosm2$limite
  superior=dadosm2$superior
  letra=dadosm2$letra
  stat_box=ggplot(datam1,aes(x=trats,y=response))+geom_boxplot()
  superior=ggplot_build(stat_box)$data[[1]]$ymax
  dadosm2$superior=superior+sup
  grafico=ggplot(datam1,aes(x=trats,y=response))
  if(fill=="trat"){grafico=grafico+geom_boxplot(aes(fill=trats))}
  else{grafico=grafico+
    geom_boxplot(aes(fill=trats),fill=fill)}
  grafico=grafico+
    geom_text(data=dadosm2,
              aes(y=superior,
                  label=letra),
              family = family)}
  grafico=grafico+
    theme+
    ylab(ylab)+
    xlab(xlab)+
    theme(text = element_text(size=textsize,color="black", family = family),
          axis.text = element_text(size=textsize,color="black", family = family),
          axis.title = element_text(size=textsize,color="black", family = family),
          legend.position = "none")
  if(angle !=0){grafico=grafico+
    theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
  if(CV==TRUE){grafico=grafico+
    labs(caption=paste("p-value = ", if(a$`Pr(>F)`[1]<0.0001){paste("<", 0.0001)}
                                                  else{paste("=", round(a$`Pr(>F)`[1],4))},"; CV = ",
                                                  round(abs(sqrt(a$`Mean Sq`[2])/mean(resp))*100,2),"%"))}
  grafico=as.list(grafico)
  }
  if(quali==FALSE){
  trat=as.numeric(as.character(trat))
  if(grau==1){graph=polynomial(trat,response, grau = 1,textsize=textsize,xlab=xlab,ylab=ylab, family=family,posi=posi)}
  if(grau==2){graph=polynomial(trat,response, grau = 2,textsize=textsize,xlab=xlab,ylab=ylab, family=family,posi=posi)}
  if(grau==3){graph=polynomial(trat,response, grau = 3,textsize=textsize,xlab=xlab,ylab=ylab, family=family,posi=posi)}
  grafico=graph[[1]]
  }}
  if(test=="noparametric"){
    krusk=kruskal(response,trat,p.adj = p.adj)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Statistics")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(krusk$statistics)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Parameters")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(krusk$parameters)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Multiple Comparison Test")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    saida=cbind(krusk$means[,c(1,3)],krusk$groups[rownames(krusk$means),])
    colnames(saida)=c("Mean","SD","Rank","Groups")
    print(saida)
    dadosm=data.frame(krusk$means,krusk$groups[rownames(krusk$means),])
    dadosm$trats=factor(rownames(dadosm),levels = unique(trat))
    dadosm$media=tapply(response,trat,mean, na.rm=TRUE)[rownames(krusk$means)]
    dadosm$std=tapply(response,trat,sd, na.rm=TRUE)[rownames(krusk$means)]
    if(addmean==TRUE){dadosm$letra=paste(format(dadosm$response,digits = dec),dadosm$groups)}
    if(addmean==FALSE){dadosm$letra=dadosm$groups}
    trats=dadosm$trats
    limite=dadosm$limite
    media=dadosm$media
    std=dadosm$std
    letra=dadosm$letra
    if(geom=="bar"){grafico=ggplot(dadosm,
                                   aes(x=trats,y=response))
      if(fill=="trat"){grafico=grafico+
        geom_col(aes(fill=trats),color=1)}
      else{grafico=grafico+
        geom_col(aes(fill=trats),fill=fill,color=1)}
    if(errorbar==TRUE){grafico=grafico+
      geom_text(aes(y=media+sup+if(sup<0){-std}else{std},
                    label=letra),family=family)}
    if(errorbar==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,label=letra),family=family)}
    if(errorbar==TRUE){grafico=grafico+
      geom_errorbar(data=dadosm,aes(ymin=response-std,
                                    ymax=response+std,
                                    color=1),
                    color="black",width=0.3)}}
    if(geom=="point"){grafico=ggplot(dadosm,
                                     aes(x=trats,
                                         y=response))
    if(errorbar==TRUE){grafico=grafico+
      geom_text(aes(y=media+sup+if(sup<0){-std}else{std},
                    label=letra),
                family=family)}
    if(errorbar==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,
                    label=letra),
                family=family)}
    if(errorbar==TRUE){grafico=grafico+
      geom_errorbar(data=dadosm,
                    aes(ymin=response-std,
                        ymax=response+std,
                        color=1),
                    color="black",width=0.3)}
    if(fill=="trat"){grafico=grafico+
      geom_point(aes(color=trats),size=5)}
    else{grafico=grafico+
      geom_point(aes(color=trats),
                 color="black",
                 fill=fill,shape=21,size=5)}}
    if(geom=="box"){
    datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),response)
    dadosm2=data.frame(krusk$means)
    dadosm2$trats=rownames(dadosm2)
    dadosm2$limite=dadosm2$response+dadosm2$std
    dadosm2$letra=paste(format(dadosm2$response,digits = dec),
                        dadosm$groups)
    dadosm2=dadosm2[unique(as.character(trat)),]
    trats=dadosm2$trats
    limite=dadosm2$limite
    letra=dadosm2$letra
    stat_box=ggplot(datam1,aes(x=trats,y=response))+geom_boxplot()
    superior=ggplot_build(stat_box)$data[[1]]$ymax
    dadosm2$superior=superior+sup

    grafico=ggplot(datam1,
                   aes(x=trats,
                       y=response))
      if(fill=="trat"){grafico=grafico+
        geom_boxplot(aes(fill=1))}
    else{grafico=grafico+
      geom_boxplot(aes(fill=trats),fill=fill)}
    grafico=grafico+
      geom_text(data=dadosm2,
                aes(y=superior,
                    label=letra),
                family = family)}
    grafico=grafico+theme+
      ylab(ylab)+
      xlab(xlab)+
      theme(text = element_text(size=textsize,color="black", family = family),
            axis.title = element_text(size=textsize,color="black", family = family),
            axis.text = element_text(size=textsize,color="black", family = family),
            legend.position = "none")
    if(angle !=0){grafico=grafico+theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
}
  print(grafico)
  graficos=list(grafico)#[[1]]
}
