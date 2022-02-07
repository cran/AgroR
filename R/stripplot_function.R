#' Analysis: DBC experiments in strip-plot
#' @description Analysis of an experiment conducted in a block randomized design in a strit-plot scheme using fixed effects analysis of variance.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with plot levels
#' @param f2 Numeric or complex vector with subplot levels
#' @param block Numeric or complex vector with blocks
#' @param response Numeric vector with responses
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param constant Add a constant for transformation (enter value)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param textsize Font size (\emph{default} is 12)
#' @param labelsize Label size (\emph{default} is 4)
#' @import ggplot2
#' @importFrom crayon green
#' @importFrom crayon bold
#' @importFrom crayon italic
#' @importFrom crayon red
#' @importFrom crayon blue
#' @import stats
#' @keywords DBC
#' @export
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
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett). The function also returns a standardized residual plot.
#' @examples
#'
#' #===================================
#' # Example tomate
#' #===================================
#' # Obs. Consider that the "tomato" experiment is a block randomized design in strip-plot.
#' library(AgroR)
#' data(tomate)
#' with(tomate, STRIPLOT(parc, subp, bloco, resp))

STRIPLOT=function(f1,
                  f2,
                  block,
                  response,
                  norm="sw",
                  alpha.f=0.05,
                  transf=1,
                  textsize=12,
                  labelsize=4,
                  constant=0){
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")

  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  fator1=f1
  fator2=f2
  fator1a=fator1
  fator2a=fator2
  bloco=block
  fac = c("F1", "F2")
  cont <- c(1, 3)

  Fator1 <- factor(fator1, levels = unique(fator1))
  Fator2 <- factor(fator2, levels = unique(fator2))
  bloco <- factor(block)
  nv1 <- length(summary(Fator1))
  nv2 <- length(summary(Fator2))
  lf1 <- levels(Fator1)
  lf2 <- levels(Fator2)
  num=function(x){as.numeric(x)}
  graph=data.frame(Fator1,Fator2,resp)
  mod=aov(resp~f1*f2+f2*block+block:f1)
  anava=anova(mod)
  anava=anava[c(3,1,6,2,5,4,7),]
  anava$`F value`[1]=anava$`Mean Sq`[1]/(anava$`Mean Sq`[3]+anava$`Mean Sq`[5]-anava$`Mean Sq`[7])
  anava$`F value`[2]=anava$`Mean Sq`[2]/anava$`Mean Sq`[3]
  anava$`F value`[4]=anava$`Mean Sq`[4]/anava$`Mean Sq`[5]
  anava[c(3,5),4:5]=NA
  anava$`Pr(>F)`[1]=1-pf(anava$`F value`[1],anava$Df[1],anava$Df[3])
  anava$`Pr(>F)`[2]=1-pf(anava$`F value`[2],anava$Df[2],anava$Df[3])
  anava$`Pr(>F)`[4]=1-pf(anava$`F value`[4],anava$Df[4],anava$Df[5])
  rownames(anava)=c("Block","F1","Error A","F2","Error B","F1:F2","Residuals")
  resids=residuals(mod,scaled=TRUE)
  Ids=ifelse(resids>3 | resids<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(resids,Ids),
                   aes(y=resids,x=1:length(resids)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(resids),label=1:length(resids),
              color=Ids,size=labelsize)+
    scale_x_continuous(breaks=1:length(resids))+
    theme_classic()+theme(axis.text.y = element_text(size=textsize),
                          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)
  # Normalidade dos erros
  if(norm=="sw"){norm1 = shapiro.test(resid(mod))}
  if(norm=="li"){norm1=lillie.test(resid(mod))}
  if(norm=="ad"){norm1=ad.test(resid(mod))}
  if(norm=="cvm"){norm1=cvm.test(resid(mod))}
  if(norm=="pearson"){norm1=pearson.test(resid(mod))}
  if(norm=="sf"){norm1=sf.test(resid(mod))}
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

  homog1=bartlett.test(resid(mod)~Fator1)
  homog2=bartlett.test(resid(mod)~Fator2)
  homog3=bartlett.test(resid(mod)~paste(Fator1,Fator2))
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Homogeneity of Variances")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Interaction\n")))
  statistic3=homog3$statistic
  phomog3=homog3$p.value
  method3=paste("Bartlett test","(",names(statistic3),")",sep="")
  homoge3=data.frame(Method=method3,
                     Statistic=statistic3,
                     "p-value"=phomog3)
  rownames(homoge3)=""
  print(homoge3)
  cat("\n")
  message(if(homog3$p.value[1]>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV1 (%) = ",round(sqrt(anava$`Mean Sq`[3])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nCV2 (%) = ",round(sqrt(anava$`Mean Sq`[5])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nCV3 (%) = ",round(sqrt(anava$`Mean Sq`[7])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nMean = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian = ",round(median(response,na.rm=TRUE),4)))
  #cat("\nPossible outliers = ", out)
  cat("\n")

  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  anava=data.frame(anava)
  colnames(anava)=c("Df","Sum Sq ","Mean Sq","F value","Pr(>F)")
  print(as.matrix(anava),na.print="",quote = FALSE)
  if(transf==1 && norm1$p.value<0.05 |  transf==1 &&homog1$p.value<0.05){
    message("\n \nYour analysis is not valid, suggests using a try to transform the data\n")}else{}

  message(if(transf !=1){blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})}
