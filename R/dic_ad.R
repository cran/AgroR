#' Analysis: Completely randomized design with an additional treatment for quantitative factor
#'
#' @description Statistical analysis of experiments conducted in a completely randomized with an additional treatment and balanced design with a factor considering the fixed model.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param response Numerical vector containing the response of the experiment.
#' @param responsead Numerical vector with additional treatment responses
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param family Font family
#' @param pointsize Point size
#' @param linesize line size (Trendline and Error Bar)
#' @param width.bar width of the error bars of a regression graph.
#' @param posi Legend position
#' @param point Defines whether to plot mean ("mean"), mean with standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se"). For quali=FALSE or quali=TRUE.
#' @note In some experiments, the researcher may study a quantitative factor, such as fertilizer doses, and present a control, such as a reference fertilizer, treated as a qualitative control. In these cases, there is a difference between considering only the residue in the unfolding of the polynomial, removing or not the qualitative treatment, or since a treatment is excluded from the analysis. In this approach, the residue used is also considering the qualitative treatment, a method similar to the factorial scheme with additional control.
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk ("sw"), Lilliefors ("li"), Anderson-Darling ("ad"), Cramer-von Mises ("cvm"), Pearson ("pearson") and Shapiro-Francia ("sf")), the test of homogeneity of variances (Bartlett ("bt") or Levene ("levene")), the test of independence of Durbin-Watson errors, adjustment of regression models up to grade 3 polynomial. The function also returns a standardized residual plot.
#' @keywords DIC
#' @keywords additional treatment
#' @export
#' @examples
#' doses = c(rep(c(1:5),e=3))
#' resp = c(3, 4, 3, 5, 5, 6, 7, 7, 8, 4, 4, 5, 2, 2, 3)
#' dic.ad(doses, resp, rnorm(3,6,0.1),grau=2)

dic.ad=function(trat,
                response,
                responsead,
                grau = 1,
                norm="sw",
                homog="bt",
                alpha.f=0.05,
                theme=theme_classic(),
                ylab="response",
                xlab="independent",
                family="sans",
                posi="top",
                pointsize=4.5,
                linesize=0.8,
                width.bar=NA,
                point="mean_sd"){
  if(is.na(width.bar)==TRUE){width.bar=0.1*mean(trat)}
  if(is.na(grau)==TRUE){grau=1}
  trat1=as.factor(trat)
  mod=aov(resp~trat1)
  an=anova(mod)
  trati=as.factor(c(trat,rep("Controle",length(responsead))))
  mod1=aov(c(resp,responsead)~trati)
  an1=anova(mod1)
  anava1=rbind(an[1,],an1[1,],an1[2,])
  anava1$Df[2]=1
  anava1$`Sum Sq`[2]=anava1$`Sum Sq`[2]-sum(anava1$`Sum Sq`[1])
  anava1$`Mean Sq`[2]=anava1$`Sum Sq`[2]/anava1$Df[2]
  anava1$`F value`[1:2]=anava1$`Mean Sq`[1:2]/anava1$`Mean Sq`[3]
  rownames(anava1)[1:2]=c("Factor","Factor vs control")
  for(i in 1:nrow(anava1)-1){
    anava1$`Pr(>F)`[i]=1-pf(anava1$`F value`[i],anava1$Df[i],anava1$Df[3])
  }
  respad=mod1$residuals/sqrt(anava1$`Mean Sq`[3])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}
  if(norm=="sw"){norm1 = shapiro.test(mod1$res)}
  if(norm=="li"){norm1=nortest::lillie.test(mod1$residuals)}
  if(norm=="ad"){norm1=nortest::ad.test(mod1$residuals)}
  if(norm=="cvm"){norm1=nortest::cvm.test(mod1$residuals)}
  if(norm=="pearson"){norm1=nortest::pearson.test(mod1$residuals)}
  if(norm=="sf"){norm1=nortest::sf.test(mod1$residuals)}
  if(homog=="bt"){
    homog1 = bartlett.test(mod1$res ~ trati)
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = levenehomog(mod1$res~trati)[1,]
    statistic=homog1$`F value`[1]
    phomog=homog1$`Pr(>F)`[1]
    method="Levene's Test (center = median)(F)"
    names(homog1)=c("Df", "statistic","p.value")}
  indep = dwtest(mod1)
  Ids=ifelse(respad>3 | respad<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(respad,Ids),aes(y=respad,x=1:length(respad)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(respad),label=1:length(respad),color=Ids,size=4)+
    scale_x_continuous(breaks=1:length(respad))+
    theme_classic()+theme(axis.text.y = element_text(size=12),
                          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)
  print(residplot)
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
  cat(green(bold("\n-----------------------------------------------------------------\n")))
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
  message(if(indep$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected.Therefore, errors are not independent"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV (%) = ",round(sqrt(anava1$`Mean Sq`[3])/mean(response,na.rm=TRUE)*100,2)))
  cat(paste("\nMean = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian = ",round(median(response,na.rm=TRUE),4)))
  cat("\nPossible outliers = ", out)
  cat("\n")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(anava1)
  print(if(anava1$`Pr(>F)`[1]>alpha.f){"H0 is not rejected"})
  a=AgroR::polynomial(trat,response,DFres = anava1$Df[3],
                      SSq = anava1$`Sum Sq`[3],
                      ylab = ylab,
                      xlab=xlab,
                      theme = theme,
                      point = point,
                      grau = grau,
                      posi = posi,
                      family = family,
                      pointsize = pointsize,
                      linesize = linesize,
                      width.bar = width.bar)}


