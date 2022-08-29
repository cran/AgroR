#' Analysis: Joint analysis of experiments in randomized block design in scheme factorial double
#'
#' @description Function of the AgroR package for joint analysis of experiments conducted in a randomized factorial double in block design with balanced data. The function generates the joint analysis through two models. Model 1: F-test of the effects of Factor 1, Factor 2 and F1 x F2 interaction are used in reference to the mean square of the interaction with the year. Model 2: F-test of the Factor 1, Factor 2 and F1 x F2 interaction effects are used in reference to the mean square of the residual.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param block Numerical or complex vector with blocks
#' @param experiment Numeric or complex vector with locations or times
#' @param response Numerical vector containing the response of the experiment.
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param model Define model of the analysis of variance
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @return Returns the assumptions of the analysis of variance, the assumption of the joint analysis by means of a QMres ratio matrix and analysis of variance
#'
#' @note The function is still limited to analysis of variance and assumptions only.
#'
#' @references
#'
#' Ferreira, P. V. Estatistica experimental aplicada a agronomia. Edufal, 2018.
#'
#' Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Practical Nonparametrics Statistics. W.J. Conover, 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' @keywords Double factorial
#' @keywords Joint Analysis
#' @export
#' @examples
#' library(AgroR)
#' ano=factor(rep(c(2018,2019,2020),e=48))
#' f1=rep(rep(c("A","B","C"),e=16),3)
#' f2=rep(rep(rep(c("a1","a2","a3","a4"),e=4),3),3)
#' resp=rnorm(48*3,10,1)
#' bloco=rep(c("b1","b2","b3","b4"),36)
#' dados=data.frame(ano,f1,f2,resp,bloco)
#' conjfat2dbc(f1,f2,bloco,ano,resp, model=1)

conjfat2dbc=function(f1,
                     f2,
                     block,
                     experiment,
                     response,
                     transf = 1,
                     model=1,
                     norm = "sw",
                     homog = "bt",
                     alpha.f = 0.05,
                     alpha.t = 0.05) {
  sup=0.2*mean(response, na.rm=TRUE)
  requireNamespace("crayon")
  requireNamespace("ggplot2")

  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  f1=factor(f1,unique(f1))
  f2=factor(f2,unique(f2))
  bloco=as.factor(block)
  local=as.factor(experiment)

  #============================================================
  dados=data.frame(f1,f2,bloco,local,resp,response)
  anova=c()
  nlocal=length(levels(local))
  qmres=data.frame(QM=1:nlocal)

  #============================================================
  for(i in 1:length(levels(local))){
    anova[[i]]=anova(aov(resp~f1*f2+bloco,
                         data=dados[dados$local==levels(dados$local)[i],]))
    qm=anova[[i]]$`Mean Sq`[5]
    qmres[i,1]=c(qm)
    names(anova)[i]=levels(local)[i]}
  matriza=matrix(rep(qmres[[1]],e=length(qmres[[1]])),
                 ncol=length(qmres[[1]]))
  matrizb=matrix(rep(qmres[[1]],length(qmres[[1]])),
                 ncol=length(qmres[[1]]))
  ratio=matriza/matrizb
  rownames(ratio)=levels(local)
  colnames(ratio)=levels(local)
  razao=data.frame(resp1=c(ratio),
                   var1=rep(rownames(ratio),e=length(rownames(ratio))),
                   var2=rep(colnames(ratio),length(colnames(ratio))))
  var1=razao$var1
  var2=razao$var2
  resp1=razao$resp1

  ratioplot=ggplot(razao,
                   aes(x=var2,
                       y=var1,
                       fill=resp1))+
    geom_tile(color="gray50",size=1)+
    scale_x_discrete(position = "top")+
    scale_fill_distiller(palette = "RdBu",direction = 1)+
    ylab("Numerator")+
    xlab("Denominator")+
    geom_label(aes(label=format(resp1,digits=2)),fill="white")+
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
  QMRES=as.vector(qmres$QM)
  qmresmedio=max(QMRES)/min(QMRES)

  # modelo 1
  # modelo tratando f1, f2 e interacao em relação a interacao (f1 x f2 x ano)
  if(model==1){interacao=anova(aov(resp~f1+f2+f1:f2+ano:bloco+f1:f2:ano))
  interacao$`Mean Sq`[6]=mean(QMRES)
  interacao$`Sum Sq`[6]=interacao$`Mean Sq`[6]*interacao$Df[6]
  interacao$`F value`[5]=interacao$`Mean Sq`[5]/interacao$`Mean Sq`[6]
  interacao$`Pr(>F)`[5]=1-pf(interacao$`F value`[5],interacao$Df[6],interacao$Df[6])

  interacao$`F value`[1:4]=interacao$`Mean Sq`[1:4]/interacao$`Mean Sq`[5]
  interacao$`Pr(>F)`[1:4]=1-pf(interacao$`F value`[1:4],interacao$Df[1:4],interacao$Df[5])
  pfint=interacao$`Pr(>F)`[5]}

  # modelo 2
  # todos fixos
  if(model==2){
  interacao <- anova(aov(resp ~ f1*f2+ano:bloco + f1:f2:ano))
  pfint=interacao$`Pr(>F)`[5]}

  #=======================================================================
  d=aov(resp~f1+f2+f1:f2+ano:bloco+f1:f2:ano)
  if(norm=="sw"){norm1 = shapiro.test(d$res)}
  if(norm=="li"){norm1=nortest::lillie.test(d$residuals)}
  if(norm=="ad"){norm1=nortest::ad.test(d$residuals)}
  if(norm=="cvm"){norm1=nortest::cvm.test(d$residuals)}
  if(norm=="pearson"){norm1=nortest::pearson.test(d$residuals)}
  if(norm=="sf"){norm1=nortest::sf.test(d$residuals)}

  if(homog=="bt"){
    homog1 = bartlett.test(d$res ~ paste(f1,f2))
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = levenehomog(d$res~paste(f1,f2))
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
  cat(if(norm1$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
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
  cat(if(homog1$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
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
  cat(black(if(indep$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
            else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors are not independent"}))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Test Homogeneity of experiments")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(qmresmedio)
  cat("\nBased on the analysis of variance and homogeneity of experiments, it can be concluded that: ")
  if(qmresmedio<7 && pfint[1]>alpha.f){
    message(black("The experiments can be analyzed together"))}else{
      message("Experiments cannot be analyzed together (Separate by experiment)")}
  cat("\n\n")
  modres=anova(d)
  respad=d$res/sqrt(modres$`Mean Sq`[6])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}
  resids=d$res/sqrt(modres$`Mean Sq`[6])
  Ids=ifelse(resids>3 | resids<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(resids,Ids),aes(y=resids,x=1:length(resids)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(resids),label=1:length(resids),color=Ids,size=4)+
    scale_x_continuous(breaks=1:length(resids))+
    theme_classic()+theme(axis.text.y = element_text(size=12),
                          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)
  print(residplot)

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(as.matrix(interacao),na.print = "")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
}


