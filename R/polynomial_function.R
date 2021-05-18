#' Analysis: Linear regression graph
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Linear regression analysis of an experiment with a quantitative factor or isolated effect of a quantitative factor
#' @param resp Numerical vector containing the response of the experiment.
#' @param trat Numerical vector with treatments (Declare as numeric)
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param grau Degree of the polynomial (1,2 or 3)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param color Graph color (\emph{default} is gray80)
#' @param posi Legend position
#' @param textsize Font size
#' @param family Font family
#' @param pointsize Point size
#' @param decimal Decimal separate
#' @param ylim y-axis scale
#' @param se Adds confidence interval (\emph{default} is FALSE)
#' @return Returns linear, quadratic or cubic regression analysis.
#' @keywords Regression
#' @keywords Experimental
#' @seealso \link{polynomial2}
#' @export
#' @examples
#' data("phao")
#' attach(phao)
#' polynomial(dose,comp, grau = 2)


polynomial=function(trat,
               resp,
               ylab="Response",
               xlab="Independent",
               grau=NA,
               theme=theme_bw(),
               color="gray80",
               posi="top",
               textsize=12,
               se=FALSE,
               ylim=NA,
               family="sans",
               pointsize=4.5,
               decimal=".")
{requireNamespace("ggplot2")
  if(is.na(grau)==TRUE){grau=1}
  # ================================
  # vetores
  # ================================
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean, na.rm=TRUE)
  mod=c()
  mod1=c()
  mod2=c()
  modm=c()
  mod1m=c()
  mod2m=c()
  text1=c()
  text2=c()
  text3=c()
  mods=c()
  mod1s=c()
  mod2s=c()
  fparcial1=c()
  fparcial2=c()
  fparcial3=c()
  media=tapply(resp, trat, mean, na.rm=TRUE)
  desvio=tapply(resp, trat, sd, na.rm=TRUE)
  dose=tapply(trat, trat, mean, na.rm=TRUE)
  moda=lm(resp~trat)
  mod1a=lm(resp~trat+I(trat^2))
  mod2a=lm(resp~trat+I(trat^2)+I(trat^3))
  mods=summary(moda)$coefficients
  mod1s=summary(mod1a)$coefficients
  mod2s=summary(mod2a)$coefficients
  modm=lm(media~dose)
  mod1m=lm(media~dose+I(dose^2))
  mod2m=lm(media~dose+I(dose^2)+I(dose^3))

  modquali=lm(resp~as.factor(trat))
  fparcial1=anova(modquali,moda)
  fparcial2=anova(modquali,mod1a)
  fparcial3=anova(modquali,mod2a)

  fparcial1=data.frame(abs(fparcial1[2,c(3,4,5,6)]))
  fparcial2=data.frame(abs(fparcial2[2,c(3,4,5,6)]))
  fparcial3=data.frame(abs(fparcial3[2,c(3,4,5,6)]))

  colnames(fparcial1)=c("GL","SQ","F","p-value");rownames(fparcial1)=""
  colnames(fparcial2)=c("GL","SQ","F","p-value");rownames(fparcial2)=""
  colnames(fparcial3)=c("GL","SQ","F","p-value");rownames(fparcial3)=""

  if(grau=="1"){r2=round(summary(modm)$r.squared, 2)}
  if(grau=="2"){r2=round(summary(mod1m)$r.squared, 2)}
  if(grau=="3"){r2=round(summary(mod2m)$r.squared, 2)}
  if(grau=="1"){s1=s <- sprintf("y == %e %s %e*x ~~~~~ italic(R^2) == %0.2f",
                                coef(moda)[1],
                                ifelse(coef(moda)[2] >= 0, "+", "-"),
                                abs(coef(moda)[2]),
                                r2)}
  if(grau=="2"){s2=s <- sprintf("y == %e %s %e * x %s %e * x^2 ~~~~~ italic(R^2) ==  %0.2f",
                                coef(mod1a)[1],
                                ifelse(coef(mod1a)[2] >= 0, "+", "-"),
                                abs(coef(mod1a)[2]),
                                ifelse(coef(mod1a)[3] >= 0, "+", "-"),
                                abs(coef(mod1a)[3]),
                                r2)}
  if(grau=="3"){s3=s <- sprintf("y == %e %s %e * x %s %e * x^2 %s %0.e * x^3 ~~~~~ italic(R^2) == %0.2f",
                                coef(mod2a)[1],
                                ifelse(coef(mod2a)[2] >= 0, "+", "-"),
                                abs(coef(mod2a)[2]),
                                ifelse(coef(mod2a)[3] >= 0, "+", "-"),
                                abs(coef(mod2a)[3]),
                                ifelse(coef(mod2a)[4] >= 0, "+", "-"),
                                abs(coef(mod2a)[4]),
                                r2)}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=as.numeric(as.character(names(media))),
                   resp=media,
                   desvio)
  grafico=ggplot(data1,aes(x=trat,y=resp))+
    geom_errorbar(aes(ymin=resp-desvio,
                      ymax=resp+desvio),width=0.15)+
    geom_point(aes(fill=as.factor(rep(1,length(resp)))),na.rm=TRUE,
               size=pointsize,shape=21,
               color="black")+
    theme+ylab(ylab)+xlab(xlab)
  if(is.na(ylim)==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}

  if(grau=="1"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x,size=1,color="black")}
  if(grau=="2"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x+I(x^2),size=1,color="black")}
  if(grau=="3"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x+I(x^2)+I(x^3),size=1,color="black")}
  if(grau=="1"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s1)),name="")}
  if(grau=="2"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s2)),name="")}
  if(grau=="3"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s3)),name="")}

  if(color=="gray"){if(grau=="1"){grafico=grafico+
    scale_fill_manual(values="black",label=c(parse(text=s1)),name="")}
    if(grau=="2"){grafico=grafico+
      scale_fill_manual(values="black",label=c(parse(text=s2)),name="")}
    if(grau=="3"){grafico=grafico+
      scale_fill_manual(values="black",label=c(parse(text=s3)),name="")}
  }

  grafico=grafico+
    theme(text = element_text(size=textsize,color="black",family=family),
          axis.text = element_text(size=textsize,color="black",family=family),
          axis.title = element_text(size=textsize,color="black",family=family),
          legend.position = posi,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)

  print(grafico)
  if(grau==1){
    print(mods)
    cat("\n----------------------------------------------------\n")
    cat("Deviations from regression")
    cat("\n----------------------------------------------------\n")
    print(fparcial1)
  }
  if(grau==2){
    print(mod1s)
    cat("\n----------------------------------------------------\n")
    cat("Deviations from regression")
    cat("\n----------------------------------------------------\n")
    print(fparcial2)
  }
  if(grau==3){
    print(mod2s)
    cat("\n----------------------------------------------------\n")
    cat("Deviations from regression")
    cat("\n----------------------------------------------------\n")
    print(fparcial3)
  }
  graficos=list(grafico)
}
