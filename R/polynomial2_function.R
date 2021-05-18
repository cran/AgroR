#' Analysis: Linear regression graph in double factorial
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Linear regression analysis for significant interaction of an experiment with two factors, one quantitative and one qualitative
#' @param fator1 Numeric or complex vector with factor 1 levels
#' @param resp Numerical vector containing the response of the experiment.
#' @param fator2 Numeric or complex vector with factor 2 levels
#' @param color Graph color (\emph{default} is NA)
#' @param grau Degree of the polynomial (1,2 or 3)
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param se Adds confidence interval (\emph{default} is FALSE)
#' @param legend.title Title legend
#' @param textsize Font size (\emph{default} is 12)
#' @param family Font family (\emph{default} is sans)
#' @param point Defines whether to plot all points ("all"), mean ("mean") or mean with standard error (\emph{default} - "mean_se").
#' @param ylim y-axis scale
#' @param posi Legend position
#' @keywords regression
#' @keywords Experimental
#' @seealso \link{polynomial}
#' @return Returns two or more linear, quadratic or cubic regression analyzes.
#' @export
#' @examples
#' dose=rep(c(0,0,0,2,2,2,4,4,4,6,6,6),3)
#' resp=c(8,7,5,23,24,25,30,34,36,80,90,80,
#' 12,14,15,23,24,25,50,54,56,80,90,40,
#' 12,14,15,3,4,5,50,54,56,80,90,40)
#' trat=rep(c("A","B","C"),e=12)
#' polynomial2(dose, resp, trat, grau=c(1,2,3))

polynomial2=function(fator1,
                     resp,
                     fator2,
                     color=NA,
                     grau=NA,
                     ylab="Response",
                     xlab="Independent",
                     theme=theme_bw(),
                     se=FALSE,
                     point="mean_se",
                     legend.title="Treatments",
                     posi="top",
                     textsize=12,
                     ylim=NA,
                     family="sans"){
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("gridExtra")
  Fator2=as.factor(fator2)
  if(is.na(color)[1]==TRUE){color=1:length(levels(Fator2))}
  if(is.na(grau)[1]==TRUE){grau=rep(1,length(levels(Fator2)))}
  curvas=c()
  texto=c()
  #desvios=c()
  data=data.frame(fator1,fator2,resp)
  grafico=ggplot(data,aes(y=resp,x=fator1))+
    theme+
    ylab(ylab)+
    xlab(xlab)

  if(point=="mean_sd"){grafico=grafico+
    stat_summary(aes(shape=fator2),fun="mean",  geom="point", size=3, na.rm=TRUE)+
    stat_summary(aes(group=fator2),fun = mean,
                 geom = "errorbar",na.rm=TRUE,
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),width=0.3)}
  if(point=="mean_se"){grafico=grafico+
    stat_summary(aes(shape=fator2,group=fator2),fun="mean",  geom="point",na.rm=TRUE, size=3)+
    stat_summary(aes(group=fator2),fun.data=mean_se, geom="errorbar",na.rm=TRUE,width=0.3)}
  if(point=="mean"){grafico=grafico+
    stat_summary(aes(shape=fator2,group=fator2),fun="mean",  geom="point",na.rm=TRUE, size=3)}
  if(point=="all"){grafico=grafico+
    geom_point(aes(shape=fator2),size=3)}

  if(is.na(ylim)==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}
  for(i in 1:length(levels(Fator2))){
    y=resp[Fator2==levels(Fator2)[i]]
    x=fator1[Fator2==levels(Fator2)[i]]
    f2=fator2[Fator2==levels(Fator2)[i]]
    d1=data.frame(y,x,f2)
    adj=grau[i]
    numero=order(levels(Fator2))[levels(Fator2)==levels(Fator2)[i]]
    if(adj==0){mod="ns"}
    if(adj==1){mod=lm(y~x)}
    if(adj==2){mod=lm(y~x+I(x^2))}
    if(adj==3){mod=lm(y~x+I(x^2)+I(x^3))}
    if(adj==1 | adj==2 | adj==3){
      ajuste=aov(y~as.factor(x))
      curvas[[i]]=summary(mod)$coefficients
      names(curvas)[i]=levels(Fator2)[i]
      }
    fats=as.character(unique(Fator2)[i])
    if(adj==1){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2),lty=unique(fator2)),
                                           method="lm", formula = y~x,na.rm=TRUE, se=se,size=0.6, show.legend=FALSE)}
    if(adj==2){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2),lty=unique(fator2)),
                                           method="lm", formula = y~x+I(x^2),na.rm=TRUE, se=se,size=0.6, show.legend=FALSE)}
    if(adj==3){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2),lty=unique(fator2)),
                                           method="lm", formula = y~x+I(x^2)+I(x^3),na.rm=TRUE, se=se,size=0.6, show.legend=FALSE)}
    m1=tapply(y,x,mean, na.rm=TRUE); x1=tapply(x,x,mean, na.rm=TRUE)
    if(adj==0){r2=0}
    if(adj==1){mod1=lm(m1~x1)
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==2){mod1=lm(m1~x1+I(x1^2))
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==3){mod1=lm(m1~x1+I(x1^2)+I(x1^3))
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==0){text=sprintf("ns")}
    if(adj==1){text=sprintf("y == %0.3e %s %0.3e*x ~~~~~ italic(R^2) == %0.2f",coef(mod)[1],
                            ifelse(coef(mod)[2] >= 0, "+", "-"),abs(coef(mod)[2]),
                            r2)}
    if(adj==2){text=sprintf("y == %0.3e %s %0.3e * x %s %0.3e * x^2 ~~~~~ italic(R^2) ==  %0.2f",
                            coef(mod)[1],
                            ifelse(coef(mod)[2] >= 0, "+", "-"),
                            abs(coef(mod)[2]),
                            ifelse(coef(mod)[3] >= 0, "+", "-"),
                            abs(coef(mod)[3]),
                            r2)}
    if(adj==3){text=sprintf("y == %0.3e %s %0.3e * x %s %0.3e * x^2 %s %0.3e * x^3 ~~~~~~ italic(R^2) == %0.2f",
                            coef(mod)[1],
                            ifelse(coef(mod)[2] >= 0, "+", "-"),
                            abs(coef(mod)[2]),
                            ifelse(coef(mod)[3] >= 0, "+", "-"),
                            abs(coef(mod)[3]),
                            ifelse(coef(mod)[4] >= 0, "+", "-"),
                            abs(coef(mod)[4]),
                            r2)}
    texto[[i]]=text
  }
  cat("\n----------------------------------------------------\n")
  cat("Regression Models")
  cat("\n----------------------------------------------------\n")
  print(curvas)
  grafico=grafico+
    scale_linetype_manual(name=legend.title,values = color,drop=FALSE,
                                     label=parse(text=paste("(\"",levels(Fator2),"\")~",unlist(texto))))+
    scale_shape_discrete(label=parse(text=paste("(\"",levels(Fator2),"\")~",unlist(texto))))+
    theme(text = element_text(size=textsize,color="black", family = family),
          axis.text = element_text(size=textsize,color="black", family = family),
          axis.title = element_text(size=textsize,color="black", family = family),
          legend.position = posi,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    labs(color=legend.title, shape=legend.title, lty=legend.title)+
    scale_color_grey(name=legend.title,
                     start = 0.12, end = 0.1,
                     label=parse(text=paste("(\"",levels(Fator2),"\")~",unlist(texto))))
  print(grafico)

  grafico=as.list(grafico)
}

polynomial2_color=function(fator1,
                           resp,
                           fator2,
                           color=NA,
                           grau=NA,
                           ylab="Response",
                           xlab="independent",
                           theme=theme_bw(),
                           se=FALSE,
                           point="mean_se",
                           legend.title="Tratamentos",
                           posi="top",
                           textsize=12,
                           ylim=NA,
                           family="sans"){
  requireNamespace("ggplot2")
  requireNamespace("gridExtra")
  Fator2=as.factor(fator2)
  if(is.na(color)[1]==TRUE){color=1:length(levels(Fator2))}
  if(is.na(grau)[1]==TRUE){grau=rep(1,length(levels(Fator2)))}
  curvas=c()
  texto=c()
  # desvios=c()
  data=data.frame(fator1,fator2,resp)
  grafico=ggplot(data,aes(y=resp,x=fator1))+
    theme+
    ylab(ylab)+
    xlab(xlab)

  if(point=="mean_sd"){grafico=grafico+
    stat_summary(aes(color=fator2),fun = "mean",  geom="point", size=3)+
    stat_summary(aes(color=fator2),fun = mean,
                 geom = "errorbar",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),width=0.2)}
  if(point=="mean_se"){grafico=grafico+
    stat_summary(aes(color=fator2),fun="mean",  geom="point", size=3)+
    stat_summary(aes(color=fator2),fun.data=mean_se, geom="errorbar",width=0.2)}
  if(point=="mean"){grafico=grafico+
    stat_summary(aes(color=fator2),fun="mean",  geom="point")}
  if(point=="all"){grafico=grafico+
    geom_point(aes(color=fator2))}

  if(is.na(ylim)==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}
  for(i in 1:length(levels(Fator2))){
    y=resp[Fator2==levels(Fator2)[i]]
    x=fator1[Fator2==levels(Fator2)[i]]
    f2=fator2[Fator2==levels(Fator2)[i]]
    d1=data.frame(y,x,f2)
    adj=grau[i]
    numero=order(levels(Fator2))[levels(Fator2)==levels(Fator2)[i]]
    if(adj==0){mod="ns"}
    if(adj==1){mod=lm(y~x)}
    if(adj==2){mod=lm(y~x+I(x^2))}
    if(adj==3){mod=lm(y~x+I(x^2)+I(x^3))}
    if(adj==1 | adj==2 | adj==3){
      ajuste=aov(y~as.factor(x))
      curvas[[i]]=summary(mod)$coefficients
      names(curvas)[i]=levels(Fator2)[i]
      }
    fats=as.character(unique(Fator2)[i])
    if(adj==1){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2)),
                                           method="lm", formula = y~x,size=0.8, se=se)}
    if(adj==2){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2)),
                                           method="lm", formula = y~x+I(x^2), size=0.8,se=se)}
    if(adj==3){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2)),
                                           method="lm", formula = y~x+I(x^2)+I(x^3), size=0.8,se=se)}
    m1=tapply(y,x,mean, na.rm=TRUE); x1=tapply(x,x,mean, na.rm=TRUE)

    if(adj==0){r2=0}
    if(adj==1){mod1=lm(m1~x1)
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==2){mod1=lm(m1~x1+I(x1^2))
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==3){mod1=lm(m1~x1+I(x1^2)+I(x1^3))
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==0){text=sprintf("ns")}
    if(adj==1){text=sprintf("y == %e %s %e*x~~~~~ italic(R^2) == %0.2f",coef(mod)[1],
                            ifelse(coef(mod)[2] >= 0, "+", "-"),abs(coef(mod)[2]),
                            r2)}
    if(adj==2){text=sprintf("y == %e %s %e * x %s %e * x^2~~~~~ italic(R^2) ==  %0.2f",
                            coef(mod)[1],
                            ifelse(coef(mod)[2] >= 0, "+", "-"),
                            abs(coef(mod)[2]),
                            ifelse(coef(mod)[3] >= 0, "+", "-"),
                            abs(coef(mod)[3]),
                            r2)}
    if(adj==3){text=sprintf("y == %e %s %e * x %s %e * x^2 %s %0.9e * x^3~~~~~~ italic(R^2) == %0.2f",
                            coef(mod)[1],
                            ifelse(coef(mod)[2] >= 0, "+", "-"),
                            abs(coef(mod)[2]),
                            ifelse(coef(mod)[3] >= 0, "+", "-"),
                            abs(coef(mod)[3]),
                            ifelse(coef(mod)[4] >= 0, "+", "-"),
                            abs(coef(mod)[4]),
                            r2)}
    texto[[i]]=text
  }
  cat("\n----------------------------------------------------\n")
  cat("Regression Models")
  cat("\n----------------------------------------------------\n")
  print(curvas)
  grafico=grafico+scale_colour_discrete(label=parse(text=paste("(\"",levels(Fator2),"\")~",
                                                               unlist(texto),sep="")))+
    theme(text = element_text(size=textsize,color="black", family = family),
          axis.text = element_text(size=textsize,color="black", family = family),
          axis.title = element_text(size=textsize,color="black", family = family),
          legend.position = posi,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    labs(color=legend.title, shape=legend.title, lty=legend.title)
  print(grafico)

  grafico=as.list(grafico)
}
