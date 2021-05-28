#' Analysis: Randomized block design evaluated over time
#'
#' @description Function of the AgroR package for analysis of experiments conducted in a balanced qualitative, single-factorial randomized block design with multiple assessments over time, however without considering time as a factor.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Gonçalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param block Numerical or complex vector with blocks
#' @param time Numerical or complex vector with times
#' @param response Numerical vector containing the response of the experiment.
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD ("lsd"), Scott-Knott ("sk"), Duncan ("duncan") and Friedman ("fd"))
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param error Add error bar
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param textsize Font size of the texts and titles of the axes
#' @param labelsize Font size of the labels
#' @param family Font family
#' @param dec Number of cells
#' @param geom Graph type (columns - "bar" or segments "point")
#' @param legend Legend title
#' @param posi Legend position
#' @param ylim y-axis scale
#' @param xnumeric Declare x as numeric (\emph{default} is FALSE)
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @keywords dbct
#' @keywords Experimental
#' @seealso \link{DBC}, \link{DICT}, \link{DQLT}
#' @references
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
#' @details The p-value of the analysis of variance, the normality test for Shapiro-Wilk errors, the Bartlett homogeneity test of variances, the independence of Durbin-Watson errors and the multiple comparison test (Tukey, Scott-Knott, LSD or Duncan).
#' @export
#' @return The function returns the p-value of Anova, the assumptions of normality of errors, homogeneity of variances and independence of errors, multiple comparison test, as well as a line graph
#' @examples
#' data(simulate2)
#' attach(simulate2)
#'
#' #===================================
#' # default
#' #===================================
#' DBCT(trat, bloco, tempo, resp)
#'
#' #===================================
#' # segment chart
#' #===================================
#' DBCT(trat, bloco, tempo, resp, geom="point")

DBCT=function(trat,
              block,
              time,
              response,
              alpha.f=0.05,
              alpha.t=0.05,
              theme=theme_classic(),
              geom="bar",
              fill="gray",
              ylab="Response",
              xlab="Independent",
              mcomp="tukey",
              textsize=12,
              labelsize=5,
              error=TRUE,
              family="sans",
              sup=0,
              addmean=FALSE,
              posi=c(0.1,0.8),
              legend="Legend",
              ylim=NA,
              dec=3,
              xnumeric=FALSE){
  requireNamespace("ScottKnott")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")
  requireNamespace("lmtest")
  requireNamespace("ggrepel")
  trat=as.factor(trat)
  resp=response
  block=as.factor(block)
  time=as.factor(time)
  dados=data.frame(resp,trat,block,time)
  if(mcomp=="tukey"){
    tukeyg=c()
    ordem=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      tukey=HSD.test(mod,"trat",alpha = alpha.t)
      tukey$groups=tukey$groups[unique(as.character(trat)),2]
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){tukey$groups=c("ns",rep(" ",length(unique(trat))-1))}
      tukeyg[[i]]=as.character(tukey$groups)
      ordem[[i]]=rownames(tukey$groups)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(tukeyg)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}

  if(mcomp=="lsd"){
    lsdg=c()
    ordem=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      lsd=LSD.test(mod,"trat",alpha = alpha.t)
      lsd$groups=lsd$groups[unique(as.character(trat)),2]
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){lsd$groups=c("ns",rep(" ",length(unique(trat))-1))}
      lsdg[[i]]=as.character(lsd$groups)
      ordem[[i]]=rownames(lsd$groups)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(lsdg)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}

  if(mcomp=="duncan"){
    duncang=c()
    ordem=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      duncan=duncan.test(mod,"trat",alpha = alpha.t)
      duncan$groups=duncan$groups[unique(as.character(trat)),2]
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){duncan$groups=c("ns",rep(" ",length(unique(trat))-1))}
      duncang[[i]]=as.character(duncan$groups)
      ordem[[i]]=rownames(duncan$groups)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(duncang)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}


  if(mcomp=="sk"){
    scott=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      letra=SK(mod,"trat",sig.level = alpha.t)
      data=data.frame(sk=letters[letra$groups])
      rownames(data)=rownames(letra$m.inf)
      data=data[unique(as.character(trat)),]
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){data=c("ns",rep(" ",length(unique(trat))-1))}
      scott[[i]]=as.character(data)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(scott)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}

  if(mcomp=="fd"){
    fdg=c()
    ordem=c()
    anovag=c()
    for(i in 1:length(levels(time))){
      fd=friedman(block,trat,resp,alpha=alpha.t)
      anovag[[i]]=mod$statistics$p.chisq
      fd$groups=fd$groups[unique(as.character(trat)),2]
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){fd$groups=c("ns",rep(" ",length(unique(trat))-1))}
      fdg[[i]]=as.character(duncan$groups)
      ordem[[i]]=rownames(fd$groups)
    }
    m=unlist(fdg)
    an=unlist(anovag)
    press=data.frame(fd);colnames(press)=c("p-value Chisq Friedman")}


  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("ANOVA and assumptions")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(press)

  dadosm=data.frame(time=as.character(rep(unique(time),e=length(unique(as.character(trat))))),
                    trat=rep(unique(as.character(trat)),length(unique(time))),
                    media=c(tapply(resp,list(trat,time),mean, na.rm=TRUE)[unique(as.character(trat)),]),
                    desvio=c(tapply(resp,list(trat,time),sd, na.rm=TRUE)[unique(as.character(trat)),]),
                    letra=m)
  if(xnumeric==TRUE){dadosm$time=as.numeric(as.character(dadosm$time))}
  time=dadosm$time
  trat=dadosm$trat
  media=dadosm$media
  desvio=dadosm$desvio
  letra=dadosm$letra

  if(geom=="point"){
    grafico=ggplot(dadosm,aes(y=media,
                              x=time))+
      geom_point(aes(shape=factor(trat, levels=unique(as.character(trat))),
                     group=factor(trat, levels=unique(as.character(trat)))),size=3)+
      geom_line(aes(lty=factor(trat, levels=unique(as.character(trat))),
                    group=factor(trat, levels=unique(as.character(trat)))),size=0.8)+
      ylab(ylab)+
      xlab(xlab)+theme+
      theme(text = element_text(size=textsize,color="black", family = family),
            axis.title = element_text(size=textsize,color="black", family = family),
            axis.text = element_text(size=textsize,color="black", family = family),
            legend.position = posi,
            legend.text = element_text(size = textsize))+labs(shape=legend, lty=legend)
    if(error==TRUE){grafico=grafico+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio), width=0.1)}
    if(addmean==FALSE && error==FALSE){grafico=grafico+
      geom_text_repel(aes(y=media+sup,label=letra),family=family,size=labelsize)}
    if(addmean==TRUE && error==FALSE){grafico=grafico+
      geom_text_repel(aes(y=media+sup,
                          label=paste(format(media,digits = dec),
                                      letra)),family=family,size=labelsize)}
    if(addmean==FALSE && error==TRUE){grafico=grafico+
      geom_text_repel(aes(y=desvio+media+sup,
                          label=letra),family=family,size=labelsize)}
    if(addmean==TRUE && error==TRUE){grafico=grafico+
      geom_text_repel(aes(y=desvio+media+sup,
                          label=paste(format(media,digits = dec),
                                      letra)),family=family,size=labelsize)}
  }

  if(geom=="bar"){
    if(sup==0){sup=0.1*mean(dadosm$media)}
    grafico=ggplot(dadosm,aes(y=media,
                              x=as.factor(time),
                              fill=factor(trat,levels = unique(trat))))+
      geom_col(position = "dodge",color="black")+
      ylab(ylab)+
      xlab(xlab)+theme+
      theme(text = element_text(size=textsize,color="black", family = family),
            axis.title = element_text(size=textsize,color="black", family = family),
            axis.text = element_text(size=textsize,color="black", family = family),
            legend.position = posi,
            legend.text = element_text(size = textsize))+labs(fill=legend)
    if(error==TRUE){grafico=grafico+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),
                    width=0.3, position = position_dodge(width=0.9))}
    if(addmean==FALSE && error==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,label=letra),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
    if(addmean==TRUE && error==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,
                    label=paste(format(media,digits = dec),letra)),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
    if(addmean==FALSE && error==TRUE){grafico=grafico+
      geom_text(aes(y=desvio+media+sup,label=letra),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
    if(addmean==TRUE && error==TRUE){grafico=grafico+
      geom_text(aes(y=desvio+media+sup,
                    label=paste(format(media,digits = dec),letra)),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
  }
  if(fill=="gray"){grafico=grafico+scale_fill_grey(start = 1, end = 0.1)}
  if(is.na(ylim)==FALSE){grafico=grafico+scale_y_continuous(breaks = ylim)}
  graficos=as.list(grafico)
  print(grafico)
}
