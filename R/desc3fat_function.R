#' Descriptive: Descriptive analysis (Three factors)
#'
#' @description Performs the descriptive graphical analysis of an experiment with three factors of interest.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param f3 Numeric or complex vector with factor 3 levels
#' @param response Numerical vector containing the response of the experiment.
#' @param legend.title Legend title
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab x name (Accepts the \emph{expression}() function)
#' @param theme ggplot theme
#' @param plot "interaction" or "box"
#' @keywords Descriptive
#' @keywords Experimental
#' @return The function returns a triple interaction graph.
#' @export
#' @examples
#' library(AgroR)
#' data(enxofre)
#' with(enxofre, desc3fat(f1, f2, f3, resp))

######################################################################################
## Analise descritiva
######################################################################################

desc3fat=function(f1,
                  f2,
                  f3,
                  response,
                  legend.title="Legend",
                  xlab="",
                  ylab="Response",
                  theme=theme_classic(),
                  plot="interaction"){
  f1=as.factor(f1)
  f2=as.factor(f2)
  f3=as.factor(f3)
  requireNamespace("ggplot2")

  #===========================
  # Fator 1 x Fator 2
  #===========================
  dados=data.frame(f1,f2,f3, response)

  #===========================
  # Geral
  #===========================

  Media = mean(response, na.rm=TRUE)
  Mediana = median(response, na.rm=TRUE)
  Minimo = min(response, na.rm=TRUE)
  Maximo = max(response, na.rm=TRUE)
  Variancia = var(response, na.rm=TRUE)
  Desvio = sd(response, na.rm=TRUE)
  CV = Desvio / Media * 100
  juntos=cbind(Media,
               Mediana,
               Minimo,
               Maximo,
               Variancia,
               Desvio,
               CV)
  colnames(juntos)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")
  rownames(juntos)=""
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(italic("General description")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos)
  #===========================
  # Fator 1
  #===========================
  Media = tapply(response, f1, mean, na.rm=TRUE)
  Mediana = tapply(response, f1, median, na.rm=TRUE)
  Minimo = tapply(response, f1, min, na.rm=TRUE)
  Maximo = tapply(response, f1, max, na.rm=TRUE)
  Variancia = tapply(response, f1, var, na.rm=TRUE)
  Desvio = tapply(response, f1, sd, na.rm=TRUE)
  CV = Desvio / Media * 100
  juntos2 = cbind(Media,
                  Mediana,
                  Minimo,
                  Maximo,
                  Variancia,
                  Desvio,
                  CV)
  colnames(juntos2)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(italic("F1")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos2)

  #===========================
  # Fator 2
  #===========================
  Media = tapply(response, f2, mean, na.rm=TRUE)
  Mediana = tapply(response, f2, median, na.rm=TRUE)
  Minimo = tapply(response, f2, min, na.rm=TRUE)
  Maximo = tapply(response, f2, max, na.rm=TRUE)
  Variancia = tapply(response, f2, var, na.rm=TRUE)
  Desvio = tapply(response, f2, sd, na.rm=TRUE)
  CV = Desvio / Media * 100
  juntos3=cbind(Media,
                Mediana,
                Minimo,
                Maximo,
                Variancia,
                Desvio,
                CV)
  colnames(juntos3)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(italic("F2")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos3)

  #===========================
  # Fator 3
  #===========================
  Media = tapply(response, f3, mean, na.rm=TRUE)
  Mediana = tapply(response, f3, median, na.rm=TRUE)
  Minimo = tapply(response, f3, min, na.rm=TRUE)
  Maximo = tapply(response, f3, max, na.rm=TRUE)
  Variancia = tapply(response, f3, var, na.rm=TRUE)
  Desvio = tapply(response, f3, sd, na.rm=TRUE)
  CV = Desvio / Media * 100
  juntos3=cbind(Media,
                Mediana,
                Minimo,
                Maximo,
                Variancia,
                Desvio,
                CV)
  colnames(juntos3)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(italic("F3")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos3)

  #===========================
  # inter
  #===========================
  Media = tapply(response, paste(f1,f2,f3), mean, na.rm=TRUE)
  Mediana = tapply(response, paste(f1,f2,f3), median, na.rm=TRUE)
  Minimo = tapply(response, paste(f1,f2,f3), min, na.rm=TRUE)
  Maximo = tapply(response, paste(f1,f2,f3), max, na.rm=TRUE)
  Variancia = tapply(response, paste(f1,f2,f3), var, na.rm=TRUE)
  Desvio = tapply(response, paste(f1,f2,f3), sd, na.rm=TRUE)
  CV = Desvio / Media * 100
  juntos4=cbind(Media,
                Mediana,
                Minimo,
                Maximo,
                Variancia,
                Desvio,
                CV)
  colnames(juntos4)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(italic("Interaction")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos4)

  if(plot=="box"){
  grafico=ggplot(dados,aes(x=f1,y=response, fill=f2))+
    stat_boxplot(geom='errorbar', linetype=1,
                 position = position_dodge(width = 0.75),width=0.5)+
    geom_boxplot()+xlab(xlab)+labs(fill=legend.title)+
    ylab(ylab)+theme+facet_wrap(~f3)
  grafico=grafico+
    theme(text = element_text(size=12,color="black"),
          axis.title = element_text(size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          strip.text = element_text(size=13))
  }
  #===========================
  # Interacao
  #===========================
  if(plot=="interaction"){
  grafico=ggplot(dados,aes(x=f1,y=response, color=f2))+
    stat_summary(fun.data = "mean_se")+
    stat_summary(aes(color=f2, group=f2),
                 geom="line", fun.data = "mean_se")+
    ylab(ylab)+xlab(xlab)+labs(fill=legend.title)+theme+facet_wrap(~f3)
  grafico=grafico+
    theme(text = element_text(size=12,color="black"),
          axis.title = element_text(size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          strip.text = element_text(size=13))}
  print(grafico)
  grafico=as.list(grafico)
}

