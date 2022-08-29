#' Analysis: t test to compare means with a reference value
#'
#' @description Sometimes the researcher wants to test whether the treatment mean is greater than/equal to or less than a reference value. For example, I want to know if the average productivity of my treatment is higher than the average productivity of a given country. For this, this function allows comparing the means with a reference value using the t test.
#' @author Gabriel Danilo Shimizu
#' @param response Numerical vector containing the response of the experiment.
#' @param trat Numerical or complex vector with treatments
#' @param mu A number indicating the true value of the mean
#' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param conf.level confidence level of the interval.
#' @importFrom gridExtra grid.arrange
#' @export
#' @return returns a list with the mean per treatment, maximum, minimum, sample standard deviation, confidence interval, t-test statistic and its p-value.
#' @note No treatment can have zero variability. Otherwise the function will result in an error.
#' @examples
#' library(AgroR)
#' data("pomegranate")
#' tonetest(resp=pomegranate$WL,
#' trat=pomegranate$trat,
#' mu=2,
#' alternative = "greater")

tonetest=function(response,
                    trat,
                    mu=0,
                    alternative="two.sided",
                    conf.level=0.95){
  trat=factor(trat,unique(trat))
  nl=nlevels(trat)
  teste=c(1:nl)
  inferior=c(1:nl)
  superior=c(1:nl)
  estatistica=c(1:nl)
  minimo=c(1:nl)
  maximo=c(1:nl)
  media=c(1:nl)
  std=c(1:nl)
  for(i in 1:nl){
    test=t.test(response[trat==levels(trat)[i]],
                    mu=mu,
                    alternative = alternative,conf.level=conf.level)
    media[i]=mean(response[trat==levels(trat)[i]])
    std[i]=sd(response[trat==levels(trat)[i]])
    minimo[i]=min(response[trat==levels(trat)[i]])
    maximo[i]=max(response[trat==levels(trat)[i]])
    teste[i]=test$p.value
    inferior[i]=test$conf.int[1]
    superior[i]=test$conf.int[2]
    estatistica[i]=test$statistic}
  names(teste)=levels(trat)
  names(inferior)=levels(trat)
  names(superior)=levels(trat)
  names(estatistica)=levels(trat)
  names(media)=levels(trat)
  names(std)=levels(trat)
  names(minimo)=levels(trat)
  names(maximo)=levels(trat)
  list("mean.ref"=mu,
       "mean"=media,
       "sd"=std,
       "min"=minimo,
       "max"=maximo,
       "Statistic"=estatistica,
       "p.value"=teste,
       "conf.inf"=inferior,
       "conf.sup"=superior)}

