#' Utils: Data transformation (Box-Cox, 1964)
#'
#' @description Estimates the lambda value for data transformation
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param response Numerical vector containing the response of the experiment.
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param f3 Numeric or complex vector with factor 3 levels
#' @param block Numerical or complex vector with blocks
#' @param line Numerical or complex vector with lines
#' @param column Numerical or complex vector with columns
#' @keywords Transformation
#' @keywords Experimental
#' @export
#' @return Returns the value of lambda and/or data transformation approximation, according to Box-Cox (1964)
#' @references
#'
#' Box, G. E., Cox, D. R. (1964). An analysis of transformations. Journal of the Royal Statistical Society: Series B (Methodological), 26(2), 211-243.
#' @examples
#'
#' #================================================================
#' # Completely randomized design
#' #================================================================
#' data("pomegranate")
#' with(pomegranate, transf(WL,f1=trat))
#'
#' #================================================================
#' # Randomized block design
#' #================================================================
#' data(soybean)
#' with(soybean, transf(prod, f1=cult, block=bloc))
#'
#' #================================================================
#' # Completely randomized design in double factorial
#' #================================================================
#' data(cloro)
#' with(cloro, transf(resp, f1=f1, f2=f2))
#'
#' #================================================================
#' # Randomized block design in double factorial
#' #================================================================
#' data(cloro)
#' with(cloro, transf(resp, f1=f1, f2=f2, block=bloco))


transf=function(response,
                f1,
                f2=NA,
                f3=NA,
                block=NA,
                line=NA,
                column=NA){
  # DIC simples
  requireNamespace("MASS")

  if(is.na(f2[1])==TRUE && is.na(f3[1])==TRUE && is.na(block[1])==TRUE &&
     is.na(line[1])==TRUE && is.na(column[1])==TRUE){
    vero=MASS::boxcox(response~f1)}

  # DBC simples
  if(is.na(f2[1])==TRUE && is.na(f3[1])==TRUE && is.na(block[1])==FALSE &&
     is.na(line[1])==TRUE && is.na(column[1])==TRUE){
    vero=MASS::boxcox(response~f1+block)}

  # DQL
  if(is.na(f2[1])==TRUE && is.na(f3[1])==TRUE && is.na(block[1])==TRUE &&
     is.na(line[1])==FALSE && is.na(column[1])==FALSE){
    vero=MASS::boxcox(response~f1+column+line)}

  # fat2.dic
  if(is.na(f2[1])==FALSE && is.na(f3[1])==TRUE && is.na(block[1])==TRUE &&
     is.na(line[1])==TRUE && is.na(column[1])==TRUE){
    vero=MASS::boxcox(response~f1*f2)}

  #fat2dbc
  if(is.na(f2[1])==FALSE && is.na(f3[1])==TRUE && is.na(block[1])==FALSE &&
     is.na(line[1])==TRUE && is.na(column[1])==TRUE){
    vero=MASS::boxcox(response~f1*f2+block)}

  # fat3dic
  if(is.na(f2[1])==FALSE && is.na(f3[1])==FALSE && is.na(block[1])==TRUE &&
     is.na(line[1])==TRUE && is.na(column[1])==TRUE){
    vero=MASS::boxcox(response~f1*f2*f3)}

  # fat3dic
  if(is.na(f2[1])==FALSE && is.na(f3[1])==FALSE && is.na(block[1])==FALSE &&
     is.na(line[1])==TRUE && is.na(column[1])==TRUE){
    vero=MASS::boxcox(response~f1*f2*f3+block)}

  maxvero = vero$x[which.max(vero$y)]
  cat("\n------------------------------------------------\n")
  cat("Box-Cox Transformation (1964)\n")
  cat("\nLambda=")
  cat(lambda=maxvero)
  cat("\n------------------------------------------------\n")
  cat("Suggestion:\n")
  cat(if(round(maxvero,2)>0.75 & round(maxvero,2)<1.25){"Do not transform data"}
      else{if(round(maxvero,2)>0.25 & round(maxvero,2)<=0.75){"square root (sqrt(Y))"}
        else{if(round(maxvero,2)>-0.75 & round(maxvero,2)<=-0.25){"1/sqrt(Y)"}
          else{if(round(maxvero,2)>-1.25 & round(maxvero,2)<=-0.75){"1/Y"}
            else{if(round(maxvero,2)>=-0.25 & round(maxvero,2)<0.25){"log(Y)"}}}}}
  )
  cat("\n\n")
  cat("ou:")
  cat("Yt=(Y^lambda-1)/lambda")
}
