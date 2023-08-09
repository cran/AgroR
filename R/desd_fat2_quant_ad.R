#' Analysis: Regression analysis by orthogonal polynomials for double factorial scheme with additional control
#'
#' @description Regression analysis by orthogonal polynomials for double factorial scheme with additional control. Cases in which the additional belongs to the regression curve, being common to the qualitative levels. In these cases, the additional (usually dose 0/control treatment) is not part of the factor arrangement. One option addressed by this function is to analyze a priori as a double factorial scheme with an additional one and correct the information a posteriore using information from the initial analysis, such as the degree of freedom and the sum of squares of the residue.
#'
#' @param output Output from a FAT2DIC.ad or FAT2DBC.ad function
#' @param ad.value Additional treatment quantitative factor level
#' @param design Type of experimental project (FAT2DIC.ad or FAT2DBC.ad)
#' @param grau  Degree of the polynomial (only for the isolated effect of the quantitative factor)
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @export
#' @examples
#'
#' #==================================================
#' # Data set
#' trat=rep(c("A","B","C"),e=12)
#' dose=rep(rep(c(200,400,600,800),e=3),3)
#' d0=c(40,45,48)
#' respo=c(60,55,56, 60,65,66, 70,75,76,
#'         80,85,86, 50,55,56, 70,75,76,
#'         60,65,66, 50,45,46, 50,45,46,
#'         50,55,66, 70,75,76, 80,85,86)
#' repe=rep(c("R1","R2","R3"),12)
#' #==================================================
#' # Analysis FAT2DIC.ad
#' resu=FAT2DIC.ad(trat,dose,repe = repe,respo,responseAd = d0,quali = c(TRUE,FALSE),grau21 = c(1,2,1))
#'
#' #==================================================
#' # Regression analysis
#' desd_fat2_quant_ad(resu,ad.value=0,design="FAT2DIC.ad")
#'
#'
#' # Data set
#' trat=rep(c("A","B"),e=12)
#' dose=rep(rep(c(200,400,600,800),e=3),2)
#' d0=c(40,45,48)
#' respo=c(60,55,56,60,65,66,70,75,76,80,85,86,50,45,46,50,55,66,70,75,76,80,85,86)
#' repe=rep(c("R1","R2","R3"),8)
#' #==================================================
#' # Analysis FAT2DIC.ad
#' resu=FAT2DIC.ad(trat,dose,repe = repe,respo,responseAd = d0,quali = c(TRUE,FALSE))
#' #==================================================
#' # Regression analysis
#' desd_fat2_quant_ad(resu,ad.value=0,design="FAT2DIC.ad",grau=1)


desd_fat2_quant_ad=function(output,
                            ad.value=0,
                            design="FAT2DIC.ad",
                            grau=1){
  alpha.f=output[[1]]$plot$alpha.f
  alpha.t=output[[1]]$plot$alpha.t
  grau21=output[[1]]$plot$grau21
  ylab=output[[1]]$plot$ylab
  xlab=parse(text=output[[1]]$plot$xlab.factor[2])
  posi=output[[1]]$plot$posi
  theme=output[[1]]$plot$theme
  textsize=output[[1]]$plot$textsize
  point=output[[1]]$plot$point
  family=output[[1]]$plot$family

  dados=output[[1]]$plot$ordempadronizado
  ana=output[[1]]$plot$anava
  nni=length(unique(dados$f1))
  respad=output[[1]]$plot$respAd
  dose0=rep(respad,nni)
  trat0=rep(unique(dados$f1),e=length(respad))
  doses0=rep(ad.value,length(dose0))
  f1=c(dados$f1,trat0)
  f2=c(dados$f2,doses0)
  resp=c(dados$resp,dose0)
  if(design=="FAT2DIC.ad"){
    if(ana$`Pr(>F)`[3]>alpha.f & ana$`Pr(>F)`[2]<alpha.f){print(polynomial(f2,
                                                                     resp,
                                                                     grau = grau,
                                                                     ylab = ylab,
                                                                     xlab=xlab,
                                                                     posi=posi,
                                                                     theme=theme,
                                                                     textsize=textsize,
                                                                     point=point,
                                                                     family=family,
                                                                     SSq = ana$`Sum Sq`[5],
                                                                     DFres = ana$Df[5])[[1]])}
    if(ana$`Pr(>F)`[3]<alpha.f){saida=polynomial2_color(f2,resp,f1,grau = grau21,SSq = ana$`Sum Sq`[5],DFres = ana$Df[5])}}
  if(design=="FAT2DBC.ad"){
    if(ana$`Pr(>F)`[4]>alpha.f & ana$`Pr(>F)`[2]<alpha.f){print(polynomial(f2,
                                                                     resp,
                                                                     grau = grau,
                                                                     ylab = ylab,
                                                                     xlab=xlab,
                                                                     posi=posi,
                                                                     theme=theme,
                                                                     textsize=textsize,
                                                                     point=point,
                                                                     family=family,
                                                                     SSq = ana$`Sum Sq`[6],
                                                                     DFres = ana$Df[6])[[1]])}
    if(ana$`Pr(>F)`[4]<alpha.f){saida=polynomial2_color(f2,resp,f1,grau = grau21,SSq = ana$`Sum Sq`[6],DFres = ana$Df[6])}}}







