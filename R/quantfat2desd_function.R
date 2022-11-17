#' Analysis: Polynomial splitting for double factorial in DIC and DBC
#' @description Splitting in polynomials for double factorial in DIC and DBC. Note that f1 must always be qualitative and f2 must always be quantitative. This function is an easier way to visualize trends for dual factor schemes with a quantitative and a qualitative factor.
#' @param factors  Define f1 and f2 and/or block factors in list form. Please note that in the list it is necessary to write `f1`, `f2` and `block`. See example.
#' @param response response variable
#' @param dec Number of cells
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @return Returns the coefficients of the linear, quadratic and cubic models, the p-values of the t test for each coefficient (p.value.test) and the p-values for the linear, quadratic, cubic model splits and the regression deviations.
#' @keywords Experimental
#' @seealso \link{FAT2DIC}, \link{FAT2DBC}
#' @export
#' @examples
#' library(AgroR)
#' data(cloro)
#' quant.fat2.desd(factors = list(f1=cloro$f1,
#' f2=rep(c(1:4),e=5,2), block=cloro$bloco),
#' response=cloro$resp)

quant.fat2.desd=function(factors=list(f1,f2,block),
                         response, dec=3){
  anlinear=function(trat,resp,sq,df){
    a=anova(lm(resp~trat))
    b=anova(lm(resp~as.factor(trat)))
    alin=rbind(a[1,],b[1,],"resid"=c(df,sq,sq/df,NA,NA))
    alin$`Sum Sq`[2]=alin$`Sum Sq`[2]-alin$`Sum Sq`[1]
    alin$Df[2]=alin$Df[2]-alin$Df[1]
    alin$`Mean Sq`=alin$`Sum Sq`/alin$Df
    alin$`F value`=c(alin$`Mean Sq`[1:2]/alin$`Mean Sq`[3],NA)
    alin$`Pr(>F)`=c(1-pf(alin$`F value`[1:2],alin$Df[1:2],alin$Df[3]),NA)
    rownames(alin)=c("Linear","Deviation","Residuals")
    alin}
  anquad=function(trat,resp,sq,df){
    a=anova(lm(resp~trat+I(trat^2)))
    b=anova(lm(resp~as.factor(trat)))
    alin=rbind(a[1:2,],b[1,],"resid"=c(df,sq,sq/df,NA,NA))
    alin$`Sum Sq`[3]=alin$`Sum Sq`[3]-sum(alin$`Sum Sq`[1:2])
    alin$Df[3]=alin$Df[3]-2
    alin$`Mean Sq`=alin$`Sum Sq`/alin$Df
    alin$`F value`=c(alin$`Mean Sq`[1:3]/alin$`Mean Sq`[4],NA)
    alin$`Pr(>F)`=c(1-pf(alin$`F value`[1:3],alin$Df[1:3],alin$Df[4]),NA)
    rownames(alin)=c("Linear","Quadratic","Deviation","Residuals")
    alin}
  ancubic=function(trat,resp,sq,df){
    a=anova(lm(resp~trat+I(trat^2)+I(trat^3)))
    b=anova(lm(resp~as.factor(trat)))
    alin=rbind(a[1:3,],b[1,],"resid"=c(df,sq,sq/df,NA,NA))
    alin$`Sum Sq`[4]=alin$`Sum Sq`[4]-sum(alin$`Sum Sq`[1:3])
    alin$Df[4]=alin$Df[4]-3
    alin$`Mean Sq`=alin$`Sum Sq`/alin$Df
    alin$`F value`=c(alin$`Mean Sq`[1:4]/alin$`Mean Sq`[5],NA)
    alin$`Pr(>F)`=c(1-pf(alin$`F value`[1:4],alin$Df[1:4],alin$Df[5]),NA)
    rownames(alin)=c("Linear","Quadratic","Cubic","Deviation","Residuals")
    alin}
  requireNamespace("crayon")
  if(length(factors)==2){
    f1=factors$f1; f1=factor(f1,unique(f1))
    f2=factors$f2; f2=factor(f2,unique(f2))
    f1a=factors$f1
    f2a=factors$f2
    nv1 <- length(summary(f1))
    nv2 <- length(summary(f2))
    lf1 <- levels(f1)
    lf2 <- levels(f2)
    mod=aov(response~f1*f2)
    #==================================================================
    mod1=aov(response~f1/f2)
    l1<-vector('list',nv1)
    names(l1)<-names(summary(f1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
      l1[[j]]<-v
      v<-numeric(0)}
    des1.tab<-summary(mod1,split=list('f1:f2'=l1))[[1]][-c(1,2),]

    cat(green("=========================================================="))
    cat(green("\nAnalyzing the interaction\n"))
    cat(green("==========================================================\n"))
    print(des1.tab)

    cat(green("\n\n==========================================================\n"))

    cate=as.list(1:nv1)
    modelo=anova(mod)
    for(i in 1:nv1){
      resp=response[f1==lf1[i]]
      trat=f2a[f1==lf1[i]]
      linear=lm(resp~trat)
      quadratico=lm(resp~trat+I(trat^2))
      cubico=lm(resp~trat+I(trat^2)+I(trat^3))
      plinear=data.frame(summary(linear)$coefficients)$Pr...t..
      pquad=data.frame(summary(quadratico)$coefficients)$Pr...t..
      pcubic=data.frame(summary(cubico)$coefficients)$Pr...t..
      # plinear=ifelse(plinear<0.0001,"p<0.0001",round(plinear,4))
      # pquad=ifelse(pquad<0.0001,"p<0.0001",round(pquad,4))
      # pcubic=ifelse(pcubic<0.0001,"p<0.0001",round(pcubic,4))
      a=anlinear(trat,resp,sq=modelo$`Sum Sq`[4],df=modelo$Df[4])
      b=anquad(trat,resp,sq=modelo$`Sum Sq`[4],df=modelo$Df[4])
      c=ancubic(trat,resp,sq=modelo$`Sum Sq`[4],df=modelo$Df[4])

      cate[[i]]=t(data.frame(intercept=c(round(coef(linear)[1],dec), round(coef(quadratico)[1],dec), round(coef(cubico)[1],dec)),
                             beta1=c(round(coef(linear)[2],dec), round(coef(quadratico)[2],dec), round(coef(cubico)[2],dec)),
                             beta2=c(NA,round(coef(quadratico)[3],dec), round(coef(cubico)[3],dec)),
                             beta3=c(NA, NA, round(coef(cubico)[4],dec)),
                             "p.value.test"=c(NA,NA,NA),
                             pb0=c(round(plinear[1],dec),round(pquad[1],dec),round(pcubic[1],dec)),
                             pb1=c(round(plinear[2],dec),round(pquad[2],dec),round(pcubic[2],dec)),
                             pb2=c(NA,round(pquad[3],dec),round(pcubic[3],dec)),
                             pb3=c(NA,NA,round(pcubic[4],dec)),
                             "p.value.mod"=c(NA,NA,NA),
                             linear=c(round(a$`Pr(>F)`[1],dec),round(b$`Pr(>F)`[1],dec),round(c$`Pr(>F)`[1],dec)),
                             quadratic=c(NA,round(b$`Pr(>F)`[2],dec),round(c$`Pr(>F)`[2],dec)),
                             cubic=c(NA,NA,round(c$`Pr(>F)`[3],dec)),
                             deviation=c(round(a$`Pr(>F)`[2],dec),round(b$`Pr(>F)`[3],dec),round(c$`Pr(>F)`[4],dec))))
      colnames(cate[[i]])=c("Linear","Quadratic","Cubic")
      }
    names(cate)=lf1
    print(cate,na.print = "")
    }
  if(length(factors)>2){
    f1=factors$f1; f1=factor(f1,unique(f1))
    f2=factors$f2; f2=factor(f2,unique(f2))
    block=factors$block
    bloco=factors$block; bloco=factor(bloco,unique(bloco))
    f1a=factors$f1
    f2a=factors$f2
    nv1 <- length(summary(f1))
    nv2 <- length(summary(f2))
    lf1 <- levels(f1)
    lf2 <- levels(f2)
    mod=aov(response~f1*f2+bloco)
    mod1=aov(response~f1/f2+bloco)
    l1<-vector('list',nv1)
    names(l1)<-names(summary(f1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
      l1[[j]]<-v
      v<-numeric(0)}
    des1.tab<-summary(mod1,split=list('f1:f2'=l1))[[1]][-c(1,2,3),]

    cat(green("=========================================================="))
    cat(green("\nDesdobramento\n"))
    cat(green("==========================================================\n"))
    print(des1.tab)

    cat(green("\n\n==========================================================\n"))

    cate=as.list(1:nv1)
    modelo=anova(mod)
    for(i in 1:nv1){
      resp=response[f1==lf1[i]]
      trat=f2a[f1==lf1[i]]
      linear=lm(resp~trat)
      quadratico=lm(resp~trat+I(trat^2))
      cubico=lm(resp~trat+I(trat^2)+I(trat^3))
      plinear=data.frame(summary(linear)$coefficients)$Pr...t..
      pquad=data.frame(summary(quadratico)$coefficients)$Pr...t..
      pcubic=data.frame(summary(cubico)$coefficients)$Pr...t..
      # plinear=ifelse(plinear<0.0001,"p<0.0001",round(plinear,4))
      # pquad=ifelse(pquad<0.0001,"p<0.0001",round(pquad,4))
      # pcubic=ifelse(pcubic<0.0001,"p<0.0001",round(pcubic,4))
      a=anlinear(trat,resp,sq=modelo$`Sum Sq`[5],df=modelo$Df[5])
      b=anquad(trat,resp,sq=modelo$`Sum Sq`[5],df=modelo$Df[5])
      c=ancubic(trat,resp,sq=modelo$`Sum Sq`[5],df=modelo$Df[5])
      cate[[i]]=t(data.frame(intercept=c(round(coef(linear)[1],dec), round(coef(quadratico)[1],dec), round(coef(cubico)[1],dec)),
                             beta1=c(round(coef(linear)[2],dec), round(coef(quadratico)[2],dec), round(coef(cubico)[2],dec)),
                             beta2=c(NA,round(coef(quadratico)[3],dec), round(coef(cubico)[3],dec)),
                             beta3=c(NA, NA, round(coef(cubico)[4],dec)),
                             "p.value.test"=c(NA,NA,NA),
                             pb0=c(round(plinear[1],dec),round(pquad[1],dec),round(pcubic[1],dec)),
                             pb1=c(round(plinear[2],dec),round(pquad[2],dec),round(pcubic[2],dec)),
                             pb2=c(NA,round(pquad[3],dec),round(pcubic[3],dec)),
                             pb3=c(NA,NA,round(pcubic[4],dec)),
                             "p.value.mod"=c(NA,NA,NA),
                             linear=c(round(a$`Pr(>F)`[1],dec),round(b$`Pr(>F)`[1],dec),round(c$`Pr(>F)`[1],dec)),
                             quadratic=c(NA,round(b$`Pr(>F)`[2],dec),round(c$`Pr(>F)`[2],dec)),
                             cubic=c(NA,NA,round(c$`Pr(>F)`[3],dec)),
                             deviation=c(round(a$`Pr(>F)`[2],dec),round(b$`Pr(>F)`[3],dec),round(c$`Pr(>F)`[4],dec))))
      colnames(cate[[i]])=c("Linear","Quadratic","Cubic")}
    names(cate)=lf1
    print(cate,na.print = "")}
}
