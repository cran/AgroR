#' Analysis: DBC experiments in split-split-plot
#' @description Analysis of an experiment conducted in a randomized block design in a split-split-plot scheme using analysis of variance of fixed effects.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with plot levels
#' @param f2 Numeric or complex vector with splitplot levels
#' @param f3 Numeric or complex vector with splitsplitplot levels
#' @param block Numeric or complex vector with blocks
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD and Duncan)
#' @param response Numeric vector with responses
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param dec Number of cells
#' @note The PSUBSUBDBC function does not present residual analysis, interaction breakdown, graphs and implementations of various multiple comparison or regression tests. The function only returns the analysis of variance and multiple comparison test of Tukey, LSD or Duncan.
#' @return Analysis of variance of fixed effects and multiple comparison test of Tukey, LSD or Duncan.
#' @keywords DBC
#' @export
#' @examples
#' library(AgroR)
#' data(enxofre)
#' attach(enxofre)
#' PSUBSUBDBC(f1, f2, f3, bloco, resp)

PSUBSUBDBC=function(f1,
                    f2,
                    f3,
                    block,
                    response,
                    alpha.f=0.05,
                    alpha.t=0.05,
                    dec=3,
                    mcomp="tukey"){
  fac.names=c("F1","F2","F3")
  fator1=as.factor(f1)
  fator2=as.factor(f2)
  fator3=as.factor(f3)
  bloco=as.factor(block)
  fatores<-data.frame(fator1,fator2,fator3)
  Fator1<-factor(fator1,levels=unique(fator1))
  Fator2<-factor(fator2,levels=unique(fator2))
  Fator3<-factor(fator3,levels=unique(fator3))
  nv1<-length(summary(Fator1))
  nv2<-length(summary(Fator2))
  nv3<-length(summary(Fator3))
  J<-(length(response))/(nv1*nv2*nv3)
  lf1<-levels(Fator1); lf2<-levels(Fator2); lf3<-levels(Fator3)
  mod=aov(response~Fator1*Fator2*Fator3+
            Error(bloco/Fator1/paste(Fator1,Fator2)))
  a=summary(mod)
  anava=rbind(data.frame(a$`Error: bloco:Fator1`[[1]]),
              data.frame(a$`Error: bloco:Fator1:paste(Fator1, Fator2)`[[1]]),
              data.frame(a$`Error: Within`[[1]]))
  anava$F.value=ifelse(is.na(anava$F.value)==TRUE,"",round(anava$F.value,5))
  anava$Pr..F.=ifelse(is.na(anava$Pr..F.)==TRUE,"",round(anava$Pr..F.,5))
  colnames(anava)=c("Df","Sum Sq","Mean Sq","F value","Pr(>F)")
  print(anava)

  fatores<-data.frame('fator 1'=fator1,
                      'fator 2' = fator2,
                      'fator 3' = fator3)
  qmres=c(as.numeric(anava[2,3]),
          as.numeric(anava[5,3]),
          as.numeric(anava[10,3]))
  GL=c(as.numeric(anava[2,1]),
       as.numeric(anava[5,1]),
       as.numeric(anava[10,1]))
  pvalor=c(as.numeric(anava[1,5]),
           as.numeric(anava[3,5]),
           as.numeric(anava[6,5]))

  ################################################################################################
  # Efeitos simples
  ################################################################################################
  if(as.numeric(anava[9,5])>alpha.f &&
     as.numeric(anava[8,5])>alpha.f &&
     as.numeric(anava[7,5])>alpha.f &&
     as.numeric(anava[4,5])>alpha.f) {
    graficos=list(1,2,3)
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold('Non-significant interaction: analyzing the simple effects')))
    cat(green(bold("\n------------------------------------------\n")))

    for(i in 1:3){
      # Comparação múltipla
      if(pvalor[i]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))
        if(mcomp=="tukey"){
        letra=HSD.test(response,
                       fatores[,i],
                       GL[i],
                       qmres[i],alpha.t)
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        print(letra1)}
        if(mcomp=="lsd"){
          letra=LSD.test(response,
                         fatores[,i],
                         GL[i],
                         qmres[i],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          print(letra1)}
        if(mcomp=="duncan"){
          letra=duncan.test(response,
                         fatores[,i],
                         GL[i],
                         qmres[i],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          print(letra1)}
        }
      if(pvalor[i]>alpha.f) {
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))
        mean.table<-tapply.stat(response,fatores[,i],mean)
        colnames(mean.table)<-c('Levels','Mean')
        print(mean.table)
        grafico=NA
        cat(green(bold("\n------------------------------------------")))}
      cat('\n')
    }
  }

  #####################################################################
  #Interacao Fator1*Fator2      +     Fator3
  #####################################################################
  if(as.numeric(anava[9,5])>alpha.f &&
     as.numeric(anava[4,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[2],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

    if(mcomp=="tukey"){
    tukeygrafico=c()
    ordem=c()
    for (i in 1:nv2) {
      trati=fatores[, 1][Fator2 == lf2[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator2 == lf2[i]]
      tukey=HSD.test(respi,trati,GL[2],qmres[2],alpha.t)
      tukeygrafico[[i]]=tukey$groups[levels(trati),2]
      ordem[[i]]=rownames(tukey$groups[levels(trati),])
    }
    letra=unlist(tukeygrafico)
    datag=data.frame(letra,ordem=unlist(ordem))
    datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
    datag=datag[order(datag$ordem),]
    letra=datag$letra

    tukeygrafico1=c()
    for (i in 1:nv1) {
      trati=fatores[, 2][Fator1 == lf1[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator1 == lf1[i]]
      tukey=HSD.test(respi,trati,GL[2],qmres[2],alpha.t)
      tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
    }
    letra1=unlist(tukeygrafico1)
    letra1=toupper(letra1)}
    if(mcomp=="lsd"){
      lsdgrafico=c()
      ordem=c()
      for (i in 1:nv2) {
        trati=fatores[, 1][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        lsd=LSD.test(respi,trati,GL[2],qmres[2],alpha.t)
        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
        ordem[[i]]=rownames(lsd$groups[levels(trati),])
      }
      letra=unlist(lsdgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      lsdgrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 2][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        lsd=LSD.test(respi,trati,GL[2],qmres[2],alpha.t)
        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]
      }
      letra1=unlist(lsdgrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="tukey"){
      duncangrafico=c()
      ordem=c()
      for (i in 1:nv2) {
        trati=fatores[, 1][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        duncan=duncan.test(respi,trati,GL[2],qmres[2],alpha.t)
        duncangrafico[[i]]=duncan$groups[levels(trati),2]
        ordem[[i]]=rownames(duncan$groups[levels(trati),])
      }
      letra=unlist(duncangrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      duncangrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 2][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        duncan=duncan.test(respi,trati,GL[2],qmres[2],alpha.t)
        duncangrafico1[[i]]=duncan$groups[levels(trati),2]
      }
      letra1=unlist(duncangrafico1)
      letra1=toupper(letra1)}

    f1=rep(levels(Fator1),e=length(levels(Fator2)))
    f2=rep(unique(as.character(Fator2)),length(levels(Fator2)))
    media=tapply(response,paste(Fator1,Fator2), mean, na.rm=TRUE)[unique(paste(f1,f2))]
    desvio=tapply(response,paste(Fator1,Fator2), sd, na.rm=TRUE)[unique(paste(f1,f2))]
    f1=factor(f1,levels = unique(f1))
    f2=factor(f2,levels = unique(f2))

    graph=data.frame(f1=f1,
                     f2=f2,
                     media,
                     desvio,
                     letra,letra1,
                     numero=format(media,digits = dec))
    numero=graph$numero
    letras=paste(graph$letra, graph$letra1, sep="")
    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),
                               ncol = length(levels(Fator1)))))
    rownames(matriz)=levels(Fator1)
    colnames(matriz)=levels(Fator2)
    print(matriz)
    message(black("\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the", mcomp, "(p<",alpha.t,")\n"))


    #Checar o Fator3
    if(as.numeric(anava[7,5])>alpha.f && as.numeric(anava[8,5])>alpha.f) {
      if(pvalor[3]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[3])))
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        if(mcomp=="tukey"){letra=HSD.test(response,fatores[,i],GL[3],qmres[3],alpha.t)}
        if(mcomp=="lsd"){letra=LSD.test(response,fatores[,i],GL[3],qmres[3],alpha.t)}
        if(mcomp=="duncan"){letra=duncan.test(response,fatores[,i],GL[3],qmres[3],alpha.t)}
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        print(letra1)
        cat(green(bold("\n-----------------------------------------------------------------")))
      }
    }
  }

  #####################################################################################################
  #Interacao Fator1*Fator3       + fator2
  #####################################################################################################
  if(as.numeric(anava[9,5])>alpha.f &&
     as.numeric(anava[7,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("\nInteraction",paste(fac.names[1],'*',fac.names[3],sep='')," significant: unfolding the interaction\n")))
    cat(green(bold("\n------------------------------------------\n")))

    if(mcomp=="tukey"){
    tukeygrafico=c()
    ordem=c()
    for (i in 1:nv3) {
      trati=fatores[, 1][Fator3 == lf3[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator3 == lf3[i]]
      tukey=HSD.test(respi,trati,GL[3],qmres[3],alpha.t)
      tukeygrafico[[i]]=tukey$groups[levels(trati),2]
      ordem[[i]]=rownames(tukey$groups[levels(trati),])}
    letra=unlist(tukeygrafico)
    datag=data.frame(letra,ordem=unlist(ordem))
    datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
    datag=datag[order(datag$ordem),]
    letra=datag$letra

    tukeygrafico1=c()
    for (i in 1:nv1) {
      trati=fatores[, 3][Fator1 == lf1[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator1 == lf1[i]]
      tukey=HSD.test(respi,trati,GL[3],qmres[3],alpha.t)
      tukeygrafico1[[i]]=tukey$groups[levels(trati),2]}
    letra1=unlist(tukeygrafico1)
    letra1=toupper(letra1)}
    if(mcomp=="lsd"){
      lsdgrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 1][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        lsd=LSD.test(respi,trati,GL[3],qmres[3],alpha.t)
        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
        ordem[[i]]=rownames(lsd$groups[levels(trati),])}
      letra=unlist(lsdgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      lsdgrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 3][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        lsd=LSD.test(respi,trati,GL[3],qmres[3],alpha.t)
        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]}
      letra1=unlist(lsdgrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="duncan"){
      duncangrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 1][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        duncan=duncan.test(respi,trati,GL[3],qmres[3],alpha.t)
        duncangrafico[[i]]=duncan$groups[levels(trati),2]
        ordem[[i]]=rownames(duncan$groups[levels(trati),])}
      letra=unlist(duncangrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      duncangrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 3][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        duncan=duncan.test(respi,trati,GL[3],qmres[3],alpha.t)
        duncangrafico1[[i]]=duncan$groups[levels(trati),2]}
      letra1=unlist(duncangrafico1)
      letra1=toupper(letra1)}

    f1=rep(levels(Fator1),e=length(levels(Fator3)))
    f3=rep(unique(as.character(Fator3)),length(levels(Fator1)))
    media=tapply(response,paste(Fator1,Fator3), mean, na.rm=TRUE)[unique(paste(f1,f3))]
    desvio=tapply(response,paste(Fator1,Fator3), sd, na.rm=TRUE)[unique(paste(f1,f3))]
    f1=factor(f1,levels = unique(f1))
    f3=factor(f3,levels = unique(f3))
    graph=data.frame(f1=f1, f3=f3,
                     media, desvio,
                     letra,letra1,
                     numero=format(media,digits = dec))
    numero=graph$numero
    letras=paste(graph$letra,graph$letra1,sep="")
    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
    rownames(matriz)=levels(Fator1)
    colnames(matriz)=levels(Fator3)
    print(matriz)
    message(black("\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the", mcomp, "(p<",alpha.t,")\n"))

    #Checar o Fator2
    if(as.numeric(anava[4,5])>alpha.f && as.numeric(anava[6,5])>alpha.f) {
      i=2
      cat(green(bold("\n------------------------------------------\n")))
      cat(green(italic('Analyzing the simple effects of the factor ',fac.names[2])))
      cat(green(bold("\n------------------------------------------\n")))
      cat(fac.names[i])
      if(mcomp=="tukey"){letra=HSD.test(response,fatores[,i], GL[2], qmres[2],alpha.t)}
      if(mcomp=="lsd"){letra=LSD.test(response,fatores[,i], GL[2], qmres[2],alpha.t)}
      if(mcomp=="duncan"){letra=duncan.test(response,fatores[,i], GL[2], qmres[2],alpha.t)}
      letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
    print(letra1)
    cat(green(bold("\n-----------------------------------------------------------------")))
  }

  ######################################################################################################################
  #Interacao Fator2*Fator3     + fator1
  ######################################################################################################################
  if(as.numeric(anava[9,5])>alpha.f &&
     as.numeric(anava[8,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

    if(mcomp=="tukey"){
    tukeygrafico=c()
    ordem=c()
    for (i in 1:nv3) {
      trati=fatores[, 2][Fator3 == lf3[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator3 == lf3[i]]
      mod=aov(respi~trati)
      tukey=HSD.test(respi,trati,GL[3],qmres[3],alpha.t)
      tukeygrafico[[i]]=tukey$groups[levels(trati),2]
      ordem[[i]]=rownames(tukey$groups[levels(trati),])}
    letra=unlist(tukeygrafico)
    datag=data.frame(letra,ordem=unlist(ordem))
    datag=data.frame(letra,ordem=unlist(ordem))
    datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
    datag=datag[order(datag$ordem),]
    letra=datag$letra
    tukeygrafico1=c()
    for (i in 1:nv2) {
      trati=fatores[, 3][Fator2 == lf2[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator2 == lf2[i]]
      mod=aov(respi~trati)
      tukey=HSD.test(respi,trati,GL[3],qmres[3],alpha.t)
      tukeygrafico1[[i]]=tukey$groups[levels(trati),2]}
    letra1=unlist(tukeygrafico1)
    letra1=toupper(letra1)}
    if(mcomp=="lsd"){
      lsdgrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 2][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        mod=aov(respi~trati)
        lsd=LSD.test(respi,trati,GL[3],qmres[3],alpha.t)
        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
        ordem[[i]]=rownames(lsd$groups[levels(trati),])}
      letra=unlist(lsdgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra
      lsdgrafico1=c()
      for (i in 1:nv2) {
        trati=fatores[, 3][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        mod=aov(respi~trati)
        lsd=LSD.test(respi,trati,GL[3],qmres[3],alpha.t)
        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]}
      letra1=unlist(lsdgrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="duncan"){
      duncangrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 2][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        mod=aov(respi~trati)
        duncan=duncan.test(respi,trati,GL[3],qmres[3],alpha.t)
        duncangrafico[[i]]=duncan$groups[levels(trati),2]
        ordem[[i]]=rownames(duncan$groups[levels(trati),])}
      letra=unlist(duncangrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra
      duncangrafico1=c()
      for (i in 1:nv2) {
        trati=fatores[, 3][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        mod=aov(respi~trati)
        duncan=duncan.test(respi,trati,GL[3],qmres[3],alpha.t)
        duncangrafico1[[i]]=duncan$groups[levels(trati),2]}
      letra1=unlist(duncangrafico1)
      letra1=toupper(letra1)}

    f2=rep(levels(Fator2),e=length(levels(Fator3)))
    f3=rep(unique(as.character(Fator3)),length(levels(Fator2)))
    media=tapply(response,paste(Fator2,Fator3), mean, na.rm=TRUE)[unique(paste(f2,f3))]
    desvio=tapply(response,paste(Fator2,Fator3), sd, na.rm=TRUE)[unique(paste(f2,f3))]
    f2=factor(f2,levels = unique(f2))
    f3=factor(f3,levels = unique(f3))
    graph=data.frame(f2=f2,
                     f3=f3,
                     media,
                     desvio,
                     letra,letra1,
                     numero=format(media,digits = dec))
    numero=graph$numero
    letras=paste(graph$letra,graph$letra1,sep="")
    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator2)))))
    rownames(matriz)=levels(Fator2)
    colnames(matriz)=levels(Fator3)
    print(matriz)
    message(black("\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the Tukey (p<",alpha.t,")\n"))
    #Checar o Fator1
    if(as.numeric(anava[4,5])>alpha.f && as.numeric(anava[7,5])>alpha.f) {
      i<-1
      if(pvalor[i]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[2])))
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        if(mcomp=="tukey"){letra=HSD.test(response,fatores[,i],GL[1],qmres[1],alpha.t)}
        if(mcomp=="lsd"){letra=LSD.test(response,fatores[,i],GL[1],qmres[1],alpha.t)}
        if(mcomp=="duncan"){letra=duncan.test(response,fatores[,i],GL[1],qmres[1],alpha.t)}
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        print(letra1)
        cat(green(bold("\n------------------------------------------\n")))
      }
    }
  }

  #########################################################################################################################
  #Para interacao tripla significativa, desdobramento
  #########################################################################################################################
  if(as.numeric(anava[9,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' inside of each level of ', fac.names[2], 'and',fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))
    ii<-0
    for(i in 1:nv2) {
      for(j in 1:nv3) {
        ii<-ii+1
        cat('\n',fac.names[1],' within the combination of levels ',lf2[i],' of  ',fac.names[2],' and ',lf3[j],' of  ',fac.names[3],"\n")
        if(mcomp=="tukey"){
        tukey=HSD.test(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                       fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                       GL[3],
                       GL[3],
                       alpha.t)
        tukey=tukey$groups;colnames(tukey)=c("Mean","letters")
        print(tukey)}
        if(mcomp=="lsd"){
          lsd=LSD.test(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                         fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                         GL[3],
                         GL[3],
                         alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("Mean","letters")
          print(lsd)}
        if(mcomp=="duncan"){
          duncan=duncan.test(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                         fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                         GL[3],
                         GL[3],
                         alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("Mean","letters")
          print(duncan)}

        }
    }

    cat('\n\n')

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' inside of each level of ', fac.names[1], 'and',fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))

    ii<-0
    for(k in 1:nv1) {
      for(j in 1:nv3) {
        ii<-ii+1
        cat('\n\n',fac.names[2],' within the combination of levels ',lf1[k],' of  ',fac.names[1],' and ',lf3[j],' of  ',fac.names[3],'\n')
        if(mcomp=="tukey"){
        tukey=HSD.test(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                       fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                       GL[3],
                       qmres[3],
                       alpha.t)
        tukey=tukey$groups;colnames(tukey)=c("Mean","letters")
        print(tukey)}
        if(mcomp=="lsd"){
          lsd=LSD.test(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                         fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                         GL[3],
                         qmres[3],
                         alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("Mean","letters")
          print(lsd)}
        if(mcomp=="duncan"){
          duncan=duncan.test(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                         fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                         GL[3],
                         qmres[3],
                         alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("Mean","letters")
          print(duncan)}

      }
    }

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' inside of each level of ', fac.names[1], 'and',fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))
    ii<-0
    for(k in 1:nv1) {
      for(i in 1:nv2) {
        ii<-ii+1
        cat('\n\n',fac.names[3],' within the combination of levels ',lf1[k],' of ',fac.names[1],' and ',lf2[i],' of  ',fac.names[2],'\n')
        if(mcomp=="tukey"){
        tukey=HSD.test(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                       fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                       GL[3],
                       qmres[3],
                       alpha.t)
        tukey=tukey$groups;colnames(tukey)=c("Mean","letters")
        print(tukey)}
        if(mcomp=="lsd"){
          lsd=LSD.test(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         GL[3],
                         qmres[3],
                         alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("Mean","letters")
          print(lsd)}
        if(mcomp=="duncan"){
          duncan=duncan.test(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         GL[3],
                         qmres[3],
                         alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("Mean","letters")
          print(duncan)}
        }
    }
  }
}
