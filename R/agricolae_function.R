#=====================================================================
# tapply.stat
#=====================================================================

tapply.stat <-
  function (y, x, stat = "mean")
  {
    k<-0
    numerico<- NULL
    if(is.null(ncol(x))){
      if(is.numeric(x)){
        k<-1
        numerico[1]<-1
      }
    }
    else {
      ncolx<-ncol(x)
      for (i in 1:ncolx) {
        if(is.numeric(x[,i])){
          k<-k+1
          numerico[k]<-i
        }}}
    cx <- deparse(substitute(x))
    cy <- deparse(substitute(y))
    x <- data.frame(c1 = 1, x)
    y <- data.frame(v1 = 1, y)
    nx <- ncol(x)
    ny <- ncol(y)
    namex <- names(x)
    namey <- names(y)
    if (nx == 2)
      namex <- c("c1", cx)
    if (ny == 2)
      namey <- c("v1", cy)
    namexy <- c(namex, namey)
    for (i in 1:nx) {
      x[, i] <- as.character(x[, i])
    }
    z <- NULL
    for (i in 1:nx) {
      z <- paste(z, x[, i], sep = "&")
    }
    w <- NULL
    for (i in 1:ny) {
      m <- tapply(y[, i], z, stat)
      m <- as.matrix(m)
      w <- cbind(w, m)
    }
    nw <- nrow(w)
    c <- rownames(w)
    v <- rep("", nw * nx)
    dim(v) <- c(nw, nx)
    for (i in 1:nw) {
      for (j in 1:nx) {
        v[i, j] <- strsplit(c[i], "&")[[1]][j + 1]
      }
    }
    rownames(w) <- NULL
    junto <- data.frame(v[, -1], w)
    junto <- junto[, -nx]
    names(junto) <- namexy[c(-1, -(nx + 1))]
    if(k==1 & nx==2) {
      junto[,numerico[1]]<-as.character(junto[,numerico[1]])
      junto[,numerico[1]]<-as.numeric(junto[,numerico[1]])
      junto<-junto[order(junto[,1]),]
    }
    if (k>0 & nx > 2) {
      for (i in 1:k){
        junto[,numerico[i]]<-as.character(junto[,numerico[i]])
        junto[,numerico[i]]<-as.numeric(junto[,numerico[i]])
      }
      junto<-junto[do.call("order", c(junto[,1:(nx-1)])),]
    }
    rownames(junto)<-1:(nrow(junto))
    return(junto)
  }

#=====================================================================
# orderPvalue
#=====================================================================

orderPvalue <-
  function (treatment, means, alpha, pvalue, console)
  {
    n <- length(means)
    z <- data.frame(treatment, means)
    letras<-c(letters[1:26],LETTERS[1:26],1:9,c(".","+","-","*","/","#","$",
                                                "%","&","^","[","]",":","@",";","_","?","!","=","#",rep(" ",2000)))
    w <- z[order(z[, 2], decreasing = TRUE), ]
    M<-rep("",n)
    k<-1
    k1<-0
    j<-1
    i<-1
    cambio<-n
    cambio1<-0
    chequeo=0
    M[1]<-letras[k]
    q <- as.numeric(rownames(w)) #Check
    while(j<n) {
      chequeo<-chequeo+1
      if (chequeo > n) break
      for(i in j:n) {
        s<-pvalue[q[i],q[j]]>alpha
        if(s) {
          if(lastC(M[i]) != letras[k])M[i]<-paste(M[i],letras[k],sep="")
        }
        else {
          k<-k+1
          cambio<-i
          cambio1<-0
          ja<-j
          for(jj in cambio:n) M[jj]<-paste(M[jj],"",sep="") # El espacio
          M[cambio]<-paste(M[cambio],letras[k],sep="")
          for( v in ja:cambio) {
            if(pvalue[q[v],q[cambio]]<=alpha) {j<-j+1
            cambio1<-1
            }
            else break
          }
          break
        }
      }
      if (cambio1 ==0 )j<-j+1
    }
    #-----------
    w<-data.frame(w,stat=M)
    trt <- as.character(w$treatment)
    means <- as.numeric(w$means)
    output <- data.frame(means, groups=M)
    rownames(output)<-trt
    if(k>81)
      cat("\n",k,"groups are estimated.The number of groups exceeded the maximum of 81 labels. change to group=FALSE.\n")
    invisible(output)
  }

#=====================================================================
# lastC
#=====================================================================

lastC <-
  function(x) {
    y<-sub(" +$", "",x)
    p1<-nchar(y)
    cc<-substr(y,p1,p1)
    return(cc)
  }

#=====================================================================
# duncan.test
#=====================================================================

duncan.test <-
  function (y, trt, DFerror, MSerror, alpha=0.05, group=TRUE,main = NULL,console=FALSE)
  {
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    clase<-c("aov","lm")
    if("aov"%in%class(y) | "lm"%in%class(y)){
      if(is.null(main))main<-y$call
      A<-y$model
      DFerror<-df.residual(y)
      MSerror<-deviance(y)/DFerror
      y<-A[,1]
      ipch<-pmatch(trt,names(A))
      nipch<- length(ipch)
      for(i in 1:nipch){
        if (is.na(ipch[i]))
          return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))
      }
      name.t<- names(A)[ipch][1]
      trt <- A[, ipch]
      if (nipch > 1){
        trt <- A[, ipch[1]]
        for(i in 2:nipch){
          name.t <- paste(name.t,names(A)[ipch][i],sep=":")
          trt <- paste(trt,A[,ipch[i]],sep=":")
        }}
      name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    Mean<-mean(junto[,1])
    CV<-sqrt(MSerror)*100/Mean
    medians<-tapply.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- tapply.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
    sds <-   tapply.stat(junto[,1],junto[,2],stat="sd") #change
    nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
    means<-data.frame(means,std=sds[,2],r=nn[,2],medians)
    names(means)[1:2]<-c(name.t,name.y)
    ntr<-nrow(means)
    Tprob<-NULL
    k<-0
    for(i in 2:ntr){
      k<-k+1
      x <- suppressWarnings(warning(qtukey((1-alpha)^(i-1), i, DFerror)))
      if(x=="NaN")break
      else Tprob[k]<-x
    }
    if(k<(ntr-1)){
      for(i in k:(ntr-1)){
        f <- Vectorize(function(x)ptukey(x,i+1,DFerror)-(1-alpha)^i)
        Tprob[i]<-uniroot(f, c(0,100))$root
      }
    }
    Tprob<-as.numeric(Tprob)
    nr <- unique(nn[,2])
    # Critical Value of Studentized Range
    if(console){
      cat("\nStudy:", main)
      cat("\n\nDuncan's new multiple range test\nfor",name.y,"\n")
      cat("\nMean Square Error: ",MSerror,"\n\n")
      cat(paste(name.t,",",sep="")," means\n\n")
      print(data.frame(row.names = means[,1], means[,2:6]))
    }
    if(length(nr) == 1 ) sdtdif <- sqrt(MSerror/nr)
    else {
      nr1 <-  1/mean(1/nn[,2])
      sdtdif <- sqrt(MSerror/nr1)
    }
    DUNCAN <- Tprob * sdtdif
    names(DUNCAN)<-2:ntr
    duncan<-data.frame(Table=Tprob,CriticalRange=DUNCAN)
    if ( group & length(nr) == 1 & console){
      cat("\nAlpha:",alpha,"; DF Error:",DFerror,"\n")
      cat("\nCritical Range\n")
      print(DUNCAN)
    }
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")
    if ( length(nr) != 1) duncan<-NULL
    Omeans<-order(means[,2],decreasing = TRUE) #correccion 2019, 1 abril.
    Ordindex<-order(Omeans)
    comb <-utils::combn(ntr,2)
    nn<-ncol(comb)
    dif<-rep(0,nn)
    DIF<-dif
    LCL<-dif
    UCL<-dif
    pvalue<-dif
    odif<-dif
    sig<-NULL
    for (k in 1:nn) {
      i<-comb[1,k]
      j<-comb[2,k]
      dif[k]<-means[i,2]-means[j,2]
      DIF[k]<-abs(dif[k])
      nx<-abs(i-j)+1
      odif[k] <- abs(Ordindex[i]- Ordindex[j])+1
      pvalue[k]<- round(1-ptukey(DIF[k]/sdtdif,odif[k],DFerror)^(1/(odif[k]-1)),4)
      LCL[k] <- dif[k] - DUNCAN[odif[k]-1]
      UCL[k] <- dif[k] + DUNCAN[odif[k]-1]
      sig[k]<-" "
      if (pvalue[k] <= 0.001) sig[k]<-"***"
      else  if (pvalue[k] <= 0.01) sig[k]<-"**"
      else  if (pvalue[k] <= 0.05) sig[k]<-"*"
      else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
    if(!group){
      tr.i <- means[comb[1, ],1]
      tr.j <- means[comb[2, ],1]
      comparison<-data.frame("difference" = dif, pvalue=pvalue,"signif."=sig,LCL,UCL)
      rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
      if(console){cat("\nComparison between treatments means\n\n")
        print(comparison)}
      groups=NULL
    }
    if (group) {
      comparison=NULL
      # The probability matrix
      Q<-matrix(1,ncol=ntr,nrow=ntr)
      p<-pvalue
      k<-0
      for(i in 1:(ntr-1)){
        for(j in (i+1):ntr){
          k<-k+1
          Q[i,j]<-p[k]
          Q[j,i]<-p[k]
        }
      }
      groups <- orderPvalue(means[, 1], means[, 2],alpha, Q,console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nMeans with the same letter are not significantly different.\n\n")
        print(groups)
      }
    }
    parameters<-data.frame(test="Duncan",name.t=name.t,ntr = ntr,alpha=alpha)
    statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    rownames(means)<-means[,1]
    means<-means[,-1]
    output<-list(statistics=statistics,parameters=parameters, duncan=duncan,
                 means=means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }

#=====================================================================
# HSD.test
#=====================================================================

HSD.test <-
  function (y, trt, DFerror, MSerror, alpha=0.05, group=TRUE,main = NULL,unbalanced=FALSE,console=FALSE)
  {
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    clase<-c("aov","lm")
    if("aov"%in%class(y) | "lm"%in%class(y)){
      if(is.null(main))main<-y$call
      A<-y$model
      DFerror<-df.residual(y)
      MSerror<-deviance(y)/DFerror
      y<-A[,1]
      ipch<-pmatch(trt,names(A))
      nipch<- length(ipch)
      for(i in 1:nipch){
        if (is.na(ipch[i]))
          return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))
      }
      name.t<- names(A)[ipch][1]
      trt <- A[, ipch]
      if (nipch > 1){
        trt <- A[, ipch[1]]
        for(i in 2:nipch){
          name.t <- paste(name.t,names(A)[ipch][i],sep=":")
          trt <- paste(trt,A[,ipch[i]],sep=":")
        }}
      name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    Mean<-mean(junto[,1])
    CV<-sqrt(MSerror)*100/Mean
    medians<-tapply.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- tapply.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
    sds <-   tapply.stat(junto[,1],junto[,2],stat="sd") #change
    nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
    means<-data.frame(means,std=sds[,2],r=nn[,2],medians)
    names(means)[1:2]<-c(name.t,name.y)
    #   row.names(means)<-means[,1]
    ntr<-nrow(means)
    Tprob <- qtukey(1-alpha,ntr, DFerror)
    nr<-unique(nn[, 2])
    nr1<-1/mean(1/nn[,2])
    if(console){
      cat("\nStudy:", main)
      cat("\n\nHSD Test for",name.y,"\n")
      cat("\nMean Square Error: ",MSerror,"\n\n")
      cat(paste(name.t,",",sep="")," means\n\n")
      print(data.frame(row.names = means[,1], means[,2:6]))
      cat("\nAlpha:",alpha,"; DF Error:",DFerror,"\n")
      cat("Critical Value of Studentized Range:", Tprob,"\n")
    }
    HSD <- Tprob * sqrt(MSerror/nr)
    statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV,MSD=HSD)
    if ( group & length(nr) == 1 & console) cat("\nMinimun Significant Difference:",HSD,"\n")
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")
    if ( length(nr) != 1) statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
    comb <-utils::combn(ntr,2)
    nn<-ncol(comb)
    dif<-rep(0,nn)
    sig<-NULL
    LCL<-dif
    UCL<-dif
    pvalue<-rep(0,nn)
    for (k in 1:nn) {
      i<-comb[1,k]
      j<-comb[2,k]
      #if (means[i, 2] < means[j, 2]){
      #comb[1, k]<-j
      #comb[2, k]<-i
      #}
      dif[k]<-means[i,2]-means[j,2]
      sdtdif<-sqrt(MSerror * 0.5*(1/means[i,4] + 1/means[j,4]))
      if(unbalanced)sdtdif<-sqrt(MSerror /nr1)
      pvalue[k]<- round(1-ptukey(abs(dif[k])/sdtdif,ntr,DFerror),4)
      LCL[k] <- dif[k] - Tprob*sdtdif
      UCL[k] <- dif[k] + Tprob*sdtdif
      sig[k]<-" "
      if (pvalue[k] <= 0.001) sig[k]<-"***"
      else  if (pvalue[k] <= 0.01) sig[k]<-"**"
      else  if (pvalue[k] <= 0.05) sig[k]<-"*"
      else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
    if(!group){
      tr.i <- means[comb[1, ],1]
      tr.j <- means[comb[2, ],1]
      comparison<-data.frame("difference" = dif, pvalue=pvalue,"signif."=sig,LCL,UCL)
      rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
      if(console){cat("\nComparison between treatments means\n\n")
        print(comparison)}
      groups=NULL
      #groups<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,4],std.err=means[,3])
    }
    if (group) {
      comparison=NULL
      # Matriz de probabilidades
      Q<-matrix(1,ncol=ntr,nrow=ntr)
      p<-pvalue
      k<-0
      for(i in 1:(ntr-1)){
        for(j in (i+1):ntr){
          k<-k+1
          Q[i,j]<-p[k]
          Q[j,i]<-p[k]
        }
      }
      groups <- orderPvalue(means[, 1], means[, 2],alpha, Q,console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nTreatments with the same letter are not significantly different.\n\n")
        print(groups)
      }
    }
    parameters<-data.frame(test="Tukey",name.t=name.t,ntr = ntr, StudentizedRange=Tprob,alpha=alpha)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    rownames(means)<-means[,1]
    means<-means[,-1]
    output<-list(statistics=statistics,parameters=parameters,
                 means=means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }

#=====================================================================
# kruskal
#=====================================================================

kruskal <-
  function (y, trt, alpha = 0.05, p.adj = c("none", "holm", "hommel", "hochberg",
                                            "bonferroni", "BH", "BY", "fdr"), group = TRUE, main = NULL,console=FALSE)
  {
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    p.adj <- match.arg(p.adj)
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    N <- nrow(junto)
    medians<-tapply.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- tapply.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    Means <- tapply.stat(junto[,1],junto[,2],stat="mean")  #change
    sds <-   tapply.stat(junto[,1],junto[,2], stat="sd")  #change
    nn <-   tapply.stat(junto[,1],junto[,2],stat="length") #change
    Means<-data.frame(Means,std=sds[,2],r=nn[,2],medians)
    rownames(Means)<-Means[,1]
    Means<-Means[,-1]
    names(Means)[1]<-name.y
    junto[, 1] <- rank(junto[, 1])
    means <- tapply.stat(junto[, 1], junto[, 2], stat = "sum")
    sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
    nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
    means <- data.frame(means, r = nn[, 2])
    names(means)[1:2] <- c(name.t, name.y)
    ntr <- nrow(means)
    nk <- choose(ntr, 2)
    DFerror <- N - ntr
    rs <- 0
    U <- 0
    for (i in 1:ntr) {
      rs <- rs + means[i, 2]^2/means[i, 3]
      U <- U + 1/means[i, 3]
    }
    S <- (sum(junto[, 1]^2) - (N * (N + 1)^2)/4)/(N - 1)
    H <- (rs - (N * (N + 1)^2)/4)/S
    p.chisq <- 1 - pchisq(H, ntr - 1)
    if(console){
      cat("\nStudy:", main)
      cat("\nKruskal-Wallis test's\nTies or no Ties\n")
      cat("\nCritical Value:", H)
      cat("\nDegrees of freedom:", ntr - 1)
      cat("\nPvalue Chisq  :", p.chisq, "\n\n")
    }
    DFerror <- N - ntr
    Tprob <- qt(1 - alpha/2, DFerror)
    MSerror <- S * ((N - 1 - H)/(N - ntr))
    means[, 2] <- means[, 2]/means[, 3]
    if(console){cat(paste(name.t, ",", sep = ""), " means of the ranks\n\n")
      print(data.frame(row.names = means[, 1], means[, -1]))
      cat("\nPost Hoc Analysis\n")
    }
    if (p.adj != "none") {
      if(console)cat("\nP value adjustment method:", p.adj)
      a <- 1e-06
      b <- 1
      for (i in 1:100) {
        x <- (b + a)/2
        xr <- rep(x, nk)
        d <- p.adjust(xr, p.adj)[1] - alpha
        ar <- rep(a, nk)
        fa <- p.adjust(ar, p.adj)[1] - alpha
        if (d * fa < 0)
          b <- x
        if (d * fa > 0)
          a <- x
      }
      Tprob <- qt(1 - x/2, DFerror)
    }
    nr <- unique(means[, 3])

    if (group & console){
      cat("\nt-Student:", Tprob)
      cat("\nAlpha    :", alpha)}

    if (length(nr) == 1)  LSD <- Tprob * sqrt(2 * MSerror/nr)
    statistics<-data.frame(Chisq=H,Df=ntr-1,p.chisq=p.chisq)
    if ( group & length(nr) == 1 & console) cat("\nMinimum Significant Difference:",LSD,"\n")
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of treatment differences and alpha level.\n")
    if ( length(nr) == 1) statistics<-data.frame(statistics,t.value=Tprob,MSD=LSD)

    comb <- utils::combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    LCL <- dif
    UCL <- dif
    pvalue <- dif
    sdtdif <- dif
    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      dif[k] <- means[i, 2] - means[j, 2]
      # S * ((N - 1 - H)/(N - ntr))
      sdtdif[k] <- sqrt(MSerror * (1/means[i,3] + 1/means[j, 3]))
      pvalue[k] <- 2*(1 - pt(abs(dif[k])/sdtdif[k],DFerror))
    }
    if (p.adj != "none") pvalue <- p.adjust(pvalue, p.adj)
    pvalue <- round(pvalue,4)
    sig <- rep(" ", nn)
    for (k in 1:nn) {
      if (pvalue[k] <= 0.001)
        sig[k] <- "***"
      else if (pvalue[k] <= 0.01)
        sig[k] <- "**"
      else if (pvalue[k] <= 0.05)
        sig[k] <- "*"
      else if (pvalue[k] <= 0.1)
        sig[k] <- "."
    }
    tr.i <- means[comb[1, ], 1]
    tr.j <- means[comb[2, ], 1]
    LCL <- dif - Tprob * sdtdif
    UCL <- dif + Tprob * sdtdif
    comparison <- data.frame(Difference = dif, pvalue = pvalue, "Signif."=sig, LCL, UCL)
    if (p.adj !="bonferroni" & p.adj !="none"){
      comparison<-comparison[,1:3]
      statistics<-data.frame(Chisq=H,p.chisq=p.chisq)
    }
    rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
    if (!group) {
      groups<-NULL
      if(console){
        cat("\nComparison between treatments mean of the ranks.\n\n")
        print(comparison)
      }
    }
    if (group) {
      comparison=NULL
      # Matriz de probabilidades
      Q<-matrix(1,ncol=ntr,nrow=ntr)
      p<-pvalue
      k<-0
      for(i in 1:(ntr-1)){
        for(j in (i+1):ntr){
          k<-k+1
          Q[i,j]<-p[k]
          Q[j,i]<-p[k]
        }
      }
      groups <- orderPvalue(means[, 1], means[, 2],alpha, Q, console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nTreatments with the same letter are not significantly different.\n\n")
        print(groups)
      }
    }
    ranks=means
    Means<-data.frame(rank=ranks[,2],Means)
    Means<-Means[,c(2,1,3:9)]
    parameters<-data.frame(test="Kruskal-Wallis",p.ajusted=p.adj,name.t=name.t,ntr = ntr,alpha=alpha)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    output<-list(statistics=statistics,parameters=parameters,
                 means=Means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }

#=====================================================================
# friedman
#=====================================================================

friedman <-
  function(judge,trt,evaluation,alpha=0.05,group=TRUE,main=NULL,console=FALSE){
    name.x <- paste(deparse(substitute(judge)))
    name.y <- paste(deparse(substitute(evaluation)))
    name.t <- paste(deparse(substitute(trt)))
    name.j <- paste(deparse(substitute(judge)))
    if(is.null(main))main<-paste(name.y,"~", name.j,"+",name.t)
    datos <- data.frame(judge, trt, evaluation)
    matriz <- by(datos[,3], datos[,1:2], function(x) mean(x,na.rm=TRUE))
    matriz <-as.data.frame(matriz[,])
    #matriz <-as.matrix(evaluation)
    name<-as.character(colnames(matriz))
    ntr <-length(name)
    m<-dim(matriz)
    v<-array(0,m)
    for (i in 1:m[1]){
      v[i,]<-rank(matriz[i,])
    }
    vv<-as.numeric(v)
    junto <- data.frame(evaluation, trt)
    medians<-tapply.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- tapply.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    Means <- tapply.stat(junto[,1],junto[,2],stat="mean")  # change
    sds <-   tapply.stat(junto[,1],junto[,2],stat="sd")    # change
    nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
    nr<-unique(nn[,2])
    s<-array(0,m[2])
    # Suma de rangos por tratamiento
    for (j in 1:m[2]){
      s[j]<-sum(v[,j])
    }
    Means<-data.frame(Means,std=sds[,2],r=nn[,2],medians)
    names(Means)[1:2]<-c(name.t,name.y)
    means<-Means[,c(1:2,4)]
    rownames(Means)<-Means[,1]
    Means<-Means[,-1]
    means[,2]<-s
    # row.names(means)<-means[,1]
    rs<-array(0,m[2])
    rs<-s-m[1]*(m[2]+1)/2
    T1<-12*t(rs)%*%rs/(m[1]*m[2]*(m[2]+1))
    T2<-(m[1]-1)*T1/(m[1]*(m[2]-1)-T1)
    # Impresion de resultados
    if(console){
      cat("\nStudy:",main,"\n\n")
      cat(paste(name.t,",",sep="")," Sum of the ranks\n\n")
      print(data.frame(row.names = means[,1], means[,-1]))
      cat("\nFriedman's Test")
      cat("\n===============")
    }
    A1<-0
    for (i in 1:m[1]) A1 <- A1 + t(v[i,])%*%v[i,]
    DFerror <-(m[1]-1)*(m[2]-1)
    Tprob<-qt(1-alpha/2,DFerror)
    #
    LSD<-as.numeric(Tprob*sqrt(2*(m[1]*A1-t(s)%*%s)/DFerror))

    C1 <-m[1]*m[2]*(m[2]+1)^2/4
    T1.aj <-(m[2]-1)*(t(s)%*%s-m[1]*C1)/(A1-C1)
    T2.aj <-(m[1]-1)*T1.aj/(m[1]*(m[2]-1)-T1.aj)
    p.value<-1-pchisq(T1.aj,m[2]-1)
    p.noadj<-1-pchisq(T1,m[2]-1)
    PF<-1-pf(T2.aj, ntr-1, (ntr-1)*(nr-1) )
    if(console){
      cat("\nAdjusted for ties")
      cat("\nCritical Value:",T1.aj)
      cat("\nP.Value Chisq:",p.value)
      cat("\nF Value:",T2.aj)
      cat("\nP.Value F:",PF,"\n")
      cat("\nPost Hoc Analysis\n")
    }
    #...............
    statistics<-data.frame(Chisq=T1.aj,Df=ntr-1,p.chisq=p.value,F=T2.aj,DFerror=DFerror,p.F=PF,t.value=Tprob,LSD)
    if ( group & length(nr) == 1 & console){
      cat("\nAlpha:",alpha,"; DF Error:",DFerror)
      cat("\nt-Student:",Tprob)
      cat("\nLSD:", LSD,"\n")
    }
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of treatment differences and alpha level(",alpha,")\n")
    if ( length(nr) != 1) statistics<-data.frame(Chisq=T1.aj,Df=ntr-1,p.chisq=p.value,F=T2.aj,DFerror=DFerror,p.F=PF)
    comb <-utils::combn(ntr,2)
    nn<-ncol(comb)
    dif<-rep(0,nn)
    pvalue<-rep(0,nn)
    LCL<-dif
    UCL<-dif
    sig<-NULL
    LSD<-rep(0,nn)
    stat<-rep("ns",nn)
    for (k in 1:nn) {
      i<-comb[1,k]
      j<-comb[2,k]
      dif[k]<-s[comb[1,k]]-s[comb[2,k]]
      sdtdif<- sqrt(2*(m[1]*A1-t(s)%*%s)/DFerror)
      pvalue[k]<- round(2*(1-pt(abs(dif[k])/sdtdif,DFerror)),4)
      LSD[k]<-round(Tprob*sdtdif,2)
      LCL[k] <- dif[k] - LSD[k]
      UCL[k] <- dif[k] + LSD[k]
      sig[k]<-" "
      if (pvalue[k] <= 0.001) sig[k]<-"***"
      else  if (pvalue[k] <= 0.01) sig[k]<-"**"
      else  if (pvalue[k] <= 0.05) sig[k]<-"*"
      else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
    if(!group){
      tr.i <- means[comb[1, ],1]
      tr.j <- means[comb[2, ],1]
      comparison<-data.frame("difference" = dif, pvalue=pvalue,"signif."=sig,LCL,UCL)
      rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
      if(console){cat("\nComparison between treatments\nSum of the ranks\n\n")
        print(comparison)}
      groups=NULL
    }
    if (group) {
      # Matriz de probabilidades
      Q<-matrix(1,ncol=ntr,nrow=ntr)
      p<-pvalue
      k<-0
      for(i in 1:(ntr-1)){
        for(j in (i+1):ntr){
          k<-k+1
          Q[i,j]<-p[k]
          Q[j,i]<-p[k]
        }
      }
      groups <- orderPvalue(means[, 1], means[, 2],alpha, Q,console)
      names(groups)[1]<-"Sum of ranks"
      if(console) {
        cat("\nTreatments with the same letter are not significantly different.\n\n")
        print(groups)
      }
      comparison<-NULL
    }
    parameters<-data.frame(test="Friedman",name.t=name.t,ntr = ntr,alpha=alpha)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    Means<-data.frame(rankSum=means[,2],Means)
    Means<-Means[,c(2,1,3:9)]
    output<-list(statistics=statistics,parameters=parameters,
                 means=Means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }

#=====================================================================
# LSD
#=====================================================================

LSD.test = function(y,
                    trt,
                    DFerror,
                    MSerror,
                    alpha = 0.05,
                    p.adj = c("none","holm","hommel","hochberg", "bonferroni", "BH", "BY", "fdr"),
                    group = TRUE,
                    main = NULL,
                    console=FALSE){
  p.adj <- match.arg(p.adj)
  clase <- c("aov", "lm")
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  if(is.null(main))main<-paste(name.y,"~", name.t)
  if ("aov" %in% class(y) | "lm" %in% class(y)) {
    if(is.null(main))main<-y$call
    A <- y$model
    DFerror <- df.residual(y)
    MSerror <- deviance(y)/DFerror
    y <- A[, 1]
    ipch <- pmatch(trt, names(A))
    nipch<- length(ipch)
    for(i in 1:nipch){
      if (is.na(ipch[i]))
        return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))}
    name.t<- names(A)[ipch][1]
    trt <- A[, ipch]
    if (nipch > 1){
      trt <- A[, ipch[1]]
      for(i in 2:nipch){
        name.t <- paste(name.t,names(A)[ipch][i],sep=":")
        trt <- paste(trt,A[,ipch[i]],sep=":")
      }}
    name.y <- names(A)[1]}
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  Mean<-mean(junto[,1])
  CV<-sqrt(MSerror)*100/Mean
  medians<-tapply.stat(junto[,1],junto[,2],stat="median")
  for(i in c(1,5,2:4)) {
    x <- tapply.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
    medians<-cbind(medians,x[,2])}
  medians<-medians[,3:7]
  names(medians)<-c("Min","Max","Q25","Q50","Q75")
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  std.err <- sqrt(MSerror)/sqrt(nn[, 2]) # change sds[,2]
  Tprob <- qt(1 - alpha/2, DFerror)
  LCL <- means[, 2] - Tprob * std.err
  UCL <- means[, 2] + Tprob * std.err
  means <- data.frame(means, std=sds[,2], r = nn[, 2], LCL, UCL,medians)
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  nk <- choose(ntr, 2)
  if (p.adj != "none") {
    a <- 1e-06
    b <- 1
    for (i in 1:100) {
      x <- (b + a)/2
      xr <- rep(x, nk)
      d <- p.adjust(xr, p.adj)[1] - alpha
      ar <- rep(a, nk)
      fa <- p.adjust(ar, p.adj)[1] - alpha
      if (d * fa < 0)
        b <- x
      if (d * fa > 0)
        a <- x}
    Tprob <- qt(1 - x/2, DFerror)
  }
  nr <- unique(nn[, 2])
  if(console){
    cat("\nStudy:", main)
    if(console)cat("\n\nLSD t Test for", name.y, "\n")
    if (p.adj != "none")cat("P value adjustment method:", p.adj, "\n")
    cat("\nMean Square Error: ", MSerror, "\n\n")
    cat(paste(name.t, ",", sep = ""), " means and individual (",
        (1 - alpha) * 100, "%) CI\n\n")
    print(data.frame(row.names = means[, 1], means[, 2:8]))
    cat("\nAlpha:", alpha, "; DF Error:", DFerror)
    cat("\nCritical Value of t:", Tprob, "\n")}
  statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
  if (length(nr) == 1)  LSD <- Tprob * sqrt(2 * MSerror/nr)
  if ( group & length(nr) == 1 & console) {
    if(p.adj=="none") cat("\nleast Significant Difference:",LSD,"\n")
    else cat("\nMinimum Significant Difference:",LSD,"\n")}
  if ( group & length(nr) != 1 & console)
    cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")

  if ( length(nr) == 1 & p.adj=="none") statistics<-data.frame(statistics, t.value=Tprob,LSD=LSD)
  if ( length(nr) == 1 & p.adj!="none") statistics<-data.frame(statistics, t.value=Tprob,MSD=LSD)
  LSD=" "
  comb <- utils::combn(ntr, 2)
  nn <- ncol(comb)
  dif <- rep(0, nn)
  pvalue <- dif
  sdtdif <- dif
  sig <- rep(" ", nn)
  for (k in 1:nn) {
    i <- comb[1, k]
    j <- comb[2, k]
    dif[k] <-means[i, 2] - means[j, 2]
    sdtdif[k] <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,4]))
    pvalue[k] <- 2 * (1 - pt(abs(dif[k])/sdtdif[k], DFerror))}
  if (p.adj != "none")
    pvalue <- p.adjust(pvalue, p.adj)
  pvalue <- round(pvalue,4)
  for (k in 1:nn) {
    if (pvalue[k] <= 0.001)
      sig[k] <- "***"
    else if (pvalue[k] <= 0.01)
      sig[k] <- "**"
    else if (pvalue[k] <= 0.05)
      sig[k] <- "*"
    else if (pvalue[k] <= 0.1)
      sig[k] <- "."}
  tr.i <- means[comb[1, ], 1]
  tr.j <- means[comb[2, ], 1]
  LCL <- dif - Tprob * sdtdif
  UCL <- dif + Tprob * sdtdif
  comparison <- data.frame(difference = dif, pvalue = pvalue, "signif."=sig, LCL, UCL)
  if (p.adj !="bonferroni" & p.adj !="none"){
    comparison<-comparison[,1:3]
    #    statistics<-statistics[,1:4]
  }
  rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
  if (!group) {
    if(console){
      cat("\nComparison between treatments means\n\n")
      print(comparison)
    }
    groups <- NULL
  }
  if (group) {
    comparison=NULL
    # Matriz de probabilidades
    Q<-matrix(1,ncol=ntr,nrow=ntr)
    p<-pvalue
    k<-0
    for(i in 1:(ntr-1)){
      for(j in (i+1):ntr){
        k<-k+1
        Q[i,j]<-p[k]
        Q[j,i]<-p[k]}
    }
    groups <- orderPvalue(means[, 1], means[, 2],alpha, Q,console)
    names(groups)[1]<-name.y
    if(console) {
      cat("\nTreatments with the same letter are not significantly different.\n\n")
      print(groups)
    }
  }
  parameters<-data.frame(test="Fisher-LSD",p.ajusted=p.adj,name.t=name.t,ntr = ntr,alpha=alpha)
  rownames(parameters)<-" "
  rownames(statistics)<-" "
  rownames(means)<-means[,1]
  means<-means[,-1]
  output<-list(statistics=statistics,parameters=parameters,
               means=means,comparison=comparison,groups=groups)
  class(output)<-"group"
  invisible(output)
}


#=====================================================================
# design.crd
#=====================================================================
design.crd <-
  function(trt,r,serie=2,seed=0,kinds="Super-Duper",randomization=TRUE)
  {
    number<-0
    if(serie>0) number<-10^serie
    junto<-data.frame(trt,r)
    junto<-junto[order(junto[,1]),]
    TR<-as.character(junto[,1])
    r<-as.numeric(junto[,2])
    y <- rep(TR[1], r[1])
    tr <- length(TR)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    parameters<-list(design="crd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
    for (i in 2:tr) y <- c(y, rep(TR[i], r[i]))
    trat<-y
    if(randomization)trat <- sample(y, length(y), replace = FALSE)
    plots <- number+1:length(trat)
    dca<-data.frame(plots, trat)
    dca[,1]<-as.numeric(dca[,1])
    xx<-dca[order(dca[,2],dca[,1]),]
    r1<-seq(1,r[1])
    for (i in 2:length(r)) {
      r1<-c(r1,seq(1,r[i]))
    }
    yy<-data.frame(plots=xx[,1],r=r1,xx[,2])
    book<-yy[order(yy[,1]),]
    rownames(book)<-rownames(yy)
    names(book)[3]<-c(paste(deparse(substitute(trt))))
    outdesign<-list(parameters=parameters,book=book)
    return(outdesign)
  }

#=====================================================================
# design.rcbd
#=====================================================================

design.rcbd <-
  function (trt, r,serie=2,seed=0,kinds="Super-Duper",first=TRUE,continue=FALSE,randomization=TRUE )
  {
    number<-10
    if(serie>0) number<-10^serie
    ntr <- length(trt)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    parameters<-list(design="rcbd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
    mtr <-trt
    if(randomization)mtr <- sample(trt, ntr, replace = FALSE)
    block <- c(rep(1, ntr))
    for (y in 2:r) {
      block <- c(block, rep(y, ntr))
      if(randomization)mtr <- c(mtr, sample(trt, ntr, replace = FALSE))
    }
    if(randomization){
      if(!first) mtr[1:ntr]<-trt
    }
    plots <- block*number+(1:ntr)
    book <- data.frame(plots, block = as.factor(block), trt = as.factor(mtr))
    names(book)[3] <- c(paste(deparse(substitute(trt))))
    names(book)[3]<-c(paste(deparse(substitute(trt))))
    if(continue){
      start0<-10^serie
      if(serie==0) start0<-0
      book$plots<-start0+1:nrow(book)
    }
    outdesign<-list(parameters=parameters,sketch=matrix(book[,3], byrow = TRUE, ncol = ntr),book=book)
    return(outdesign)
  }

#=====================================================================
# design.lsd
#=====================================================================


design.lsd <-
  function (trt,serie=2,seed=0,kinds="Super-Duper",first=TRUE,randomization=TRUE)
  {
    number<-10
    if(serie>0) number<-10^serie
    r <- length(trt)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    parameters<-list(design="lsd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
    a <- 1:(r * r)
    dim(a) <- c(r, r)
    for (i in 1:r) {
      for (j in 1:r) {
        k <- i + j - 1
        if (k > r)
          k <- i + j - r - 1
        a[i, j] <- k
      }
    }
    m<-2:r
    if(randomization)m<-sample(2:r,r-1)
    a<-a[,c(1,m)]
    if(randomization){
      if (first) {
        m<-sample(1:r,r)
        a<-a[m,]
      }}
    trat<-trt[a]
    columna <- rep(gl(r, 1), r)
    fila <- gl(r, r)
    fila <- as.character(fila)
    fila <- as.numeric(fila)
    plots <- fila*number+(1:r)
    book <- data.frame(plots, row = as.factor(fila), col = as.factor(columna),
                       trat = as.factor(trat))
    names(book)[4] <- c(paste(deparse(substitute(trt))))
    outdesign<-list(parameters=parameters,sketch=matrix(book[,4], byrow = TRUE, ncol = r),book=book)
    return(outdesign)
  }

#=====================================================================
# design.split
#=====================================================================

design.split <-
  function (trt1, trt2,r=NULL, design=c("rcbd","crd","lsd"),serie = 2, seed = 0, kinds = "Super-Duper",
            first=TRUE,randomization=TRUE )
  {
    n1<-length(trt1)
    n2<-length(trt2)
    if (seed == 0) {
      genera<-runif(1)
      seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    design <- match.arg(design)
    number<-10^serie +1
    if (design == "crd") {
      plan<-design.crd(trt1,r,serie, seed, kinds,randomization)
      k<-3
    }
    if (design == "rcbd"){
      plan<-design.rcbd(trt1,r,serie, seed, kinds, first,randomization)
      k<-3
    }
    if (design == "lsd") {
      plan<-design.lsd(trt1,serie, seed, kinds, first,randomization)
      r<-n1
      k<-4
    }
    book<-plan$book
    parameters<-plan$parameters
    names(parameters)[2]<-"trt1"
    parameters$applied<-parameters$design
    parameters$design<-"split"
    parameters$trt2<-trt2
    j<-0
    B<-list()
    for(i in c(1,7,2,8,3:6)){
      j<-j+1
      B[[j]]<-parameters[[i]]
      names(B)[j]<-names(parameters)[i]
    }
    nplot<-nrow(book)
    d<-NULL
    if(randomization){
      for(i in 1:nplot)d<-rbind(d,sample(trt2,n2))
    }
    else{
      d<-rbind(d,trt2[1:n2])
    }
    aa<-data.frame(book,trt2=d[,1])
    for(j in 2:n2) aa<-rbind(aa,data.frame(book,trt2=d[,j]))
    aa<-aa[order(aa[,1]),]
    splots<-rep(gl(n2,1),nplot)
    book <- data.frame(plots=aa[,1],splots,aa[,-1])
    rownames(book)<-1:(nrow(book))
    names(book)[k+1] <- c(paste(deparse(substitute(trt1))))
    names(book)[k+2] <- c(paste(deparse(substitute(trt2))))
    outdesign<-list(parameters=B,book=book)
    return(outdesign)
  }

#=====================================================================
# design.ab
#=====================================================================

design.ab <-
  function(trt, r=NULL,serie=2,design=c("rcbd","crd","lsd"),seed=0,kinds="Super-Duper",
           first=TRUE,randomization=TRUE ){
    design <- match.arg(design)
    if( design=="rcbd" | design=="crd") posicion <- 3
    else posicion <- 4
    serie<-serie; seed<-seed; kinds<-kinds; first<-first;

    # Process to trt to factorial
    ntr<-length(trt)
    fact<-NULL
    tr0<-1:trt[1]
    k<-0
    a<-trt[1];b<-trt[2]
    for(i in 1:a){
      for(j in 1:b){
        k<-k+1
        fact[k]<-paste(tr0[i],j)
      }
    }

    if(ntr >2) {
      for(m in 3:ntr){
        k<-0
        tr0<-fact
        fact<-NULL
        a<-a*b
        b<-trt[m]
        for(i in 1:a){
          for(j in 1:b){
            k<-k+1
            fact[k]<-paste(tr0[i],j)
          }
        }
      }
    }
    #------------------------------
    if(design=="rcbd")plan<-design.rcbd(trt=fact, r, serie, seed, kinds, first,randomization )
    if(design=="crd")plan<-design.crd(trt=fact, r, serie, seed, kinds,randomization)
    if(design=="lsd")plan<-design.lsd(trt=fact, serie, seed, kinds, first,randomization )
    parameters<-plan$parameters
    parameters$applied<-parameters$design
    parameters$design<-"factorial"
    plan<-plan$book
    trt<-as.character(plan[,posicion])
    nplan<-nrow(plan)
    A<-rep(" ",nplan*ntr)
    dim(A)<-c(nplan,ntr)
    colnames(A)<-LETTERS[1:ntr]

    for(i in 1:nplan) {
      A[i,]<-unlist(strsplit(trt[i], " "))
    }
    A<-as.data.frame(A)
    book<-data.frame(plan[,1:(posicion-1)],A)
    outdesign<-list(parameters=parameters,book=book)
    return(outdesign)
  }

##=========================================================
## Sk
##=========================================================
sk_triple<-function(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
{
  sk <- function(medias,s2,dfr,prob){
    bo <- 0
    si2 <- s2
    defr <- dfr
    parou <- 1
    np <- length(medias) - 1
    for (i in 1:np){
      g1 <- medias[1:i]
      g2 <- medias[(i+1):length(medias)]
      B0 <- sum(g1)^2/length(g1) + sum(g2)^2/length(g2) - (sum(g1) + sum(g2))^2/length(c(g1,g2))
      if (B0 > bo)
      {bo <- B0
        parou <- i}
    }

    g1 <- medias[1:parou]
    g2 <- medias[(parou+1):length(medias)]
    teste <- c(g1,g2)

    sigm2 <- (sum(teste^2) - sum(teste)^2/length(teste) + defr*si2)/(length(teste) + defr)
    lamb <- pi*bo/(2*sigm2*(pi-2))
    v0 <- length(teste)/(pi-2)
    p <- pchisq(lamb,v0,lower.tail = FALSE)

    if (p < prob) {
      for (i in 1:length(g1)){
        cat(names(g1[i]),"\n",file="skresult",append=TRUE)}
      cat("*","\n",file="skresult",append=TRUE)}

    if (length(g1)>1){sk(g1,s2,dfr,prob)}
    if (length(g2)>1){sk(g2,s2,dfr,prob)}
  }

  medias <- sort(tapply(y,trt,mean),decreasing=TRUE)
  dfr <- DFerror

  rep <- tapply(y,trt,length)
  s0 <- MSerror <-SSerror/DFerror
  s2 <- s0/rep[1]
  prob <- alpha
  sk(medias,s2,dfr,prob)
  f <- names(medias)
  names(medias) <- 1:length(medias)
  resultado <- data.frame("r"=0,"f"=f,"m"=medias)
  if (file.exists("skresult") == FALSE) {stop} else{
    xx <- read.table("skresult")
    file.remove("skresult")
    x <- xx[[1]]
    x <- as.vector(x)
    z <- 1

    for (j in 1:length(x)){
      if (x[j] == "*")	{z <- z+1}
      for (i in 1:length(resultado$f)){
        if (resultado$f[i]==x[j]){
          resultado$r[i] <- z;}
      }
    }

  }
  letras<-letters
  if(length(resultado$r)>26) {
    l<-floor(length(resultado$r)/26)
    for(i in 1:l) letras<-c(letras,paste(letters,i,sep=''))
  }
  res <- 1
  for (i in 1:(length(resultado$r)-1))
  {
    if (resultado$r[i] != resultado$r[i+1]){
      resultado$r[i] <- letras[res]
      res <- res+1
      if (i == (length(resultado$r)-1)){
        resultado$r[i+1] <- letras[res]}
    }
    else{
      resultado$r[i] <- letras[res]
      if (i == (length(resultado$r)-1)){
        resultado$r[i+1] <- letras[res]
      }
    }
  }
  names(resultado) <- c("Groups","Tratamentos","Medias")
  resultado1=resultado[,c(3,1)]
  rownames(resultado1)=resultado$Tratamentos
  final=list(resultado1)[[1]]
}
