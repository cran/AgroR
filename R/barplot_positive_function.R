#' Graph: Positive barplot
#'
#' @description Column chart with two variables that assume a positive response and represented by opposite sides, such as dry mass of the area and dry mass of the root
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @param a Object of DIC, DBC or DQL functions
#' @param b Object of DIC, DBC or DQL functions
#' @param ylab Y axis names
#' @param var_name Name of the variable
#' @param fill_color Bar fill color
#' @param legend.title Legend title
#' @param width.bar Width error bar
#' @param width.col Width Column
#' @seealso \link{radargraph}, \link{sk_graph}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @return The function returns a column chart with two positive sides
#' @note When there is only an effect of the isolated factor in the case of factorial or subdivided plots, it is possible to use the barplot_positive function.
#' @export
#' @examples
#' data("passiflora")
#' attach(passiflora)
#' a=with(passiflora, DBC(trat, bloco, MSPA))
#' b=with(passiflora, DBC(trat, bloco, MSR))
#' barplot_positive(a, b, var_name = c("DMAP","DRM"), ylab = "Dry root (g)")
#'
#' a=with(passiflora, DIC(trat, MSPA,test = "noparametric"))
#' b=with(passiflora, DIC(trat, MSR))
#' barplot_positive(a, b, var_name = c("DMAP","DRM"), ylab = "Dry root (g)")

barplot_positive=function(a,
                          b,
                          ylab="Response",
                          var_name=c("Var1","Var2"),
                          legend.title="Variable",
                          fill_color=c("darkgreen",
                                       "brown"),
                          width.col=0.9,
                          width.bar=0.2){
  requireNamespace("ggplot2")
  if(a[[1]]$plot$test=="parametric" & b[[1]]$plot$test=="parametric"){
  dataA=a[[1]]$data
  dataB=b[[1]]$data
  dataB$media=dataB$media*-1
  dataB$desvio=dataB$desvio*-1
  dataB$limite=dataB$limite*-1.1
  dataA$limite=dataA$limite*1.1
  if(colnames(dataA)[3]=="respO"){dataA=dataA[,-3]}
  if(colnames(dataB)[3]=="respO"){dataB=dataB[,-3]}
  data=rbind(dataA,
             dataB)
  data$vari=rep(c("Var1","Var2"),
                e=length(rownames(dataA)))
  data$vari=as.factor(data$vari)
  levels(data$vari)=var_name
  trats=data$trats
  media=data$media
  vari=data$vari
  desvio=data$desvio
  limite=data$limite
  letra=data$letra
  graph=ggplot(data,aes(x=trats,
                  y=media,
                  fill=vari))+
    geom_col(color="black",width = width.col)+
    geom_errorbar(aes(ymin=media-desvio,
                      ymax=media+desvio),
                  width=width.bar)+
    scale_y_continuous(breaks = pretty(media*1.5),
                       labels = abs(pretty(media*1.5)))+
    theme_classic()+xlab("")+ylab(ylab)+
    geom_text(aes(y=limite,
                  label=letra),family = a[[1]]$plot$family)+
    scale_fill_manual(values=fill_color,
                      labels = c(var_name))+
    geom_hline(yintercept=0)+
    labs(fill=legend.title)+
    theme(axis.text = element_text(size=a[[1]]$theme$axis.text$size,
                                   color = a[[1]]$theme$axis.text$colour,family = a[[1]]$plot$family),
          axis.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
          legend.text = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
          legend.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family))}
  if(a[[1]]$plot$test=="noparametric" & b[[1]]$plot$test=="noparametric"){
    dataA=a[[1]]$data
    dataB=b[[1]]$data
    dataB$media=dataB$media*-1
    dataB$std=dataB$std*-1
    dataB$limite=dataB$limite*-1.1
    dataA$limite=dataA$limite*1.1
    # if(colnames(dataA)[3]=="respO"){dataA=dataA[,-3]}
    # if(colnames(dataB)[3]=="respO"){dataB=dataB[,-3]}
    data=rbind(dataA,
               dataB)
    data$vari=rep(c("Var1","Var2"),
                  e=length(rownames(dataA)))
    data$vari=as.factor(data$vari)
    levels(data$vari)=var_name
    trats=data$trats
    media=data$media
    vari=data$vari
    desvio=data$std
    limite=data$limite
    letra=data$letra
    graph=ggplot(data,aes(x=trats,
                    y=media,
                    fill=vari))+
      geom_col(color="black",width = width.col)+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),
                    width=width.bar)+
      scale_y_continuous(breaks = pretty(media*1.5),
                         labels = abs(pretty(media*1.5)))+
      theme_classic()+xlab("")+ylab(ylab)+
      geom_text(aes(y=limite,
                    label=letra),family = a[[1]]$plot$family)+
      scale_fill_manual(values=fill_color,
                        labels = c(var_name))+
      geom_hline(yintercept=0)+
      labs(fill=legend.title)+
      theme(axis.text = element_text(size=a[[1]]$theme$axis.text$size,
                                     color = a[[1]]$theme$axis.text$colour,family = a[[1]]$plot$family),
            axis.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
            legend.text = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
            legend.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family))}
  if(a[[1]]$plot$test=="parametric" & b[[1]]$plot$test=="noparametric"){
    dataA=a[[1]]$data
    dataB=b[[1]]$data
    dataB$media=dataB$media*-1
    dataB$desvio=dataB$std*-1
    dataB$limite=dataB$limite*-1.1
    dataA$limite=dataA$limite*1.1
    dataB=dataB[,c(12,13,14,15,16)]
    if(a[[1]]$plot$transf == 1){dataA=dataA[,c(5,3,6,7,4)]}
    if(a[[1]]$plot$transf != 1){dataA=dataA[,c(6,4,7,8,5)]}
    data=rbind(dataA,
               dataB)
    data$vari=rep(c("Var1","Var2"),
                  e=length(rownames(dataA)))
    data$vari=as.factor(data$vari)
    levels(data$vari)=var_name
    trats=data$trats
    media=data$media
    vari=data$vari
    desvio=data$desvio
    limite=data$limite
    letra=data$letra
    graph=ggplot(data,aes(x=trats,
                          y=media,
                          fill=vari))+
      geom_col(color="black",width = width.col)+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),
                    width=width.bar)+
      scale_y_continuous(breaks = pretty(media*1.5),
                         labels = abs(pretty(media*1.5)))+
      theme_classic()+xlab("")+ylab(ylab)+
      geom_text(aes(y=limite,
                    label=letra),family = a[[1]]$plot$family)+
      scale_fill_manual(values=fill_color,
                        labels = c(var_name))+
      geom_hline(yintercept=0)+
      labs(fill=legend.title)+
      theme(axis.text = element_text(size=a[[1]]$theme$axis.text$size,
                                     color = a[[1]]$theme$axis.text$colour,family = a[[1]]$plot$family),
            axis.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
            legend.text = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
            legend.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family))
    }
  if(a[[1]]$plot$test=="noparametric" & b[[1]]$plot$test=="parametric"){
    dataA=a[[1]]$data
    dataB=b[[1]]$data
    dataB$media=dataB$media*-1
    dataB$desvio=dataB$desvio*-1
    dataB$limite=dataB$limite*-1.1
    dataA$limite=dataA$limite*1.1
    dataA=dataA[,c(12,13,14,15,3)]
    colnames(dataA)[5]="desvio"
    if(a[[1]]$plot$transf == 1){dataB=dataB[,c(5,3,6,7,4)]}
    if(a[[1]]$plot$transf != 1){dataB=dataB[,c(6,4,7,8,5)]}
    data=rbind(dataA,
               dataB)
    data$vari=rep(c("Var1","Var2"),
                  e=length(rownames(dataA)))
    data$vari=as.factor(data$vari)
    levels(data$vari)=var_name
    trats=data$trats
    media=data$media
    vari=data$vari
    desvio=data$desvio
    limite=data$limite
    letra=data$letra
    graph=ggplot(data,aes(x=trats,
                          y=media,
                          fill=vari))+
      geom_col(color="black",width = width.col)+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),
                    width=width.bar)+
      scale_y_continuous(breaks = pretty(media*1.5),
                         labels = abs(pretty(media*1.5)))+
      theme_classic()+xlab("")+ylab(ylab)+
      geom_text(aes(y=limite,
                    label=letra),family = a[[1]]$plot$family)+
      scale_fill_manual(values=fill_color,
                        labels = c(var_name))+
      geom_hline(yintercept=0)+
      labs(fill=legend.title)+
      theme(axis.text = element_text(size=a[[1]]$theme$axis.text$size,
                                     color = a[[1]]$theme$axis.text$colour,family = a[[1]]$plot$family),
            axis.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
            legend.text = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family),
            legend.title = element_text(size=a[[1]]$plot$textsize,family = a[[1]]$plot$family))
    }

  graph}
