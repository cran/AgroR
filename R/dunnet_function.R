#' Analysis: Dunnett test
#' @export
#' @description The function performs the Dunnett test
#' @param trat Numerical or complex vector with treatments
#' @param resp Numerical vector containing the response of the experiment.
#' @param control Treatment considered control (write identical to the name in the vector)
#' @param model Experimental design (DIC, DBC or DQL)
#' @param block Numerical or complex vector with blocks
#' @param line Numerical or complex vector with lines
#' @param column Numerical or complex vector with columns
#' @param alpha.t Significance level (\emph{default} is 0.05)
#' @param label Variable label
#' @param pointsize Point size
#' @param pointshape Shape
#' @param textsize Font size
#' @param linesize Line size
#' @param labelsize Label size
#' @param errorsize Errorbar size
#' @param widthsize Width errorbar
#' @param fontfamily font family
#' @note Do not use the "-" symbol or space in treatment names
#' @return I return the Dunnett test for experiments in a completely randomized design, randomized blocks or Latin square.
#' @importFrom multcomp glht
#' @importFrom multcomp mcp
#' @examples
#'
#' #====================================================
#' # complete randomized design
#' #====================================================
#' data("pomegranate")
#' with(pomegranate,dunnett(trat=trat,resp=WL,control="T1"))
#'
#' #====================================================
#' # randomized block design in factorial double
#' #====================================================
#' library(AgroR)
#' data(cloro)
#' attach(cloro)
#' respAd=c(268, 322, 275, 350, 320)
#' a=FAT2DBC.ad(f1, f2, bloco, resp, respAd,
#'              ylab="Number of nodules",
#'              legend = "Stages",mcomp="sk")
#' data=rbind(data.frame(trat=paste(f1,f2,sep = ""),bloco=bloco,resp=resp),
#'            data.frame(trat=c("Test","Test","Test","Test","Test"),
#'                       bloco=unique(bloco),resp=respAd))
#' with(data,dunnett(trat = trat,
#'                   resp = resp,
#'                   control = "Test",
#'                   block=bloco,model = "DBC"))


dunnett=function(trat,
                 resp,
                 control,
                 model="DIC",
                 block=NA,
                 column=NA,
                 line=NA,
                 alpha.t=0.05,
                 pointsize=5,
                 pointshape=21,
                 linesize=1,
                 labelsize=4,
                 textsize=12,
                 errorsize=1,
                 widthsize=0.2,
                 label="Response",
                 fontfamily="sans"){
  trat1=factor(trat,unique(trat))
  trat=factor(trat,unique(trat))
  levels(trat1)=paste("T",1:length(levels(trat1)),sep = "")
  controle=as.character(trat1[trat==control][1])
  if(model=="DIC"){mod=aov(resp~trat1)}
  if(model=="DBC"){
    block=as.factor(block)
    mod=aov(resp~trat1+block)}
  if(model=="DQL"){
    column=as.factor(column)
    line=as.factor(line)
    mod=aov(resp~trat1+column+line)}
  requireNamespace("multcomp")
  dados=data.frame(trat1,resp)
  contras=unique(trat1)[!unique(trat1)==controle]
  a=confint(glht(mod,
               linfct = mcp(trat1=paste(contras,"-",
                                       controle,
                                       "==0",sep=""))),
          level = 1-alpha.t)
  a=summary(a)
  teste=cbind(a$confint,
        round(a$test$tstat,4),
        round(a$test$pvalues,4))
  nomes=rownames(teste)
  nomes1=t(matrix(unlist(strsplit(nomes," - ")),nrow=2))[,1]
  nomes1=factor(nomes1,unique(nomes1))
  levels(nomes1)=levels(trat)[!levels(trat)==control]
  rownames(teste)=paste(control," - ",nomes1)
  teste=data.frame(teste)
  colnames(teste)=c("Estimate","IC-lwr","IC-upr","t value","p-value")
  teste$sig=ifelse(teste$`p-value`>alpha.t,"ns",
                   ifelse(teste$`p-value`<alpha.t,"*",""))
  print(teste)
  data=data.frame(teste)
  `IC-lwr`=data$IC.lwr
  `IC-upr`=data$IC.upr
  sig=data$sig
  Estimate=data$Estimate
  graph=ggplot(data,aes(y=rownames(data),x=Estimate))+
    geom_errorbar(aes(xmin=`IC-lwr`,xmax=`IC-upr`),width=widthsize,size=errorsize)+
    geom_point(shape=pointshape,size=pointsize,color="black",fill="gray")+
    theme_classic()+
    labs(y="")+
    geom_vline(xintercept = 0,lty=2,linewidth=linesize)+
    geom_label(aes(label=paste(round(Estimate,3),
                               sig)),fill="lightyellow",size=labelsize,
               vjust=-0.5,family=fontfamily)+
    theme(axis.text = element_text(size=textsize,family = fontfamily),
          axis.title = element_text(size=textsize,family = fontfamily))
  plot(graph)
  }
