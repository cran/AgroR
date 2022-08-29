#' Graphics: Graphic for t test to compare means with a reference value
#'
#' @description Sometimes the researcher wants to test whether the treatment mean is greater than/equal to or less than a reference value. For example, I want to know if the average productivity of my treatment is higher than the average productivity of a given country. For this, this function allows comparing the means with a reference value using the t test.
#' @author Gabriel Danilo Shimizu
#' @param tonetest t.one.test object
#' @param alpha confidence level.
#'
#' @export
#' @return returns a density plot and a column plot to compare a reference value with other treatments.
#' @examples
#' library(AgroR)
#' data("pomegranate")
#' resu=tonetest(resp=pomegranate$WL, trat=pomegranate$trat, mu=2)
#' plot_tonetest(resu)

plot_tonetest=function(tonetest,alpha=0.95){
  xp=as.list(1:length(tonetest$min))
  resp=as.list(1:length(tonetest$min))
  for(i in 1:length(tonetest$min)){
    xp1=seq(tonetest$mean[i]-5*tonetest$sd[i],
            tonetest$mean[i]+5*tonetest$sd[i],length=500)
    yp=dnorm(xp1,mean = tonetest$mean[i],sd = tonetest$sd[i])
    resp[[i]]=yp
    xp[[i]]=xp1}
  names(resp)=names(tonetest$min)
  data=data.frame(xp=unlist(xp),yp=unlist(resp),trat=rep(names(resp),e=500))
  a=ggplot(data,aes(x=xp,y=yp,fill=trat))+labs(fill="",x="Response",y="Density")+
    geom_polygon(color="black",alpha=0.3)+
    theme_classic()+
    geom_vline(xintercept = tonetest$mean.ref,lty=2,size=3,alpha=0.3)+
    theme(axis.text = element_text(size=12),legend.position = c(0.9,0.8))
  media=c();trat=c()
  b=ggplot(data.frame(media=tonetest$mean,trat=names(tonetest$mean)))+
    geom_col(aes(x=trat,y=media,fill=trat),show.legend = FALSE,
             color="black")+
    theme_classic()+labs(x="",y="Response")+
    geom_hline(yintercept = tonetest$mean.ref,lty=2,size=3,alpha=0.3)+
    theme(axis.text = element_text(size=12))
  if(tonetest$conf.inf[1]!=-Inf & tonetest$conf.sup[1]!=Inf){
    b=b+geom_errorbar(aes(ymin=tonetest$conf.inf,
                          ymax=tonetest$conf.sup,x=trat),width=0.2)+
      geom_text(aes(x=trat,
                    y=tonetest$conf.sup+0.1*media,
                    label=ifelse(tonetest$p.value<1-alpha,"*","ns")),size=5)}
  if(tonetest$conf.inf[1]==-Inf){
    b=b+geom_errorbar(aes(ymin=media,
                          ymax=tonetest$conf.sup,x=trat),width=0.2)+
      geom_text(aes(x=trat,
                    y=tonetest$conf.sup+0.1*media,
                    label=ifelse(tonetest$p.value<1-alpha,"*","ns")),size=5)}

  if(tonetest$conf.sup[1]==Inf){
    b=b+geom_errorbar(aes(ymin=tonetest$conf.inf,
                          ymax=media,x=trat),width=0.2)+
      geom_text(aes(x=trat,
                    y=media+0.1*media,
                    label=ifelse(tonetest$p.value<1-alpha,"*","ns")),size=5)}
  gridExtra::grid.arrange(a,b,ncol=2)
}
