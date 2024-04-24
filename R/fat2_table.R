#' Utils: Summary of the analysis for factor arrangement with two qualitative factors.
#' @description Summarizes the output returned in the summarise_anova function in list form. The advantage is that the table, in the case of significant interaction, is returned in a format that facilitates assembly in terms of scientific publication.
#' @author Gabriel Danilo Shimizu
#' @param output Output of summarise_anova function for FAT2DIC, FAT2DIC.ad, FAT2DBC, FAT2DBC.ad, PSUBDIC and PSUBDBC design.
#' @param nf1 Number of levels of factor 1
#' @param nf2 Number of levels of factor 2
#' @param column Variable column
#' @return returns a list containing analysis output for experiments in FAT2DIC, FAT2DIC.ad, FAT2DBC, FAT2DBC.ad, PSUBDIC and PSUBDBC design.
#' @export
#' @examples
#'
#' #==============================================================
#' data(corn)
#' attach(corn)
#' a=FAT2DIC(A, B, Resp, quali=c(TRUE, TRUE))
#' output_1=summarise_anova(list(a),design="FAT2DIC",divisor = FALSE)
#' fat2_table(output_1,nf1=3,nf2=2,column=1)
#'
#' #==============================================================
#' data(cloro)
#' respAd=c(268, 322, 275, 350, 320)
#' resu=with(cloro, FAT2DIC.ad(f1, f2, bloco, resp, respAd))
#' output_2=summarise_anova(list(resu),design="FAT2DIC.ad",divisor = FALSE)
#' fat2_table(output_2,nf1=2,nf2=4,column=1)



fat2_table=function(output,nf1,nf2,column=1){
  output_f1=data.frame(output[-c((nf1+1):nrow(output)),column])
    rownames(output_f1)=rownames(output)[-c((nf1+1):nrow(output))]
    output_f2=data.frame(output[-c(1:nf1,(nf1+nf2+1):nrow(output)),column])
    rownames(output_f2)=rownames(output)[-c(1:nf1,(nf1+nf2+1):nrow(output))]
    colnames(output_f1)=colnames(output)[column]
    colnames(output_f2)=colnames(output)[column]
    output_f1f2=data.frame(matrix(output[-c(1:(nf1+nf2),(nf1+nf2+nf1*nf2+1):nrow(output)),column],ncol=nf1))
    colnames(output_f1f2)=rownames(output)[c(1:nf1)]
    rownames(output_f1f2)=rownames(output)[c((nf1+1):(nf1+nf2))]
    output_info=data.frame(output[-c(1:(nf1+nf2+nf1*nf2)),column])
    rownames(output_info)=rownames(output)[-c(1:(nf1+nf2+nf1*nf2))]
    colnames(output_info)=colnames(output)[column]
    saida=list(f1=output_f1,f2=output_f2,f1xf2=output_f1f2,"inf"=output_info)
  saida
}
