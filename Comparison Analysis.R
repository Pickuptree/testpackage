#' This function does comparison analysis based on ATAC-seq data. 
#' It assumes that the independent variable is the stimulation.
#' There are two conditions (GF and SPF).
#' Four tsv report files will be produced for the following analysis:
#' Condition GF upregulation
#' Condition GF downregulation
#' Condition SPF upregulation
#' Condition SPF downregulation
#'
#' @param metadataFileName Path to metadata file for ATAC-seq
#' @param countsFileName Path to file with read counts 
#' @param reference Stimulation (Cytokine or PBS) that will be used as the reference 
#'
#' @export
comparisonAnalysis = function(metadataFileName,countsFileName,reference){
  metadata=read.table(metadataFileName,header = T,row.names = 1,sep = '\t',check.names = F)
  Counts=read.table(countsFileName,header = T,row.names = 1,sep = '\t',check.names = F)
  
  # Iterate over all conditions
  for (condition in list("GF","SPF")){
    
    stimulantComparison.metadata=metadata[metadata$Condition==condition,]
    stimulantComparison=Counts[,rownames(stimulantComparison.metadata)]
    row_sub = apply(stimulantComparison, 1, function(row) sum(row !=0 )>=0.25*length(row))
    stimulantComparison=stimulantComparison[row_sub,]
    
    stimulantComparison=DESeqDataSetFromMatrix(stimulantComparison,colData = stimulantComparison.metadata,design = ~Stimulation)
    stimulantComparison$Condition=relevel(stimulantComparison$Stimulation,ref = reference)
    stimulantComparison=DESeq(stimulantComparison)
    stimulantComparison=results(stimulantComparison,name = resultsNames(stimulantComparison)[2])
    stimulantComparison=data.frame(stimulantComparison)
    
    stimulantComparison=stimulantComparison[!is.na(stimulantComparison$padj),]
    
    stimulantComparison.upregulated.region=stimulantComparison[stimulantComparison$padj<0.05&stimulantComparison$log2FoldChange>1,]
    stimulantComparison.downregulated.region=stimulantComparison[stimulantComparison$padj<0.05&stimulantComparison$log2FoldChange<(-1),]
    
    # Iterate over regulation changes
    i = 1
    fileNames = list(".upregulated.region.tsv",".downregulated.region.tsv")
    for (regulation in list(stimulantComparison.upregulated.region,stimulantComparison.downregulated.region)) {
      if (length(rownames(regulation))==0){
        print(paste("no significant regulation change observed for",condition, sep=" "))
      }
      else{
        regulation.annotation=stringr::str_split_fixed(rownames(regulation),'\\.',3)
        regulation.annotation=data.frame(regulation.annotation)
        colnames(regulation.annotation)=c('Chr','Start','End')
        regulation.annotation=makeGRangesFromDataFrame(regulation.annotation)
        regulation.annotation=annotatePeak(regulation.annotation,TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                           annoDb = 'org.Mm.eg.db',level = 'gene')
        fileName = paste(condition, fileNames[i], sep="")
        write.table(data.frame(regulation.annotation),fileName,quote = F,sep = '\t')
        i = i+1
      }
    }
  }
}
