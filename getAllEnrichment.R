getAllEnrichment <- function( 
  data = enriched.Up
  ,  goTYPE = "goCC"
  , dbgo = "GO_Cellular_Component_2017"
  , direction = "UP"
  , revigo = T
  , TF = F
  , TF2match = NULL 
  , pVal = 0.01) {
  
  output = list()

  Data <- paste(goTYPE,direction, sep="." )
  assign( Data, data[[dbgo]] )
  assign( Data, get(Data)[get(Data)$Adjusted.P.value < pVal , ] )
  assign( Data, get(Data)[order( get(Data)[,"Adjusted.P.value"] ),] )
  data_1<-get(Data)
  data_1$Overlap <- gsub("/.*","",data_1$Overlap)
  data_1 <- data_1[, c("Term","Overlap","Adjusted.P.value","Genes")]
  assign(Data, data_1)
  
  
  if( TF == T) {
  data_1<-get(Data)
  data_1$Term <- gsub("_.*" ,"", data_1$Term)
  data_1 <- subset(data_1, data_1$Term %in% toupper(TF2match)  ) 
  assign(Data, data_1)
  }
  
  File <- paste(Data, "LFC" , round(LFC,2), project, "txt", sep=".")
  if(!file.exists(File)) { write.table( get(Data), file=File, row.names = F , quote=F, sep = "\t" ) }
  output[["Enrichment"]] <- get(Data)
  
  if (revigo == T) {
    revData <- paste("revigo",goTYPE,direction, sep="." )
    assign( revData, getRevigoMatrix(get(Data) ) )
    revFile <- paste(revData, "LFC" , round(LFC,2), project, "txt", sep=".")
    if(!file.exists(revFile)) { write.table( get(revData), file=revFile, row.names = F, quote=F  ) }
    output[["revigo"]] <- get(revData)
  }
  output[["Enrichment"]]$Description <- gsub (" \\(.*", "", output[["Enrichment"]]$Term )
  output[["Enrichment"]]$Term <- gsub (".*\\(|)", "", output[["Enrichment"]]$Term )
  return(output)
}
