function (data = pertTF.UP, direction = "UP", genesUp = uniqueUp, genesDown = uniqueDown
                         , type = "pertTF" ) {
  
  if ( type == "pertTF" ) {
    pertAllowed <-  c("OE", "KO", "KD", "SIRNA", "ACTIVATION", "UNK", "MUT", "INACTIVATION", "DEPLETION", "DEFICIENCY", "SHRNA",       "ABLATION", "INHIBITION", "KNOCKIN", "SILENCING","EWSFUSION","KNOCKOUT","DELETION","INDUCED","INDUCTION","TRUNCATED","MUTATION", "DOMNEG", "MUTANT","HYPERACTIVE") 
    
    pertUp <- c("OE", "ACTIVATION", "KNOCKIN", "INDUCED", "INDUCTION", "HYPERACTIVE", "EWSFUSION")
    
    pertDown <- c("KO", "KD", "SIRNA", "MUT", "INACTIVATION", "DEPLETION", "DEFICIENCY", "SHRNA", "ABLATION", "INHIBITION", "SILENCING","KNOCKOUT","DELETION","TRUNCATED","MUTATION","DOMNEG", "MUTANT") 
  }
  else {
    pertDown <- c("knockdown", "knockout","druginhibition","mutant","defectivemutant")
    pertUp <- c("overexpression","activemutant","drugactivation")
    pertAllowed <- c(pertDown,pertUp)
  }
  
  
  ## LookUp perturbations and create column defining it's type (Search and add to new column)
  data$pertType <- findReplaceWildcard(findIn = data$Term, replaceWith = pertAllowed, start = "^.*_", end="_.*$")
  
  ##Separate Upregulated perturbations and define direction in new column
  data.pert.up <- subset(data, (data$pertType %in% pertUp ) )  ## Up Pertb
  data.pert.up$pertDir <- "UP"
  
  data.pert.down <- subset(data, (data$pertType %in% pertDown ) ) ##Down pertb
  data.pert.down$pertDir <- "DOWN"
  
  data <- rbind(data.pert.up,data.pert.down) ## Up and Down Pertb merged
  
  ## TF Pert: Resultant gene expression given in Term name
  if( type ==  "pertTF" ) {
    data$geneDir <- gsub(".*_","",data$Term)
  } 
  ## Kinases: We already know from database direction of resultant gene expression
  else if ( type == "kinases" ) {
    data$geneDir <- direction
  }
  
  ## Completely Isolate Kinase Symbol
  data$Term <- gsub("_.*","",data$Term)
  
  ## Special exception for UNK, where you upregulate and knockout
  data$pertType <- ifelse(data$pertType=="UNK", data$geneDir , data$pertType )
  
  ## Sort by AdjPVal
  data <- data[order( data$Adjusted.P.value) ,]
  #data <- data[!duplicated(data$Term),]
  
  ## If Term if Upregulated, and genes get upregulated(or vice versa), it's direct correlation
  if( direction == "UP" ){
  corrDirect <- data[ data$pertDir=="UP" & data$geneDir == direction ,]
  corrDirect <- subset(corrDirect, corrDirect$Term %in% toupper(genesUp)  ) 
  corrInverse <- data[ data$pertDir=="DOWN" & data$geneDir == direction,]
  corrInverse <- subset(corrInverse, corrInverse$Term %in% toupper(genesDown)  ) 
  
  ## If Term if Upregulated, and genes get Downregulated (or vice versa), it's direct correlation
  } else if ( direction == "DOWN" ) {
    corrDirect <- data[ data$pertDir=="DOWN" & data$geneDir == direction ,]
    corrDirect <- subset(corrDirect, corrDirect$Term %in% toupper(genesDown)  ) 
    corrInverse <- data[ data$pertDir=="UP" & data$geneDir == direction,]
    corrInverse <- subset(corrInverse, corrInverse$Term %in% toupper(genesUp)  ) 
  }
  
# 
  corrDirectSummary <- mergeGeneRowEnrichR(corrDirect)
  corrInverseSummary <- mergeGeneRowEnrichR(corrInverse)
  
  output <- list()
  output[["Direct"]] <- corrDirect 
  output[["Inverse"]] <- corrInverse
  output[["DirectSummary"]] <- corrDirectSummary
  output[["InverseSummary"]] <- corrInverseSummary
  return(output)
}
