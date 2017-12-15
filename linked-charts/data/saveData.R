library(tidyverse)
library(jsonlite)

datadir <- path.expand( "." )

mAP <- 
  read_csv( file.path( datadir, "LG10_SecB_AP.csv" ) ) %>% 
  filter( gene != "ddl" ) %>% 
  column_to_rownames( "gene" ) %>% as.matrix

mTT <- 
  read_csv( file.path( datadir, "LG11_SecB_TT.csv" ) ) %>% 
  filter( gene != "ddl" ) %>% 
  column_to_rownames( "gene" ) %>% as.matrix

mAP <- mAP[ rownames(mAP) %in% rownames(mTT), ]
mTT <- mTT[ rownames(mTT) %in% rownames(mAP), ]

jsData <- list()

for(gene in rownames(mAP)){
  AP <- mAP[gene, ]
  TT <- mTT[gene, ]
  maxInd <- max(which(AP + TT != 0))
  AP <- AP[1:maxInd]
  TT <- TT[1:maxInd]
  if(sum(AP + TT) > 10 && maxInd > 10)
    jsData[[gene]] <- list(AP = AP, TT = TT);
}

writeLines(toJSON(jsData), "rpf_input.json")
