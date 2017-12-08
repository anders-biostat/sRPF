library( tidyverse )
library( forcats )
library( binom )

datadir <- path.expand( "~/w/diary-simon/2017-Q3/Lena" )

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

stopifnot( identical( dim(mAP), dim(mTT) ) )
stopifnot( identical( rownames(mAP), rownames(mTT) ) )

# Parameters:
thresh <- .6
binsize <- 6

tibble(
  gene = rownames(mTT)[row(mTT)],
  pos = as.vector( col(mTT) ),
  TT = as.vector(mTT),
  AP = as.vector(mAP) ) %>%
group_by( gene ) %>%
filter( pos <= max( pos[ AP>0 | TT>0 ] ) ) %>%
mutate( bin = floor( pos / binsize ) ) %>%
group_by( gene, bin ) %>%
summarise( 
  AP = sum(AP),
  TT = sum(TT),
  leftpos = min(pos) ) %>%
do( bind_cols( ., 
   binom.agresti.coull( .$AP, .$TT + .$AP ) %>% select( lower, upper ) ) ) -> 
tbl

prob2odds <- function( p ) p / ( 1 - p )

tbl %>%
ungroup %>%
mutate( 
  lower = pmax( 0, prob2odds( lower ) ),
  upper = pmin( 1, prob2odds( upper ) ) ) %>% 
mutate( binding = 
   ifelse( lower > thresh, "+",
     ifelse( upper < thresh, "-",
        "?" ) ) ) %>%
group_by( gene ) %>%
filter( sum( binding != "?" ) > 10 ) ->
tbl2

tbl2 %>%
group_by( gene ) %>%
summarise( d = sum(binding=="+") - sum(binding=="-") ) %>%
arrange( d ) %>%
{ .$gene } ->
genes_ordered

pdf( "tallplot.pdf", width=20, height=300 )
tbl2 %>%
ungroup %>%
mutate( gene = factor( gene, genes_ordered ) ) %>%
ggplot +
geom_tile(aes( 
    x = leftpos + binsize/2,
    y = gene,
    fill = binding ),
  width = binsize, 
  height = 1 )  +
scale_fill_manual( values = c( `-` = "coral3", `+` = "deepskyblue4", `?` = "gray70" ) )
dev.off()
  
