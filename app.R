#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library( tidyverse )
library( forcats )
library( binom )

datadir <- path.expand( "~/w/diary-simon/2017-Q3/Lena" )

mAP <- 
  read_csv( file.path( datadir, "LG10_SecB_AP.csv" ), n_max=10 ) %>% 
  filter( gene != "ddl" ) %>% 
  column_to_rownames( "gene" ) %>% as.matrix

mTT <- 
  read_csv( file.path( datadir, "LG11_SecB_TT.csv" ), n_max=10 ) %>% 
  filter( gene != "ddl" ) %>% 
  column_to_rownames( "gene" ) %>% as.matrix

mAP <- mAP[ rownames(mAP) %in% rownames(mTT), ]
mTT <- mTT[ rownames(mTT) %in% rownames(mAP), ]

stopifnot( identical( dim(mAP), dim(mTT) ) )
stopifnot( identical( rownames(mAP), rownames(mTT) ) )

# Function to get data frame, with binning
get_df <- function( gene, bin=3 ) {
  df <- tibble( 
    AP = mAP[ gene, ],
    TT = mTT[ gene, ],
    pos = 0:(ncol(mAP)-1) ) %>%  
    filter( pos <= max(pos[TT>0]) ) %>%
    mutate( bin = floor( pos / bin ) ) %>%
    group_by( bin ) %>%
    summarise( AP = sum(AP), TT=sum(TT) ) %>%
    mutate( 
      lo = prob2odds( pmax( 0, binom.agresti.coull( AP, AP + TT )$lower ) ),
      hi = prob2odds( pmin( 1, binom.agresti.coull( AP, AP + TT )$upper ) ) ) %>%
    mutate( lo = ifelse( is.nan(lo), 0, lo ) ) %>%
    mutate( hi = ifelse( is.nan(hi), Inf, hi ) ) %>%
    mutate( binding = 
              ifelse( lo > thresh, "+",
                      ifelse( hi < thresh, "-",
                              "?" ) ) )
  attr( df, "gene" ) <- gene
  attr( df, "bin" ) <- bin
  df
}

make_plot <- function( df ) {
  binsize <- attr( df, "bin" )
  ggplot( df ) +
    scale_y_continuous( trans="log2", limits = 2^c(-3.8,4.8), oob=scales::squish, expand=c(0,0),
                        breaks=2^(-3:4)) +
    scale_alpha_continuous( trans="sqrt", breaks=c(3,10,30,100), limits=c(0,100), name="prec", range=c(0.02,1) ) +
    theme_bw() + theme( panel.grid.minor.y = element_blank() ) +
    geom_rect( aes( xmin=bin*binsize/3, ymin=lo/thresh, xmax=(bin+1)*binsize/3, ymax=hi/thresh, alpha=1/(1/TT+1/AP),
                    fill = binding ) ) +
    scale_fill_manual( values = c( `-` = "red", `+` = "blue", `?` = "gray20" ) ) +
    ggtitle( attr( df, "gene" ) ) + xlab( "AA" ) + ylab( "AP/TT" )
}


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("secB data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        textInput("gene", label = h3("Gene"), value = "ybcL")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("genePlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   df <- reactive( get_df( input$gene ) );
  
   observe( output$genePlot <- renderPlot(
     {
        make_plot( df() )
     },
     width = nrow(df()) * 2 + 50 ) )
   
}

# Run the application 
shinyApp(ui = ui, server = server)

