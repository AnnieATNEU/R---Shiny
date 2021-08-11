# Loading Pacman package to check if the necessary packages are
# installed and can be loaded, installs them from CRAN and loads 
# them if they aren't already present

#CREDITS TO BRIDGET>>>>

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(shiny, DBI, DT, RSQLite, dplyr, devtools, tidyverse, shinythemes,
       profvis, shinyWidgets, shinyBS, stringr, shinyjs, ggplot2, hrbrthemes,
       gridExtra, RColorBrewer)

dbname = 'vcf_df.db'
# Installing package 'radiator' separately because it requires devtools 
if(!require("radiator")) {
  devtools::install_github("thierrygosselin/radiator")
}

server <- shinyServer(
  
  function(input, output, session) {
    
    # * * * * * * * * * * Helper Functions * * * * * * * * * * * 
    is.not.null <- function(x) !is.null(x) 
    
    
    # Connection to database
    mydb <- dbConnect(RSQLite::SQLite(), dbname)
    
    tableList <- dbListTables(mydb) # Makes list of table names to use for choices later
    print(tableList)
    
    selected_data <- reactiveVal()
    
    # Requires tabnames and columns input (below) -- responds to observeEvent inputs with filtered data
    # Creates a dataset value ("filtered") based on the inputs 
    filtered_data <- eventReactive(c(input$reload,input$tabnames, input$cols, input$subset,
                                     input$CHROM, input$sourcefile, input$var_type,
                                     input$position, input$sub_type, input$DP, input$query,
                                     input$querytype, input$numquery, input$numquerytype),{
      req(input$cols)
      req(input$tabnames)
      filtered = selected_data()[, input$cols] # "filtered" is the valueExpression
      # This is the default data that will be returned without filters if there are none applied
      
      # Assuring that there are at least 2 columns to display, returns null for valueExpr if <1
      if(length(input$cols)<2) {
        
        showNotification("Please select at least 2 columns to show!")
        return(NULL)
        
      } 
      else { # Subset filters (difference)
        
        if(input$subset=="diff"){
          
          filtered = (filtered %>%
                        group_by(CHROM,var_type,POS) %>%
                        filter(n()==1)) # %>% from dplyr, shows the rows that only occur once  
          
        } 
        else if(input$subset=="sim"){ 
          num_src = length(unique(selected_data()[,'sourcefile']))
          if(length(input$sourcefile)==num_src) {
            sql <- 'select a.* from vcf a
                    join (select CHROM, var_type, POS, count(*)
                    from vcf
                    group by CHROM, var_type, POS
                    having count(*) =(
                    select count(distinct sourcefile)
                    from vcf)) b
                    on a.CHROM = b.CHROM
                    and a.var_type = b.var_type
                    and a.POS = b.POS
                    order by a.CHROM;'
          }
          else {
            sql <- 'select a.* from vcf a
                    join (select CHROM, var_type, POS, count(*)
                    from vcf
                    group by CHROM, var_type, POS
                    having count(*) > 1) b
                    on a.CHROM = b.CHROM
                    and a.var_type = b.var_type
                    and a.POS = b.POS
                    order by a.CHROM;'
          }
          # Subset filter if "similar" is checked
          # Shows the rows that have duplicate values in the database 

          filtered = dbGetQuery(mydb, sql)
        }
        
        # Numeric range filters
        if(!is.null(input$position)) {
          filtered = (filtered %>%
            filter(POS >= input$position[1],
                   POS <= input$position[2]))
        }
        
        if(!is.null(input$DP)) {
          filtered = (filtered %>%
                        filter(DP >= input$DP[1],
                               DP <= input$DP[2]))
        }
        
        # Sample choices
        if((input$sample!="") && ('sample' %in% input$cols)) { filtered=filtered[which(filtered$sample %in% input$sample),] }
        
        # CHROM, var_type, sourcefile filters 
        if('var_type' %in% input$cols) { filtered = filtered[which(filtered$var_type %in% input$var_type),] }
        if('CHROM' %in% input$cols) { filtered = filtered[which(filtered$CHROM %in% input$CHROM),] }
        if('sourcefile' %in% input$cols) { filtered = filtered[which(filtered$sourcefile %in% input$sourcefile),] }
        if('sub_type' %in% input$cols) { filtered = filtered[which(filtered$sub_type %in% input$sub_type),] }
        
        if(input$numquery!="") {
          search <- input$numquery
          term <- input$numquerytype
          queries = unlist(str_split(search, ", "))
          sql_query = ''
          
          first = word(queries[1], 1:3)
          sql_query = paste0("SELECT * FROM vcf WHERE ", first[1], " BETWEEN ", first[2], " AND ", first[3])
          from <- 2
          to <- length(queries)
          
          if(length(queries)>1) {
            for(x in from:to) {
              parts = word(queries[x], 1:3)
              this_sql <- paste0(" ", term," ", parts[1], " BETWEEN ", parts[2], " AND ", parts[3])
              
              if(x!=length(queries)) { this_sql <- paste0(this_sql,' ', term, ' ') }
              sql_query <- paste0(sql_query, this_sql)
            }
          }
          
          filtered = dbGetQuery(mydb, sql_query)
        }
        
        if(input$query!="") {
          search <- input$query
          term <- input$querytype
          queries = unlist(str_split(search, ", "))
          sql_query = ''
          
          first = word(queries[1], 1:2)
          sql_query = paste0("SELECT * FROM vcf WHERE ", first[1], " == '", first[2],"'")
          from <- 2
          to <- length(queries)
          
          if(length(queries)>1) {
            for(x in from:to) {
              parts = word(queries[x], 1:3)
              this_sql <- paste0(" ", term," ", parts[1], " == '", parts[2],"'")
              
              if(x!=length(queries)) { this_sql <- paste0(this_sql,' ', term,' ') }
              sql_query <- paste0(sql_query, this_sql)
            }
          }
          filtered = dbGetQuery(mydb, sql_query)
        }

        
        return(filtered)
        # Returns data with all applied filters (the filters can be stacked)
      }
    })
    
    
    # Stores tabnames (table names) input on server side using provided options/choices
    updateSelectizeInput(session, "tabnames", choices = tableList, server = TRUE)
    
    # General event handling for selecting table & updating column/samples choices based on 
    # selected table
    observeEvent(c(input$tabnames,input$reset),{
            req(input$tabnames)
            mydata <<- dbReadTable(mydb, input$tabnames)
            
            print(colnames(mydata))
            
            # Updates all filter choices separately based on unique values in supplied dataset
            updateSelectizeInput(session, "cols"
                                 , server=TRUE, choices = colnames(mydata),
                                 selected = colnames(mydata))
            
            updateSelectizeInput("sample",session = session, server=TRUE, choices = mydata[,'sample'],
                                  selected = mydata[,'sample']) 
            
            updateSelectizeInput("CHROM",session=session, server=TRUE, choices=mydata[,'CHROM'],
                                 selected = mydata[,'CHROM'])
    
            updateSelectizeInput("sourcefile",session=session, server=TRUE, choices=mydata[,'sourcefile'],
                                 selected=mydata[,'sourcefile'])
            
            updateSelectizeInput("var_type",session=session, server=TRUE, choices=mydata[,'var_type'],
                                 selected=mydata[,'var_type'])
            
            updateSelectizeInput("sub_type",session=session, server=TRUE, choices=mydata[,'sub_type'],
                                 selected=mydata[,'sub_type'])
            
            updateNumericRangeInput("position",session=session,
                                    value=c( 
                                      min(mydata[,'POS']),
                                      max(mydata[,'POS'])
                                    ))
            
            updateNumericRangeInput("DP",session=session,
                                    value=c( 
                                      min(mydata[,'DP']),
                                      max(mydata[,'DP'])
                                    ))
            
            updateTextInput("numquery",session=session,value="")
            
            updateTextInput("query", session=session, value="")
            
            selected_data(mydata) # Passes in the choices to the reactiveValue variable selected_data
    })
    
    # * * * * * * * * * * * * * * * PLOT OUTPUTS * * * * * * * * * * * * * * * #
    #
    #          * * * Sourcefile v. Mutations Stacked Bar Chart * * *
    # 
    # Plot generation builds a SQL query based on the number of sourcefiles
    # and variant types selected, and creates an array of the sources, vars,
    # and the count of each var type in each source. This array is used in the
    # params for the stacked bar chart.
    
    svmplot <- eventReactive(c(input$reload,input$tabnames, input$cols, input$subset,
                               input$CHROM, input$sourcefile, input$var_type,
                               input$position, input$sub_type, input$DP, input$query,
                               input$querytype, input$numquery, input$numquerytype), {
      req(input$sourcefile)
      req(input$var_type)
      
      svmarray <- filtered_data() %>% 
                    group_by(sourcefile) %>%
                    count(var_type)
      
      # Final plot creation
      plot <- ggplot(svmarray, aes(fill=var_type,y=n,x=sourcefile)) +
              geom_bar(position="stack", stat="identity") +
              scale_fill_brewer(palette = "Set3") +
              xlab("SAMPLES") +
              ylab("MUTATIONS") +
              scale_y_continuous(expand=c(0,0)) +
              ggtitle("Mutations per Variant-Type") +
              theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14))
      
      return(plot)
    })
    
    
    
    #        * * * Sourcefile v. Sub Mutations Stacked Bar Chart * * *
    # 
    # This plot generation uses the same algorithm as the Sourcefile v. Var Type
    # stacked bar chart.
    #
    
    subplot <- eventReactive(c(input$reload,input$tabnames, input$cols, input$subset,
                               input$CHROM, input$sourcefile, input$var_type,
                               input$position, input$sub_type, input$DP, input$query,
                               input$querytype, input$numquery, input$numquerytype), {
     req(input$sourcefile)
     req(input$sub_type)
     
     subarray <- filtered_data() %>% 
       group_by(sourcefile) %>%
       count(sub_type)
      
      # Final plot creation
      plot <- ggplot(subarray, aes(fill=sub_type,y=n,x=sourcefile)) +
        geom_bar(position="stack", stat="identity") +
        xlab("SAMPLES") +
        ylab("MUTATIONS") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_brewer(palette = "Set3") +
        ggtitle("Mutations per Sub Variant-Type") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14))

      
      return(plot)
    })
    
    
    
    # Output for grid of two sourcefile stacked bar charts
    output$mutplots <- renderPlot({
      p1 <- svmplot() # Sourcefile v. Mutations Plot
      p2 <- subplot() # Sourcefile v. Sub Mutation Plot
      grid.arrange(p1,p2, ncol=2, widths = c(1,1))
    })
    
    
    # Table output handling after filtering, handled on server side
    output$table <- DT::renderDataTable(
      {
        dataTableProxy(outputId = 'table') %>%
          hideCols(hide = c(7,8,9))
        filtered_data()
      }, server=TRUE, filter="top", options=list(pageLength=10,
                                                 searchHighlight=TRUE,
                                                 autoWidth=TRUE,
                                                 columnDefs = list(
                                                   list(width = '200px', targets=c(1:16)),
                                                   list(width = '5px', targets=0),
                                                   list(className = 'dt-center', targets = '_all')),
                                                 scrollX=TRUE)
      )

    

    
    session$onSessionEnded(function() { dbDisconnect(mydb) })
  })
  


ui_panel <- navbarPage("VSQV", theme=shinytheme("yeti"),
  
  tabPanel("Browse & Filter Data",
           
    sidebarLayout(
      sidebarPanel( 
     
      h4("Cohort Selection:"),
     
      # * * * * * * UI for filters in leftside panel * * * * * *  
     
      selectizeInput("tabnames",
                     tags$span("Table",
                               tags$i(
                                 id="tableinfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     choices=NULL),
  
      selectizeInput("cols",
                     tags$span("Columns",
                               tags$i(
                                 id="colsinfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     choices=NULL, multiple = T),
      
      selectizeInput("sourcefile",
                     tags$span("Sourcefile",
                               tags$i(
                                 id="srcinfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     c("Select all"=""), multiple=T),
      
      selectizeInput("sample",
                     tags$span("Samples",
                               tags$i(
                                 id="samplesinfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     c("Select all"=""), multiple=T ),  
     
      radioButtons("subset",
                   tags$span("Subset",
                             tags$i(
                               id="subsetinfo",
                               class = "glyphicon glyphicon-info-sign",
                               style = "color:#0072B2",
                               title = ""
                             )),
                   c("All"="all","Unique to one sample"="diff",
                                     "Appears in selected sourcefiles"="sim")),
     
      selectizeInput("CHROM",
                     tags$span("CHROM",
                               tags$i(
                                 id="chrominfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     c("Select all"=""),multiple=T),
     
      selectizeInput("var_type",
                     tags$span("Variant Type",
                               tags$i(
                                 id="varinfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     c("Select all"=""), multiple=T),
      
      selectizeInput("sub_type",
                     tags$span("Sub Variant Type",
                               tags$i(
                                 id="subinfo",
                                 class = "glyphicon glyphicon-info-sign",
                                 style = "color:#0072B2",
                                 title = ""
                               )),
                     c("Select all"=""), multiple=T),
      
      radioButtons("numquerytype",
                   tags$span("Range Query Type",
                             tags$i(
                               id="querytypeinfo",
                               class = "glyphicon glyphicon-info-sign",
                               style = "color:#0072B2",
                               title = ""
                             )), 
                   c('AND'='AND','OR'='OR')
                   ),
      
      textInput("numquery",
                tags$span("Range Query",
                          tags$i(
                            id="numqueryinfo",
                            class = "glyphicon glyphicon-info-sign",
                            style = "color:#0072B2",
                            title = ""
                          ))
                ),
      
      radioButtons("querytype",
                   tags$span("Text Query Type",
                             tags$i(
                               id="querytypeinfo",
                               class = "glyphicon glyphicon-info-sign",
                               style = "color:#0072B2",
                               title = ""
                             )), 
                   c('AND'='AND','OR'='OR')
      ),
      
      textInput("query",
                tags$span("Text Query",
                          tags$i(
                            id="queryinfo",
                            class = "glyphicon glyphicon-info-sign",
                            style = "color:#0072B2",
                            title = ""
                          ))
      ),
      
      
      # * * * * * * * * * Numeric Filters * * * * * * * *
      numericRangeInput(
        inputId = 'position', label = (
          tags$span("POS Range:",
                    tags$i(
                      id="posinfo",
                      class = "glyphicon glyphicon-info-sign",
                      style = "color:#0072B2",
                      title = ""
                    ))
        ),
        value = c(NULL,NULL)
      ),
      
      numericRangeInput(
        inputId = 'DP', label = (
          tags$span("DP Range:",
                    tags$i(
                      id="dpinfo",
                      class = "glyphicon glyphicon-info-sign",
                      style = "color:#0072B2",
                      title = ""
                    ))
        ),
        value = c(NULL,NULL)
      ),
      
      # * * * * * * * * * Tooltip UI * * * * * * * * * * *
      
      # Tooltips for Filters * * *
      
      bsTooltip(id="numqueryinfo",title="Query a range of values in one or more numeric (integer) columns in this format: ColumnName RangeStart RangeEnd",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="queryinfo",title="Query a range of values in one or more alphanumeric columns in this format: ColumnName Value",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="numquerytypeinfo",title="To filter based on one column only, select OR. To filter based on different columns, select WHERE.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="querytypeinfo",title="To filter based on one column only, select OR. To filter based on different columns, select WHERE.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="samplesinfo",title="Select the sample(s) for data to be shown from.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="srcinfo",title="Select which uploaded files are included in the data shown.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="colsinfo",title="Select which columns should be shown. Click anywhere in this box to reselect hidden columns.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="chrominfo",title="Records from only the scaffolds selected will be shown.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="varinfo",title="Records that match the type(s) of variants selected will be shown.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="subinfo",title="Records that match the type(s) of sub-variants selected will be shown.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="subsetinfo",title="Unique to one sample will select those records that occur only once in all the sourcefiles selected; difference will select records that occur in all sourcefiles selected.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="dpinfo", title="Enter the lower and upper values that the DP values should fall between.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="posinfo", title="Enter the lower and upper values that the Position values should fall between.",
                placement="right",
                trigger="hover"),
      
      bsTooltip(id="tableinfo", title="Choose which table to show data from.",
                placement="right",
                trigger="hover"),

      
      # * * * * * * * * * Bottom Row Buttons * * * * * * *
      splitLayout(
        # Resets the inputs
        actionButton("reset","Reset")
      )
  ), 

         # mainPanel that holds tabs Tables, Plots, and IGV tabs (Plots and IGV to be added)
         mainPanel(
          tabsetPanel( 
            
            tabPanel("Tables", 
                     
                     dataTableOutput("table") ),
            
            tabPanel("Plots",
                     fluidRow( plotOutput("mutplots") )
            )
        )
      )
    )
  ),
  tabPanel("IGV")
)

ui <- shinyUI(ui_panel)

shinyApp(ui=ui, server=server)
