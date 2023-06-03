#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Author: Belinda B. Garana; last edit: 2023-05-09

library(shiny);library(utils);library(DT);library(dplyr);

# set illegal filename characters
illegal.chars <- c("#","<","%",">","!","`","&","'","=","}","/",":","@") # source: https://www.mtu.edu/umc/services/websites/writing/characters-avoid/; would be nice to protect against " and \ too
illegal.chars.need.brackets <- c("$","+","*","|","{","?") # for these chars, gsub needs brackets but str_contains can't have brackets

# load functions
as.filename <- function(moa.name){
  # if moa.name contains illegal file name characters, replace these characters with "-" or "_"
  if(sjmisc::str_contains(moa.name, illegal.chars, logic="or")){
    for(j in 1:length(illegal.chars)){
      moa.name <- gsub(illegal.chars[j], "-", moa.name)
    }
  }
  if(sjmisc::str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
    for(j in 1:length(illegal.chars.need.brackets)){
      moa.name <- gsub(paste0("[",illegal.chars.need.brackets[j],"]"), "-", moa.name)
    }
  }
  if(sjmisc::str_contains(moa.name, " ")){
    moa.name <- gsub(" ", "_", moa.name)
  }
  return(moa.name)
}

prep.files <- function(drug.dir, moa.dir, moa = "", drug = ""){
  # get all results
  cat(file=stderr(), "About to get results from GitHub", "\n")
  allResults <- utils::read.csv(paste0(moa.dir, "MOAL_results.csv"))
  qcResults <- utils::read.csv(paste0(moa.dir, "MOAL_qc_results.csv"))
  longResults <- utils::read.csv(paste0(moa.dir, "MOAL_DMEA_results.csv"))
  allDrugResults <- utils::read.csv(paste0(drug.dir, "MOAL_drug_DMEA_results.csv"))
  misResults <- utils::read.csv(paste0(drug.dir, "MOAL_drugs_failing_qc.csv"))
  
  # get results for MOA of interest (if any) TO DO: check if works if MOA not found ####
  if(moa != ""){
    cat(file=stderr(), "About to get results for MOA of interest", "\n")
    moaResults <- longResults[longResults$moa == moa, ]  
  }else{
    moaResults <- NULL
  }
  
  # get results for drug of interest (if any) TO DO: check if works if drug not found ####
  if(drug != ""){
    cat(file=stderr(), "About to get results for drug of interest", "\n")
    drugResults <- allDrugResults[allDrugResults$Drug == drug, ]
    
    # if drug not found, try to match with synonym
    if(nrow(drugResults) == 0){
      cat(file=stderr(), "About to check for drug synonym", "\n")
      drug.synonym <- NULL
      
      # access drug synonyms
      synonyms <- utils::read.csv(
        paste0("https://raw.github.com/BelindaBGarana/DMEA/",
               "shiny-app/Inputs/PRISM_drug_synonyms.csv")
      )
      
      # check if drug name of interest is a CID for a PRISM drug
      if (drug %in% synonyms$PubChem_CID) {
        drug.synonym <-
          synonyms[synonyms$PubChem_CID == drug, ]$PRISM_drug_name
      } else {
        # else check if unannotated drug is a synonym for a PRISM drug name
        for (j in seq_len(nrow(synonyms))) {
          temp.synonyms <- strsplit(synonyms$PubChem_synonyms[j], "[|]")
          if (drug %in% temp.synonyms) {
            drug.synonym <- synonyms[j, ]$PRISM_drug_name
            break
          }
          
          # get results for drug synonym
          if(!is.null(drug.synonym)){
            cat(file=stderr(), "About to get results for drug synonym", "\n")
            drugResults <- allDrugResults[allDrugResults$Drug == drug.synonym, ] 
          }
        }
      }
    }
  }else{
    drugResults <- NULL
  }
  
  return(list(allResults = allResults, 
              qcResults = qcResults, 
              longResults = longResults, 
              allDrugResults = allDrugResults, 
              misResults = misResults, 
              moaResults = moaResults, 
              drugResults = drugResults))
}

# Define UI for application
ui <- fluidPage(
  # Application title
  titlePanel("Mechanism-of-action Landscape: Web Application"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      # bubble select: which drug screen? default: prism
      radioButtons(inputId = "screen", label = "Which drug screen would you like to query?",
                   choices = c("PRISM (1351 drugs with 85 MOAs tested in 320 adherent cancer cell lines)" = "PRISM",
                               "CTRPv2 (481 drugs with 14 MOAs tested in 406 adherent cancer cell lines)" = "CTRPv2"),
                   selected = "PRISM"),
      
      conditionalPanel(
        condition = "input.screen == 'PRISM'",
        # drop down select: tissue type?
        selectInput(inputId = "tissue", label = "Would you like to filter by tissue type?",
                    choices = c("No tissue type selected" = "none",
                                "Lung (63 cell lines)" = "lung",
                                "Skin (29 cell lines)" = "skin",
                                "CNS (25 cell lines)" = "cns",
                                "Ovary (24 cell lines)" = "ovary",
                                "Pancreas (23 cell lines)" = "pancreas",
                                "Esophagus (20 cell lines)" = "esophagus",
                                "Breast (17 cell lines)" = "breast",
                                "Urinary Tract (16 cell lines)" = "urinary",
                                "Liver (16 cell lines)" = "liver",
                                "Upper Aerodigestive (14 cell lines)" = "ua",
                                "Uterus (13 cell lines)" = "uterus",
                                "Gastric (11 cell lines)" = "gastric",
                                "Kidney (11 cell lines)" = "kidney",
                                "Colorectal (10 cell lines)" = "crc",
                                "Soft Tissue (8 cell lines)" = "soft",
                                "Bone (6 cell lines)" = "bone",
                                "Peripheral Nervous System (5 cell lines)" = "pns",
                                "Bile Duct (4 cell lines)" = "bile_duct",
                                "Thyroid (3 cell lines)" = "thyroid"))
      ),
        
      # bubble select: remove potentially misannotated drugs? default: fALSE
      radioButtons(inputId = "rm_qc_fail", label = "Remove potentially misannotated drugs from consideration?",
                   choices = c("Yes" = TRUE,
                               "No" = FALSE), selected = TRUE),
 
      # get drug moa of interest (if any)
      textInput(inputId = "moa", label = "Optional: enter a MOA of interest to view its results separately (case-sensitive; e.g., HMGCR inhibitor)"),
      
      # get drug of interest (if any)
      textInput(inputId = "drug", label = "Optional: enter a drug of interest to view its results separately (case-sensitive; e.g., pitavastatin)"),
      
      # include "Run" button
      actionButton(inputId = "run", label = "Run Query"),
      
      # including loading message
      conditionalPanel(
        condition = "input.run && !output.msg",
        textOutput(outputId = "fyi")
      ),
      
      # indicate when run is completed
      textOutput("msg")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      # display plots in tabs
      tabsetPanel(
        tabPanel("Results for all MOA pairs", DT::dataTableOutput("allResults"),
                 conditionalPanel(
                   condition = "output.msg=='Query completed'",
                   downloadButton(outputId = "allDownload", label = "Download results for all MOA pairs")
                 )),
        tabPanel("Results for MOA of interest", DT::dataTableOutput("moaResults"),
                 conditionalPanel(
                   condition = "output.msg=='Query completed'",
                   downloadButton(outputId = "moaDownload", label = "Download results for MOA pair of interest")
                 )),
        tabPanel("Results for drug of interest", DT::dataTableOutput("drugResults"),
                 conditionalPanel(
                   condition = "output.msg=='Query completed'",
                   downloadButton(outputId = "drugDownload", label = "Download results for drug of interest")
                 )),
      ),
      
      # output result files (.zip)
      conditionalPanel(
        condition = "output.msg=='Query completed'",
        downloadButton(outputId = "results", label = "Download all results")
      ),
      
      uiOutput("info"),
      textOutput("private"),
      textOutput("refresh")
    )
  )
)

# Define backend
server <- function(input, output) {
  url <- a("https://belindabgarana.github.io/MOAL", href = "https://belindabgarana.github.io/MOAL")
  output$info <- renderUI({tagList("For more information or to contact us, please visit: ", url)})
  output$refresh <- renderText({"Please note that '0' values are < 0.001 and refresh your web browser after each query to avoid errors. You will also need to refresh this webpage after 5 minutes of inactivity."})
  output$fyi <- renderText({"Running query... Please allow 1 to 5 minutes of run time"})
  observeEvent(input$run, {
    # get input parameters
    screen <- input$screen
    tissue <- input$tissue
    rm_qc_fail <- input$rm_qc_fail
    moa <- input$moa
    drug <- input$drug
    
    # get screen directory
    screen.dir <- paste0("https://raw.github.com/BelindaBGarana/MOAL/shiny-app/Results/", screen, "/")
    
    if(screen == "CTRPv2"){
      ##### Using CTRPv2 drug screen (all adherent cancer cell lines) #####
      # get directory for drug-level results
      screen.dir <- paste0(screen.dir, "all_adherent/")
      
      # get directory for MOA-level results
      if(rm_qc_fail){
        moa.dir <- paste0(screen.dir, "wo_drugs_failing_qc/")
      }else{
        moa.dir <- screen.dir
      }
      
      # get files
      files <- prep.files(screen.dir, moa.dir, moa, drug)
    }else if(screen == "PRISM"){
      ##### Using PRISM drug screen (all adherent cancer cell lines) #####
      # get directory
      if(tissue != "none"){
        screen.dir <- paste0(screen.dir, tissue, "/")
      }else{
        screen.dir <- paste0(screen.dir, "all_adherent/")
      }
      
      if(rm_qc_fail){
        moa.dir <- paste0(screen.dir, "wo_drugs_failing_qc/")
      }else{
        moa.dir <- screen.dir
      }
      
      # get files
      files <- prep.files(screen.dir, moa.dir, moa, drug)
    }
    
    # format data tables for view
    allResultsTable <- files$allResults %>% dplyr::mutate_if(is.numeric, ~round(., 3))
    
    if(moa != ""){
      moaResultsTable <- files$moaResults %>% dplyr::mutate_if(is.numeric, ~round(., 3))
    }
    
    if(drug != ""){
      drugResultsTable <- files$drugResults %>% dplyr::mutate_if(is.numeric, ~round(., 3))
    }
    
    # output data tables for view
    cat(file=stderr(), "About to output results for view", "\n")
    
    output$allResults <- DT::renderDataTable({
      DT::datatable(allResultsTable[allResultsTable$qc_pass & !allResultsTable$self, 
                                     c("pair", "NES_avg", "NES_sd", "Fisher_p", "percent_overlap")],
                    options = list(order = list(list(2, 'desc'))))
    })
    output$moaResults <- DT::renderDataTable({
      DT::datatable(moaResultsTable[ , c(3:6, 8:ncol(moaResultsTable))],
                    options = list(order = list(list(3, 'desc'))))
    })
    output$drugResults <- DT::renderDataTable({
      DT::datatable(drugResultsTable[ , c(3:6, 8:(ncol(drugResultsTable)-2))],
                    options = list(order = list(list(3, 'desc'))))
    })
    
    # output all results for download TO DO (also check on orderClasses above) ####
    cat(file=stderr(), "About to output results for download", "\n")
    output$results <- downloadHandler(
      filename = function(){
        main.name <- paste0("MOAL_results_", screen, "_")
        
        # add tissue selection (if any)
        if(tissue != "none"){
          main.name <- paste0(main.name, tissue, "_")
        }
        
        # add rm_qc_fail parameter (if used)
        if(rm_qc_fail){
          main.name <- paste0(main.name, "wo_drugs_failing_qc_")
        }
        
        # add date and file type
        paste0(main.name, Sys.Date(), ".zip")
        },
      
      content = function(file){
        all.files <- c("MOAL_results.csv",
                       "MOAL_qc_results.csv",
                       "MOAL_DMEA_results.csv",
                       "MOAL_drug_DMEA_results.csv",
                       "MOAL_drugs_failing_qc.csv")
        
        # add file name for results for MOA of interest (if any)
        if(moa != ""){
          all.files <- c(all.files,
                         paste0("MOAL_DMEA_results_", as.filename(moa), ".csv"))
        }
        
        # add file name for results for drug of interest (if any)
        if(drug != ""){
          all.files <- c(all.files,
                         paste0("MOAL_drug_DMEA_results_", as.filename(drug), ".csv"))
        }
        
        # output files
        for(i in seq_len(length(all.files))){
          utils::write.csv(files[i], all.files[i], row.names = FALSE)
        }
        zip(zipfile=file, files=all.files)}
    )
    
    output$msg <- renderText({"Query completed"})
  }
  )
}

# Run the application
shinyApp(ui = ui, server = server)

