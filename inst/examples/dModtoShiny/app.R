#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinydashboard)
library(deSolve)
library(plotly)
library(dMod)
library(plyr)


#abspath <- getwd() #"homemfehlingShinyAppsdModtoShiny" #getwd() #achtung: geht nix mit slash

models <- system("ls -d */",intern = TRUE)

#list.dirs(recursive = FALSE)

# Define UI for application that draws a histogram ----
ui <- dashboardPage(
  dashboardHeader(title = "dMod"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      selectInput(inputId = "model", label = "Model:", choices = models, selected = models[1], multiple = FALSE),
      menuItem("Plot Prediction", tabName = "model", icon = icon("line-chart")),
      menuItem("Optima Selection", tabName = "parameters", icon = icon("list-ol")),
      menuItem("Parameter Tuning", tabName = "parvalue", icon = icon("sliders")),
      menuItem("Run on Profiles", tabName = "profiles", icon = icon("hand-pointer-o")),
      menuItem("More Plots", tabName = "addplots", icon = icon("area-chart"),
               menuSubItem("Fluxes", tabName = "addplotsflux"),
               menuSubItem("Residuals", tabName = "addplotsres")
      ),
      menuItem("View Reactions", tabName = "reactions", icon = icon("exchange")),
      uiOutput("moreinfo")   
    )
  ),
  dashboardBody(
    
    tabItems(
      tabItem(tabName = "model",
              fluidRow(
                column(width = 3, 
                       box(width = 12,
                           tabBox(id = "plotselector", width = 12,
                                  tabPanel("general",
                                           sliderInput("times", "Time points:", min = 0, max = 300, value = c(0, 300)),
                                           selectInput("color", "Color:", list(), multiple = FALSE),
                                           radioButtons("facetsetup","Facets:", c("wrap" = "wrap", "grid (y~1, x~2,3,...)" =  "grid")),
                                           selectInput("facet", label = NULL, list(), multiple = TRUE),
                                           shiny::p(strong("y-axis scale:")),
                                           checkboxInput("log","Log y-axis", value = FALSE),
                                           checkboxInput("yfree","free y-axis between panels", value = FALSE),
                                           shiny::p(strong("Uncertainties:")),
                                           checkboxInput("errdata","data", value = FALSE),
                                           uiOutput("errmod")
                                           
                                  ),
                                  tabPanel("subsetting",
                                           selectInput("states", "States:", list(), multiple = TRUE),
                                           uiOutput("morecons"),
                                           selectInput("conditions", "Conditions:", list(), multiple = TRUE)
                                  )
                           ),
                           sliderInput('plotHeight', 'Height of plot (in pixels)', min = 100, max = 2000, value = 600),
                           checkboxInput(inputId = "plotlypred", "plotly rendering (interactive plot)", value = FALSE),
                           actionButton("plot", "plot", style = "color: white; background-color: #337ab7; border-color: #2e6da4;")
                       )
                ),
                column(width = 9,  
                       useShinyjs(),
                       box(width = 12, height = 650,
                           #plotOutput("distPlot")
                           #withSpinner(plotlyOutput("plotnew"))
                           uiOutput("plotnew")
                           #uiOutput("plot.ui")
                       )
                )
              )
              
      ),
      tabItem(tabName = "parameters",
              fluidRow(
                box(width = 3,
                    sliderInput("tol", "Tolerance:", min = 0.1, max = 1, value = 0.5),
                    selectInput("jump", "Steps:", list(1), selected = c(1), multiple = TRUE),
                    #p(strong("Parameter values:")),
                    radioButtons("allstepsval","Parameter values:", c("All fits" = "all", "Selected steps"= "sel"), selected = "sel"),
                    actionButton("plotvals", "plot parameter values"),
                    actionButton("exportsteppars","export Parameters")
                ),
                box(width = 9,
                    plotOutput("stepplot"),
                    withSpinner(plotOutput("parvalplot"))
                )
              )
      ),
      tabItem(tabName = "parvalue", # set individual parameter values
              fluidRow(
                column(width = 3,
                       box(width = 12,
                           actionButton("plotind","Reset Parameters"),
                           uiOutput("parssel"),
                           sliderInput('plotHeightind', 'Height of plot (in pixels)', min = 100, max = 2000, value = 600)
                       ),
                       box(width = 12,
                           title = "Hint:",
                           "for fast processing select few conditions to be shown",br(),
                           actionButton("select2","select"),
                           actionButton("selectflux","fluxes")
                       )
                ),
                column(width = 9,
                       tabBox(id = "indparplots", width = 12, height = 600,
                              
                              tabPanel(value = "indparpredTab","Prediction",
                                       uiOutput("indparpredplotui")
                              ),
                              tabPanel(value = "indparfluxTab","Fluxes",
                                       uiOutput("plotfluxind")
                              )
                       )
                       
                )
              )
      ),
      tabItem(tabName = "profiles",
              fluidRow(
                column(width = 3,
                       box(width = 12,
                           title = "Hint:",
                           "for fast processing select few conditions to be shown",br(),
                           actionButton("select","select")
                       ),
                       box(width = 12,
                           uiOutput("profssel"),
                           selectInput("modeprof", "Select contributions:",list(), multiple = TRUE),
                           actionButton("plotprofiles","Plot profiles"),                          
                           selectInput("whichprof", "Choose profile:", list(), multiple = FALSE),
                           withSpinner(plotOutput("singleprof")),
                           sliderInput("profval","", min = NA, max = NA, value = NA, ticks = FALSE),
                           actionButton("plotprof","Plot prediction"),
                           sliderInput('plotHeightprof', 'Height of plot (in pixels)', min = 100, max = 2000, value = 600),
                           actionButton("exportprofpars","export Parameters")
                       )
                ),
                column(width = 9,
                       tabBox(id = "profileplots", width = 12, height = 600,
                              
                              tabPanel(value = "profilesTab","Profiles",
                                       uiOutput("profplotui")
                                       #withSpinner(plotOutput("profplot",height = 600))
                              ),
                              tabPanel(value = "pathTab","Path",
                                       uiOutput("profpathui")
                                       #withSpinner(plotOutput("profpath",height = 600))
                              ),
                              tabPanel(value = "predictionTab", "Prediction",
                                       uiOutput("profpredui")
                                       #withSpinner(plotOutput("profpred",height = 600))
                              )
                       )
                       #                       box(width = 12,
                       #                    plotOutput("predprof")
                       #                    )
                )
              )
      ),
      tabItem(tabName = "addplotsflux",
              fluidRow(
                column(width = 3,
                       box(width = 12,
                           selectInput("whichflux", "Select fluxes:", list(), multiple = TRUE),
                           textInput("customflux", "Custom fluxes:",placeholder = "e.g. A,B or k1*A-k2*B"),
                           actionButton("plotflux","Plot", 
                                        style = "color: white; background-color: #337ab7; border-color: #2e6da4;")
                       ),
                       box(width = 12,
                           sliderInput('plotHeightFlux', 'Height of plot (in pixels)', min = 100, max = 2000, value = 600),
                           checkboxInput(inputId = "plotlyflux", "plotly rendering (interactive plot)", value = FALSE)
                       )
                ),
                column(width = 9,
                       box( width = 12, height = 650,
                            uiOutput("plotfluxesui")
#                            withSpinner(plotlyOutput("plotfluxes"))
                       )
                )
              )
      ),
      tabItem(tabName = "addplotsres",
              fluidRow(
                column(width = 3,
                       box(width = 12,
                           radioButtons("splitres", "Option:", choices=c(
                             "time resolved" = "time", 
                             "per optimum" = "jump", 
                             "per condition and name" = "conditionname",
                             "per condition" = "condition",
                             "per name" = "name"
                           )),
                           h5("uses selected optima (default: 1)"),
                           actionButton("plotresiduals","Plot")
                       )
                ),
                column(width = 9,
                       box( width = 12,height = 650,
                            withSpinner(plotOutput("plotres",height = 600))
                       )
                )
              )
      ),
      tabItem(tabName = "reactions",
              tableOutput("reactionTable")
              
      )
    )
    
  ))

####################################################################################################
# Define server logic

server <- shinyServer(function(input, output, session) {
  
  # variables where input choices are defined in model setup
  vmodel <- reactiveValues(vstates = NULL, vconditions = NULL, vgridconditions = character(0), vtimes = NULL, vwhichprof = NULL, vprofssel = 1, vfacetsetup = "wrap", vmodeprof = NULL, vparameters = NULL, vparametersInd = NULL)
  vextra <- reactiveValues(ngrid = 0)
  vmodelold <- reactiveValues(name = NULL, myso = NULL) # needed to reset plot, path
  vplots <- reactiveValues(profsplot = NULL, plotdata = NULL, pf = NULL) 
  
  # loading of a model ----
  model <- eventReactive(input$model,{
    vplots$plotdata = NULL
    vplots$profsplot = NULL
    vplots$pf = NULL
    
    if(!input$plot){
      vmodelold$name <- input$model
    }

    
    # loading libraries and input
    for(so in vmodelold$myso){
      dyn.unload(so)  # unfortunatelly this does not seem to work
    }

    load(file = paste0(input$model,"input_shiny.RData"))

    myso <- paste0(input$model, list.files(path = input$model, "*.so"))
    
    for(so in myso){
      dyn.load(paste0(so))
    }
    
    if(!exists("errmodel")){
      errmodel <- NULL
    }
    
    if(!exists("pubref")){
      pubref <- "none"
    }
    
    # setting up specific input Selectors
    conditions <- getConditions(x)
    conditionsData <- getConditions(data)
    namesData <- sort(unique(do.call(c,lapply(conditionsData, function(con){
      as.character(unique(data[[con]]$name))
    }))))
    states <- reactions$states
    updateSelectInput(session, "states", choices = unique(c(namesData,states)), selected = namesData)
    vmodel$vstates <- namesData
    updateSelectInput(session, "conditions", choices = conditions, selected = conditionsData)
    vmodel$vconditions <- conditionsData
    vmodel$vgridconditions <- setdiff(names(attr(data,"condition.grid")),"condition")
    updateRadioButtons(session = session, "facetsetup", selected = "wrap")
    vmodel$vfacetsetup <- "wrap"
    times <- sort(unique(do.call(c,lapply(conditionsData, function(c) data[[c]]$time))))
    updateSliderInput(session, "times", min = min(0,times), max = 2*max(times), value = c(min(times),max(times)))
    vmodel$vtimes <- c(min(times),max(times))
    profs <- profiles
    if(is.list(profs) & !("parframe" %in% class(profs))){
      profs <- profs[[1]]
      vmodel$vprofssel <- 1
    }
    vmodel$vparameters <- as.parvec(parameters)
    vmodel$vparametersInd <- as.parvec(parameters)
    profnames <- unique(profs$whichPar)
    updateSelectInput(session=session, "whichprof", choices = profnames, selected = profnames[1])
    vmodel$vwhichprof <- profnames[1]
    profmodes <- attr(profs, "obj.attributes")
    vmodel$vmodeprof <- profmodes
    updateSelectInput(session=session, "modeprof", choices = profmodes, selected = profmodes)
    
    # reset other reactive values
    vpred$jump <- 1
    
    # put all main objects into the model() return
    list(x = x, data = data, reactions = reactions, fixed = fixed, parameters = parameters, profiles = profiles, err = errmodel, timesData = times, ref = file.path(pubref))
  })
  
  # update of input choices ----
  observe({
    vmodel$vstates <- input$states
  })
  observe({
    vmodel$vconditions <- input$conditions
  })
  observe({
    vmodel$vtimes <- input$times
  })
  observe({
    vmodel$vwhichprof <- input$whichprof
  })
  observe({
    vmodel$vmodeprof <- input$modeprof
  })
  observe({
    vmodel$vfacetsetup <- input$facetsetup
  })
  
  
  # more inputs for subsetting coming from condition.grid in data ----
  output$morecons <- renderUI({
    mycons <- vmodel$vgridconditions
    vextra$ngrid <- length(mycons)
    if(length(mycons)<1)
      return()
    congrid <- attr(model()$data,"condition.grid")
    
    myinputs <- lapply(1:length(mycons), function(i){
      mychoices <-  unique(congrid[,i])
      selectInput(inputId = mycons[i], label = mycons[i], choices = mychoices, selected = mychoices, multiple = TRUE)
    })
  })
  
  output$errmod <- renderUI({
    errors <- model()$err
    if(is.null(errors))
      return()
    
    checkboxInput("errors","errormodel", value = FALSE)
  })
  
  output$parssel <- renderUI({
    pars <- vmodel$vparameters
    myinputs <- lapply(1:length(pars), function(pp){
      pari <- pars[pp]
      sliderInput(inputId = names(pari), label = names(pari), min = round(pari-3,digits = 2), max = round(pari+3),value = pari,step = 0.01, round = FALSE)
    })
    vmodel$vparametersInd <- pars
    return(myinputs)
  })
  
  output$profssel <- renderUI({
    profs <- model()$profiles
    if(is.list(profs) & !("parframe" %in% class(profs))){
        sI <- selectInput("profssel",label = "More than one set of profiles (select one):", choices = 1:length(profs), selected = 1)
        rB <- radioButtons(inputId = "plotallprofs", label = "Plotted profiles:", choices = c("all","selected"),inline = TRUE)
        return(list(sI,rB))
    }else
      return()
    })
  
  observe({
    mycons <- isolate(vmodel$vgridconditions)
    if(vextra$ngrid>0){
      congrid <- isolate(attr(model()$data,"condition.grid"))
      for(i in 1:length(mycons)){
        con <- mycons[i]
        change <- input[[con]]
        change <- change[1:length(change)]
        congrid <- congrid[congrid[[con]] %in% change, ]
      }
      myconditions <- rownames(congrid)
      updateSelectInput(session, "conditions", selected = myconditions)
      vmodel$vconditions <- myconditions
    }
  })
  
  observe({
    pars <- vmodel$vparameters
    #pars <- isolate(vmodel$vparametersInd)
    if(!is.null(pars)){
    for(i in 1:length(pars)){
      change <- input[[names(pars[i])]]
      change <- change[1:length(change)]
      if(!is.null(change)) pars[i] <- change
    }
    vmodel$vparametersInd <- pars
    }
  })
  
  observe({
    vmodel$vprofssel <- as.numeric(input$profssel)
  })
  
  ############### model prediction ########################
  
  vpred <- reactiveValues(jump = 1, color = "name", facet = "condition")
  
  observe({
    choicesPlot <- unique(c("condition", "name", names(attr(model()$data,"condition.grid"))))
    selC <- "name"
    selF <- "condition"
    if(length(vpred$jump) > 1){
      choicesPlot <- c(choicesPlot, "step")
      selC <- "step"
      selF <- c("name","condition")
    }
    updateSelectInput(session, "color", choices = choicesPlot, selected = selC)
    updateSelectInput(session, "facet", choices = choicesPlot, selected = selF)
    vpred$color <- selC
    vpred$facet <- selF
  })
  
  observe({
    vpred$jump <- input$jump
  })
  
  observe({
    vpred$color <- input$color
  })
  
  observe({
    vpred$facet <- input$facet
  })
  
  prediction <- reactive({
    fun <- model()$x
    steps <- vpred$jump
    parametersSteps <- model()$parameters[as.numeric(steps)]
    
    nfits <- nrow(parametersSteps)
    predFrame <- do.call(rbind, lapply(1:nfits, function(i) {
      pars <- as.parvec(parametersSteps, i)
      pred <- fun(seq(from = min(0,model()$timesData), to = 2*max(model()$timesData), length.out = 500), pars, fixed= model()$fixed, deriv = FALSE)
      pred <- as.data.frame(pred, data = model()$data, errfn = model()$err)
      cbind(pred, step = as.character(steps[i]))
    }))
    return(predFrame)
  })
  
  predictionSubset <- reactive({
    pred <- prediction()
    subset(pred, name %in% vmodel$vstates & condition %in% vmodel$vconditions & time > vmodel$vtimes[1] & time < vmodel$vtimes[2])
  })
  
  dataSubset <- reactive({
    data <- cbind(as.data.frame(model()$data), step = "1")
    subset(data, name %in% vmodel$vstates & condition %in% vmodel$vconditions & time > vmodel$vtimes[1] & time < vmodel$vtimes[2])
  })
  
  observeEvent(input$plot, {
    colorP <- vpred$color
    predP <- predictionSubset()
    dataP <- dataSubset()
    facetsP <- vpred$facet
    facetsetupP <- vmodel$vfacetsetup
    P <- ggplot(predP, aes_string(x = "time", y = "value", color = colorP, fill = colorP, ymin = "value-sigma", ymax = "value+sigma")) + 
      theme_dMod() + scale_color_dMod() + scale_fill_dMod() + 
      geom_line() + #theme(legend.position = "none") +
      geom_point(data = dataP)  
      #xlab("time") + ylab("value")
    scalesP <- "fixed"
    if(input$yfree)
      scalesP <- "free_y"
    if(!is.null(facetsP)){
      if(facetsetupP=="wrap")
        P <- P + facet_wrap(facetsP, scales = scalesP)
      else{
        P <- P + facet_grid(paste0(facetsP[1], "~" , paste0(facetsP[2:length(facetsP)], collapse = "*")), scales = scalesP)
      }
    }
    if(input$errdata)
      P  <- P + geom_errorbar(data = dataP, width = 0)
    if(!is.null(model()$err)){
      if(input$errors){
        P <- P + geom_ribbon(alpha = .3, lty = 0)
      }
    }
  
    
    if(input$log)
      P <- P + scale_y_log10()
    
    vplots$plotdata <- P

  })
  
  output$plotnewly <- renderPlotly({
    p <- ggplotly(vplots$plotdata, height = input$plotHeight, autosize=TRUE)
    p$elementId <- NULL
    p
  })
  
  output$plotnewgg <- renderPlot({
    vplots$plotdata
  })
  
  output$plotnew <- renderUI({
    if(input$plotlypred)
      withSpinner(plotlyOutput("plotnewly"))
    else
      withSpinner(plotOutput("plotnewgg",height = input$plotHeight))
  })
  
  
  pstep <- reactive({
    tol <- input$tol
    P <- plotValues(model()$parameters,tol = tol)
    jumps <- attr(P,"jumps")
    jumpsPrev <- isolate({vpred$jump})
    commonjumps <- intersect(jumpsPrev, jumps)
    if(is.null(commonjumps))
      commonjumps <- c(1)
    updateSelectInput(session, "jump", choices = jumps, selected = commonjumps)
    vpred$jump <- commonjumps
    P
  })
  
  output$stepplot <- renderPlot({
    pstep() + geom_vline(xintercept = as.numeric(vpred$jump), lty = 3, size = 2, color = 2) 
  })
  
  pvals <- eventReactive(input$plotvals, {
    tol <- input$tol
    mypar <- model()$parameters
    if(input$allstepsval == "sel")
      mypar <- mypar[vpred$jump]
    plotPars(mypar,tol = tol)
  })
  
  output$parvalplot <- renderPlot({
    print(pvals())
  })
  
  ############# Profiles ###############
  vprofs <-   reactiveValues(plotProfs = TRUE, label = "plot Prediction", profval = 0)
  
  observe({
    vprofs$profval <- input$profval
  })
  
   switchTab <- observeEvent(input$select, {
    updateTabItems(session = session, inputId = "tabs", selected = "model")
    updateTabsetPanel(session = session, "plotselector", "subsetting")
  })
   
  switchTab <- observeEvent(input$select2, { #did not find how to use the same button on to menu items -> duplicated
     updateTabItems(session = session, inputId = "tabs", selected = "model")
     updateTabsetPanel(session = session, "plotselector", "subsetting")
   })
  
  switchTab <- observeEvent(input$selectflux, { # choose fluxes
    updateTabItems(session = session, inputId = "tabs", selected = "addplotsflux")
  })
  
  profilesSelect <- reactive({
    profs <- model()$profiles
    if(is.list(profs) & !("parframe" %in% class(profs)))
      profs <- profs[[vmodel$vprofssel]]
    return(profs)
  })
  
  singleprofvalues <- reactive({
    profs <- profilesSelect()
    par <- vmodel$vwhichprof
    profs <- profs[profs$whichPar %in% par]
    zeroC <- which(profs$constraint == 0)
    lengthC <- dim(profs)[1]
    profs$counter <- as.character((1-zeroC):(lengthC-zeroC))
    #values <- profs[,par]
    #ini <- values[which(profs$constraint == 0)]
    #stepsize <- (max(values) - min(values))/nrow(values)
    updateSliderInput(session = session, "profval",min = 1-zeroC, max = lengthC-zeroC, value = 0, step = 1)
    vprofs$profval <- 0
    return(profs)
  }) 
  
  singleprofplot <- reactive({
    mymodes <- vmodel$vmodeprof
    myprof <- singleprofvalues()
    attr(myprof, "obj.attributes") <- mymodes
    plotProfile(myprof)
  })
  
  profval <- reactive({
    profs <- singleprofvalues()
    pval <- profs[profs$counter==as.character(vprofs$profval),vmodel$vwhichprof]
  })
  
  profvalmin <- reactive({
    profs <- singleprofvalues()
    pval <- profs[profs$counter=="0",vmodel$vwhichprof]
  })
  
  parametersProf <- reactive({
    profs <- singleprofvalues()
    pos <- which(profs$counter==as.character(vprofs$profval)) 
    as.parvec(profs[pos, ])
  })
  
  predictionProfPlot <- reactive({
    # null prediction to compare with
    fun <- model()$x
    profs <- singleprofvalues()
    pos <-which(profs$counter=="0")
    pars <- as.parvec(profs[pos, ])
    pred <- fun(seq(from = vmodel$vtimes[1], to = vmodel$vtimes[2], length.out = 250), pars, fixed= model()$fixed, conditions = vmodel$vconditions, deriv = FALSE)
    pred <- as.data.frame(pred, data = model()$data)
    pred <- subset(pred, name %in% vmodel$vstates)
    
    colorP <- vpred$color
    dataP <- dataSubset()
    facetsP <- setdiff(vpred$facet,"step")
    facetsetupP <- vmodel$vfacetsetup
    if(colorP == "step"){
      colorP <- facetsP[1]
    }
    
    P <- ggplot(pred, aes_string(x = "time", y = "value", color = colorP, ymin = "value-sigma", ymax = "value+sigma")) + 
      theme_dMod() + scale_color_dMod() +
      geom_line(lty = 2) + 
      geom_point(data = dataP) 
      #xlab("time") + ylab("value")
    scalesP <- "fixed"
    if(input$yfree)
      scalesP <- "free_y"
    if(!is.null(facetsP)){
      if(facetsetupP=="wrap")
        P <- P + facet_wrap(facetsP, scales = scalesP)
      else{
        P <- P + facet_grid(paste0(facetsP[1], "~" , paste0(facetsP[2:length(facetsP)], collapse = "*")), scales = scalesP)
      }
    }
    if(input$errdata)
      P  <- P + geom_errorbar(data = dataP, width = 0)
    if(input$log)
      P <- P + scale_y_log10()
    return(P)
  })
  
  predictionProf <- reactive({
    fun <- model()$x
    pars <- parametersProf()
    pred <- fun(seq(from = vmodel$vtimes[1], to = vmodel$vtimes[2], length.out = 250), pars, fixed= model()$fixed, conditions = vmodel$vconditions, deriv = FALSE)
    pred <- as.data.frame(pred, data = model()$data)
    subset(pred, name %in% vmodel$vstates)
  })

  observeEvent(vmodel$vparameters,{ # reset slider values to changed vparameters
    pars <- vmodel$vparameters
    for(pp in 1:length(pars)){
      pari <- pars[pp]
      vali <- as.numeric(pari)
     updateSliderInput(session = session, as.character(names(pari)), value = vali, min = round(vali-3, digits = 2),max = round(vali+3, digits = 2))
    }
    vmodel$vparametersInd <- pars
  })
  
  observeEvent(input$plotind,{ # reset slider values to values stored in vparameters
    pars <- vmodel$vparameters
    for(pp in 1:length(pars)){
      pari <- pars[pp]
      updateSliderInput(session = session, as.character(names(pari)), value = as.numeric(pari))
    }
    vmodel$vparametersInd <- pars
  })
  
  predictionIndPlot <- reactive({
    # null prediction to compare with
    fun <- model()$x
    pars <- vmodel$vparameters
    pred <- fun(seq(from = vmodel$vtimes[1], to = vmodel$vtimes[2], length.out = 250), pars, fixed= model()$fixed, conditions = vmodel$vconditions, deriv = FALSE)
    pred <- as.data.frame(pred, data = model()$data)
    pred <- subset(pred, name %in% vmodel$vstates)
    
    colorP <- vpred$color
    dataP <- dataSubset()
    facetsP <- setdiff(vpred$facet,"step")
    facetsetupP <- vmodel$vfacetsetup
    if(colorP == "step"){
      colorP <- facetsP[1]
    }
    
    P <- ggplot(pred, aes_string(x = "time", y = "value", color = colorP, ymin = "value-sigma", ymax = "value+sigma")) + 
      theme_dMod() + scale_color_dMod() +
      geom_line(lty = 2) + 
      geom_point(data = dataP) 
    #xlab("time") + ylab("value")
    scalesP <- "fixed"
    if(input$yfree)
      scalesP <- "free_y"
    if(!is.null(facetsP)){
      if(facetsetupP=="wrap")
        P <- P + facet_wrap(facetsP, scales = scalesP)
      else{
        P <- P + facet_grid(paste0(facetsP[1], "~" , paste0(facetsP[2:length(facetsP)], collapse = "*")), scales = scalesP)
      }
    }
    if(input$errdata)
      P  <- P + geom_errorbar(data = dataP, width = 0)
    if(input$log)
      P <- P + scale_y_log10()
    return(P)
  })
  
  predictionInd <- reactive({
    fun <- model()$x
    pars <- vmodel$vparametersInd
    pred <- fun(seq(from = vmodel$vtimes[1], to = vmodel$vtimes[2], length.out = 250), pars, fixed= model()$fixed, conditions = vmodel$vconditions, deriv = FALSE)
    pred <- as.data.frame(pred, data = model()$data)
    subset(pred, name %in% vmodel$vstates)
  })
  
  observeEvent(input$exportprofpars,{ 
    vmodel$vparameters <- parametersProf()
    updateTabItems(session = session, inputId = "tabs", selected = "parvalue")
  })

  observeEvent(input$exportsteppars,{
    vmodel$vparameters <- as.parvec(model()$parameters,as.numeric(vpred$jump)[1])
    updateTabItems(session = session, inputId = "tabs", selected = "parvalue")
  })
  
  output$singleprof <- renderPlot({
    singleprofplot() + geom_vline(xintercept = as.numeric(profval()), colour = "red")+ geom_vline(xintercept = as.numeric(profvalmin()), colour = "red", lty = 2)
  })
  
  # switchPlot <- observeEvent(input$plotprof, {
  #   if(vprofs$plotProfs){
  #     vprofs$plotProfs <- FALSE
  #     vprofs$label <- "plot Profiles"
  #     updateTabsetPanel(session = session, "profileplots", "predictionTab")
  #   }else{
  #     vprofs$plotProfs <- TRUE
  #     vprofs$label <- "plot Prediction"
  #     updateTabsetPanel(session = session, "profileplots", "profilesTab")
  #   }
  #   updateActionButton(session = session, inputId = "plotprof", label = vprofs$label)
  # })
  
  observeEvent(input$plotprof, {
      updateTabsetPanel(session = session, "profileplots", "predictionTab")
  })
  
  observeEvent(input$plotprofiles,{
    updateTabsetPanel(session = session, "profileplots", "profilesTab")
  })
  
  profilesPlot <- observeEvent(input$plotprofiles, {
    profs <- model()$profiles
    if(is.list(profs) & !("parframe" %in% class(profs))){
      if(input$plotallprofs == "selected")
        profs <- profs[vmodel$vprofssel]
    }
    vplots$profsplot <- plotProfile(profs)
  })
  
  output$profplot <- renderPlot({
    vplots$profsplot
  })
  
  output$profplotui <- renderUI({
    #if(input$plotlyflux)
     # withSpinner(plotlyOutput("plotfluxesly"))
    #else
    withSpinner(plotOutput("profplot",height = input$plotHeightprof))
  })
  
  output$profpath <- renderPlot({
    plotPaths(model()$profiles, whichPar = vmodel$vwhichprof)
  })
  
  output$profpathui <- renderUI({
    #if(input$plotlyflux)
    # withSpinner(plotlyOutput("plotfluxesly"))
    #else
    withSpinner(plotOutput("profpath",height = input$plotHeightprof))
  })
  
  output$indparpred <- renderPlot({
    predictionIndPlot()  + geom_line(data = predictionInd())
  })
  
  output$profpred <- renderPlot({
    predictionProfPlot()  + geom_line(data = predictionProf())
  })
  
  output$indparpredplotui <- renderUI({
    plotOutput("indparpred",height = input$plotHeightind)
  })
  
  output$profpredui <- renderUI({
    #if(input$plotlyflux)
    # withSpinner(plotlyOutput("plotfluxesly"))
    #else
    plotOutput("profpred",height = input$plotHeightprof)
  })
  ################ Fluxes #################
  
  observe({
    rates <- model()$reactions$rates
    updateSelectInput(session = session, inputId = "whichflux", choices = rates, selected = rates[1])
  })
  
  observeEvent(input$plotflux,{
    parametersSteps <- model()$parameters[as.numeric(vpred$jump)[1]]
    pars <- as.parvec(parametersSteps)
    fluxes <- input$whichflux
    if(!(input$customflux=="")){
      fluxcustom <- input$customflux
      if(grepl(",",fluxcustom))
        fluxcustom <- unlist(strsplit(gsub(" ","", fluxcustom),","))
      fluxes <- c(fluxes, fluxcustom)
    }
    
    fun <- function(...) (model()$x)(..., conditions = vmodel$vconditions)
    P <- plotFluxes(pouter = pars, x = fun, times = seq(from = vmodel$vtimes[1], to = vmodel$vtimes[2], length.out = 500), fixed = model()$fixed, fluxEquations = fluxes)
    vplots$pf <- P
  })
  
  output$plotfluxesly <- renderPlotly({
    pp <- ggplotly(vplots$pf, height = input$plotHeightFlux, autosize=TRUE) 
    pp$elementId <- NULL
    pp
  })
  
  output$plotfluxesgg <- renderPlot({
    vplots$pf
  })
  
  output$plotfluxesui <- renderUI({
    if(input$plotlyflux)
      withSpinner(plotlyOutput("plotfluxesly"))
    else
      withSpinner(plotOutput("plotfluxesgg",height = input$plotHeightFlux))
  })
  
  ### for individual parameters tab:
  plotFluxInd <- reactive({
    pars <- vmodel$vparametersInd
    fluxes <- input$whichflux
    if(!(input$customflux=="")){
      fluxcustom <- input$customflux
      if(grepl(",",fluxcustom))
        fluxcustom <- unlist(strsplit(gsub(" ","", fluxcustom),","))
      fluxes <- c(fluxes, fluxcustom)
    }
    
    fun <- function(...) (model()$x)(..., conditions = vmodel$vconditions)
    P <- plotFluxes(pouter = pars, x = fun, times = seq(from = vmodel$vtimes[1], to = vmodel$vtimes[2], length.out = 500), fixed = model()$fixed, fluxEquations = fluxes)
    return(P)
  })
  
  output$plotFluxIndplot <- renderPlot({
    plotFluxInd()
  })
  
  output$plotfluxind <- renderUI({
    plotOutput("plotFluxIndplot",height = input$plotHeightind)
  })
  
  ################ reactions #############################
  output$reactionTable <- renderTable({
    getReactions(model()$reactions)
  })
  
  #################residuals ############################
  residuals <- eventReactive(input$plotresiduals, {
    splits <- NULL
    steps <- vpred$jump
    if(input$splitres=="time"){
      splits <- c("time","index","condition","name")
      if(length(steps) <2 )
        splits <- c("time","condition","name")
    }
    if(input$splitres=="jump")
      splits <- c("index","name","condition")
    if(input$splitres=="condition"){
      splits <- c("condition","index")
      if(length(steps) <2 )
        splits <- c("condition")
    }
    if(input$splitres=="conditionname"){
      splits <- c("condition","name","index")
      if(length(steps) <2 )
        splits <- c("condition","name")
    }
    if(input$splitres=="name"){
      splits <- c("name","index")
      if(length(steps) <2 )
        splits <- c("name")
    }
    parametersSteps <- model()$parameters[as.numeric(steps)]
    
    plotResiduals(parframe = parametersSteps,x = model()$x, data = model()$data,errmodel = model()$err, conditions = getConditions(model()$data), split = splits)
  })
  
  output$plotres <- renderPlot({
    print(residuals())
  })
  
  output$moreinfo <- renderUI({
    if(model()$ref=="none"){return()}
    list(
      actionButton(inputId='ab1', label="", 
                   icon = icon("newspaper-o"),
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                   onclick =paste0("window.open('",model()$ref,"')")
      )
    )
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)

