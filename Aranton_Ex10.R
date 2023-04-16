#This program has a graphical user interface for both exer8 and exer9
#Aranton, Andreau O.
#CMSC 150 AB3L
#11/26/21

library(shiny)
#install.packages("shinythemes")
library(shinythemes)
#install.packages("rhandsontable")
library(rhandsontable)
#install.packages("shinyalert")
library(shinyalert)

source("ArantonEx08.R") #imports exer8
source("ArantonEx09.R") #imports exer9

initialMat <- problem1 #problem1 is imported from the exer9 file
initialObjFunc <- "10x1 + 8x2 + 6x3 + 5x4 + 4x5 + 6x6 + 5x7 + 4x8 + 3x9 + 6x10 + 3x11 + 4x12 + 5x13 + 5x14 + 9x15"
initialCons <- "x1 + x2 + x3 + x4 + x5 <= 310\nx6 + x7 + x8 + x9 + x10 <= 260\nx11 + x12 + x13 + x14 + x15 <= 280\nx1 + x6 + x11 >= 180\nx2 + x7 + x12 >= 80\nx3 + x8 + x13 >= 200\nx4 + x9 + x14 >= 160\nx5 + x10 + x15 >= 220"
ui <- fluidPage(theme = shinytheme("yeti"), #uses a shiny theme
  useShinyalert(), #uses shinyalert
  navbarPage( #navigation bar
    "Aranton_Exer10",
    tabPanel("Quadratic Spline Interpolation", #for the exer 8
             sidebarPanel(
               h2("Quadratic Spline Interpolation"), #text
               br(),
               h4("Data:"),
               rHandsontableOutput("test_table"), #outputs the table
               br(),
               numericInput("comp",h4("x:"), 5), #numeric input
               actionButton("compute", "compute") #compute button
               
             ),
             mainPanel( #main panel
               h3(textOutput("currentX")), #when x is
               br(),
               textOutput("qsiTitle"), #qsi text
               verbatimTextOutput("qsiFunctions"), #qsi functions
               br(),
               textOutput("yTitle"), #y text
               verbatimTextOutput("yAns") #y value
             )
    ),
    tabPanel("Simplex Method", #exer9
             h2("Simplex Method Solver"),
             hr(),
             h3("Inputs:"),
             mainPanel(
               textInput("objFunc", "Objective Function/Z =", value = initialObjFunc, width = "1000px"),
               textAreaInput("constraints", "Constraints", initialCons, width = "1000px", height = "300px"),
               checkboxInput("max", "Maximize (Uncheck to Minimize)", FALSE), #checkbox for isMax
               checkboxInput("shipNum", "Show shipping num (Must only be checked if the problem is the Exer09 Lab Problem)", TRUE), #checkbox for problem
               actionButton("solve", "solve"), #solve button
               br(),
               hr(),
               hr(),
               verbatimTextOutput("sol"), #solution text
               tags$head(tags$style("#sol{font-size:30px;}")), #adjust text size
               br(),
               verbatimTextOutput("finalT"), #final tableau text
               tags$head(tags$style("#finalT{font-size:20px;}")), #adjust text size
               tableOutput("fTableau"), #final tableau
               br(),
               verbatimTextOutput("basicS"), #basic solution text 
               tags$head(tags$style("#basicS{font-size:20px;}")), #adjust text size
               tableOutput("basicSol"), #basic solution matrix
               br(),
               verbatimTextOutput("optimalV"), #optimal value text
               tags$head(tags$style("#optimalV{font-size:20px;}")), #adjust text size
               br(),
               verbatimTextOutput("shippingN"), #shipping num text
               tags$head(tags$style("#shippingN{font-size:20px;}")), #adjust text size
               tableOutput("shippingNum") #shipping num matrix
             )
    )
  )
)

server <- shinyServer(function(input, output, session) {
  dataf = data.frame( #data frame for QSI
    x = c(3,4.5,7,9), #lec sample
    y = c(2.5,1,2.5,0.5) #lec sample
  )

  output$test_table <- renderRHandsontable({ #shows the data from as rhandsontable
    rhandsontable(dataf, width = 300)
  })
  
  
  values <- reactiveValues(data = NULL)
  
  observe({
    values$data <- hot_to_r(input$test_table) #values$data is the dataframe (every time rhandsontable changes, the values$data gets updated)
  })
  
  observeEvent(input$compute, { #when the compute button is clicked
    toComp <- input$comp #toComp is the input x value
    
    if(is.na(toComp) || toComp < values$data[1,1] || toComp > values$data[nrow(values$data), 1]) { #if the input x is null or if it is out of bounds
      output$qsiTitle <- renderText({ #blank
        ""
      })
      output$qsiFunctions <- renderText({ #blank
        ""
      })
      output$yTitle <- renderText({ #blank
        ""
      })
      output$yAns <- renderText({ #blank
        ""
      })
      output$currentX <- renderText({ #blank
        ""
      })
      if(is.na(toComp)) { #if the input x is null
        shinyalert("Oops!","x must not be blank!", type="error") #error alert
      } else { #if it is out of bounds
        shinyalert("Oops!","The input x is not within the possible range of x values!", type="error") #error alert
      }
    } else { #if rthere are no errors in the input x
      xVec <- values$data[,1] #stores the x values
      yVec <- values$data[,2] #stores the y values
      inputD <- list(xVec, yVec) #puts them into a liist
      ans <- poly.qsi(inputD, toComp) #apply the poly qsi method
      stringQsiF <- "" #string
      for(i in 1:length(ans[[1]])) { #for the qsi functions, paste the answers into the string
        ans[[1]][[i]] <- deparse1(ans[[1]][[i]])
        stringQsiF <- paste(stringQsiF, "\t", sep = "")
        stringQsiF <- paste(stringQsiF, i, sep = "")
        stringQsiF <- paste(stringQsiF, ". ", sep = "")
        stringQsiF <- paste(stringQsiF, ans[[1]][[i]], sep = "")
        stringQsiF <- paste(stringQsiF, "\n", sep = "")
      }
      stringcurrentX <- "When x = " #string
      stringcurrentX <- paste(stringcurrentX, toComp, sep = "") #paste the input x
      stringcurrentX <- paste(stringcurrentX, ":", sep = "")
      output$currentX <- renderText({ #show the string
        stringcurrentX
      })
      output$qsiTitle <- renderText({ #show the qsi functions
        "QSI Functions:"
      })
      output$qsiFunctions <- renderText({ #show the qsi functions
        stringQsiF
      })
      ansY <- "\t" #string
      ansY <- paste(ansY, ans[[2]], sep = "") #paste the value of y
      output$yTitle <- renderText({ #show the y
        "Value of y:"
      })
      output$yAns <- renderText({ #show the value of y
        ansY
      })
    }
    
  })

  observeEvent(input$solve, { #if solve button is clicked
    isMax <- input$max #true if maximization, false if minimization
    isShipNum <- input$shipNum #true if it is checked, false if not
    objF <- input$objFunc #objective function
    objF <- gsub(" \\+ ", " ", objF) #removes "+ "
    objFF <- strsplit(objF, " ") #splits the string
    lenObj <- length(objFF[[1]]) #to get the number of variables
    matL <- list()
    matSolve <- c()
    matCons <- input$constraints #constraints
    sepMatCons <- strsplit(matCons, "\n") #splits the string
    numVars <- lenObj #number of variables
    lenCons <- length(sepMatCons[[1]]) #length of constraint

    for(a in 1:lenCons) {
      zeroVec <- c(1:(numVars+1)) * 0 #a zero vector
      sepMatCons[[1]][a] <- gsub(" \\+ ", " ", sepMatCons[[1]][a])
      currentCons <- strsplit(sepMatCons[[1]][a], " ") #splits
      lenCurCons <- length(currentCons[[1]])
      for(b in 1:(lenCurCons-2)) {
        sepX <- strsplit(currentCons[[1]][b], "x") #splits
        if(sepX[[1]][1] == "") {
          zeroVec[as.numeric(sepX[[1]][2])] <- 1 #if there is no coefficient, it is 1
        } else {
          zeroVec[as.numeric(sepX[[1]][2])] <- as.numeric(sepX[[1]][1]) #stores it
        }
        
      }
      
      lastElem <- as.numeric(currentCons[[1]][lenCurCons]) #last element
      zeroVec[numVars+1] <- lastElem #writes the last element
      if(isMax == TRUE) { #if it is maximization
        if(currentCons[[1]][lenCurCons-1] == ">=") { #if is is >=
          zeroVec <- zeroVec * -1 #multiply the vector by -1
        }
      } else { #if it is minimization
        if(currentCons[[1]][lenCurCons-1] == "<=") { #if it is <=
          zeroVec <- zeroVec * -1 #fmultiply the vector by -1
        }
      }
      matL <- append(matL, list(zeroVec)) #stores it in a list
      
      
    }
  
    zeroVec <- c(1:(numVars+1)) * 0 #a zero vector
    for(c in 1:lenObj) {
      sepX <- strsplit(objFF[[1]][c], "x") #splits
      zeroVec[as.numeric(sepX[[1]][2])] <- as.numeric(sepX[[1]][1]) #stores
    }
    zeroVec[numVars+1] <- 1 #assign
    matL <- append(matL, list(zeroVec)) #appends the list
    for(d in 1:length(matL)) {
      matSolve <- rbind(matSolve,matL[[d]]) #binds to form a matrix
    }
    
    if(isMax == FALSE) { #if it is minimization
      matSolve <- setupTableauTestCase(matSolve) #creates the initial tableau
    } else { #if it is maximization
      matSolve <- setupMax(matSolve) #creates the initial tableau
    }
    
    simplexAns <- simplex(matSolve, isMax, isShipNum) #computes via the simplex method function
    if(simplexAns$opt.val == "No feasible solution") { #if there is no feasible solution
      output$sol <- renderText({ #show the text that rhere is no feasible solution
        "No Feasible solution"
      })
      output$finalT <- renderText({ #blank
        ""
      })
      output$fTableau <- renderTable({ #blank
        
        
      })
      output$basicS <- renderText({ #blank
        ""
      })
      output$basicSol <- renderTable({ #blank
        
      })
      output$optimalV <- renderText({ #blank
        ""
      })
      output$shippingN <- renderText({ #blank
        ""
      })
      output$shippingNum <- renderTable({ #blank
        
      })
    } else { #if there is a feasible solution
      dfMat <- as.data.frame(simplexAns$final.tableau) #converts the final tableau matrix into a data frame
      output$sol <- renderText({ #show text
        "Solution"
      })
      output$finalT <- renderText({ #show text
        "Final Tableau:"
      })
      output$fTableau <- renderTable({ #show the final tableau
        dfMat
        
      },rownames = TRUE)
      output$basicS <- renderText({ #show text
        "Basic Solution:"
      })
      bSolM <- as.data.frame(simplexAns$basic.solution) #convert the basic solution matrix into data frame
      output$basicSol <- renderTable({ #show the basic solution
        bSolM
        
      },rownames = TRUE)
      optV <- "" #blank string
      if(isMax == TRUE) { #if maximization
        optV <- "Maximum Value: " #string
      } else {
        optV <- "Minimum Value: " #string
      }
      optV <- paste(optV, simplexAns$opt.val, sep = "") #paste the optimal value
      output$optimalV <- renderText({ #show optimal value
        optV
      })
      if(isShipNum == TRUE) { #if shipnum checkbox is checked
        if(isMax == FALSE) {
          if(numVars == 15) {
            output$shippingN <- renderText({  #show text
              "Number of items shipped from plants to warehouses:"
            })
            shipM <- as.data.frame(simplexAns$shipping.num) #convert the matrix to a dataframe
            output$shippingNum <- renderTable({ #show the table
              shipM
              
            },rownames = TRUE)
          }
          
        } else {
          output$shippingN <- renderText({  #show text
            ""
          })
        }
      } else { #if it is not checked
        output$shippingN <- renderText({ #blank
          ""
        })
        output$shippingNum <- renderTable({ #blank
          
        })
      }
    }
  })
})

shinyApp(ui, server)