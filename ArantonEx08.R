#This program has a method of finding the polynomials on each interval given a data set using quadratic spline interpolation.
#Aranton, Andreau O.
#CMSC 150 AB3L
#10/31/21

source("ArantonEx04.R") #imports the gauss jordan (the augcoeff method is imported on the imported file)

poly.qsi <- function(data, x) {
  xs <- data[[1]] #vector of the x's
  ys <- data [[2]] #vector of the y's
  numDataPts <- length(xs) #number of data points
  numIntervals <- numDataPts -1 #number of intervals
  coeffV <- c(1:((3 * numIntervals)+1)) * 0 #creates a vector filled with zeroes
  coeffL <- list() #this will store the list of coefficient vectors
  matFunc <- list() #this is a list that will store the function to be used in gauss jordan
  qsiF <- list() #this is the list where the qsi function per interval will be stored
  counter <- 1 #counter1
  counter2 <- 2 #counter 2
  numVec <- 1
  for(i in 2: numIntervals) { #for the condition 1, for i=2 to n
    for(j in 1:2) {
      for(k in 1:3) {
        if(k == 1) {
          coeffV[counter] <- (xs[counter2])^2 #squares the x sub i minus 1
        }
        if(k == 2) {
          coeffV[counter] <- xs[counter2] #x sub i minus 1
        }
        if(k == 3) {
          coeffV[counter] <- 1 #coefficient of C
        }
        counter <- counter + 1 #increments the counter
      }
      coeffV[(3 * numIntervals) + 1] = ys[counter2] #this is the f(X sub i-1)
      coeffL <- append(coeffL, list(coeffV)) #stores it to the list of coefficient vectors
      coeffV <- c(1:((3 * numIntervals)+1)) * 0 #resets the vector
    }
    counter <- counter - 3 #decrements the counter by 3
    counter2 <- counter2 + 1 #increments the counter
  }
  counter <- 1 #sets the counter
  for(m in 1:3) { #for the condition 2 first function
    if(m == 1) {
      coeffV[counter] = (xs[1])^2 #square of x sub 0
      counter <- counter + 1
    }
    if(m == 2) {
      coeffV[counter] = xs[1] #x sub 0
      counter <- counter + 1
    }
    if(m == 3) {
      coeffV[counter] = 1 #coefficient of C
      counter <- counter + 1
    }
  }
  coeffV[(3 * numIntervals) + 1] <- ys[1] #this is the f(x sub 0)
  coeffL <- append(coeffL, list(coeffV)) #adds it to the list
  coeffV <- c(1:((3 * numIntervals)+1)) * 0 #resets the coefficient vector
  counter <- (3 * numIntervals) - 2 #sets the counter
  for(p in 1:3) { #for the condition 2, last function
    if(p == 1) {
      coeffV[counter] = (xs[numDataPts])^2 #square of x sub n
      counter <- counter + 1
    }
    if(p == 2) {
      coeffV[counter] = xs[numDataPts] #x sub n
      counter <- counter + 1
    }
    if(p == 3) {
      coeffV[counter] = 1 #coefficient of c
      counter <- counter + 1
    }
  }
  coeffV[(3 * numIntervals) + 1] <- ys[numDataPts] #stores the f(x sub n)
  coeffL <- append(coeffL, list(coeffV)) #adds it to the list
  counter <- 1 #sets the counter
  counter2 <- 2 #sets the counter 2
  for(q in 2: numIntervals) { #for the condition 3, i=2 to n
    coeffV <- c(1:((3 * numIntervals)+1)) * 0 #sets the vector
    coeffV[counter] <- 2 * xs[q] #2 times x sub i minus 1
    counter <- counter + 3
    coeffV[counter2] <- 1 #b sub i minus 1
    counter2 <- counter2 + 3
    coeffV[counter] <- -2 * xs[q] #negative 2 times x sub i minus 1
    coeffV[counter2] <- -1 #negative b sub i minus 1
    coeffL <- append(coeffL, list(coeffV)) #append it to the list
  }
  for(r in 1:length(coeffL)) { #for the formation of the string for the augmatrix/gauss jordan
    stringFunc <- "function (" #this is for the first sub, because a sub 0 is removed
    stringFunc <- paste(stringFunc, "b", sep = "") #b
    stringFunc <- paste(stringFunc, "1", sep = "")
    stringFunc <- paste(stringFunc, ", ", sep = "")
    stringFunc <- paste(stringFunc, "c", sep = "") #c
    stringFunc <- paste(stringFunc, "1", sep = "")
    stringFunc <- paste(stringFunc, ", ", sep = "")
    
    for(s in 2:numIntervals) { #from a sub 1
      stringFunc <- paste(stringFunc, "a", sep = "") #a
      stringFunc <- paste(stringFunc, s, sep = "")
      stringFunc <- paste(stringFunc, ", ", sep = "")
      stringFunc <- paste(stringFunc, "b", sep = "") #b
      stringFunc <- paste(stringFunc, s, sep = "")
      stringFunc <- paste(stringFunc, ", ", sep = "")
      stringFunc <- paste(stringFunc, "c", sep = "") #c
      stringFunc <- paste(stringFunc, s, sep = "")
      if(s != numIntervals) { #so that it will not place a comma at the end
        stringFunc <- paste(stringFunc, ", ", sep = "")
      }
    }
    stringFunc <- paste(stringFunc, ") ", sep = "") #)
    stringFunc <- paste(stringFunc, coeffL[[r]][2], sep = "") #initial, no a sub1
    stringFunc <- paste(stringFunc, " * ", sep = "")
    stringFunc <- paste(stringFunc, "b1", sep = "") #b1
    stringFunc <- paste(stringFunc, " + ", sep = "")
    stringFunc <- paste(stringFunc, coeffL[[r]][3], sep = "")
    stringFunc <- paste(stringFunc, " * ", sep = "")
    stringFunc <- paste(stringFunc, "c1", sep = "") #c1
    stringFunc <- paste(stringFunc, " + ", sep = "")

    counter <- 4 #sets the counter
    counter2 <- 2 #sets the counter2
    for(t in 1:(numIntervals - 1)) { #completes the string
      stringFunc <- paste(stringFunc, coeffL[[r]][counter], sep = "")
      counter <- counter + 1
      stringFunc <- paste(stringFunc, " * ", sep = "")
      stringFunc <- paste(stringFunc, "a", sep = "") #a
      stringFunc <- paste(stringFunc, counter2, sep = "")
      stringFunc <- paste(stringFunc, " + ", sep = "")
      stringFunc <- paste(stringFunc, coeffL[[r]][counter], sep = "")
      counter <- counter + 1
      stringFunc <- paste(stringFunc, " * ", sep = "")
      stringFunc <- paste(stringFunc, "b", sep = "") #b
      stringFunc <- paste(stringFunc, counter2, sep = "")
      stringFunc <- paste(stringFunc, " + ", sep = "")
      stringFunc <- paste(stringFunc, coeffL[[r]][counter], sep = "")
      counter <- counter + 1
      stringFunc <- paste(stringFunc, " * ", sep = "")
      stringFunc <- paste(stringFunc, "c", sep = "") #c
      stringFunc <- paste(stringFunc, counter2, sep = "")
      stringFunc <- paste(stringFunc, " + ", sep = "")
      counter2 <- counter2 +1
    }
    stringFunc <- paste(stringFunc, " -", sep = "") #paste - for the to be rhs
    stringFunc <- paste(stringFunc, coeffL[[r]][counter], sep = "") #paste the last value
    evalMatFunc <- eval(parse(text = stringFunc)) #makes it a function
    matFunc <- append(matFunc, evalMatFunc) #adds iut to the list of functions
  }
  getValues <- GaussJordan(matFunc) #apply gauss jordan to know the values
  variableValues <-getValues$solutionSet #get the solution set
  qsiStrFunc <- "function (x) " #string to be used for the function making
  qsiStrFunc <- paste(qsiStrFunc, variableValues[1], sep = "") #paste the b sub 1
  qsiStrFunc <- paste(qsiStrFunc, " * x + ", sep = "")
  qsiStrFunc <- paste(qsiStrFunc, variableValues[2], sep = "") #paste c
  qsiEvalFunc <- eval(parse(text = qsiStrFunc)) #the first func is done and is converted to a function
  qsiF <- append(qsiF, qsiEvalFunc) #adds it to the list
  counter <- 3 #sets the counter
  for(u in 2:numIntervals) { #for the creation of eq for other interval
    qsiStrFunc <- "function (x) " #resets the string
    qsiStrFunc <- paste(qsiStrFunc, variableValues[counter], sep = "") #paste the a
    counter <- counter + 1
    qsiStrFunc <- paste(qsiStrFunc, " * ", sep = "")
    qsiStrFunc <- paste(qsiStrFunc, "x^2", sep = "")
    qsiStrFunc <- paste(qsiStrFunc, " + ", sep = "")
    qsiStrFunc <- paste(qsiStrFunc, variableValues[counter], sep = "") #paste the b
    counter <- counter + 1
    qsiStrFunc <- paste(qsiStrFunc, " * ", sep = "")
    qsiStrFunc <- paste(qsiStrFunc, "x", sep = "")
    qsiStrFunc <- paste(qsiStrFunc, " + ", sep = "")
    qsiStrFunc <- paste(qsiStrFunc, variableValues[counter], sep = "") #paste the c
    counter <- counter + 1
    qsiEvalFunc <- eval(parse(text = qsiStrFunc)) #converts it to functioj
    qsiF <- append(qsiF, qsiEvalFunc) #add it to the list
  }
  eqToUse <- 0 #variable to hold which eq is to be used
  for(v in 1:(numDataPts - 1)) { #for the comparing
    if(xs[v] <= x && x <= xs[v+1]) { #if x is between the 2 data points
      eqToUse <- v #v is the equation to be used
      break
    }
  }
  qsiAns <- qsiF[[eqToUse]](x) #computes for the value of y given the x
  finalList <- list(qsi.fxns = qsiF, y = qsiAns) #sets up the list to be returned
  return(finalList) #returns the list
}

xval <- c(3, 4.5, 7, 9) #x values
yval <- c(2.5, 1, 2.5, 0.5) #corresponding y values
data <- list(xval, yval) #stores the vectors to the list
try <- poly.qsi(data, 5) #performs the qsi function
try #displays the result