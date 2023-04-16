#This program has two functions named GaussianElimination and GaussJordanElimination which both take a list that contains mathematical functions as input (similar to Exercise No. 3), and perform the Gaussian and Gauss-Jordan Elimination algorithms correspondingly.
#Aranton, Andreau O.
#CMSC 150 AB3L
#10/1/21
source("ArantonEx03.R") #imports exercise 3

GaussianElimination <- function(system) { #gaussian elimination
  result1 <- AugCoeffMatrix(system) #uses the augcoeffmatrix function from exer 3 to have an aug matrix
  x <- c() #vector to store ther solution set
  AugMat <- result1$augcoeffmatrix #the aug matrix
  varUsed <- result1$variables #the variables used
  n <- length(system) #also the number of equations
  
  for(i in 1:(n-1)) {
    pivotRow <- AugMat[i,] #initially sets the pivot row to i
    tempPivEle <- pivotRow[i] #initially sets the pivot element
    rowNumber <- i #sets the row number to i

    for(j in (i+1):n) {
      if(abs(AugMat[j,i]) > abs(tempPivEle)) { #compares the absolute values
        pivotRow <- AugMat[j,] #gets the real pivot row
        rowNumber <- j #gets the row number of the pivot row
      }
    }
    
    if((AugMat[rowNumber, i]) == 0) { #if there is no solution
      print("No solution exists")
      noSol <- list(solutionSet = NA, Variables = varUsed, matrix = NA)
      return(noSol)
      break
    } else {
      #swap if there's a need to
      temp <- AugMat[i,];
      AugMat[i,] <- AugMat[rowNumber,]
      AugMat[rowNumber,] <- temp
      #print(AugMat)
      
      for(k in (i+1):n) {
        pivotElement <- AugMat[i,i] #sets the pivot element
        multiplier <- AugMat[k,i] / pivotElement #sets the multiplier
        tempVec <- AugMat[i,] * multiplier #temporary vector
        resultingVec <- AugMat[k,] - tempVec #resulting vector
        AugMat[k,] <- resultingVec #updates the matrix
      }
    }
  }
  rhs <- AugMat[,n+1] #gets the rhs
  x[n] = rhs[n] / AugMat[n,n] #computes for the value of the last variable used
  
  for(m in (n-1):1) {
    sum <- 0
    tempA <- 0
    for(o in (m+1):n) {
      tempA <- x[o] * AugMat[m, o] #multiplies the coefficient/s to the value of the variable/s
      sum <- sum + tempA #accumulate the sum
    }
    x[m] = (rhs[m] - sum) / AugMat[m,m] #computes for the value of the variable
  }
  finalList <- list(solutionSet = x, Variables = varUsed, matrix = AugMat) #final list
  return(finalList)
}

GaussJordan <- function(system) { #gauss jordan
  result1 <- AugCoeffMatrix(system) #uses the auegcoeffmatrix function from exercise 3
  AugMat <- result1$augcoeffmatrix #aug coeff matrix
  varUsed <- result1$variables #variables used
  n <- length(system) #also the number of equations
  solSet <- c() #vector for the solution set
  
  for(i in 1:n) {
    if(i != n) {
      #rowNumber <- 1
      for(k in i:(n-1)) {
        pivotRow <- AugMat[k,] #temporary sets the pivot row
        tempPivEle <- pivotRow[k] #temporary sets the pivot element
        rowNumber <- k #makes k the row number(temporarily)

        for(j in (i+1):n) {
          if(abs(AugMat[j,k]) > abs(tempPivEle)) { #compare absolute values
            pivotRow <- AugMat[j,] #gets the pivot row
            rowNumber <- j #gets the row number of the pivot row
          }
        }

        if((AugMat[rowNumber, k]) == 0) { #if there are no solution
          print("No solution exists")
          noSol <- list(solutionSet = NA, Variables = varUsed, matrix = NA)
          return(noSol)
          break
        } else {
          #swaps/reiterate the row number/pivot row
          temp <- AugMat[k,]; 
          AugMat[k,] <- AugMat[rowNumber,]
          AugMat[rowNumber,] <- temp
        }
      }
    }
    tempNorm <- AugMat[i,] / AugMat[i,i] #normalizing the pivot row
    AugMat[i,] = tempNorm
    
    for(o in 1:n) {
      if(i == o) { #whenever it is the diagonal
        next
      }
      multiplier <- AugMat[o,i] #sets the ,multiplier
      tempVec <- AugMat[i,] * multiplier #computes for the temporary vector
      resultingVec <- AugMat[o,] - tempVec #computes for the resutltng vector
      AugMat[o,] <- resultingVec #updates the matrix
    }
  }
  for(p in 1:n) { #gets the solution set
    solSet[p] = AugMat[p,n+1]
  }
  finalList <- list(solutionSet = solSet, Variables = varUsed, matrix = AugMat) #creates the final list
  return(finalList)
}
#given
E1 <- function (x1, x2, x3) 0.3 * x1 + -0.2 * x2 + 10 * x3 + -71.4
E2 <- function (x1, x2, x3) 3 * x1 + -0.2 * x3 + -0.1 * x2 + -7.85
E3 <- function (x1, x2, x3) 0.1 * x1 + 7 * x2 + -0.3 * x3 + 19.3
