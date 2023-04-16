#This program has an R function named AugCoeffMatrix which takes a list which contains mathematical functions as input. It creates the augmented coefficient matrix for the functions and return it in a single list variable.
#Aranton, Andreau O.
#CMSC 150 AB3L
#09/21/21

AugCoeffMatrix <- function(system) {
  numVars <- 0 #for the number of variables used
  varUsed <- c() #vector of the variables
  rhs <- c() #stores the rhs
  counter1 <- 1 #counters
  counter2 <- 1
  tempVec <- c() #temporary vector
  allVar <- c() #all variables that will be used in swapping
  coeff <- c() #vector that would be used in storing the coefficients
  tempVec2 <- c() #temp vector
  rowNms <- c() #vector for the row names
  colNms <- c() #vector for the col names
  for(i in 1:(length(system))) { #for every eq
    n <- length(system) #in is the number of eqs
    tempVec3 <- c() #temp vector
    tempVec4 <- c() #temp vector
    
    A <- deparse1(system[[i]]) #deparses the eq
    #print(A)
    #A <- A[1] #stores the variables together with other stuff
    A <- gsub("function ", "", A) #deletes some unnecesarry stuff
    A <- gsub("\\)", "", A)
    A <- gsub("\\(", "", A)
    #A <- gsub(" ", "", A)
    A <- strsplit(A, "  ") #separates the variabables used
    #print(A)
    B <- A[[1]][1]
    #print(B)
    B <- gsub(",", "", B)
    #print(B)
    B <- strsplit(B, " ")
    #print(length(B[[1]]))
    if(i == 1) { #in the first iteration
      numVars <- length(B[[1]]) #counts the number of variables
      for(z in 1:numVars) {
        varUsed[z] <- B[[1]][z] #gets the variables used
      }
    } else {
      if(numVars != length(B[[1]])) { #if the number of variables on an eq is not equal to the num of variables of another eq
        return(NA) #return NA
      }
    }
    C <- A[[1]][2] #the other part of the eq
    #print(C)
    C <- strsplit(C, " \\+ ")
    for(j in 1:(n+1)) {
      C[[1]][j] <- gsub(" ", "", C[[1]][j])
    }
    #print(C)
    rhs[counter1] <- C[[1]][n+1] #stores the rhs
    #print(rhs)
    counter1 <- counter1 + 1
    for(k in 1:n) {
      tempVec[k] <- C[[1]][k]
      tempVec[k] <- gsub("\\*", " ", tempVec[k])
      tempVec2[k] <- strsplit(tempVec[k], " ")
      tempVec3[k] <- tempVec2[[k]][1] #coeff
      tempVec4[k] <- tempVec2[[k]][2] #var
      
    }
    for(i in 1:n) {
      coeff[counter2] <- tempVec3[i] #stores all the coeff into an array
      counter2 <- counter2 + 1
    }
  }
  coeff <- as.numeric(coeff) #makes the coefficients numerical
  matA = matrix(coeff, nrow = n, byrow = TRUE) #creates a matrix of coefficients
  #print(rhs)
  rhs <- as.numeric(rhs) * (-1) #makes the rhs numerical and flips the sign
  matRhs <- matrix(rhs) #converts the rhs to a matrix
  matAug = cbind(matA, matRhs) #binds the rhs to the coef matrix
  for(y in 1:numVars) { #for the rownames
    rowNms[y] = y
  }
  rownames(matAug) <- rowNms #updates the matrix's row names
  colNms <- append(varUsed, "RHS") #for the col names
  colnames(matAug) <- colNms #updates the colnames
  finalList <- list(variables = varUsed, augcoeffmatrix = matAug) #creates a list which contains a vector of the var used and the Augmented matrix
  return(finalList)
}

E1 <- function (x1, x2, x3) 0.3 * x1 + -0.2 * x2 + 10 * x3 + -71.4
E2 <- function (x1, x2, x3) 3 * x1 + -0.2 * x3 + -0.1 * x2 + -7.85
E3 <- function (x1, x2, x3) 0.1 * x1 + 7 * x2 + -0.3 * x3 + 19.3

G1 <- function (a0, a1, a2, a3) 31 * a0 + 32275 * a1 + 40758893 * a2 + 57390772057 * a3 + -6010
G2 <- function (a0, a1, a2, a3) 32275 * a0 + 40758893 * a1 + 57390772057 * a2 + 86610992822897 * a3 + -6897861
G3 <- function (a0, a1, a2, a3) 40758893 * a0 + 57390772057 * a1 + 86610992822897 * a2 + 137211629709233072 * a3 + -9005235271
G4 <- function (a0, a1, a2, a3) 57390772057 * a0 + 86610992822897 * a1 + 137211629709233072 * a2 + 2.25332881428815e+20 * a3 + -12805253783307
