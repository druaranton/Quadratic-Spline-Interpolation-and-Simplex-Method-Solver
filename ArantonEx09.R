#This program has generic function to implement the simplex method for maximization and minimization. Note that the input for the simplex method should already be ready for maximization
#Aranton, Andreau O.
#CMSC 150 AB3L
#11/20/21

setupMax <- function(testcase) {
  maxCase <- testcase
  numS <- ncol(maxCase)-1
  lastcol <- maxCase[,ncol(maxCase)]
  maxCase <- maxCase[1:nrow(maxCase),1:(ncol(maxCase)-1)]
  #print(maxCase)
  for(j in 1:(nrow(maxCase)-1)) {
    zVec <- c(1:nrow(maxCase)) * 0
    zVec[j] <- 1
    maxCase <- cbind(maxCase,zVec)
  }
  #print(maxCase)
  Z <- c(1:nrow(maxCase)) * 0
  Z[nrow(maxCase)] <- -1
  maxCase <- cbind(maxCase,Z)
  maxCase <- cbind(maxCase,lastcol)
  #print(maxCase)
  maxCase[nrow(maxCase),] <- maxCase[nrow(maxCase),] * -1
  cNms <- c() #for the column names
  counter <- 1 #counter
  #print(maxCase)
  for(k in 1:ncol(maxCase)) { #for the column names
    stringS <- "S"
    stringX <- "X"
    if(k <= numS) {
      stringX <- paste(stringX,k,sep = "")
      cNms <- append(cNms,stringX)
    } else {
      if(k == ncol(maxCase)-1) {
        cNms <- append(cNms,"Z")
      } else {
        if(k == ncol(maxCase)) {
          cNms <- append(cNms,"Solution")
        } else {
          stringS <- paste(stringS,counter,sep = "")
          cNms <- append(cNms,stringS)
          counter <- counter + 1
        }
      }
    }
  }
  maxCase[nrow(maxCase),ncol(maxCase)] <- 0 #assigns the value
  colnames(maxCase) <- cNms #sets the column names of the matrix
  rowNms <- c() #for the row names
  for(w in 1:(nrow(maxCase)-1)) {
    rowNms <- append(rowNms, paste("constraint", w, sep = " "))
  }
  rowNms <- append(rowNms, "Objective function")
  rownames(maxCase) <- rowNms
  return(maxCase)
}

setupTableauTestCase <- function(testcase) { #this function is for converting a matrix into the initial tableau of the dual problem. Note that this is only intended for the exercise problem. Hence, it may not work for other problem/s
  tmat <- t(testcase) #transposes the matrix
  numS <- ncol(tmat)-1 #counts the S (the variables of the dual problem)
  lastCol <- tmat[,ncol(tmat)] #gets the last column
  tmat <- tmat[1:nrow(tmat), 1:(ncol(tmat)-1)] #removes the last col of the matrix
  for(j in 1:(nrow(tmat)-1)) { #for the slack variables of the dual problem
    zVec <- c(1:nrow(tmat)) *0
    zVec[j] <- 1
    tmat <- cbind(tmat,zVec) #binds it to the matrix
  }
  Z <- c(1:nrow(tmat)) *0 #creates a vector for the Z
  Z[nrow(tmat)] <- -1 #places a value
  tmat <- cbind(tmat,Z) #bind
  tmat <- cbind(tmat,lastCol) #bind
  tmat[nrow(tmat),] <- tmat[nrow(tmat),] * -1 #flips the sign of the last row of the matrix
  cNms <- c() #for the column names
  counter <- 1 #counter
  for(k in 1:ncol(tmat)) { #for the column names
    stringS <- "S"
    stringX <- "X"
    if(k <= numS) {
      stringS <- paste(stringS,k,sep = "")
      cNms <- append(cNms,stringS)
    } else {
      if(k == ncol(tmat)-1) {
        cNms <- append(cNms,"Z")
      } else {
        if(k == ncol(tmat)) {
          cNms <- append(cNms,"Solution")
        } else {
          stringX <- paste(stringX,counter,sep = "")
          cNms <- append(cNms,stringX)
          counter <- counter + 1
        }
      }
    }
  }
  tmat[nrow(tmat),ncol(tmat)] <- 0 #assigns the value
  colnames(tmat) <- cNms #sets the column names of the matrix
  rowNms <- c() #for the row names
  for(w in 1:(nrow(tmat)-1)) {
    rowNms <- append(rowNms, paste("constraint", w, sep = " "))
  }
  rowNms <- append(rowNms, "Objective function")
  rownames(tmat) <- rowNms
  
  return(tmat) #returns the matrix which is ready for maximization
}

checkLROW <- function(tableau) { #function to check if the last row of the matrix still has negative nalue/s
  isExit <- TRUE #initially sets it to true
  for(i in 1:(ncol(tableau)-1)) {
    if(tableau[nrow(tableau),i] < 0) { #if there is a negative value
      isExit <- FALSE #sets to false
      break
    }
  }
  return(isExit) #return
}

simplex <- function(tableau, isMax, problem) { #function of simplex method
  exit <- checkLROW(tableau) #checks if it is already done
  
  while(exit == FALSE) { #while there is a negative value in the last row
    pivotCol <- 1 #initializes the temporary pivot column
    pivotElement <- tableau[1,pivotCol] #temporary pivot element
    
    for(i in 1:(ncol(tableau)-1)) { #this will find the pivot column
      if(tableau[nrow(tableau), i] < tableau[nrow(tableau), pivotCol]) { #if the ith column is lesser than the currently assigned pivot col
        pivotCol <- i #updates the pivot column
      }
    }
    pivotRow <- 1 #temporary pivot row
    pivotElement <- tableau[pivotRow,pivotCol] #tmporary pivot element
    for(z in 1:(nrow(tableau))) {#this will try to find a pivot element which is a positive number (temporarily)
      if(tableau[z,pivotCol] > 0) {
        pivotElement <- tableau[z,pivotCol] #updates the temprary pivot element
        pivotRow <- z #updates the temporary pivot row
        break
      }
    }
    for(j in 1:(nrow(tableau)-1)) { #this will find the true pivot element in this iteration
      if((tableau[j,ncol(tableau)] == 0) && (tableau[j,pivotCol] == 0)) { #if candidate element is zero, skip
        next
      }
      if((tableau[j,ncol(tableau)]/tableau[j,pivotCol]) <= 0) { #if the test ratio is less than or equal to zero, skip
        next
      }
      if(((tableau[j,ncol(tableau)]/tableau[j,pivotCol]) < (tableau[pivotRow,ncol(tableau)]/tableau[pivotRow,pivotCol]))) { #if the tr is lower than the currently assigned's tr
        pivotRow <- j #updates the pivot row
        pivotElement <- tableau[pivotRow,pivotCol] #updates the pivot element
      }
    }
    pivotRowElements <- tableau[pivotRow,] #this is the pivot row elements
    if((tableau[pivotRow,ncol(tableau)]/tableau[pivotRow,pivotCol]) <= 0) { #if the test ratio is zero or negative,, happens when the initial pivot row has not changed and that it's tr is zero or negative
      if(isMax == TRUE) {
        finalList <- list(final.tableau = "No feasible solution", basic.solution = "No feasible solution", opt.val = "No feasible solution")
        return(finalList) #return no feasible solution
      } else {
        if(problem == FALSE) {
          finalList <- list(final.tableau = "No feasible solution", basic.solution = "No feasible solution", opt.val = "No feasible solution")
          return(finalList) #return no feasible solution
        } else {
          finalList <- list(final.tableau = "No feasible solution", basic.solution = "No feasible solution", opt.val = "No feasible solution", shipping.num = "No feasible solution")
          return(finalList) #return no feasible solution
        }
      }
    }
    normalizedPivotRow <- pivotRowElements/pivotElement #gets the normalized pivot row
    tableau[pivotRow,] <- normalizedPivotRow #updates the matrix
    for(k in 1:nrow(tableau)) {
      if(k == pivotRow) { #if the kth row is the pivot row, skip
        next
      }
      currentR <- tableau[k,] #row elements of the row
      currentE <- currentR[pivotCol] #the C element of the row
      currentR <- currentR - (normalizedPivotRow * currentE) #gets the new row values for the matrix
      tableau[k,] <- currentR #updates the matrix
      
    }
    exit <- checkLROW(tableau) #checks if there is still negative value/s in the last row
  }
  if(isMax == TRUE) { #if it is maximization
    solVec <- c() # the solutions vector
    for(t in 1:(ncol(tableau)-1)) { #for the calculation of the basic sol values
      countZ <- 0 #counter
      locRow <- 0 #location of row
      locCol <- 0 #location of column
      for(v in 1:nrow(tableau)) {
        if(tableau[v,t] != 0) { #if the element is not zero
          countZ <- countZ + 1 #increments the counter
          locRow <- v #gets the row
          locCol <- t #gets the col
        }
      }
      if(countZ == 1) { #if the non zero is only one
        val <- tableau[locRow,ncol(tableau)] / tableau[locRow,locCol] #divides the value in its last col to the element
        solVec <- append(solVec,val) #stores it to the solution vector
      } else { #if there is more than one non-zero
        solVec <- append(solVec, 0) #puts zero
      }
    }
    solMat <- matrix(solVec, nrow = 1) #creates a matrix of the soluion
    colnames(solMat) <- colnames(tableau)[1:(ncol(tableau)-1)] #sets the column names
    rownames(solMat) <- "Solution" #sets row name
    finalList <- list(final.tableau = tableau, basic.solution = solMat, opt.val = tableau[nrow(tableau),ncol(tableau)]) #list to be returned
    return(finalList) #returns the list
    
  } else { #if it is minimization
    solMat <- matrix(tableau[nrow(tableau),1:(ncol(tableau)-1)], nrow = 1) #creates a solution matrix
    solMat[nrow(solMat),ncol(solMat)] <- tableau[nrow(tableau),ncol(tableau)] #sets the value
    colnames(solMat) <- colnames(tableau)[1:(ncol(tableau)-1)] #sets the column names
    rownames(solMat) <- "Solution" #sets the row name
    if(problem == FALSE) { #if the problem parameter is false
      finalList <- list(final.tableau = tableau, basic.solution = solMat, opt.val = tableau[nrow(tableau),ncol(tableau)]) #sets the final list to be returned
      return(finalList) #returns the list
    } else { #if the problem parameter is true
      shipVec <- c() #creates a vector for the values of x's
      if(ncol(tableau) == 25) { #to make sure that it is the exer9 problem
        for(s in 9:(ncol(tableau)-2)) {
          shipVec <- append(shipVec, tableau[nrow(tableau),s]) #append the vector with the values of x
        }
        shipMat <- matrix(shipVec, nrow = 3, ncol = 5,byrow = TRUE) #creates a matrix for the ship num
        colnames(shipMat) <- c("SAC","SL","ALB", "CHI", "NYC") #sets the col names
        rownames(shipMat) <- c("Denver", "Phoenix", "Dallas") #sets the row names
        finalList <- list(final.tableau = tableau, basic.solution = solMat, opt.val = tableau[nrow(tableau),ncol(tableau)], shipping.num = shipMat) #final list to be returnd
        return(finalList) #returns the final list
      } else {
        finalList <- list(final.tableau = tableau, basic.solution = solMat, opt.val = tableau[nrow(tableau),ncol(tableau)]) #list to be returned
        return(finalList) #returns the list
      }
      
    }
  }
}
#This is the problem to be solved in Exer 9
problem1 <- setupTableauTestCase(matrix(c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,-310,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,-260,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-280,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,180,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,80,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,200,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,160,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,220,10,8,6,5,4,6,5,4,3,6,3,4,5,5,9,1),nrow = 9,ncol = 16,byrow = TRUE))
problemSolve <- simplex(problem1,FALSE,TRUE)
problemSolve

#testcases for the exer 9
testcase2 <- setupTableauTestCase(matrix(c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,-200,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,-200,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-200,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,100,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,100,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,100,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,100,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,100,5,6,7,8,9,6,7,8,9,10,3,5,7,11,13,1),nrow = 9,ncol = 16,byrow = TRUE))
testSolve2 <- simplex(testcase2,FALSE,TRUE)
testSolve2

testcase3 <- setupTableauTestCase(matrix(c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,-1400,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,-400,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-200,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,431,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,332,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,350,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,450,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,400,30,29,31,35,33,26,24,23,25,27,11,13,15,20,17,1),nrow = 9,ncol = 16,byrow = TRUE))
testSolve3 <- simplex(testcase3,FALSE,TRUE)
testSolve3

testcase4 <- setupTableauTestCase(matrix(c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,-100,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,-100,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-100,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,20,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,20,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,20,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,20,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,20,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,1),nrow = 9,ncol = 16,byrow = TRUE))
testSolve4 <- simplex(testcase4,FALSE,TRUE)
testSolve4

testcase5 <- setupTableauTestCase(matrix(c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,-50,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,-50,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-50,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,20,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,25,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,90,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,60,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,70,30,29,31,35,33,26,24,23,25,27,11,13,15,20,17,1),nrow = 9,ncol = 16,byrow = TRUE))
testSolve5 <- simplex(testcase5,FALSE,TRUE)
testSolve5

#this is a sample problem for minimization from the simplex method handout
testt <- matrix(c(1,7,1,0,0,14,2,6,0,1,0,20,-4,-20,0,0,1,0), nrow = 3, ncol = 6, byrow = TRUE, dimname = list(c("constraint 1","constraint 2", "Objective function"), c("S1","S2","X1","X2","Z","Solution")))
answerr <- simplex(testt, FALSE, FALSE)
answerr

#ths is a sample problem for maximization from the simplex method handout
maxHandout <- matrix(c(7,11,1,0,0,0,0,77,10,8,0,1,0,0,0,80,1,0,0,0,1,0,0,9,0,1,0,0,0,1,0,6,-150,-175,0,0,0,0,1,0), nrow = 5, ncol = 8, byrow = TRUE, dimname = list(c("constraint 1","constraint 2","constraint 3","constraint 4", "Objective function"),c("X1","X2","S1","S2","S3","S4","Z","Solution")))
ans2 <- simplex(maxHandout,TRUE,FALSE)
ans2
