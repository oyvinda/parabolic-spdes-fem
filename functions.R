##############################################################################
# functions.R                                                                #
#    Functions for nested FEM in 1D                                          #
##############################################################################

# Load libraries
library(Matrix)

# Function to create mass matrix of FEM basis
#    Assumes D = (0, xMax) and regular grid with nGrid basis functions
#    Computes C_{ij} = <\phi_i, \phi_j>
makeMassMatrix = function(xMax, nGrid){
  # Grid resolution
  h = xMax/(nGrid+1)
  
  # Diagonal entries
  iVec = 1:nGrid
  jVec = 1:nGrid
  vVec = h*rep(2/3, nGrid)
  
  # Append entries below diagonal
  iVec = c(iVec, 2:nGrid)
  jVec = c(jVec, 1:(nGrid-1))
  vVec = c(vVec, rep(h/6, nGrid-1))
  
  # Append entries above diagonal
  iVec = c(iVec, 1:(nGrid-1))
  jVec = c(jVec, 2:nGrid)
  vVec = c(vVec, rep(h/6, nGrid-1))
  
  # Assemble matrix
  Cmat = sparseMatrix(i = iVec, 
                      j = jVec, 
                      x = vVec)  
  
  return(Cmat)
}

# Function to create stiffness matrix of FEM basis
#    Assumes D = (0, xMax) and regular grid with nGrid basis functions
#    Computes G_{ij} = <\nabla \phi_i, \nabla \phi_j>
makeStiffMatrix = function(xMax, nGrid){
  # Resolution
  h = xMax/(nGrid+1)
  
  # Diagonal entries
  iVec = 1:nGrid
  jVec = 1:nGrid
  vVec = rep(2, nGrid)/h
  
  # Append subdiagonal entries
  iVec = c(iVec, 2:nGrid)
  jVec = c(jVec, 1:(nGrid-1))
  vVec = c(vVec, rep(-1/h, nGrid-1))
  
  # Append superdiagonal entries
  iVec = c(iVec, 1:(nGrid-1))
  jVec = c(jVec, 2:nGrid)
  vVec = c(vVec, rep(-1/h, nGrid-1))
  
  # Assemble matrix
  Gmat = sparseMatrix(i = iVec,
                      j = jVec, 
                      x = vVec)
  
  return(Gmat)
}

# make (\kappa^2-D_xx) with Dirichlet = 0 on boundary so need to remove it
makeSpatialOperator = function(C, G, kappa){
  A = kappa^2*C+G
  return(A)
}

# This functions computes L2 norms exactly between fine res and coarser resolutions
computeL2norm = function(ySol, xMax, Cmass){
  yFine = list()
  for(kRes in 1:length(ySol)){
    yFine[[kRes]] = ySol[[kRes]]
  }
  for(kRes in 1:(length(ySol)-1)){
    for(i in 1:(length(ySol)-kRes)){
      yTmp = rep(0, length(ySol[[kRes+1]]))
      yOld = yFine[[i]]
      for(j in 1:length(yOld)){
        yTmp[2*j] = yOld[j]
      }
      #  yTmp[2] = (yOld[1]+yOld[2])/2
      #  yTmp[length(yTmp)-1] = (yOld[length(yOld)]+yOld[length(yOld)-1])/2
      yOld = c(0, yOld, 0)
      for(j in 1:(length(yOld)-1)){
        yTmp[2*j-1] = (yOld[j]+yOld[j+1])/2
      }
      yFine[[i]] = yTmp
    }
  }
  
  nVecAbs = rep(0, length(ySol))
  nVecRel = nVecAbs
  idxBest = length(ySol)
  nBest = sqrt(as.vector(yFine[[idxBest]])%*%Cmass%*%as.vector(yFine[[idxBest]]))
  dx = rep(0, length(ySol))
  for(kRes in 1:length(ySol)){
    nVecAbs[kRes] = sqrt(as.vector(yFine[[kRes]]-yFine[[idxBest]])%*%Cmass%*%as.vector(yFine[[kRes]]-yFine[[idxBest]]))
    nVecRel[kRes] = nVecAbs[kRes]/nBest
    dx[kRes] = xMax/(length(ySol[[kRes]])+1)
  }
  return(list(L2abs = nVecAbs,
              L2rel = nVecRel,
              L2best = nBest,
              dx = dx))
}

# This functions computes L2 norms exactly between fine res and coarser resolutions
# This NEW function considers different resolutions in space and in time
computeL2normNEW = function(ySol, xMax, Cmass){
  
  # Resolution
  tRes = length(ySol)
  xRes = length(ySol[[1]])
  
  # Store all solutions at finest spatial resolution
  yFine = list()
  for(tIdx in 1:tRes){
    yFine[[tIdx]] = list()
    for(kRes in 1:xRes){
      yFine[[tIdx]][[kRes]] = ySol[[tIdx]][[kRes]]
    }
    for(kRes in 1:(xRes-1)){
      for(i in 1:(xRes-kRes)){
        yTmp = rep(0, length(ySol[[tIdx]][[kRes+1]]))
        yOld = yFine[[tIdx]][[i]]
        for(j in 1:length(yOld)){
          yTmp[2*j] = yOld[j]
        }
        #  yTmp[2] = (yOld[1]+yOld[2])/2
        #  yTmp[length(yTmp)-1] = (yOld[length(yOld)]+yOld[length(yOld)-1])/2
        yOld = c(0, yOld, 0)
        for(j in 1:(length(yOld)-1)){
          yTmp[2*j-1] = (yOld[j]+yOld[j+1])/2
        }
        yFine[[tIdx]][[i]] = yTmp
      }
    }
  }
  
  nVecAbs = matrix(0, nrow = tRes, ncol = xRes)
  nVecRel = nVecAbs
  nBest = as.numeric(sqrt(as.vector(yFine[[tRes]][[xRes]])%*%Cmass%*%as.vector(yFine[[tRes]][[xRes]])))
  dx = rep(0, xRes)
  for(tIdx in 1:tRes){
    for(kRes in 1:xRes){
      nVecAbs[tIdx, kRes] = as.numeric(sqrt(as.vector(yFine[[tIdx]][[kRes]]-yFine[[tRes]][[xRes]])%*%Cmass%*%as.vector(yFine[[tIdx]][[kRes]]-yFine[[tRes]][[xRes]])))
      nVecRel[tIdx, kRes] = nVecAbs[tIdx, kRes]/nBest
      dx[kRes] = xMax/(length(ySol[[tIdx]][[kRes]])+1)
    }
  }
  return(list(L2abs = nVecAbs,
              L2rel = nVecRel,
              L2best = nBest,
              dx = dx))
}
