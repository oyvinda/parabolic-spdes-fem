##############################################################################
# Example.R                                                                  #
#    This script demonstrates how to solve                                   #
#       D_t y = -(\kappa_1^2-D_{xx})y + (\kappa_2^2-D_{xx})^{-\gamma} D_t W  #
#    where                                                                   #
#       W: Generalized Wiener process with covariance operator I             #
#       \kappa_1, \kappa_2 > 0: Spatial scaling parameters for the two       #
#                               operators                                    #
#       0 < \gamma <= 1: Spatial smoothness parameter                        #
#                                                                            #
#    The script solves for y(.,T), where T is terminal time, for multiple    #
#    resolutions simultaneously.                                             #
##############################################################################

# Load libraries
library(doParallel)
registerDoParallel(cores=5)
library(foreach)

## Get functions
source('functions.R')

## Configure problem
  # Parameters
  kappa1 = 0
  kappa2 = 1
  gamma = 0.5
  
  # Temporal settings
  Tend = 1      # Domain [0, Tend]
  kMax = 22
  kMin = 13
  K = 2^kMax      # Regular time discretization 0 = t_0 < ... < t_K = Tend
  
  # Spatial settings
  xMax = 1     # Domain (0, xMax)
  mMin = 3      # Minimum 2^mMin+1 basis elements
  mMax = 11     # Maximum 2^mMax+1 basis elements
  
  # Quadrature resolution
  kQ = 0.5
  
## Initialize storage for solutions and spatial grids
  ySol = list()
  dt = list()
  for(i in 1:(kMax-kMin+1)){
    ySol[[i]] = list()
    dt[[i]] = Tend/2^(kMin+i-1)
  }
  xGrid = list()
  for(kRes in mMin:mMax){
    idx = kRes-mMin+1
    for(i in 1:(kMax-kMin+1)){
      ySol[[i]][[idx]] = matrix(NA, nrow = 2^kRes-1, ncol = 1)
      
      # Initial condition
      ySol[[i]][[idx]][,1] = 0
    }
    
    # Grid
    nGrid = nrow(ySol[[1]][[idx]])
    xGrid[[idx]] = seq(xMax/(nGrid+1), xMax-xMax/(nGrid+1), length.out = nGrid)
  }
  
## Preparation of quadrature
  # Compute quadrature variables
  NQ = ceiling(pi^2/(2*gamma*kQ^2))
  MQ = ceiling(pi^2/(2*(1-gamma)*kQ^2))
  CQ = kQ*sin(pi*gamma)/pi
  if(gamma > 1-1e-8){
    gamma = 1
    NQ = 0
    MQ = 0
  }
  if(gamma < 1e-8){
    gamma = 0
    NQ = 0
    MQ = 0
  }

## Prepare spatial matrices
  # Mass matrices and stiffness matrices
  Cmat = list()
  Gmat = list()
  for(i in 1:length(ySol[[1]])){
    nGrid = length(ySol[[1]][[i]])
    Cmat[[i]] = makeMassMatrix(xMax, nGrid)
    Gmat[[i]] = makeStiffMatrix(xMax, nGrid)
  }
  
  # Evolution and noise operators
  A1mat = list()
  A2mat = list()
  for(i in 1:length(ySol[[1]])){
    A1mat[[i]] = makeSpatialOperator(Cmat[[i]], Gmat[[i]], kappa1)
    A2mat[[i]] = makeSpatialOperator(Cmat[[i]], Gmat[[i]], kappa2)
  }
  
  # Operators needed in quadrature
  Q = list() # list of each resolution
  for(i in 1:length(ySol[[1]])){
    Q[[i]] = list() # List of matrices needed in quadrature for a given resolution
    for(j in (-MQ):NQ){
      idx = j + MQ + 1
      if(gamma > 1-1e-8){
        Q[[i]][[idx]] = A2mat[[i]]
        break
      }
      if(gamma < 1e-8){
        Q[[i]][[idx]] = Cmat[[i]]
        break
      }
      yj = j*kQ
      Q[[i]][[idx]] = (exp(yj)*Cmat[[i]] + A2mat[[i]])/(CQ*exp((1-gamma)*yj))
    }
  }

## Precompute Cholesky factorizations
  Lc = list()
  Llhs = list()
  for(i in 1:(kMax-kMin+1)){
    Llhs[[i]] = list()
  }
  Lq = list()
  for(i in 1:length(ySol[[1]])){
    # Precompute Cholesky factor of mass matrix
    Lc[[i]] = chol(Cmat[[i]], pivot = F)
    Lc[[i]] = t(Lc[[i]])
    
    # Precompute Cholesky factorization of LHS
    for(j in 1:(kMax-kMin+1)){
      Alhs = Cmat[[i]] + dt[[j]]*A1mat[[i]]
      Llhs[[j]][[i]] = Cholesky(Alhs)
    }
    
    # Precompute Cholesky factorization of matrices for quadrature
    Lq[[i]] = list()
    for(j in 1:length(Q[[1]])){
      Lq[[i]][[j]] = Cholesky(Q[[i]][[j]])
    }
  }
  
## Prepare to map noise to lower resolutions
  Pcoarse = list()
  for(kRes in 1:(length(ySol[[1]])-1)){
    n2 = length(ySol[[1]][[kRes]])
    n1 = length(ySol[[1]][[kRes+1]])
    
    # Distance in finest grid
    iVec = rep(1:n2, each = 3)
    jVec = as.vector(rbind(seq(1, n1-2, length.out = n2),
                           seq(2, n1-1, length.out = n2),
                           seq(3, n1, length.out = n2)))
    vVec = rep(c(1/2,1,1/2), n2)
    Pcoarse[[kRes]] = sparseMatrix(i = iVec,
                                   j = jVec,
                                   x = vVec)
  }
  
## Function for solving SPDE
  ## Move forward in time
  solveSPDE = function(rSeed){
    # Set seed
    set.seed(rSeed)
    
    ySolTmp = list()
    for(j in 1:(kMax-kMin+1)){
      ySolTmp[[j]] = list()
      for(i in 1:length(ySol[[1]])){
        ySolTmp[[j]][[i]] = rep(0, length(ySol[[1]][[i]]))
      }
    }
    # Compute solution in batches so all solutions are moved forward at the
    # coarsest solution
    for(bIdx in 1:(2^kMin)){
      print(bIdx)
      
      # Simulate white noise in V (highest spatial resolution)
      idxMax = length(ySol[[1]])
      nGrid = nrow(ySol[[1]][[idxMax]])
      simNoiseFine = Lc[[idxMax]]%*%matrix(rnorm(nGrid*2^(kMax-kMin)), nrow = nGrid)
      simNoise = list()
      simNoise[[idxMax]] = as.matrix(simNoiseFine)
      
      # Project noise to coarser spatial resolutions in V
      for(kRes in (idxMax-1):1){
        simNoise[[kRes]] = as.matrix(Pcoarse[[kRes]]%*%simNoise[[kRes+1]])
      }
      
      rhs = list()
      # Smooth driving noise
      for(kRes in 1:length(ySol[[1]])){
        nGrid = length(ySol[[1]][[kRes]])
        
        # Use quadrature 
        rhs[[kRes]] = 0*simNoise[[kRes]]
        for(i in 1:length(Lq[[kRes]])){
          rhs[[kRes]] = rhs[[kRes]] + as.matrix(solve(Lq[[kRes]][[i]], simNoise[[kRes]], system = "A"))
        }
      }
      
      # Move solutions forward at different resolutions within coarse resolution
      for(idxTimeRes in 1:(kMax-kMin+1)){
        # Move solution forward in time
        for(kRes in 1:length(ySol[[1]])){
          for(idxTime in 1:2^(idxTimeRes-1)){
            uPrev = ySolTmp[[idxTimeRes]][[kRes]]
            
            # All operations are only acting in V_0
            idxSum = (idxTime-1)*2^(kMax-kMin+1-idxTimeRes) + (1:2^(kMax-kMin+1-idxTimeRes))
            uTmp1 = Cmat[[kRes]]%*%(uPrev+sqrt(dt[[kMax-kMin+1]])*rowSums(rhs[[kRes]][,idxSum, drop=FALSE]))
            uNext = solve(Llhs[[idxTimeRes]][[kRes]], uTmp1, system = "A")
            
            # Store solution
            ySolTmp[[idxTimeRes]][[kRes]] = as.vector(uNext)
          }
        }
      }
    }
    
    return(ySolTmp)
  }
  
## Solve SPDE multiple times
  # Number of simulations
  nSim = 50
  set.seed(14133)
  rSeeds = sample(1:1000000, 
                  size = nSim,
                  replace = T)
  
  ySolAll = foreach(i = 1:nSim) %dopar% solveSPDE(rSeeds[i])
  save(file="Des2_ST_gamma_05_50.RData", ySolAll, rSeeds, gamma, ySol, xMax, Cmat, dt)
