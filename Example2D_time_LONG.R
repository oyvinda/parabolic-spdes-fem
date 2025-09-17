##############################################################################
# Example.R                                                                  #
#    This script demonstrates how to solve                                   #
#       D_t y = -(\kappa_1^2-\Delta)y + (\kappa_2^2-\Delta)^{-\gamma} D_t W  #
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
source('mesh2d_experiment.R')


## Configure problem
  # Parameters
  kappa1 = 0
  kappa2 = 1
  gamma = 0.00
  
  # Temporal settings
  Tend = 1      # Domain [0, Tend]
  kMax = 20
  kMin = 12
  K = 2^kMax      # Regular time discretization 0 = t_0 < ... < t_K = Tend
  
  # Spatial settings
  meshList = list()
  mMin = 4
  mMax = 9
  for(i in mMin:mMax){
    meshList[[i-mMin+1]] = makeMesh(i)
  }
  
  # Quadrature resolution
  kQ = 0.5
  
## Initialize storage for solutions and spatial grids
  ySol = list()
  dt = list()
  for(i in 1:(kMax-kMin+1)){
    ySol[[i]] = list()
    dt[[i]] = Tend/2^(kMin+i-1)
  }
  for(kRes in mMin:mMax){
    idx = kRes-mMin+1
    for(i in 1:(kMax-kMin+1)){
      ySol[[i]][[idx]] = matrix(NA, nrow = length(meshList[[idx]]$int), ncol = 1)
      
      # Initial condition
      ySol[[i]][[idx]][,1] = 0
    }
    
    # Grid
    nGrid = nrow(ySol[[1]][[idx]])
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
    Cmat[[i]] = makeMassMatrixMesh(meshList[[i]])
    Gmat[[i]] = makeStiffMatrixMesh(meshList[[i]])
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
    tmpA = fm_evaluator_mesh_2d(meshList[[kRes]], meshList[[kRes+1]]$loc[meshList[[kRes+1]]$int,])$A
    tmpA = t(tmpA[,meshList[[kRes]]$int])
    Pcoarse[[kRes]] = tmpA
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
  nSim = 5
  set.seed(14133)
  rSeeds = sample(1:1000000, 
                  size = nSim,
                  replace = T)
  
  ySolAll = foreach(i = 1:nSim) %dopar% solveSPDE(rSeeds[i])
  save(file="Des10_ST_2D_gamma_000_TMP_LONG.RData", ySolAll, rSeeds, gamma, ySol, Cmat, dt, meshList, Pcoarse)
  
  
  print("solved and saved for gamma = 0.0. Computing other gammas")
  
# Solve for different gamma
ySolAllGamma = list()
allGamma = c(0.25, 0.5, 0.75, 1.0)
for(idxG in 1:length(allGamma)){
  
  gamma = allGamma[idxG]
  print("Gamma =")
  print(gamma)
  
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
  
  Lq = list()
  for(i in 1:length(ySol[[1]])){
    # Precompute Cholesky factorization of matrices for quadrature
    Lq[[i]] = list()
    for(j in 1:length(Q[[1]])){
      print(paste("kRes =", i, "  x  Q =", j, "of", MQ+NQ+1))
      Lq[[i]][[j]] = Cholesky(Q[[i]][[j]])
    }
  }
  
  ySolAllNEW = ySolAll
  for(i in 1:length(ySolAll)){
    for(idxTimeRes in 1:(kMax-kMin+1)){
      for(kRes in 1:length(ySol[[1]])){
        print(paste("i =", i, "  x  idxTimeRes =", idxTimeRes, "  x  kRes = ", kRes))
        uNextAlt = Cmat[[kRes]]%*%as.vector(ySolAll[[i]][[idxTimeRes]][[kRes]])
        uNextAltS = 0*uNextAlt
        for(j in 1:length(Lq[[kRes]])){
          uNextAltS = uNextAltS + as.vector(solve(Lq[[kRes]][[j]], uNextAlt, system = "A"))
        }
        ySolAllNEW[[i]][[idxTimeRes]][[kRes]] = as.vector(uNextAltS)
      }
    }
  }
  ySolAllGamma[[idxG]] = ySolAllNEW
  
#  save(file="Des10_ST_2D_5.RData", ySolAll, ySolAllGamma, allGamma, rSeeds, ySol, Cmat, dt, meshList, Pcoarse)
}

save(file="Des10_ST_2D_5_LONG.RData", ySolAll, ySolAllGamma, allGamma, rSeeds, ySol, Cmat, dt, meshList, Pcoarse)

  
  # xRes = seq(0, 1, length.out = 300)
  # yRes = seq(0, 1, length.out = 300)
  # xs = rep(xRes, each = 300)
  # ys = rep(yRes, 300)
  # 
  # kRes = 8
  # tRes = 7
  # Atmp = fm_evaluator_mesh_2d(mesh = meshList[[kRes]], as.matrix(cbind(xs, ys)))$A
  # Atmp = Atmp[,meshList[[kRes]]$int]
  # zs = Atmp%*%ySolAllGamma[[4]][[1]][[8]][[kRes]]
  # library(fields)
  # image.plot(x = xRes, y = yRes, z = matrix(zs, ncol = 300))
  # 
  #            