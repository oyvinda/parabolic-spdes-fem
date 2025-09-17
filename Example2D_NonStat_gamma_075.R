################################################################################
# Statistics inspired example                                                  #
################################################################################

# Libraries
library(doParallel)
registerDoParallel(cores=3)
library(foreach)

## Get functions
source('functions.R')
source('mesh2d_experiment.R')


## Lindgren et al.
  # Smoothness 0.5 in space
  gamma = 0.75
  
  # Spatial range rs = 0.25
  #  rs = 0.25
  #  gammaS = sqrt(8*2*gamma)/rs
  gammaS = 14
  kappa1Sq = gammaS^2
  
  # Temporal range rT = 0.3125
  #  rT = 0.3
  #  scaleTime = rT*gammaS^2/sqrt(4)
  scaleTime = 30
  
  # Variance \sigma^2 = 1
  #alphT = 1
  #alph = 2*gamma+1
  #Ct = gamma(alphT-0.5)/gamma(alphT)/sqrt(4*pi)
  #Cs = gamma(alph-1)/gamma(alph)/(4*pi)
  #noiseScale = sqrt(Ct*Cs/scaleTime/gammaS^(2*alph-2))
  noiseScale = 5.6e-4
  
  # Diffusion
  a0 = 1


# Configure problem
  # Dampening
  kappaSqFun = function(loc){
    return(rep(kappa1Sq, nrow(loc)))
  }
  
  # Diffusion
  aFun = function(loc){
    return(rep(a0,nrow(loc)))
  }
  vxFun = function(loc){
    return(-8*(loc[,2]-0.5))
  }
  vyFun = function(loc){
    return(8*(loc[,1]-0.5))
  }
  
  # Advection
  wxFun = function(loc){
    return(rep(200, nrow(loc)))
  }
  wyFun = function(loc){
    return(rep(0, nrow(loc)))
  }

## Discretization
  # Temporal settings
  Tend = 1      # Domain [0, Tend]
  kMax = 12
  kMin = 12
  K = 2^kMax      # Regular time discretization 0 = t_0 < ... < t_K = Tend
  
  # Spatial settings
  meshList = list()
  mMin = 6
  mMax = 7
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
  # Mass matrix, A1mat, and A2mat
  Cmat = list()
  A1mat = list()
  A2mat = list()
  zeroFun = function(loc){
    return(0*loc[,1])
  }
  for(i in 1:length(ySol[[1]])){
    Cmat[[i]] = makeMassMatrixMesh(meshList[[i]])
    A1mat[[i]] = makeOpMat(mesh = meshList[[i]],
                           kappa2 = kappaSqFun,
                           a = aFun,
                           vx = vxFun,
                           vy = vyFun, 
                           wx = wxFun, 
                           wy = wyFun)
    A2mat[[i]] = makeOpMat(mesh = meshList[[i]],
                           kappa2 = kappaSqFun,
                           a = aFun,
                           vx = vxFun,
                           vy = vyFun, 
                           wx = zeroFun, 
                           wy = zeroFun)
    # Scale time derivative
    A1mat[[i]] = A1mat[[i]]/scaleTime
    A2mat[[i]] = A2mat[[i]]*scaleTime^(1/gamma)
    
    # Scale noise
    A2mat[[i]] = A2mat[[i]]*noiseScale^{1/gamma}
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
      # NB: NO LONGER SYMMETRIC -- NO CHOLESKY FACTORIZATION
      Llhs[[j]][[i]] = Alhs
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
  solveSPDE = function(rSeed, tSave){
    # Set seed
    set.seed(rSeed)
    
    ySolTmp = list()
    for(j in 1:(kMax-kMin+1)){
      ySolTmp[[j]] = list()
      for(i in 1:length(ySol[[1]])){
        ySolTmp[[j]][[i]] = rep(0, length(ySol[[1]][[i]]))
      }
    }
    
    # Store at desired resolution
    ySolRes = list()
    for(j in 1:(kMax-kMin+1)){
      ySolRes[[j]] = list()
      for(i in 1:length(ySol[[1]])){
        ySolRes[[j]][[i]] = matrix(0, nrow = length(ySol[[1]][[i]]), ncol = 2^tSave+1)
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
            uNext = solve(Llhs[[idxTimeRes]][[kRes]], uTmp1)
            
            # Store solution
            ySolTmp[[idxTimeRes]][[kRes]] = as.vector(uNext)
          }
        }
      }
      
      if(bIdx%%2^(kMin-tSave) == 0){
        for(j in 1:(kMax-kMin+1)){
          for(i in 1:length(ySol[[1]])){
            ySolRes[[j]][[i]][,bIdx/2^(kMin-tSave)+1] = ySolTmp[[j]][[i]]
          }
        }
      }
    }
    
    return(ySolRes)
  }
  
## Solve SPDE multiple times
  # Number of simulations
  nSim = 3
  set.seed(14133)
  rSeeds = sample(1:1000000, 
                  size = nSim,
                  replace = T)
  
  tSave= 10
  ySolAll = foreach(i = 1:nSim) %dopar% solveSPDE(rSeeds[i], tSave)
  save(file="Des12_2D_statistics_g075.RData", ySolAll, rSeeds, gamma, ySol, Cmat, dt, meshList, Pcoarse)
  
  