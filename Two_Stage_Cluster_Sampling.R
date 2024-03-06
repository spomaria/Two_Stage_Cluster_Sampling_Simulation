# This script is meant to run simulation for 
# Two stage Cluster Sampling in the presence of 
# Sampling Non-response

# @param N is the number of clusters in the population
# @param n is the number of clusters to be sampled
# @param M is the number of units in each cluster
# @param r is the number of non-response (could vary across clusters or remain constant)
# @param m is the number of units to be sampled in each selected cluster
# @param mu is a vector of means for the variables y, x, and z
# @param Sigma is a symmetric matrix. It is the variance-covariance matrix of the variables
# @param Case takes only two possible values: "A" or "B"

TwoStageClusterSampling <- function(
  M, N, m, n, r, mPrime, nPrime, mu, Sigma, Case  
  ){
  # Check that the number of clusters to be sampled does not 
  # exceed the total number of clusters in the population
 if ( n > N | n == N) return("Error: n should be less than N")
  # Check that the number of units to be sampled per cluster 
  # does not exceed the cluster size
  else if (m > M | m == M) return("Error: m should be less than M")
  else if (m > mPrime | m == mPrime) return("Error: m should be less than mPrime")
  else if (n > nPrime | n == nPrime) return("Error: n should be less than nPrime")
  # Check that the number of non-response per cluster 
  # does not exceed the sample size per cluster
  else if (r < 0 | r > m-2) return("Error: r should range between 0 and m-2")
  else if (!(length(mu) == ncol(Sigma)) | !(length(mu) == nrow(Sigma))) return("length of mu and dim of Sigma inconsistent")
  else {
    # Defining some of the derived constants needed in the various computations
    p = r/m
    q = 1 - p
    f = 1/n - 1/N
    fPrime = 1/nPrime - 1/N
    fmPrime = 1/mPrime - 1/M
    f1 = 1/n - 1/nPrime
    fm = 1/m - 1/M
    fmr = 1/(m * q + 2 * p) - 1/M
    fmPrimeR = 1/(mPrime * q + 2 * p) - 1/M
    
    
    totalSim = M * N
    # Simulating the entire population
    # Load the 'MASS' library to enable us generate 
    # random numbers from a multivariate normal distribution
    library(MASS)
    population = mvrnorm(totalSim, mu = mu, Sigma = Sigma)
    # detach the library after the simulation
    detach("package:MASS", unload = TRUE)
    
    # Demarcation of the various clusters in the population
    # We split the population into clusters by converting
    # the population into a matrix where each column 
    # represents a cluster
    # This is carried out for all three variables
    clustersY = matrix(population[,1], ncol = N)
    clustersX = matrix(population[,2], ncol = N)
    clustersZ = matrix(population[,3], ncol = N)
    
    # @clusterSpace is a sequence of numbers representing the index
    # of each of the clusters.
    clusterSpace = c(1:N)
    # @fsu1Clusters is a sequence of numbers corresponding to 
    # clusters to be sampled across all variables
    fsu1Clusters = sample(clusterSpace, nPrime)
    if (Case = "A") fsu2Clusters = sample(fsu1Clusters, n)
    else if (Case = "B") fsu2Clusters = sample(clusterSpace, n)
    # @sampledClusters is a matrix whose columns represent the sampled clusters
    # This refers to the first stage units (fsu) or primary stage units (psu)
    sampledClustersY = matrix(NA, ncol = n, nrow = M)
    sampledClustersX = matrix(NA, ncol = n, nrow = M)
    sampledClustersZ = matrix(NA, ncol = n, nrow = M)
    for (i in 1:n){
      # The main variable Y
      sampledClustersY[,i] = clustersY[,fsu2Clusters[i]]
      # The auxiliary variable X
      sampledClustersX[,i] = clustersX[,fsu2Clusters[i]]
      # The auxiliary variable Z
      sampledClustersZ[,i] = clustersZ[,fsu2Clusters[i]]
    }
    
    # @sampleSpace is a sequence of numbers representing the index
    # of units that make up each cluster.
    sampleSpace = c(1:M)
    # @indexSampleSpace is a sequence of numbers corresponding to 
    # units to be sampled in each cluster across all variables
    indexOfSampleSpace = sample(sampleSpace, m)
    # non-response space
    # To make this exercise more practical, we shall vary the index
    # of non-responsive units across the different clusters
    
    # response space
    # to obtain the response space, we write a function that computes the 
    # relative complement of any given set (or vector)
    # Note that in the function below,
    # @param a should be the main set while 
    # @param b should be the subset i.e. bearing in mind that the  
    # non-response set is a subset of the of sample space
    relcomp <- function(a, b) {
      # initiate an empty vector where elements of the complement set
      # will be included
      comp <- vector()
      
      for (i in a) {
        # Recall that we intend to find the complement of set 'a'
        # relative to set 'b'
        # The below condition checks for elements in set 'a' that are
        # not contained in set 'b'
        if (i %in% a && !(i %in% b)) {
          # The append function below includes the element 'i'
          # to the set (vector) 'comp' from the rear 
          # i.e. after the last element of 'comp'
          comp <- append(comp, i)
        }
      }
      
      return(comp)
    }
    
    
    # Y, X and Z matrices whose columns represent the sampled units
    # per cluster i.e. the second stage units (ssu)
    # The main variable Y
    Y = matrix(NA, ncol = n, nrow = m)
    # The auxiliary variable X
    X = matrix(NA, ncol = n, nrow = m)
    # The auxiliary variable Z
    Z = matrix(NA, ncol = n, nrow = m)
    
    # sampling m units from a total of M units in each cluster
    
    # This looping approach ensures that for each unit selected,
    # the three variables are observed accordingly so that
    # there is no mixing of variables across selected units
    
    if (Case == A){
      # In case A, the second sample is drawn as a subset of the first sample 
      for (i in 1:n){
        # simulate a different non-response space for each cluster
        nonresponseSpace = sample(indexOfSampleSpace, r)
        
        # vector (or set) of units that responded during sampling
        response <- relcomp(indexOfSampleSpace, nonresponseSpace)
        # since only 'n-r' responded to the exercise, we re-sample r units 
        # from the 'n-r' that responded so that we end up with 'n' responses
        # as planned i.e. (n-r) + r = n
        indexOfFinalSample <- c(response, sample(indexOfSampleSpace, r))
        
        for (j in 1:m){
          
          # The main variable Y
          Y[j,i] = sampledClustersY[indexOfFinalSample[j],i]
          # The auxiliary variable X
          X[j,i] = sampledClustersX[indexOfFinalSample[j],i]
          # The auxiliary variable Z
          Z[j,i] = sampledClustersZ[indexOfFinalSample[j],i]
        }
      }
    } else{
      # In case B, the second sample is drawn independent of the first sample
      for (i in 1:n){
        # simulate a different non-response space for each cluster
        nonresponseSpace = sample(indexOfSampleSpace, r)
        
        # vector (or set) of units that responded during sampling
        response <- relcomp(indexOfSampleSpace, nonresponseSpace)
        # since only 'n-r' responded to the exercise, we re-sample r units 
        # from the 'n-r' that responded so that we end up with 'n' responses
        # as planned i.e. (n-r) + r = n
        
        # where the second sample is independent of the first sample
        indexOfFinalSample <- c(response, sample(sampleSpace, r))
        
        for (j in 1:m){
          
          # The main variable Y
          Y[j,i] = sampledClustersY[indexOfFinalSample[j],i]
          # The auxiliary variable X
          X[j,i] = sampledClustersX[indexOfFinalSample[j],i]
          # The auxiliary variable Z
          Z[j,i] = sampledClustersZ[indexOfFinalSample[j],i]
        }
      }
    }
    
    # Computing the several components of the estimators
    # Components relating to variable Y
    S0i_2 <- c()
    for (i in 1:n){
      S0i_2[i] = (Y[,i] - mean(Y[,i]))^2/(m-1)
    }
    S0bar_2 <- mean(S0i_2)
    
    Sibar <- c()
    for (i in 1:n){
      Sibar[i] = mean(Y[,i])
    }
    S0_2 <- sum(Sibar - mean(Y))^2/(n-1)
    
    # Components relating to variable X
    S1i_2 <- c()
    for (i in 1:n){
      S1i_2[i] = (X[,i] - mean(X[,i]))^2/(m-1)
    }
    S1bar_2 <- mean(S1i_2)
    
    Sxibar <- c()
    for (i in 1:n){
      Sxibar[i] = mean(X[,i])
    }
    S1_2 <- sum(Sxibar - mean(X))^2/(n-1)
    
    # Components relating to variable Z
    S2i_2 <- c()
    for (i in 1:n){
      S2i_2[i] = (Z[,i] - mean(Z[,i]))^2/(m-1)
    }
    S2bar_2 <- mean(S2i_2)
    
    Szibar <- c()
    for (i in 1:n){
      Szibar[i] = mean(Z[,i])
    }
    S2_2 <- sum(Szibar - mean(Z))^2/(n-1)
    
    # Components relating to variables Y and X
    S01i <- c()
    for (i in 1:n){
      S01i[i] = sum(prod((Y[,i] - mean(Y[,i])), (X[,i] - mean(X[,i]))))/(m-1)
    }
    S01bar <- mean(S01i)
    
    S01istar <- c()
    for (i in 1:n){
      S01istar[i] = prod((mean(Y[,i]) - mean(Y)), (mean(X[,i]) - mean(X)))
    }
    S01star <- sum(S01istar)/(n-1)
    
    # Components relating to variables Y and Z
    S02i <- c()
    for (i in 1:n){
      S02i[i] = sum(prod((Y[,i] - mean(Y[,i])), (Z[,i] - mean(Z[,i]))))/(m-1)
    }
    S02bar <- mean(S02i)
    
    S02istar <- c()
    for (i in 1:n){
      S02istar[i] = prod((mean(Y[,i]) - mean(Y)), (mean(Z[,i]) - mean(Z)))
    }
    S02star <- sum(S02istar)/(n-1)
    
    # Components relating to variables X and Z
    S12i <- c()
    for (i in 1:n){
      S12i[i] = sum(prod((X[,i] - mean(X[,i])), (Z[,i] - mean(Z[,i]))))/(m-1)
    }
    S12bar <- mean(S12i)
    
    S12istar <- c()
    for (i in 1:n){
      S12istar[i] = prod((mean(X[,i]) - mean(X)), (mean(Z[,i]) - mean(Z)))
    }
    S12star <- sum(S12istar)/(n-1)
    
    # Calculating the basic variance which is unbiased
    # Calculating mean square between cluster means
    sbsquareUnitsY = c()
    for (i in 1:n){
      sbsquareUnitsY[i] = (mean(Y[,i]) - mean(Y))^2
    }
    sbsquare = 1/(n - 1) *sum(sbsquareUnitsY)
    
    # Calculating mean square within the clusters
    swsquareUnitsY = c()
    for (i in 1:n){
      swsquareUnitsY[i] = sum(Y[,i] - mean(Y[,i]))^2
    }
    swsquare = 1/n * 1/(m - 1) *sum(swsquareUnitsY)
    
    vYbarnm = f * sbsquare + 1/n *fm*swsquare
    
    return (list(clustersY = clustersY, indexOfSampleSpace = indexOfSampleSpace, 
                indexOfFinalSample = indexOfFinalSample, sampledClustersY = sampledClustersY, Y = Y,
                 vYbarnm = vYbarnm))
  }
}

Sigma <- matrix(c(50, 2, 3, 2, 100, 7, 3, 7 , 50), 3,3)

a <- TwoStageClusterSampling(M = 15, N = 20, m = 7, n = 7, r = 2, mPrime = 5,
                             nPrime = 6, p = 0.5, mu = c(10, 50, 20),
                             Sigma = Sigma)
a
a$fsu2Clusters
a$clusters
a$sampledClusters
