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
    M, N, m, n, r, mPrime, nPrime, mu, Sigma, Case = "A", Procedure, seed_num=4113  
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
    #set.seed(seed_num)
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
    
    if (Case == "A") {
      fsu1Clusters = sample(clusterSpace, nPrime)
      fsuClusters = sample(fsu1Clusters, n)
    }
    else if (Case == "B") fsuClusters = sample(clusterSpace, n)
    # @sampledClusters is a matrix whose columns represent the sampled clusters
    # This refers to the first stage units (fsu) or primary stage units (psu)
    sampledClustersY = matrix(NA, ncol = n, nrow = M)
    sampledClustersX = matrix(NA, ncol = n, nrow = M)
    sampledClustersZ = matrix(NA, ncol = n, nrow = M)
    for (i in 1:n){
      # The main variable Y
      sampledClustersY[,i] = clustersY[,fsuClusters[i]]
      # The auxiliary variable X
      sampledClustersX[,i] = clustersX[,fsuClusters[i]]
      # The auxiliary variable Z
      sampledClustersZ[,i] = clustersZ[,fsuClusters[i]]
    }
    
    # @sampleSpace is a sequence of numbers representing the index
    # of units that make up each cluster.
    sampleSpace = c(1:M)
    # @indexSampleSpace is a sequence of numbers corresponding to 
    # units to be sampled in each cluster across all variables
    if (Procedure ==1) indexOfSampleSpace = sample(sampleSpace, m)
    else if (Procedure ==2){
      ssu1Space = sample(sampleSpace, mPrime)
      indexOfSampleSpace = sample(ssu1Space, m)
    }
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
    
    for (i in 1:n){
      # simulate a different non-response space for each cluster
      nonresponseSpace = sample(indexOfSampleSpace, r)
      
      # vector (or set) of units that responded during sampling
      response <- relcomp(indexOfSampleSpace, nonresponseSpace)
      # since only 'n-r' responded to the exercise, we re-sample r units 
      # from the 'n-r' that responded so that we end up with 'n' responses
      # as planned i.e. (n-r) + r = n
      indexOfFinalSample <- c(response, sample(response, r))
      
      for (j in 1:m){
        
        # The main variable Y
        Y[j,i] = sampledClustersY[indexOfFinalSample[j],i]
        # The auxiliary variable X
        X[j,i] = sampledClustersX[indexOfFinalSample[j],i]
        # The auxiliary variable Z
        Z[j,i] = sampledClustersZ[indexOfFinalSample[j],i]
      }
    }
    
    # Computing the several components of the estimators
    # Components relating to variable Y
    S0i_2 <- c()
    for (i in 1:n){
      S0i_2[i] = sum((Y[,i] - mean(Y[,i]))^2)/(m-1)
    }
    S0bar_2 <- mean(S0i_2)
    
    Sibar <- c()
    for (i in 1:n){
      Sibar[i] = mean(Y[,i])
    }
    S0Star_2 <- sum((Sibar - mean(Y))^2)/(n-1)
    
    # Components relating to variable X
    S1i_2 <- c()
    for (i in 1:n){
      S1i_2[i] = sum((X[,i] - mean(X[,i]))^2)/(m-1)
    }
    S1bar_2 <- mean(S1i_2)
    
    Sxibar <- c()
    for (i in 1:n){
      Sxibar[i] = mean(X[,i])
    }
    S1Star_2 <- sum((Sxibar - mean(X))^2)/(n-1)
    
    # Components relating to variable Z
    S2i_2 <- c()
    for (i in 1:n){
      S2i_2[i] = sum((Z[,i] - mean(Z[,i]))^2)/(m-1)
    }
    S2bar_2 <- mean(S2i_2)
    
    Szibar <- c()
    for (i in 1:n){
      Szibar[i] = mean(Z[,i])
    }
    S2Star_2 <- sum((Szibar - mean(Z))^2)/(n-1)
    
    # Components relating to variables Y and X
    S01i <- c()
    for (i in 1:n){
      #elems <- c()
      #for (j in 1:m){
      #  elems[j] = prod((Y[j,i] - mean(Y[,i])), (X[j,i] - mean(X[,i])))
      #}
      S01i[i] = sum(drop(crossprod((Y[j,i] - mean(Y[,i])), (X[j,i] - mean(X[,i])))))/(m-1)
    }
    S01bar <- mean(S01i)
    
    S01istar <- c()
    for (i in 1:n){
      S01istar[i] = drop(crossprod((mean(Y[,i]) - mean(Y)), (mean(X[,i]) - mean(X))))
    }
    S01star <- sum(S01istar)/(n-1)
    
    # Components relating to variables Y and Z
    S02i <- c()
    for (i in 1:m){
      S02i[i] = sum(drop(crossprod((Y[,i] - mean(Y[,i])), (Z[,i] - mean(Z[,i])))))/(m-1)
    }
    S02bar <- mean(S02i)
    
    S02istar <- c()
    for (i in 1:m){
      S02istar[i] = drop(crossprod((mean(Y[,i]) - mean(Y)), (mean(Z[,i]) - mean(Z))))
    }
    S02star <- sum(S02istar)/(n-1)
    
    # Components relating to variables X and Z
    S12i <- c()
    for (i in 1:m){
      S12i[i] = sum(drop(crossprod((X[,i] - mean(X[,i])), (Z[,i] - mean(Z[,i])))))/(m-1)
    }
    S12bar <- mean(S12i)
    
    S12istar <- c()
    for (i in 1:m){
      S12istar[i] = drop(crossprod((mean(X[,i]) - mean(X)), (mean(Z[,i]) - mean(Z))))
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
    
    if (Procedure == 1){
      if (Case == "A"){
        mt1_bopt_num = 2*(f1*S01star + 1/n*fmr*S01bar) - 
          mean(X)/mean(Z)*(f*S02star + 1/n*fm*S02bar) -
          mean(Y)/mean(Z)*(f1*S12star + 1/n*fm*S12bar) +
          mean(Y)*mean(X)/(2*mean(Z)^2)*(f*S2Star_2 + 1/n*fm*S2bar_2)
        
        mt1_bopt_den = 2*(f1*S1Star_2 + 1/n*fmr*S1bar_2) +
          1/2*(mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) -
          2*mean(X)/mean(Z)*(f1*S12star + 1/n*fm*S12bar)
        
        mt1_bopt = mt1_bopt_num / mt1_bopt_den
        
        MT1opt = (f*S0Star_2 + 1/n*fmr*S0bar_2) + mt1_bopt^2 *(f1*S1Star_2 + 1/n*fmr*S1bar_2) +
          1/4*((mean(Y) - mt1_bopt*mean(X))/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) - 
          2*mt1_bopt*(f1*S01star + 1/n*fmr*S01bar) +
          ((mt1_bopt*mean(X) - mean(Y))/mean(Z))*(f*S02star + 1/n*fm*S02bar) +
          ((mt1_bopt*mean(Y) - mt1_bopt^2*mean(X))/mean(Z))*(f1*S12star + 1/n*fm*S12bar)
        
      } else if (Case == "B"){
        mt1_bopt_num = 2*(f*S01star + 1/n*fmr*S01bar) - 
          mean(X)/mean(Z)*(f*S02star + 1/n*fm*S02bar) -
          mean(Y)/mean(Z)*(f*S12star + 1/n*fm*S12bar) +
          mean(Y)*mean(X)/(2*mean(Z)^2)*(f*S2Star_2 + 1/n*fm*S2bar_2)
        
        mt1_bopt_den = 2*((f + fPrime)*S1Star_2 + 1/n*fmr*S1bar_2) +
          1/2*(mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) -
          2*mean(X)/mean(Z)*(f*S12star + 1/n*fm*S12bar)
        
        mt1_bopt = mt1_bopt_num / mt1_bopt_den
        
        MT1opt = (f*S0Star_2 + 1/n*fmr*S0bar_2) + mt1_bopt^2 *((f+ fPrime)*S1Star_2 + 1/n*fmr*S1bar_2) +
          1/4*((mean(Y) - mt1_bopt*mean(X))/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) - 
          2*mt1_bopt*(f*S01star + 1/n*fmr*S01bar) +
          ((mt1_bopt*mean(X) - mean(Y))/mean(Z))*(f*S02star + 1/n*fm*S02bar) +
          ((mt1_bopt*mean(Y) - mt1_bopt^2*mean(X))/mean(Z))*(f*S12star + 1/n*fm*S12bar)
        
      }
      
    }else if (Procedure == 2){
      mt1_bopt_num = mean(X)/mean(Z)*(f*S02star + 1/n*fmPrime*S02bar) -
        2/n*(fmPrimeR - fm)*S01bar
      
      mt1_bopt_den = 1/2*(mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fmPrime*S2bar_2) -
        2/n*(fmPrimeR - fm)*S1bar_2
      
      mt1_bopt = mt1_bopt_num / mt1_bopt_den
      
      MT1opt = (f*S0Star_2 + 1/n*fm*S0bar_2) - mt1_bopt^2/n *(fmPrimeR - fm)*S1bar_2 +
        1/4*(mt1_bopt*mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fmPrime*S2bar_2) + 
        2*mt1_bopt/n*(fmPrimeR -fm)*S01bar -
        mt1_bopt*mean(X)/mean(Z)*(f*S02star + 1/n*fmPrime*S02bar) 
    }
    
    
    return (list(clustersY = clustersY, indexOfSampleSpace = indexOfSampleSpace, 
                 indexOfFinalSample = indexOfFinalSample, sampledClustersY = sampledClustersY, Y = Y,
                 vYbarnm = vYbarnm, MT1opt = MT1opt))
    
  }
}

# This function performs the task multiple times and takes the average
replicateSampling <- function(nRep, M, N, m, n, r, mPrime,
                              nPrime, mu, Sigma,
                              Case, Procedure, seed_num =4113){
  vYbarnm_list = c()
  MT1opt_list = c()
  
  set.seed(seed_num)
  for (i in 1:nRep){
    rep_i = TwoStageClusterSampling(M, N, m, n, r, mPrime,
                                    nPrime, mu, Sigma,
                                    Case, Procedure, seed_num)
    vYbarnm_list[i] = rep_i$vYbarnm
    MT1opt_list[i] = rep_i$MT1opt
    
  }
  return(c(list(vYbarnm = mean(vYbarnm_list), MT1opt = mean(MT1opt_list),
                vYbarnm_list = vYbarnm_list[1:10], 
                MT1opt_list = MT1opt_list[1:10], anyNegative = any(MT1opt_list<0))))
}

Sigma <- matrix(c(50, 2, 3, 2, 100, 7, 3, 7 , 50), 3,3)

aa <- replicateSampling(nRep = 100, M = 15, N = 20, m = 5, 
                        n = 6, r = 2, mPrime = 7,
                        nPrime = 7, mu = c(20, 50, 40), Sigma = Sigma,
                        Case = "B", Procedure = 2, seed_num = 543)

aa




a <- TwoStageClusterSampling(M = 15, N = 20, m = 5, n = 6, r = 2, mPrime = 7,
                             nPrime = 7, mu = c(20, 50, 40), Sigma = Sigma,
                             Case = "B", Procedure = 1)
a
