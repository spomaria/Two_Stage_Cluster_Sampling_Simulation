# This script is meant to run simulation for 
# Two stage Cluster Sampling in the presence of 
# Sampling Non-response

# @param N is the number of clusters in the pop
# @param n is the number of clusters to be sampled
# @param M is the number of units in each cluster
# @param p is the probability of non-response (could vary across clusters or remain constant)
# @param m is the number of units to be sampled in each selected cluster
# @param mu is a vector of means for the variables y, x, and z
# @param Sigma is a symmetric matrix. It is the variance-covariance matrix of the variables
# @param Case takes only two possible values: "A" or "B"

TwoStageClusterSampling <- function(
    M, N, m, n, p, mPrime, nPrime, mu, Sigma, Case = "A", Procedure, seed_num=4113  
){
  # Check that the number of clusters to be sampled does not 
  # exceed the total number of clusters in the pop
  if ( n > N | n == N) return("Error: n should be less than N")
  # Check that the number of units to be sampled per cluster 
  # does not exceed the cluster size
  else if (m > M | m == M) return("Error: m should be less than M")
  else if (m > mPrime | m == mPrime) return("Error: m should be less than mPrime")
  else if (n > nPrime | n == nPrime) return("Error: n should be less than nPrime")
  # Check that the number of non-response per cluster 
  # does not exceed the sample size per cluster
  else if (p <= 0 | p >= 1) return("Error: p should range exclusively between 0 and 1")
  else if (!(length(mu) == ncol(Sigma)) | !(length(mu) == nrow(Sigma))) return("length of mu and dim of Sigma inconsistent")
  else {
    # Defining some of the derived constants needed in the various computations
    #p = r/m
    q = 1 - p
    f = 1/n - 1/N
    fPrime = 1/nPrime - 1/N
    fmPrime = 1/mPrime - 1/M
    f1 = 1/n - 1/nPrime
    fm = 1/m - 1/M
    fmr = 1/(m * q + 2 * p) - 1/M
    fmPrimeR = 1/(mPrime * q + 2 * p) - 1/M
    
    totalSim = M * N
    # Simulating the entire pop
    # Load the 'MASS' library to enable us generate 
    # random numbers from a multivariate normal distribution
    library(MASS)
    #set.seed(seed_num)
    pop = mvrnorm(totalSim, mu = mu, Sigma = Sigma)
    # detach the library after the simulation
    detach("package:MASS", unload = TRUE)
    
    # Demarcation of the various clusters in the pop
    # We split the pop into clusters by converting
    # the pop into a matrix where each column 
    # represents a cluster
    # This is carried out for all three variables
    
    mu_y = mu[1]; mu_x = mu[2]; mu_z = mu[3]
    sig_y = sqrt(Sigma[1,1]); sig_x = sqrt(Sigma[2,2]); sig_z = sqrt(Sigma[3,3]);
    rho_xy = 0.7; rho_xz = 0.8
    
    clustersY = matrix(mu_y + sig_y *(rho_xy*pop[,2] + sqrt(1-rho_xy^2)*pop[,1]), ncol = N)
    clustersX = matrix(mu_x + sig_x *pop[,2], ncol = N)
    clustersZ = matrix(mu_z + sig_z *(rho_xz*pop[,2] + sqrt(1-rho_xz^2)*pop[,3]), ncol = N)
    
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
    if (Procedure ==1) indexOfFinalSample = sample(sampleSpace, m, replace = TRUE)
    else if (Procedure ==2){
      ssu1Space = sample(sampleSpace, mPrime)
      # indexOfFinalSample = sample(ssu1Space, m)
      # Since non-response is accounted for by sampling again, we allow the 
      # second sampling to be done with replacement
      indexOfFinalSample = sample(ssu1Space, m, replace = TRUE)
    }
    # non-response space
    # To make this exercise more practical, we shall vary the index
    # of non-responsive units across the different clusters
    
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
      Yi = Y[,i] - mean(Y[,i]); Xi = X[,i] - mean(X[,i])
      # the 'crossprod' function computes the cross product of two vectors and returns the 
      # result as a 1 by 1 matrix
      # the 'drop' function converts the 1 by 1 matrix to a scalar
      S01i[i] = sum(drop(crossprod(Yi, Xi)))/(m-1)
    }
    S01bar <- mean(S01i)
    
    S01istar <- c()
    for (i in 1:n){
      S01istar[i] = drop(crossprod((mean(Y[,i]) - mean(Y)), (mean(X[,i]) - mean(X))))
    }
    S01star <- sum(S01istar)/(n-1)
    
    # Components relating to variables Y and Z
    S02i <- c()
    for (i in 1:n){
      Yi = Y[,i] - mean(Y[,i]); Zi = Z[,i] - mean(Z[,i])
      S02i[i] = sum(drop(crossprod(Yi, Zi)))/(m-1)
    }
    S02bar <- mean(S02i)
    
    S02istar <- c()
    for (i in 1:n){
      S02istar[i] = drop(crossprod((mean(Y[,i]) - mean(Y)), (mean(Z[,i]) - mean(Z))))
    }
    S02star <- sum(S02istar)/(n-1)
    
    # Components relating to variables X and Z
    S12i <- c()
    for (i in 1:n){
      Xi = X[,i] - mean(X[,i]); Zi = Z[,i] - mean(Z[,i])
      S12i[i] = sum(drop(crossprod(Xi, Zi)))/(m-1)
    }
    S12bar <- mean(S12i)
    
    S12istar <- c()
    for (i in 1:n){
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
      swsquareUnitsY[i] = sum((Y[,i] - mean(Y[,i]))^2)
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
    
    
    return (list(clustersY = clustersY, indexOfFinalSample = indexOfFinalSample, 
                 indexOfFinalSample = indexOfFinalSample, sampledClustersY = sampledClustersY, Y = Y,
                 vYbarnm = vYbarnm, MT1opt = MT1opt))
    
  }
}

# This function performs the task multiple times and takes the average
replicateSampling <- function(nRep, M, N, m, n, p, mPrime,
                              nPrime, mu, Sigma,
                              Case, Procedure, seed_num =4113){
  vYbarnm_list = c()
  MT1opt_list = c()
  
  set.seed(seed_num)
  for (i in 1:nRep){
    rep_i = TwoStageClusterSampling(M, N, m, n, p, mPrime,
                                    nPrime, mu, Sigma,
                                    Case, Procedure, seed_num)
    vYbarnm_list[i] = rep_i$vYbarnm
    MT1opt_list[i] = rep_i$MT1opt
    
  }
  PRE = mean(vYbarnm_list)/mean(MT1opt_list)*100
  LOSS = (mean(MT1opt_list) - mean(vYbarnm_list))/mean(MT1opt_list)*100
  
  return(c(list(vYbarnm = mean(vYbarnm_list), MT1opt = mean(MT1opt_list),
                vYbarnm_list = vYbarnm_list[1:10], 
                MT1opt_list = MT1opt_list[1:10], anyNegative = any(MT1opt_list<0),
                PRE = PRE, LOSS = LOSS)))
}

Sigma <- matrix(c(20, 0, 0, 0, 60, 0, 0, 0 , 10), 3,3)

aa <- replicateSampling(nRep = 100, M = 15, N = 20, m = 7, 
                        n = 5, p = .05, mPrime = 8,
                        nPrime = 7, mu = c(20, 50, 40), Sigma = Sigma,
                        Case = "B", Procedure = 1, seed_num = 5431)

aa


