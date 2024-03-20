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
# @param m_error is a boolean variable that indicates if measurement error is present or not
# @param aux_param_option is a value that specifies the preferred option for auxiliary parameters

# To use this function, ensure you have the 'moments' package installed on your local machine


TwoStageClusterSampling <- function(
    M, N, m, n, p, mPrime, nPrime, mu, Sigma, Case = "A", Procedure, seed_num=4113,
    m_error = TRUE, aux_param_option = 1, data_matrix = NA
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
  else if (any(p <= 0) | any(p >= 1)) return("Error: p should range exclusively between 0 and 1")
  else if (!(length(p) == 1 | length(p) == n)) return("Error: p should either contain 1 or n entries")
  else if (!(length(mu) == ncol(Sigma)) | !(length(mu) == nrow(Sigma))) return("length of mu and dim of Sigma inconsistent")
  #else if (!(!is.matrix(data_matrix) | ncol(data_matrix) == 3)) return("Error: data_matrix should have 3 columns if it exists")
  else {
    # --- Defining some of the derived constants needed in the various computations
    #p = r/m
    q = 1 - p
    f = 1/n - 1/N
    fPrime = 1/nPrime - 1/N
    fmPrime = 1/mPrime - 1/M
    f1 = 1/n - 1/nPrime
    fm = 1/m - 1/M
    fmr = 1/(m * q + 2 * p) - 1/M
    fmPrimeR = 1/(mPrime * q + 2 * p) - 1/M
    
    if (!is.matrix(data_matrix)){
      # Where no real life data is provided, the computer simulates data from the 
      # normal distribution
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
      
      clustersY = matrix(mu_y + sig_y *(rho_xy*pop[,2] + sqrt(1-rho_xy^2)*pop[,1]) + pop[,4], ncol = N)
      clustersX = matrix(mu_x + sig_x *pop[,2] + pop[,5], ncol = N)
      clustersZ = matrix(mu_z + sig_z *(rho_xz*pop[,2] + sqrt(1-rho_xz^2)*pop[,3]) + pop[,6], ncol = N)
      # creating the errors relating to the variables
      clustersU = matrix(pop[,4], ncol = N)
      clustersV = matrix(pop[,5], ncol = N)
      clustersW = matrix(pop[,6], ncol = N)
      
    } else {
      # Using real life data
      clustersY = matrix(data_matrix[,1], ncol = N, byrow = FALSE)
      clustersX = matrix(data_matrix[,2], ncol = N, byrow = FALSE)
      clustersZ = matrix(data_matrix[,3], ncol = N, byrow = FALSE)
      
    }
    
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
    
    if (!is.matrix(data_matrix) & m_error){
      # The errors of the variables
      sampledClustersU = matrix(NA, ncol = n, nrow = M)
      sampledClustersV = matrix(NA, ncol = n, nrow = M)
      sampledClustersW = matrix(NA, ncol = n, nrow = M)
    }
    for (i in 1:n){
      # The main variable Y
      sampledClustersY[,i] = clustersY[,fsuClusters[i]]
      # The auxiliary variable X
      sampledClustersX[,i] = clustersX[,fsuClusters[i]]
      # The auxiliary variable Z
      sampledClustersZ[,i] = clustersZ[,fsuClusters[i]]
      
      if (!is.matrix(data_matrix) & m_error){
        # The error of the main variable Y
        sampledClustersU[,i] = clustersU[,fsuClusters[i]]
        # The error of the auxiliary variable X
        sampledClustersV[,i] = clustersV[,fsuClusters[i]]
        # The error of the auxiliary variable Z
        sampledClustersW[,i] = clustersW[,fsuClusters[i]]
      }
      
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
    
    if (!is.matrix(data_matrix) & m_error){
      # The error of the main variable Y
      U = matrix(NA, ncol = n, nrow = m)
      # The error of the auxiliary variable X
      V = matrix(NA, ncol = n, nrow = m)
      # The error of the auxiliary variable Z
      W = matrix(NA, ncol = n, nrow = m)
    }
    
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
        
        if (!is.matrix(data_matrix) & m_error){
          # The error of the main variable Y
          U[j,i] = sampledClustersU[indexOfFinalSample[j],i]
          # The error of the auxiliary variable X
          V[j,i] = sampledClustersV[indexOfFinalSample[j],i]
          # The error of the auxiliary variable Z
          W[j,i] = sampledClustersW[indexOfFinalSample[j],i]
        }
        
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
    
    
    if (!is.matrix(data_matrix) & m_error){
      
      # Considering the variance component that incorporates the error term
      sbsquareUnitsU = c()
      for (i in 1:n){
        sbsquareUnitsU[i] = (mean(U[,i]) - mean(U))^2
      }
      sbsquare_e = 1/(n - 1) *sum(sbsquareUnitsU)
      
      # Calculating mean square within the clusters
      swsquareUnitsU = c()
      for (i in 1:n){
        swsquareUnitsU[i] = sum((U[,i] - mean(U[,i]))^2)
      }
      swsquare_e = 1/n * 1/(m - 1) *sum(swsquareUnitsU)
      
      
      #--- Incorporating the error components
      # Computing the several components of the estimators
      # Components relating to error variable U
      Sui_2 <- c()
      for (i in 1:n){
        Sui_2[i] = sum((U[,i] - mean(U[,i]))^2)/(m-1)
      }
      Subar_2 <- mean(Sui_2)
      
      Suibar <- c()
      for (i in 1:n){
        Suibar[i] = mean(U[,i])
      }
      SuStar_2 <- sum((Suibar - mean(U))^2)/(n-1)
      
      # Components relating to error variable V
      Svi_2 <- c()
      for (i in 1:n){
        Svi_2[i] = sum((V[,i] - mean(V[,i]))^2)/(m-1)
      }
      Svbar_2 <- mean(Svi_2)
      
      Svibar <- c()
      for (i in 1:n){
        Svibar[i] = mean(V[,i])
      }
      SvStar_2 <- sum((Svibar - mean(V))^2)/(n-1)
      
      # Components relating to error variable W
      Swi_2 <- c()
      for (i in 1:n){
        Swi_2[i] = sum((W[,i] - mean(W[,i]))^2)/(m-1)
      }
      Swbar_2 <- mean(Swi_2)
      
      Swibar <- c()
      for (i in 1:n){
        Swibar[i] = mean(W[,i])
      }
      SwStar_2 <- sum((Swibar - mean(W))^2)/(n-1)
      
      # Components relating to error variables U and V
      Suvi <- c()
      for (i in 1:n){
        Ui = U[,i] - mean(U[,i]); Vi = V[,i] - mean(V[,i])
        # the 'crossprod' function computes the cross product of two vectors and returns the 
        # result as a 1 by 1 matrix
        # the 'drop' function converts the 1 by 1 matrix to a scalar
        Suvi[i] = sum(drop(crossprod(Ui, Vi)))/(m-1)
      }
      Suvbar <- mean(Suvi)
      
      Suvistar <- c()
      for (i in 1:n){
        Suvistar[i] = drop(crossprod((mean(U[,i]) - mean(U)), (mean(V[,i]) - mean(V))))
      }
      Suvstar <- sum(Suvistar)/(n-1)
      
      # Components relating to error variables U and W
      Suwi <- c()
      for (i in 1:n){
        Ui = U[,i] - mean(U[,i]); Wi = W[,i] - mean(W[,i])
        Suwi[i] = sum(drop(crossprod(Ui, Wi)))/(m-1)
      }
      Suwbar <- mean(Suwi)
      
      Suwistar <- c()
      for (i in 1:n){
        Suwistar[i] = drop(crossprod((mean(U[,i]) - mean(U)), (mean(W[,i]) - mean(W))))
      }
      Suwstar <- sum(Suwistar)/(n-1)
      
      # Components relating to error variables V and W
      Svwi <- c()
      for (i in 1:n){
        Vi = V[,i] - mean(V[,i]); Wi = W[,i] - mean(W[,i])
        Svwi[i] = sum(drop(crossprod(Vi, Wi)))/(m-1)
      }
      Svwbar <- mean(Svwi)
      
      Svwistar <- c()
      for (i in 1:n){
        Svwistar[i] = drop(crossprod((mean(V[,i]) - mean(V)), (mean(W[,i]) - mean(W))))
      }
      Svwstar <- sum(Svwistar)/(n-1)  
      
    }else {
      # Where the error term is not incorporated
      sbsquare_e = swsquare_e = 0
      
      Subar_2 = SuStar_2 = Svbar_2 = SvStar_2 = Swbar_2 = 0
      SwStar_2 = Suvbar = Suvstar = Suwbar = Suwstar = 0 
      Svwbar = Svwstar = 0
      Sui_2 = Svi_2 = Swi_2 = Suvi_2 = Suwi_2 = Svwi_2 = 0
      Suvi = Suwi = Svwi = 0
    }
    
    # The variance to be used as baseline for comparison with other MSEs
    vYbarnm_e = f * (sbsquare + sbsquare_e) + 1/n *fm*(swsquare + swsquare_e)
    
    # Compute other constants relating to Skewness and Kurtosis
    # We load the moments library to enable us calculate the skewness and kurtosis
    library(moments)
    # Coefficient of Skewness of variable X
    B1_X2 = skewness(X)
    # Coefficient of Kurtosis of variable X
    B2_X2 = kurtosis(X)
    
    detach("package:moments", unload = TRUE)
    # Standard deviation of variable X
    S_X2 = sd(X)
    # Coefficient of variation of variable X
    C_X2 = S_X2/mean(X)
    
    # Setting the Auxiliary Parameter Options
    if (aux_param_option == 1){
      A_X2 = 1; B_X2 = 0
    } else if (aux_param_option == 2){
      A_X2 = 1; B_X2 = B1_X2
    } else if (aux_param_option == 3){
      A_X2 = 1; B_X2 = B2_X2
    } else if (aux_param_option == 4){
      A_X2 = 1; B_X2 = C_X2
    } else if (aux_param_option == 5){
      A_X2 = 1; B_X2 = S_X2
    } else if (aux_param_option == 6){
      A_X2 = B1_X2; B_X2 = B2_X2
    } else if (aux_param_option == 7){
      A_X2 = B1_X2; B_X2 = C_X2
    } else if (aux_param_option == 8){
      A_X2 = B1_X2; B_X2 = S_X2
    } else if (aux_param_option == 9){
      A_X2 = B2_X2; B_X2 = B1_X2
    } else if (aux_param_option == 10){
      A_X2 = B2_X2; B_X2 = C_X2
    } else if (aux_param_option == 11){
      A_X2 = B2_X2; B_X2 = S_X2
    } else if (aux_param_option == 12){
      A_X2 = C_X2; B_X2 = B1_X2
    } else if (aux_param_option == 13){
      A_X2 = C_X2; B_X2 = B2_X2
    } else if (aux_param_option == 14){
      A_X2 = C_X2; B_X2 = S_X2
    } else if (aux_param_option == 15){
      A_X2 = S_X2; B_X2 = B1_X2
    } else if (aux_param_option == 16){
      A_X2 = S_X2; B_X2 = B2_X2
    } else if (aux_param_option == 17){
      A_X2 = S_X2; B_X2 = C_X2
    } else{ return("Auxiliary Parameter Option should lie between 1 - 17")}
    
    theta_X = A_X2*mean(Z)/(A_X2*mean(Z) + B_X2)
    
    if (length(p) == 1){
      
      # Computations for Ibrahim's Estimator
      c0e = (f*(S0Star_2 + SuStar_2) + 1/n*1/n*fmr*(S0bar_2 + Subar_2))/mean(Y)^2
      
      c1e = (f*(S1Star_2 + SvStar_2) + 1/n*1/n*fmr*(S1bar_2 + Svbar_2))/mean(X)^2
      
      c2e = fPrime*(S1Star_2 + SvStar_2)/mean(X)^2
      
      c3e = (f*(S2Star_2 + SwStar_2) + 1/n*fm*(S2bar_2 + Swbar_2))/mean(Z)^2
      
      c4e = (f*(S01star + Suvstar) + 1/n*1/n*fmr*(S01bar + Suvbar))/(mean(Y)*mean(X))
      
      c5e = fPrime*(S01star + Suvstar)/(mean(Y)*mean(X))
      
      c6e = (f*(S02star + Suwstar) + 1/n*fm*(S02bar + Suwbar))/(mean(Y)*mean(Z))
      
      c7e = (f*(S12star + Svwstar) + 1/n*fm*(S12bar + Svwbar))/(mean(X)*mean(Z))
      
      c8e = fPrime*(S12star + Svwstar)/(mean(X)*mean(Z))
      
      c9e = (f*(S0Star_2 + SuStar_2) + 1/n*fm*(S0bar_2 + Subar_2))/mean(Y)^2
      
      c10e = (f*(S1Star_2 + SvStar_2) + 1/n*fmPrimeR*(S1bar_2 + Svbar_2))/mean(X)^2
      
      c11e = (f*(S2Star_2 + SwStar_2) + 1/n*fmPrime*(S2bar_2 + Swbar_2))/mean(Z)^2
      
      c12e = (f*(S1Star_2 + SvStar_2) + 1/n*fm*(S1bar_2 + Svbar_2))/mean(X)^2
      
      c13e = (f*(S01star + Suvstar) + 1/n*fmPrimeR*(S01bar + Suvbar))/(mean(Y)*mean(X))
      
      c14e = (f*(S02star + Suwstar) + 1/n*fmPrime*(S02bar + Suwbar))/(mean(Y)*mean(Z))
      
      c15e = (f*(S01star + Suvstar) + 1/n*fm*(S01bar + Suvbar))/(mean(Y)*mean(X))
      
      c16e = (f*(S12star + Svwstar) + 1/n*fmPrime*(S12bar + Svwbar))/(mean(X)*mean(Z))
      
      c17e = (f*(S1Star_2 + SvStar_2) + 1/n*fmPrimeR*(S1bar_2 + Svbar_2))/mean(X)^2
      
      c18e = c16e
      
      if (Procedure == 1){
        if (Case == "A"){
          #-- Note that
          # 'mt1_... stands for Maji while
          # 'mt2_... stands for Ibrahim
          
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
          
          # Now the estimator proposed by Ibrahim
          mt2_bopt_num = mean(Y)*(
            p^2*(theta_X + 1)^2*c3e + c4e - c5e -p*(theta_X + 1)* (c6e + c7e - c8e)  
          )
          
          mt2_bopt_den = mean(X)*(
            c1e - c2e + p^2*(theta_X + 1)^2*c3e - 2*p*(theta_X + 1)*(c7e - c8e)
          )
          
          mt2_bopt = mt2_bopt_num / mt2_bopt_den
          
          MT2opt = mean(Y)^2*(
            c0e + p^2*(theta_X + 1)^2*c3e - 2*p*(theta_X + 1)*c6e
          ) - 2*mt2_bopt*mean(Y)*mean(X)*(
            p^2*(theta_X + 1)^2*c3e + c4e - c5e - p*(theta_X + 1)* (c6e + c7e - c8e)
          ) + mt2_bopt^2*mean(X)^2*(
            c1e - c2e + p^2*(theta_X + 1)^2*c3e - 2*p*(theta_X + 1)*(c7e - c8e)
          )
          
          
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
          
          # Computation for Ibrahim's Estimator
          mt2_bopt_num = mean(Y)*(
            p^2*(theta_X + 1)^2*c3e + c4e -(theta_X + 1)* p*(c6e + c7e)  
          )
          
          mt2_bopt_den = mean(X)*(
            c1e + c2e + p^2*(theta_X + 1)^2*c3e - 2*p*(theta_X + 1)*c7e
          )
          
          mt2_bopt = mt2_bopt_num / mt2_bopt_den
          
          MT2opt = mean(Y)^2*(
            c0e + p^2*(theta_X + 1)^2*c3e - 2*p*(theta_X + 1)*c6e
          ) - 2*mt2_bopt*mean(Y)*mean(X)*(
            p^2*(theta_X + 1)^2*c3e + c4e - (theta_X + 1)* p*(c6e + c7e)
          ) + mt2_bopt^2*mean(X)^2*(
            c1e + c2e + p^2*(theta_X + 1)^2*c3e - 2*p*(theta_X + 1)*c7e
          )
          
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
        
        # Computation of Ibrahim's estimator
        mt2_bopt_num = mean(Y)*(
          (theta_X + 1)* p*c14e - c13e + c15e  
        )
        
        mt2_bopt_den = mean(X)*(
          c10e + p^2*(theta_X + 1)^2*c11e + c12e - 2*c17e
        )
        
        mt2_bopt = mt2_bopt_num / mt2_bopt_den
        
        MT2opt = mean(Y)^2*c9e + mt2_bopt^2*mean(X)^2*(
          c10e + p^2*(theta_X + 1)^2*c11e + c12e - 2*c17e
        ) + 2*mt2_bopt*mean(Y)*mean(X)*(
          c13e - (theta_X + 1)* p*c14e - c15e
        )
        
      }
      
    } else{
      
      # Computations for Ibrahim's Estimator
      c0e = (f*(S0Star_2 + SuStar_2) + 
               1/n*1/n*sum(drop(crossprod(fmr,S0i_2 + Sui_2))))/mean(Y)^2
      
      c1e = (f*(S1Star_2 + SvStar_2) + 
               1/n*1/n*sum(drop(crossprod(fmr,S1i_2 + Svi_2))))/mean(X)^2
      
      c2e = fPrime*(S1Star_2 + SvStar_2)/mean(X)^2
      
      c3e = (f*(S2Star_2 + SwStar_2) + 1/n*fm*(S2bar_2 + Swbar_2))/mean(Z)^2
      
      c4e = (f*(S01star + Suvstar) + 
               1/n*1/n*sum(drop(crossprod(fmr,S01i + Suvi))))/(mean(Y)*mean(X))
      
      c5e = fPrime*(S01star + Suvstar)/(mean(Y)*mean(X))
      
      c6e = (f*(S02star + Suwstar) + 1/n*fm*(S02bar + Suwbar))/(mean(Y)*mean(Z))
      
      c7e = (f*(S12star + Svwstar) + 1/n*fm*(S12bar + Svwbar))/(mean(X)*mean(Z))
      
      c8e = fPrime*(S12star + Svwstar)/(mean(X)*mean(Z))
      
      c9e = (f*(S0Star_2 + SuStar_2) + 1/n*fm*(S0bar_2 + Subar_2))/mean(Y)^2
      
      c10e = (f*(S1Star_2 + SvStar_2) + 
               1/n*1/n*sum(drop(crossprod(fmPrimeR,S1i_2 + Svi_2))))/mean(X)^2
      
      c11e = (f*(S2Star_2 + SwStar_2) + 1/n*fmPrime*(S2bar_2 + Swbar_2))/mean(Z)^2
      
      c12e = (f*(S1Star_2 + SvStar_2) + 1/n*fm*(S1bar_2 + Svbar_2))/mean(X)^2
      
      c13e = (f*(S01star + Suvstar) + 
                1/n*1/n*sum(drop(crossprod(fmPrimeR,S01i + Suvi))))/(mean(Y)*mean(X))
      
      c14e = (f*(S02star + Suwstar) + 1/n*fmPrime*(S02bar + Suwbar))/(mean(Y)*mean(Z))
      
      c15e = (f*(S01star + Suvstar) + 1/n*fm*(S01bar + Suvbar))/(mean(Y)*mean(X))
      
      c16e = (f*(S12star + Svwstar) + 1/n*fmPrime*(S12bar + Svwbar))/(mean(X)*mean(Z))
      
      c17e = c10e
      
      
      if (Procedure == 1){
        if (Case == "A"){
          mt1_bopt_num = 2*(f1*S01star + 1/n*1/n*sum(drop(crossprod(fmr,S01i)))) - 
            mean(X)/mean(Z)*(f*S02star + 1/n*fm*S02bar) -
            mean(Y)/mean(Z)*(f1*S12star + 1/n*fm*S12bar) +
            mean(Y)*mean(X)/(2*mean(Z)^2)*(f*S2Star_2 + 1/n*fm*S2bar_2)
          
          mt1_bopt_den = 2*(f1*S1Star_2 + 1/n*1/n*sum(drop(crossprod(fmr,S1i_2)))) +
            1/2*(mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) -
            2*mean(X)/mean(Z)*(f1*S12star + 1/n*fm*S12bar)
          
          mt1_bopt = mt1_bopt_num / mt1_bopt_den
          
          MT1opt = (f*S0Star_2 + 1/n*1/n*sum(drop(crossprod(fmr,S0i_2)))) + 
            mt1_bopt^2 *(f1*S1Star_2 + 1/n*1/n*sum(drop(crossprod(fmr,S1i_2)))) +
            1/4*((mean(Y) - mt1_bopt*mean(X))/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) - 
            2*mt1_bopt*(f1*S01star + 1/n*1/n*sum(drop(crossprod(fmr,S01i)))) +
            ((mt1_bopt*mean(X) - mean(Y))/mean(Z))*(f*S02star + 1/n*fm*S02bar) +
            ((mt1_bopt*mean(Y) - mt1_bopt^2*mean(X))/mean(Z))*(f1*S12star + 1/n*fm*S12bar)
          
          
          # Computations for Ibrahim's Estimator
          
          mt2_bopt_num = mean(Y)*(
            (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e + c4e - 
              c5e -sum(p)/(n*m)*(theta_X + 1)* (c6e + c7e - c8e)  
          )
          
          mt2_bopt_den = mean(X)*(
            c1e - c2e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e - 
              2*sum(p)/(n*m)*(theta_X + 1)* (c7e - c8e)
          )
          
          mt2_bopt = mt2_bopt_num / mt2_bopt_den
          
          MT2opt = mean(Y)^2*(
            c0e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2* c3e - 
              2*sum(p)/(n*m)*(theta_X + 1)* c6e
          ) - 2*mt2_bopt*mean(Y)*mean(X)*(
            (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e + c4e - c5e - 
              sum(p)/(n*m)*(theta_X + 1)*(c6e + c7e - c8e)
          ) + mt2_bopt^2*mean(X)^2*(
             c1e - c2e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e - 
               2*sum(p)/(n*m)*(theta_X + 1)* (c7e - c8e)
          )
          
        } else if (Case == "B"){
          mt1_bopt_num = 2*(f*S01star + 1/n*1/n*sum(drop(crossprod(fmr,S01i)))) - 
            mean(X)/mean(Z)*(f*S02star + 1/n*fm*S02bar) -
            mean(Y)/mean(Z)*(f*S12star + 1/n*fm*S12bar) +
            mean(Y)*mean(X)/(2*mean(Z)^2)*(f*S2Star_2 + 1/n*fm*S2bar_2)
          
          mt1_bopt_den = 2*((f + fPrime)*S1Star_2 + 1/n*1/n*sum(drop(crossprod(fmr,S1i_2)))) +
            1/2*(mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) -
            2*mean(X)/mean(Z)*(f*S12star + 1/n*fm*S12bar)
          
          mt1_bopt = mt1_bopt_num / mt1_bopt_den
          
          MT1opt = (f*S0Star_2 + 1/n*1/n*sum(drop(crossprod(fmr,S0i_2)))) + 
            mt1_bopt^2 *((f+ fPrime)*S1Star_2 + 1/n*1/n*sum(drop(crossprod(fmr,S1i_2)))) +
            1/4*((mean(Y) - mt1_bopt*mean(X))/mean(Z))^2*(f*S2Star_2 + 1/n*fm*S2bar_2) - 
            2*mt1_bopt*(f*S01star + 1/n*1/n*sum(drop(crossprod(fmr,S01i)))) +
            ((mt1_bopt*mean(X) - mean(Y))/mean(Z))*(f*S02star + 1/n*fm*S02bar) +
            ((mt1_bopt*mean(Y) - mt1_bopt^2*mean(X))/mean(Z))*(f*S12star + 1/n*fm*S12bar)
          
          # Computations for Ibrahim's Estimator
          
          mt2_bopt_num = mean(Y)*(
            (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e + c4e -
              sum(p)/(n*m)*(theta_X + 1)* (c6e + c7e)  
          )
          
          mt2_bopt_den = mean(X)*(
            c1e + c2e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e - 
              2*sum(p)/(n*m)*(theta_X + 1)*c7e
          )
          
          mt2_bopt = mt2_bopt_num / mt2_bopt_den
          
          MT2opt = mean(Y)^2*(
            c0e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2* c3e - 
              2*sum(p)/(n*m)*(theta_X + 1)* c6e
          ) - 2*mt2_bopt*mean(Y)*mean(X)*(
            (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e + c4e - 
              sum(p)/(n*m)*(theta_X + 1)*(c6e + c7e)
          ) + mt2_bopt^2*mean(X)^2*(
            c1e + c2e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c3e - 
              2*sum(p)/(n*m)*(theta_X + 1)*c7e
          )
          
        }
        
      }else if (Procedure == 2){
        mt1_bopt_num = mean(X)/mean(Z)*(f*S02star + 1/n*fmPrime*S02bar) -
          2/n*(1/n*sum(drop(crossprod(fmPrimeR, S01i))) - fm*S01bar)
        
        mt1_bopt_den = 1/2*(mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fmPrime*S2bar_2) -
          2/n*(1/n*sum(drop(crossprod(fmPrimeR, S1i_2))) - fm*S1bar_2)
        
        mt1_bopt = mt1_bopt_num / mt1_bopt_den
        
        MT1opt = (f*S0Star_2 + 1/n*fm*S0bar_2) - 
          mt1_bopt^2/n *(1/n*sum(drop(crossprod(fmPrimeR, S1i_2))) - fm*S1bar_2) +
          1/4*(mt1_bopt*mean(X)/mean(Z))^2*(f*S2Star_2 + 1/n*fmPrime*S2bar_2) + 
          2*mt1_bopt/n*(1/n*sum(drop(crossprod(fmPrimeR, S01i))) -fm*S01bar) -
          mt1_bopt*mean(X)/mean(Z)*(f*S02star + 1/n*fmPrime*S02bar) 
        
        # Computations for Ibrahim's Estimator
        
        mt2_bopt_num = mean(Y)*(
          sum(p)/(n*m)*(theta_X + 1)*c14e - c13e + c15e  
        )
        
        mt2_bopt_den = mean(X)*(
          c10e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c11e + c12e - 2*c17e
        )
        
        mt2_bopt = mt2_bopt_num / mt2_bopt_den
        
        MT2opt = mean(Y)^2*c9e + mt2_bopt^2*mean(X)^2*(
          c10e + (m*var(p) + 1/n*sum(p))*(theta_X + 1)^2*c11e + c12e - 2*c17e
        ) + 2*mt2_bopt*mean(Y)*mean(X)*(
          c13e - sum(p)/(n*m)*(theta_X + 1)*c14e - c15e
        )
        
      }
      
    }
  
  return (list(clustersY = clustersY, indexOfFinalSample = indexOfFinalSample, 
               indexOfFinalSample = indexOfFinalSample, sampledClustersY = sampledClustersY, Y = Y,
               vYbarnm_e = vYbarnm_e, MT1opt = MT1opt, MT2opt = MT2opt))
  }
}

# This function performs the sampling task multiple times and takes the average
replicateSampling <- function(nRep, M, N, m, n, p, mPrime,
                              nPrime, mu, Sigma,
                              Case, Procedure, seed_num =4113, m_error = TRUE,
                              aux_param_option = 1, data_matrix = NA){
  vYbarnm_elist = c()
  MT1opt_list = c()
  MT2opt_list = c()
  
  set.seed(seed_num)
  for (i in 1:nRep){
    rep_i = TwoStageClusterSampling(M, N, m, n, p, mPrime,
                                    nPrime, mu, Sigma,
                                    Case, Procedure, seed_num, m_error,
                                    aux_param_option, data_matrix)

    vYbarnm_elist[i] = rep_i$vYbarnm_e
    MT1opt_list[i] = rep_i$MT1opt
    MT2opt_list[i] = rep_i$MT2opt
    
  }
  # PRE = mean(vYbarnm_list)/mean(MT1opt_list)*100
  # LOSS = (mean(MT1opt_list) - mean(vYbarnm_list))/mean(MT1opt_list)*100
  PRE_Maji = mean(vYbarnm_elist)/mean(MT1opt_list)*100
  LOSS_Maji = (mean(MT1opt_list) - mean(vYbarnm_elist))/mean(MT1opt_list)*100
  
  PRE_Ibro = mean(vYbarnm_elist)/mean(MT2opt_list)*100
  LOSS_Ibro = (mean(MT2opt_list) - mean(vYbarnm_elist))/mean(MT2opt_list)*100
  
  Rel_Performance = matrix(c(PRE_Maji, PRE_Ibro, LOSS_Maji, LOSS_Ibro), nrow = 1, 
                           dimnames = list(c(), c("PRE_Maji", "PRE_Ibro", 
                                                  "LOSS_Maji", "LOSS_Ibro")))
  
  return(c(list(vYbarnm_elist = vYbarnm_elist[1:10], 
                MT1opt_list = MT1opt_list[1:10], anyNegative_Maji = any(MT1opt_list<0),
                anyNegative_Ibro = any(MT2opt_list<0), Rel_Performance = Rel_Performance)))
}

Sigma <- matrix(c(2, 0, 0, 0, 0, 0, 
                  0, 6, 0, 0, 0, 0, 
                  0, 0, 3, 0, 0, 0,
                  0, 0, 0, 5, 0, 0,
                  0, 0, 0, 0, 3, 0,
                  0, 0, 0, 0, 0, 7), 6,6)

aa <- replicateSampling(nRep = 100, M = 100, N = 10, m = 50, 
                        n = 5, p = c(rep(0.15,4), 0.05), mPrime = 60,
                        nPrime = 7, mu = c(20, 50, 40, 0, 0, 0), Sigma = Sigma,
                        Case = "B", Procedure = 2, seed_num = 531, m_error = TRUE,
                        aux_param_option = 4, data_matrix = NA)

ab <- replicateSampling(nRep = 100, M = 100, N = 10, m = 50, 
                        n = 5, p = c(rep(0.15,4), 0.05), mPrime = 60,
                        nPrime = 7, mu = c(20, 50, 40, 0, 0, 0), Sigma = Sigma,
                        Case = "B", Procedure = 2, seed_num = 531, m_error = FALSE,
                        aux_param_option = 4)

# aux_param_options that result in warning messages are: 
# 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16

aa$Rel_Performance
ab$Rel_Performance



# Using the real life data used by Maji
# Read in the data saved as a csv file
realData = read.csv("https://raw.githubusercontent.com/spomaria/Two_Stage_Cluster_Sampling_Simulation/main/CommHospitalsData.csv")
# realData = read.csv(file = "https://github.com/spomaria/Two_Stage_Cluster_Sampling_Simulation/blob/main/CommHospitalsData.csv")

# convert the csv file into a matrix
realData
real_data = matrix(nrow = 50, ncol = 3)
real_data
real_data[,1] = realData$Y
real_data[,2] = realData$X
# remove the commas from the numbers in the data frame
real_data[,3] = as.numeric(gsub(",", "", realData$Z))

# using the real data set
aa <- replicateSampling(nRep = 100, M = 10, N = 5, m = 6, 
                        n = 3, p = 0.15, mPrime = 8,
                        nPrime = 4, mu = c(20, 50, 40, 0, 0, 0), Sigma = Sigma,
                        Case = "A", Procedure = 1, seed_num = 531, m_error = TRUE,
                        aux_param_option = 4, data_matrix = real_data)

a_e <- replicateSampling(nRep = 100, M = 5, N = 10, m = 3, 
                        n = 5, p = 0.05, mPrime = 4,
                        nPrime = 7, mu = c(20, 50, 40, 0, 0, 0), Sigma = Sigma,
                        Case = "A", Procedure = 1, seed_num = 531, m_error = FALSE,
                        aux_param_option = 4, data_matrix = real_data)

aa$Rel_Performance
a_e$Rel_Performance
