
#' Etimation of species richness 
#' 
#' \code{ChaoSpecies} is a function to estimate species richness in one community.
#' @param data a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param k cut-off point. It is a value that separates frequency counts into abundant and rare groups
#' @param conf a confidence coefficient between 0 and 1.
#' @return a vector of species richness estimator and its confidence interval.
#' @author Anne Chao, K. H. Ma and T. C. Hsieh
#' @examples
#' data(SpecDemoAbu)
#' ChaoSpecies(SpecDemoAbu, datatype="abundance", k = 10, conf = 0.95)
#' data(SpecDemoInci)
#' ChaoSpecies(SpecDemoInci, datatype="incidence", k = 10, conf = 0.95)
#' 
#' @export


ChaoSpecies <- function(data, datatype = c("abundance", "incidence"), k = 10, conf = 0.95){
    method <- "all"
    if (k != round(k) || k < 0) 
      stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) 
      stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
    
    if (is.matrix(data) == T || is.data.frame(data) == T){
      if (ncol(data) != 1 & nrow(data) != 1)
        stop("Error: The data format is wrong.")
      if (ncol(data) == 1){
        data <- data[, 1]
      } else {
        data <- data[1, ]
      }
    }
    #  data <- as.numeric(round(data))
    
    if (datatype == "abundance"){
      f <- function(i, data){length(data[which(data == i)])}
      if (f(1, data) == sum(data)){
        stop("Error: The information of data is not enough.")}
      z <- (list(Basic.Data.Information = basicAbun(data, k)[[1]], Rare.Species.Group = RareSpeciesGroup(data, k), 
                 Species.Table = round(SpecAbunOut(data, method, k, conf), 3)))
    } else {
      dat <- data[-1]; 
      Q <- function(i, data){length(data[which(data == i)])}
      if (Q(1, dat) == sum(dat)){
        stop("Error: The information of data is not enough.")}
      z <- (list(Basic.Data.Information = basicInci(data, k)[[1]], Infreq.Species.Group = InfreqSpeciesGroup(data, k),
                 Species.Table = round(SpecInciOut(data, method, k, conf),3)))
    }
    class(z) <- c("species")
    z
  }


#' Estimation of the number of shared species in two communities
#' 
#' \code{ChaoShared} is a function to provide a shared species estimator in tow communities based on the following two sampling schemes. \cr
#' \enumerate{
#'   \item{In each community, a random sample of individuals is taken and species frequencies or abundances are recored;}
#'   \item{Each communtiy is sample several times or the whole area is divided into several quadrats and species presence/absence data for multiple samples/quadrats are recorded.
#'   }
#' }
#' @param data a numerical matrix or a data frame with two columns. Each column represent the species abundances or incidence frequencies of each community.
#' If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequencies in each column (See examples).
#' @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or presence/absence sample-base incidence data (\code{datatype = "incidence"}).
#' @param se calculate bootstrap standard error and show confidence interval; default is \code{TRUE}.
#' @param nboot the number of bootstrap resampling times, default is \code{200}.
#' @param conf a confidence coefficient between 0 and 1.
#' @return returns a list of two. First is the basic data informaiton; the other is a table of various estimators, their standard error and \code{100*conf} \% confidence interval.
#' @author Anne Chao, K. H. Ma and T. C. Hsieh
#' @examples
#' # load the individual-base (abundance) data
#' data(SharedSpecDemoAbu)
#' # Estimation of shared species
#' ChaoShared(SharedSpecDemoAbu, datatype="abundance", se=TRUE, nboot=200, conf=0.95)
#' # load the presence/absence sample-base (incidence) data
#' data(SharedSpecDemoInci)
#' # Estimation of shared species
#' ChaoShared(SharedSpecDemoInci, datatype="incidence", se=TRUE, nboot=200, conf=0.95)
#' 
#' @export
ChaoShared <-
  function(data, datatype = c("abundance", "incidence"), 
           se = TRUE, nboot = 200, conf = 0.95) {
    
    method <- "all"
    if (se == TRUE) {
      if (nboot < 1)
        nboot <- 1
      if (nboot == 1)
        cat("Warning: When \"nboot\" =" ,nboot, ", the bootstrap s.e. and confidence interval can't be calculated.", 
            "\n\n")  
    }
    
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
      cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
          "\n")
      cat("          We use \"conf\" = 0.95 to calculate!", 
          "\n\n")
      conf <- 0.95
    }
    
    datatype <- match.arg(datatype)
    if (datatype == "abundance") {
      x1 <- data[, 1]
      x2 <- data[, 2]
      Basic <- BasicFun(x1, x2, nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Ind(x1, x2, method, nboot, conf, se)
      colnames(output) <- c("Estimator", "Est_s.e.", paste(conf*100,"% Lower Bound",sep=""), paste(conf*100,"% Upper Bound",sep=""))
    }
    if (datatype == "incidence") {
      y1 <- data[, 1]
      y2 <- data[, 2]
      Basic <- BasicFun(y1, y2, B=nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Sam(y1, y2, method, conf, se)
      colnames(output) <- c("Estimator", "Est_s.e.", paste(conf*100,"% Lower Bound",sep=""), paste(conf*100,"% Upper Bound",sep=""))
    }
    out <- list(BASIC_DATA_INFORMATION=Basic, 
                ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES=output)
    class(out) <- c("ChaoShared")
    return(out)
  }


#' Etimation of species diversity 
#' 
#' \code{Diversity} This part features various diversity indices including the Shannon??s index and its effective number of species (diversity of order 1, or Shannon diversity), the Simpson??s index and its effective number of species (diversity order 2, or Simpson diversity), species richness (diversity of order 0).
#' @param X a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a list of species diversity estimator with order q and its confidence interval.
#' @author Anne Chao, K. H. Ma and T. C. Hsieh
#' @examples
#' \dontrun{
#' data(DivDemoAbu)
#' Diversity(DivDemoAbu,datatype="abundance")
#' }
#' @export

Diversity=function(X, datatype=c("abundance","incidence"))
{
  X=X[,1]
  
  BASIC.DATA <- matrix(round(c(sum(X), sum(X>0), 1-sum(X==1)/sum(X), CV.Ind(X)),3), ncol = 1)
  nickname <- matrix(c("n", "D", "C", "CV"), ncol = 1)
  BASIC.DATA <- cbind(nickname, BASIC.DATA)
  
  colnames(BASIC.DATA) <- c("Variable", "Value")
  rownames(BASIC.DATA) <- c("Number of observed individuals", "Number of observed species",
                            "Estimated sample coverage",
                            "Estimated CV")
  BASIC.DATA <- data.frame(BASIC.DATA)
  
  table0 <- matrix(0,4,4)
  table0[1,]=c(Chao1(X)[-5])
  table0[2,]=c(Chao1_bc(X))
  table0[3,]=round(c(SpecAbunAce(X)),1)
  table0[4,]=round(c(SpecAbunAce1(X)),1)
  colnames(table0) <- c("Estimator", "Est_s.e.", paste(Chao1(X)[5]*100,"% Lower Bound", sep=""), paste(Chao1(X)[5]*100,"% Upper Bound", sep=""))
  rownames(table0) <- c(" Chao1 (Chao, 1984)"," Chao1-bc "," ACE (Chao & Lee, 1992)",
                        " ACE-1 (Chao & Lee, 1992)")
  
  SHANNON=Shannon_index(X)
  table1=round(SHANNON[c(1:5),],3)
  colnames(table1) <- c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
  rownames(table1) <- c(" MLE"," MLE_bc"," Jackknife",
                        " Chao & Shen"," Chao (2013)")
  
  table1_exp=round(SHANNON[c(6:10),],3)
  colnames(table1_exp) <- c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
  rownames(table1_exp) <- c(" MLE"," MLE_bc"," Jackknife",
                            " Chao & Shen"," Chao (2013)")
  
  table2=round(Simpson_index(X)[c(1:2),],5)
  colnames(table2) <- c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
  rownames(table2) <- c(" MVUE"," MLE")
  
  table2_recip=round(Simpson_index(X)[c(3:4),],5)
  colnames(table2_recip) <- c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
  rownames(table2_recip) <- c(" MVUE"," MLE")
  
  Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", from=0, to=3, interval=0.25, B=50, conf=0.95))
  Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
  Hill<-round(Hill,3)
  Hill <- data.frame(Hill)
  colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
  
  z <- list("BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0, 
            "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
            "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
            "HILL.NUMBERS"= Hill)
  class(z) <- c("spadeDiv")
  return(z) 
}







#' Etimation of two-community similarity index 
#' 
#' \code{SimilarityPair} This part features various similarity indices for comparing data from two assemblages based on either multiple-sample incidence data or abundance data. The incidence-based indices include the classic Jaccard, Sorensen and Lennon et al. (2001) indices, and the abundance-based indices include the Bray-Curtis, Morisita-Horn, Horn and two newly developed abundance-based Jaccard and Sorensen indices.
#' @param X a matrix of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st row must be the number of sampling unit) for two communities.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param nboot the number of bootstrap resampling times, default is 200.
#' @return include two communities data summary and the various similarity indice estimation
#' @author Anne Chao, K. H. Ma and T. C. Hsieh
#' @examples
#' data(TwoComSimDemoAbu)
#' SimilarityPair(TwoComSimDemoAbu, datatype="abundance")
#' data(TwoComSimDemoInci)
#' SimilarityPair(TwoComSimDemoInci,datatype="incidence")
#' @export

SimilarityPair=function(X, datatype = c("abundance","incidence"),nboot=200)
{
  if(datatype=="abundance")
  {
    type="abundance"
    info1 <- c("S.total"=sum(rowSums(X)>0), "n1"=sum(X[,1]), "n2"=sum(X[,2]), 
               "D1"=sum(X[,1]>0), "D2"=sum(X[,2]>0), "D12"=sum(X[,1]>0 & X[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("f[11]"=sum(X[,1]==1 & X[,2]==1), 
               "f[1+]"=sum(X[,1]==1 & X[,2]>0), "f[+1]"=sum(X[,1]>0 & X[,2]==1),
               "f[2+]"=sum(X[,1]==2 & X[,2]>0), "f[+2]"=sum(X[,1]>0 & X[,2]==2))
    
    temp <- list()
    temp[[1]] <- Jaccard_Sorensen_Abundance_equ(datatype,X[,1],X[,2],nboot)
    temp[[1]] <- matrix(sapply(c(temp[[1]]), function(x) ifelse(x==0,"",x)),10)
    temp[[1]] <- rbind(temp[[1]][1:5,],
                       #c(round(C1n_equ(method="relative",X[,c(1,2)],nboot),4)[1:2],"","","",""),
                       NA,
                       c(round(C1n_equ(method="absolute",X[,c(1,2)],nboot),4)[1:2],"","","",""),
                       temp[[1]][6:10,])
    
    colnames(temp[[1]]) <- c("Estimate", "Bootstrap se.", "U_hat*", "U_hat* se.", "V_hat**", "V_hat** se.")
    rownames(temp[[1]]) <- c("Jaccard incidence", "Sorensen incidence", "Lennon et al. (2001)",
                             "Bray-Curtis", "Morisita-Horn", "Morisita Original",
                             "Horn (relative)", "Horn (absolute)","Jaccard Abundance (unadjusted)",
                             "Jaccard Abundance (adjusted)", "Sorensen Abundance(unadjusted)", "Sorensen Abundance(adjusted)")
    temp[[1]] <- as.data.frame(temp[[1]])
    
    temp[[2]] <- rbind(Chao1_equ(X[,1],conf=0.95)[-5],Chao1_bc_equ(X[,1],conf=0.95),
                       SpecAbunAce_equ(X[,1], k=10, conf=0.95), SpecAbunAce1_equ(X[,1] ,k=10, conf=0.95))
    colnames(temp[[2]]) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
    rownames(temp[[2]]) <- c("Chao1", "Chao1-bc", "ACE", "ACE-1")
    temp[[2]] <- as.data.frame(temp[[2]])
    
    temp[[3]] <- rbind(Chao1_equ(X[,2],conf=0.95)[-5],Chao1_bc_equ(X[,2],conf=0.95),
                       SpecAbunAce_equ(X[,2], k=10, conf=0.95), SpecAbunAce1_equ(X[,2] ,k=10, conf=0.95))
    colnames(temp[[3]]) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
    rownames(temp[[3]]) <- c("Chao1", "Chao1-bc", "ACE", "ACE-1")
    temp[[3]] <- as.data.frame(temp[[3]])
    z <- list("datatype"=type,"info1"=info1, "info2"=info2, "similarity"=temp[[1]], "assemblage1"=temp[[2]], "assemblage2"=temp[[3]])
  }      
  ##---------------------------------------------------------------
  if(datatype=="incidence")
  {
    no.assemblage=length(X[1,])
    Y=X[-1,]  
    "type"="incidence"
    info1 <- c("S.total"=sum(rowSums(Y)>0), "w"=X[1,1], "z"=sum(X[1,2]), 
               "D1"=sum(Y[,1]>0), "D2"=sum(Y[,2]>0), "D12"=sum(Y[,1]>0 & Y[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("Q[11]"=sum(Y[,1]==1 & Y[,2]==1), 
               "Q[1+]"=sum(Y[,1]==1 & Y[,2]>0), "Q[+1]"=sum(Y[,1]>0 & Y[,2]==1),
               "Q[2+]"=sum(Y[,1]==2 & Y[,2]>0), "Q[+2]"=sum(Y[,1]>0 & Y[,2]==2))
    
    temp <- list()
    temp[[1]] <- Jaccard_Sorensen_Abundance_equ(datatype,X[,1],X[,2],nboot)
    temp[[1]] <- matrix(sapply(c(temp[[1]]), function(x) ifelse(x==0,"",x)),10)
    colnames(temp[[1]]) <- c("Estimate", "Bootstrap se.", "U_hat*", "U_hat* se.", "V_hat**", "V_hat** se.")
    rownames(temp[[1]]) <- c("Jaccard incidence", "Sorensen incidence", "Lennon et al. (2001)",
                             "Bray-Curtis", "Morisita-Horn", "Morisita Original",
                             "Incidence-based Jaccard (unadjusted)", "Incidence-based Jaccard (adjusted)", 
                             "Incidence-based Sorensen(unadjusted)", "Incidence-based Sorensen(adjusted)")
    temp[[1]] <- temp[[1]]
    
    temp[[2]] <- rbind(SpecInciChao2(X[,1],conf=0.95)[-5],SpecInciChao2bc(X[,1],conf=0.95),
                       SpecInciModelh(X[,1], k=10, conf=0.95), SpecInciModelh1(X[,1] ,k=10, conf=0.95))
    colnames(temp[[2]]) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
    rownames(temp[[2]]) <- c("Chao2", "Chao2-bc", "ICE", "ICE-1")
    temp[[2]] <- as.data.frame(temp[[2]])
    
    temp[[3]] <- rbind(SpecInciChao2(X[,2],conf=0.95)[-5],SpecInciChao2bc(X[,2],conf=0.95),
                       SpecInciModelh(X[,2], k=10, conf=0.95), SpecInciModelh1(X[,2] ,k=10, conf=0.95))
    colnames(temp[[3]]) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
    rownames(temp[[3]]) <- c("Chao2", "Chao2-bc", "ICE", "ICE-1")
    temp[[3]] <- as.data.frame(temp[[3]])
    z <- list("datatype"=type, "info1"=info1, "info2"=info2, "similarity"=temp[[1]], "assemblage1"=temp[[2]], "assemblage2"=temp[[3]])
  }  
  class(z) <- c("spadeTwo")
  return(z)   
}


#' Etimation of multiple-community similarity measure 
#' 
#' \code{SimilarityMult} This part computes the generalized Sorensen, Horn similarity, and Morisita similarity/dissimilarity indices for comparing frequency or abundance data from more than two communities.
#' @param X a matrix of species sample frequencies (for abundance data), row is species number, and column is community number.
#' @param q For q=0, the SpadeR computes the estimates of Sorensen index for paired-wise communities. For q=1, the SpadeR computes the estimates of 1-Horn index for paired-wise communities. For q=2, the SpadeR computes the estimates of Morisita index for paired-wise communities.
#' @param nboot the number of bootstrap resampling times, default is 200.
#' @return include every community data summary, overlap measure in all communities, and pair-wise similarity measure depended on which order you choose in two communities. 
#' @author Anne Chao, K. H. Ma and T. C. Hsieh
#' @examples
#' data(MulComSimDemoAbu)
#' SimilarityMult(MulComSimDemoAbu,q=2,nboot=200)
#' @export


SimilarityMult=function(X,q=2,nboot=200)
{ 
  type <- "abundance"
  method<-"absolute"
  N <- no.community <- ncol(X)
  temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
  n <- apply(X,2,sum)
  D <- apply(X,2,function(x)sum(x>0))
  
  if(N > 2){
    temp1 <- temp2 <- rep(0, N*(N-1)/2)
    k <- 1
    for(i in 1:(N-1)){     
      for(j in (i+1):N){
        temp1[k] <- paste('D',i,j,sep="")
        temp2[k] <- sum(X[,i]>0 & X[,j]>0)
        k <- k + 1
      }
    }
  }
  names(temp2) <- temp1
  names(n) <- paste('n',1:N, sep="")
  names(D) <- paste('D',1:N, sep="")
  info <- c(temp, n, D, temp2)
  if(N == 3) info <- c(temp, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
  info <- c(info, nboot=nboot)
  
  
  Cqn=rbind(Cqn_se_equ(X,q=0,nboot),
            #C1n_equ(method="relative",X,nboot),
            NA,
            C1n_equ(method="absolute",X,nboot), 
            Cqn_se_equ(X,q=2,nboot)[1:4])
  if(N==3){Cqn <- rbind(Cqn, C33_se_equ(X,nboot)[1:4])}
  colnames(Cqn) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
  rownames(Cqn) <- c(paste("C0",N,sep=""),paste("C1",N,sep=""),paste("C1",N,"*",sep=""),paste("C2",N,sep=""),if(N==3) "C33")
  
  
  if(q==0 || q==1){Cqn_PC=matrix(0,choose(no.community,2),4)}
  if(q==2)        {Cqn_PC=matrix(0,choose(no.community,2),6)}
  k=1
  temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
  for(i in 1:(N-1)){  
    for(j in (i+1):N){
      Cqn_PC[k,] <- Cqn_se_equ(X[,c(i,j)],q,nboot,method)
      temp_PC[k] <- paste("C22(",i,",",j,")", sep="")
      temp_PD[k] <- paste("1-C22(",i,",",j,")", sep="")
      k <- k+1
    }
  }
  if(q==0 || q==1){
    colnames(Cqn_PC) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
    rownames(Cqn_PC) <- temp_PC
  }
  if(q==2){
    colnames(Cqn_PC) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL", "D.95%.LCL", "D.95%.UCL")
    rownames(Cqn_PC) <- temp_PC
  }
  
  C_SM=matrix(1,N,N)
  k <- 1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      C_SM[i,j] <- C_SM[j,i] <- Cqn_PC[k,1]
      k <- k+1
    }
  }
  
  z <- list("datatype"=type, "info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  class(z) <- c("spadeMult")
  z
}