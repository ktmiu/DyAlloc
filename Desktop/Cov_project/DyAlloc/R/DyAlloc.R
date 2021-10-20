#' Creates a marginal imbalance score for each treatment arm of N treatment groups.
#' This is a low level function to be used in subsequent functions and will most
#' likely not be used directly by a practitioner.
#' @export
#' @param df data frame. The profiles of existing subjects in the study.
#' @param N numeric variable. The number of treatment arms in the study.
#' @return A vector of marginal imbalance scores for each treatment arm is produced (a total of N scores)
MarImb <- function(df,N) {
  #INPUT: a data frame and number of treatment arms
  #OUTPUT: marginal imbalance score for each treatment arm (a total of N scores)

  imbs_arm <- vector()
  final_imbs <- vector()

  #for each treatment i, calculate the imbalance of each arm if the subject were to
  #assign to treatment j, then sum all possible j assignments
  trt <- vector()
  #checking to see if it's an empty/NA data frame
  if (is.na(df$treatment[1]) || is.null(df$treatment[1]) ){
    assigned <- 0
    trt <- c(0,0,0,0)
  } else{
    #number of assigned subjects already in the study
    assigned <- length(df$treatment)
    #number of subjects in each treatment arm
    for (t in 1:4){
      trt[t] <- sum(ifelse(df$treatment==t,1,0))
    }
  }
  for (i in 1:N){
    for (j in 1:N){
      imbs_arm[j] <- abs((trt[j]+ifelse(j==i,1,0))/(assigned+1)-(1/N))
    }
    final_imbs[i]<- sum(imbs_arm)
  }
  return(final_imbs)
}

#' Creates a weighted imbalance score for each treatment arm
#' *NOTE: the default for site here is FALSE, but one can specify site=TRUE in the argument
#' to account for site imbalance. Note that the weight for site has to be placed in the
#' third element of the 'weight' vector (and can be omitted when site=FALSE)*
#' @export
#' @param df data frame.The profiles of existing subjects in the study.
#' @param N numeric variable. The number of treatment arms in the study.
#' @param covars a vector. All the names of covariates in the study.
#' @param weight a vector. The weights for overall study, with-in stratum, (site),
#' and factors/covariates, in this particular order
#' @param newsub a vector. The factor profile of the new subject entering the study.
#' @param site a boolean value. Whether or not to account for site imbalance (default:FALSE)
#' @return A vector of weighted imbalance scores for each treatment arm is produced (a total of N scores)

WeiImb <- function(df,N,covars,weight,newsub,site=FALSE) {
  ##INPUT: df = data frame, covars = a vector of the name of covariates in the study,
  #weight = a vector of weights for overall study, with-in stratum, (site), and factors
  #newsub = a vector of the factor profile for the new subject
  #N = number of treatment arms, site = whether or not to account for site imbalance
  ##OUTPUT = weighted imbalance score for each treatment arm (a total of N scores)
  if (site==TRUE){
    stratum_profile <- paste(newsub[-1],collapse=".")
  } else{
    stratum_profile <- paste(newsub,collapse=".")
  }
  #overall data frame
  study <- df

  #stratum data frame
  dflist <- list(df[,1])
  if (is.na(dflist)){
    stratum_wimb <- c(rep(1.5,N))
  } else{
    for (i in 2:length(covars)){
      dflist <- c(dflist,list(df[,i]))
    }
    stratum_data <- split(df,dflist)
    stratum <- stratum_data[[which(names(stratum_data)==stratum_profile)]]

    stratum_wimb <- weight[2]*MarImb(stratum,N)
  }

  #factor data frame
  if (is.na(df[1,1])){
    factor_wimb <- c(rep(1.5,N))
  } else{
    factor_data <- list()
    factor <- list()
    for (j in 1:length(covars)){
      factor_data[[j]]<- split(df,df[,j])
    }
    if (site==TRUE){
      factor[[j]] <- factor_data[[j]][[which(names(factor_data[[j]])==newsub[j+1])]]
    } else{
      factor[[j]] <- factor_data[[j]][[which(names(factor_data[[j]])==newsub[j])]]
    }
    factor_wimb <- c(rep(0,N))
    for (i in 1:length(covars)){
      if (site==TRUE){
        factor_wimb <- factor_wimb + weight[i+3]*MarImb(factor[[i]],N)
      } else{
        factor_wimb <- factor_wimb + weight[i+2]*MarImb(factor[[i]],N)
      }
    }
  }

  ##accounting for site
  if (site==TRUE){
    site_data <- split(df,df$site)
    site <- site_data[[which(names(site_data)==newsub[1])]]
    site_wimb <- weight[3]*MarImb(site,N)
  } else{
    site_wimb <- c(rep(0,N))
  }
  wimb <- weight[1]*MarImb(study,N) + stratum_wimb  + site_wimb + factor_wimb
  return(wimb)
}

#' Creates a vector of probabilities of being assigned to each treatment group.
#' This function uses the imbalance scores created from the function WeiImb as input.
#' *NOTE: if the imbalance scores are the same, then they have equal probability of
#' being assigned to. E.g. if we get two equal least weighted imbalance scores, they each
#' have a probability of 0.5 for a subject to be assigned to.*
#' @export
#' @param imbalances a vector. Imbalance scores across the N treatment groups
#' (generally, the scores we get from WeiImb would be plugged in here)
#' @param alpha a numeric number. The "second best probability" parameter, generally alpha is set to 0.2.
#' This means that the treat group a subject assign to gets the probability of 0.8.
#' @return A vector of probabilities of being assigned to each treatment group.

Trt_Prob <- function(imbalances,alpha){
  #INPUT: the imbalance scores across treatment groups and the "second best probability" parameter
  #OUTPUT: a vector of probabilities of being assigned to each treatment group

  #put the least balance score(s) in FC
  FC <- imbalances[imbalances==min(imbalances)]
  #put the second least balance score(s) in FC
  if (length(imbalances[imbalances!=min(imbalances)])==0){
    SC <- vector()
  } else{
    SC <- imbalances[imbalances==min(imbalances[imbalances!=min(imbalances)])]
  }
  prob <- vector()
  ##assigning probabilities
  for (i in 1:length(imbalances)){
    if (i %in% which(imbalances==min(imbalances)) & length(FC)>=2){
      prob[i]=1/length(FC)
    } else if (i %in% which(imbalances==min(imbalances)) & length(FC)==1){
      prob[i]=1-alpha
    }else if (i %in% which(imbalances==min(imbalances[imbalances!=min(imbalances)]))
             & length(FC)==1){
      prob[i]=alpha/length(SC)
    } else{
      prob[i]=0
    }
  }
  return(prob)
}

#' Provides the final treatment group allocation for a new subject entering the study.
#'  Incorporates all the previous functions and a practitioner can use this function solely
#'  to obtain treatment assignments for the study subjects. The site default is FALSE.
#' @export
#' @param df data frame. The profiles of existing subjects in the study.
#' @param N numeric variable. The number of treatment arms
#' @param covars a vector. All the names of covariates in the study.
#' @param weight a vector. The weights for overall study, with-in stratum, (site),
#' and factors/covariates, in this particular order
#' @param newsub a vector. The factor profile of the new subject entering the study.
#' @param site a boolean value. Whether or not to account for site imbalance (default:FALSE)
#' @return A treatment group allocation

DyAlloc <- function(df, N, covars, weight, newsub, site=FALSE) {
  imb <- WeiImb(df, N, covars, weight, newsub, site)
  prob <- Trt_Prob(imb,0.2)
  final_trt <- sample(c("1","2","3","4"),1,FALSE,prob)
  return(final_trt)
}
