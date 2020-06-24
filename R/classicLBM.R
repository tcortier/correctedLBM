#----- Fichier 'R/classicLBM.R'

#' Main script for a LBM
#'
#' @param probability_distribution string
#' can be 'bernoulli' or 'poisson'
#' @param adj Matrix
#' @param initMethod string
#' can be 'hierarchical_clust' 'kmeans_clust' or 'spectral_clust'
#' @param VEMparam list
#' exemple = list(maxIter=50,fixPointIter=5,threshold=1e-3,trace=TRUE,cores=1)
#' @param MAINparam list
#' exemple = list(maxExplo=1.5,maxGroups=15,reinitialize=FALSE)
#'
#'
#' @return List of models for different number of groups
#' @export
#'
#' @examples
#' npc <- c(50,40) # nodes per class
#' Q <- c(2,3) # classes
#' n <- npc * Q # nodes
#' Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
#' Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
#' P<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2]) 
#' M<-1*(matrix(runif(n[1]*n[2]),n[1],n[2])<Z1%*%P%*%t(Z2)) ## adjacency matrix
#' classicLBM('bernoulli', P)


classicLBM<-function(probability_distribution='bernoulli',adj,initMethod="kmeans_clust",VEMparam=list(maxIter=50,fixPointIter=5,threshold=1e-3,trace=TRUE,cores=1),MAINparam=list(maxExplo=1.5,maxGroups=15,reinitialize=FALSE,cores=1)){
  k<-c(1,1)
  models<-list()
  max=-1e20
  whmax=c(0,0)
  gr=1
  cond=TRUE
  while(cond){
    if (gr==1){
      name=unlist(lapply(1:k[2],function(k2){paste(as.character(k[1]),as.character(k2),sep="&")}))
      mods=mclapply(1:k[2],function(k2){
        print(paste('k={',k[1],',',k2,'}'))
        cl0=clustinit_LBM(adj,k[1],k2,initMethod)
        tau1=clustering_indicator(cl0[[1]])
        tau2=clustering_indicator(cl0[[2]])
        model<-LBM_VEM(probability_distribution,adj,k[1],k2,tau1,tau2,VEMparam)
        return(model)
      },mc.cores = MAINparam$cores)
      k2max=which.max(sapply(mods, function(mod) mod$ICL))
      print(k2max)
      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k[1]),as.character(k2max),sep="&")]]$ICL>max){
        whmax=c(k[1],k2max)
        max=models[[paste(as.character(k[1]),as.character(k2max),sep="&")]]$ICL
      }
      LBM_plot(models)
    }
    else if (gr==2){
      name=unlist(lapply(1:k[1],function(k1){paste(as.character(k1),as.character(k[2]),sep="&")}))
      mods=mclapply(1:k[1],function(k1){
        print(paste('k={',k1,',',k[2],'}'))
        cl0=clustinit_LBM(adj,k1,k[2],initMethod)
        tau1=clustering_indicator(cl0[[1]])
        tau2=clustering_indicator(cl0[[2]])
        model<-LBM_VEM(probability_distribution,adj,k1,k[2],tau1,tau2,VEMparam)
        return(model)
      },mc.cores = MAINparam$cores)
      k1max=which.max(sapply(mods, function(mod) mod$ICL))
      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k1max),as.character(k[2]),sep="&")]]$ICL>max){
        whmax=c(k1max,k[2])
        max=models[[paste(as.character(k1max),as.character(k[2]),sep="&")]]$ICL
      }
      LBM_plot(models)
    }
    cond=((k[1]<4|k[1]<round((maxExplo*whmax[1])+0.1)|k[2]<4|k[2]<round((maxExplo*whmax[2])+0.1))&(k[1]<maxGroups)&(k[2]<maxGroups))
    if ((k[1]<max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]<max(4,round((maxExplo*whmax[2])+0.1)))){
      gr=which.min(k)
      k[which.min(k)]<-k[which.min(k)]+1
    }
    else if((k[1]>=max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]<max(4,round((maxExplo*whmax[2])+0.1)))){
      k[2]<-k[2]+1
      gr=2
    }
    else if((k[1]<max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]>=max(4,round((maxExplo*whmax[2])+0.1)))){
      k[1]<-k[1]+1
      gr=1
    }
  }
  if (reinitialize==TRUE){
    max2=max
    it<-0
    cond<-TRUE
    print(k)
    while(cond&(it<3)){
      it<-it+1
      for (k1 in 1:(k[1]-1)){
        for (k2 in 1:(k[2]-1)){
          print(paste('forward','k={',k1,',',k2,'}'))
          model<-forward_explo(models,k1,k2,probability_distribution,adj,VEMparam)
          if (model[[1]]$ICL>models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL){
            models[[paste(as.character(k1+1),as.character(k2),sep="&")]]=model[[1]]
            if (models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL
            }
          }

          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2+1),sep="&")]]$ICL){
            models[[paste(as.character(k1),as.character(k2+1),sep="&")]]=model[[2]]
            if (models[[paste(as.character(k1),as.character(k2+1),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2+1),sep="&")]]$ICL
            }
          }
          LBM_plot(models)
        }
      }
      for (k1 in c(k[1]:3)){
        for (k2 in c(k[2]:3)){
          print(paste('backward','k={',k1,',',k2,'}'))
          model<-backward_explo(models,k1,k2,probability_distribution,adj,VEMparam)
          if (model[[1]]$ICL>models[[paste(as.character(k1-1),as.character(k2),sep="&")]]$ICL){
            models[[paste(as.character(k1-1),as.character(k2),sep="&")]]=model[[1]]
            if (models[[paste(as.character(k1-1),as.character(k2),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1-1),as.character(k2),sep="&")]]$ICL
            }
          }
          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2-1),sep="&")]]$ICL){
            models[[paste(as.character(k1),as.character(k2-1),sep="&")]]=model[[2]]
            if (models[[paste(as.character(k1),as.character(k2-1),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2-1),sep="&")]]$ICL
            }
          }
          LBM_plot(models)
        }
      }
      if (max2>max){
        max=max2
      }
      else{
        cond=FALSE
      }
    }
  }
  return(models)
}
