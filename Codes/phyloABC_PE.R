#Simulate the datasets on a random phylogeny with different parameter combinations of OU and competition models

library(TESS)
library(RPANDA)
library(ape)
library(igraph)
library(phytools)
library(ggplot2)

geomat<-function(geo){
  mat<-matrix(nrow=nrow(geo),ncol=nrow(geo))
  if(nrow(geo)<2){
    mat<-matrix(2,nrow=1,ncol=1)
  }else{
    for(i in 1:nrow(geo)){
      for(j in 1:nrow(geo)){
        mat[i,j]<-sum(geo[i,]*geo[j,])
      }
    }
    mat<-(mat/ncol(geo))
    diag(mat)<-1
  }
  return(as.matrix(mat))
}


#Simulation functions
simspat<-function(trait,geo,coln,alpha2){
  newspat<-matrix(0,nrow=nrow(geo),ncol=ncol(geo))
  for(i in 1:nrow(geo)){
    for(j in 1:ncol(geo)){
      prox<-(abs(trait[i]-trait))*geo[,j]
      if(all(prox==0)){
        prob<-coln[j]
      }else{
        prox<-min(prox[prox>0])
        prob<-coln[j]*((1+exp(-prox))^(-alpha2))
      }
      rand<-runif(1)
      if(rand<prob){
        newspat[i,j]<-1
      }else{newspat[i,j]<-0}
    }
  }
  return(newspat)
}






#simulation scheme
#Framework: Species disperse, evolve and compete on 5 separate "islands"

#model: z(i)(t+1)=z(i)(t)+ sigma*dW + psi*(theta-z(i)(t))+ 
#         alpha*(sum((geographic overlap(ij)(t)*(z(i)(t)-z(j)(t))*exp(-abs(z(i)-z(j))))))
# For geo_distribution shift probabilty 
#P( L(t+1)=1 | L(t)=x )=col*((1+(min(geographic overlap*exp(-abs(z(i)-z(j))))))^alpha2)

#parameters:sigma (BM variance),theta (OU optima),psi(OU coefficient),
#alpha(competition coefficient),col(colonization/extirpation rate),
#alpha2(competition coefficient for colonization/extirpation)
#Ancestral state is always assumed to be zero.

evolnull<-function(param,coln,br){
  
  sigma  <-param$sigma
  theta  <-param$theta
  psi    <-param$psi
  alpha  <-param$alpha
  
  ntip<-br[1,1]+1
  
  tr<-as.matrix(c(0,rnorm(1,0,sigma)))
  tr<-cbind(which(br[,3]==0),tr)
  
  #randomize the initial geographic distribution for MRCA
  geo<-matrix(rep(0,2*length(coln)),2,length(coln))
  geo[1,sample(1:length(coln),1)]<-1
  geo[2,sample(1:length(coln),1)]<-1
  geo<-cbind(which(br[,3]==0),geo)
  
  tips<-1:nrow(tr)
  
  for(i in 1:100){
    #If there is a branching event, add trait and geo characteristics  
    if(i%in%br[,3]){  
      
      parent<-tr[which(tr[,1]%in%which(br[,4]==i)),1]
      
      for(i1 in parent){
        tips<-c(tips,nrow(tr)+1,nrow(tr)+2)
        
        new.ind<-which(br[,1]==br[i1,2])
        p.tr<-tr[which(tr[,1]==i1),2]
        tr<-rbind(tr,cbind(new.ind,c(p.tr,p.tr+rnorm(1,0,sigma))))
        
        tips<-tips[!tips%in%which(tr[,1]==i1)]
        
        geo<-rbind(geo,c(new.ind[1],geo[which(geo[,1]==i1),-1]))
        
        newgeo<-rep(0,length(coln))
        if(sum(geo[nrow(geo),-1])==length(coln)){
          newgeo[sample((1:length(coln)),1)]<-1
        }else{
          newgeo[sample(which(geo[nrow(geo),-1]==0),1)]<-1
        }
        geo<-rbind(geo,c(new.ind[2],newgeo))
        
      }
    }
    newtr<-tr[tips,-1]
    geo1<-geo[tips,-1]
    
    newtr<-newtr + psi*(theta-newtr) + rnorm(length(newtr),0,sigma)
    for(i2 in 1:length(newtr)){
      newtr[i2]<-newtr[i2]+alpha*sum((geomat(geo1)[i2,])*(newtr[i2]-newtr)*exp(-abs(newtr[i2]-newtr)))
    }
    
    geo1<-simspat(newtr,geo1,coln,0)
    
    #ensure that there is no total extincion for any species
    if(sum(rowSums(geo1)==0)>0){
      geo1[which(rowSums(geo1)==0),sample((1:ncol(geo1)),1)]<-1
    }
    
    tr[tips,-1]<-newtr
    geo[tips,-1]<-geo1
  }
  
  
  #Outputs: traits and geo data
  ind<-which(br[,2]%in%(1:ntip))
  tr<-tr[which(tr[,1]%in%ind),]
  geo<-geo[which(tr[,1]%in%ind),]
  br1<-br[tr[,1],1:2]
  trait<-cbind(br1[,2],tr[,2])
  trait<-trait[order(trait[,1]),2]
  geo<-cbind(br1[,2],geo[,-1])
  geo<-geo[order(geo[,1]),-1]
  
  return(list(trait=trait,geo=geo))
}


evol<-function(param,coln,br){
  
  sigma  <-param$sigma
  theta  <-param$theta
  psi    <-param$psi
  alpha  <-param$alpha
  alpha2 <-param$alpha2
  
  ntip<-br[1,1]+1
  
  tr<-as.matrix(c(0,rnorm(1,0,sigma)))
  tr<-cbind(which(br[,3]==0),tr)
  
  #randomize the initial geographic distribution for MRCA
  geo<-matrix(rep(0,2*length(coln)),2,length(coln))
  geo[1,sample(1:length(coln),1)]<-1
  geo[2,sample(1:length(coln),1)]<-1
  geo<-cbind(which(br[,3]==0),geo)
  
  tips<-1:nrow(tr)
  
  for(i in 1:100){
    #If there is a branching event, add trait and geo characteristics  
    if(i%in%br[,3]){  
      
      parent<-tr[which(tr[,1]%in%which(br[,4]==i)),1]
      
      for(i1 in parent){
        tips<-c(tips,nrow(tr)+1,nrow(tr)+2)
        
        new.ind<-which(br[,1]==br[i1,2])
        p.tr<-tr[which(tr[,1]==i1),2]
        tr<-rbind(tr,cbind(new.ind,c(p.tr,p.tr+rnorm(1,0,sigma))))
        
        tips<-tips[!tips%in%which(tr[,1]==i1)]
        
        geo<-rbind(geo,c(new.ind[1],geo[which(geo[,1]==i1),-1]))
        
        newgeo<-rep(0,length(coln))
        if(sum(geo[nrow(geo),-1])==length(coln)){
          newgeo[sample((1:length(coln)),1)]<-1
        }else{
          newgeo[sample(which(geo[nrow(geo),-1]==0),1)]<-1
        }
        geo<-rbind(geo,c(new.ind[2],newgeo))
        
      }
    }
    newtr<-tr[tips,-1]
    geo1<-geo[tips,-1]
    
    newtr<-newtr + psi*(theta-newtr) + rnorm(length(newtr),0,sigma)
    for(i2 in 1:length(newtr)){
      newtr[i2]<-newtr[i2]+alpha*sum((geomat(geo1)[i2,])*(newtr[i2]-newtr)*exp(-abs(newtr[i2]-newtr)))
    }
    
    geo1<-simspat(newtr,geo1,coln,alpha2)
    
    #ensure that there is no total extincion for any species
    if(sum(rowSums(geo1)==0)>0){
      geo1[which(rowSums(geo1)==0),sample((1:ncol(geo1)),1)]<-1
    }
    
    tr[tips,-1]<-newtr
    geo[tips,-1]<-geo1
  }
  
  
  #Outputs: traits and geo data
  ind<-which(br[,2]%in%(1:ntip))
  tr<-tr[which(tr[,1]%in%ind),]
  geo<-geo[which(tr[,1]%in%ind),]
  br1<-br[tr[,1],1:2]
  trait<-cbind(br1[,2],tr[,2])
  trait<-trait[order(trait[,1]),2]
  geo<-cbind(br1[,2],geo[,-1])
  geo<-geo[order(geo[,1]),-1]
  
  return(list(trait=trait,geo=geo))
}


#########################################################################
#Simulate trees with different parameters


#Parameters to vary across trees:
sigmas<-c(0.1)
thetas<-0
psis<-c(0,0.05)
alphas<-c(0,0.05)
alpha2s<-c(0,1)
coln<-rep(0.5,5)

paramdat<-as.matrix(expand.grid(sigmas,thetas,psis,alphas,alpha2s))
colnames(paramdat)<-c("sigma","theta","psi","alpha","alpha2")


times<-100
tree <- read.tree("finches.phy")
ed<-tree$edge
bt1<-nodeHeights(tree)
bt1<-floor(bt1*times/max(bt1))
br<-as.matrix(cbind(ed,bt1))
while(sum(br[,3]==br[,4])!=0){
  ind<-which(br[,3]==br[,4])
  for(j in ind){
    br[j,4]<-br[j,4]+1
    br[which(br[,1]==br[j,2]),3]<-br[which(br[,1]==br[j,2]),3]+1
  }
}

trdat<-c()
geodat<-list()

for(i in 1:nrow(paramdat)){
  param<-as.list(paramdat[i,])
  res<-evol(param,coln,br)
  trdat<-rbind(trdat,res$trait)
  geodat<-c(geodat,list(res$geo))
}


trdat1<-t(apply(trdat,1,function(x) (x-min(x))/(max(x)+abs(min(x)))))

#########################################################################

#Set up an ABC-SMC alogorithm for parameter estimation
#Vary sigma,psi and alpha and alpha2.
#Inputs:
#Priors for the parameters:
#psi: U(0,1)
#alpha: U(0,1)
#alpha2: U(0,5)
theta<-0
coln<-rep(0.5,5)
sigma<-0.1



#Purturbation kernels for parameters: P(y|x)~N(x,0.1*x)

#Tolerences wrt summary statistics: 
#For euclidean distance between raw trait data: for n tips
tol.tr<-c(1.25,1)
#For comparing spatial data, convert it to overlap matrix and compute Frobenius
#norm of absolute difference between observed and predicted overlap matrices.

tol.geo<-c(3.15,3)

N<-50


#Start the cluster and split the loop on #cores
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

numcores<-detectCores()
clust<-makeCluster(3)
clusterExport(clust,c("geomat","simspat","evol","paramdat","times","tree","br","trdat","geodat","theta","coln","tol.geo","tol.tr","N"))
registerDoParallel(clust)


dat<-foreach(j1=1:nrow(trdat),.packages=c("RPANDA","TESS","ape","igraph","phytools"))%dopar%
{
  tr.o<-trdat1[j1,]
  
  geo.o<-geomat(geodat[[j1]])
  
  
  #Results to save:
  points1<-matrix(NA,1,5)
  

  j<-0 #sample particle
  jdash<-0
  #Start the filtering for iteration i. This keeps going until we collect N "particles"   
  
  while(j<2){
    
    #sample parameters from priors
    sigma  <-runif(1,0.01,0.2)
    psi    <-runif(1,0,0.2)
    alpha  <-runif(1,0,0.2)
    alpha2 <-runif(1,0,3)
    param<-list(sigma=sigma,theta=theta,psi=psi,alpha=alpha,alpha2=alpha2)
    
    
    #simulate the stochastic model 10 times with the sampled particle
 
    bt<-0
    
    for(k in 1:10){
      res<-evol(param,coln,br)
      tr.p<-res$trait
      tr.p<-(tr.p-min(tr.p))/(max(tr.p)+abs(min(tr.p)))
      geo.p<-geomat(res$geo)
      euc<-c(euc,sqrt(sum((tr.o-tr.p)^2)))
      f<-c(f,norm(abs(geo.o-geo.p),type="F"))
      f.dist<-norm(abs(geo.o-geo.p),type="F")
      
      if(euc<tol.tr[1] & f.dist<tol.geo[1]){
        bt<-bt+1
      }
    }
    
    if(bt>0){
      points1<-rbind(points1,c(sigma,psi,alpha,alpha2,bt))
      j<-j+1
    }
    jdash<-jdash+1
  }
  
  
  if(nrow(points1)>1){
    points1<-points1[-1,]
    points1<-points1[rep((1:nrow(points1)),points1[,5]),]
    
    for(j2 in 1:nrow(points1)){
      points1<-rbind(points1,cbind(rnorm(5,points1[j2,1],0.005),
                                       rnorm(5,points1[j2,2],0.005),
                                       rnorm(5,points1[j2,3],0.005),
                                       rnorm(5,points1[j2,4],0.05)))
    }
    
    points2<-matrix(NA,1,5)
    
    
    for(j3 in 1:nrow(points1)){
      sigma<-  points1[j3,1]
      psi<-    points1[j3,2]
      alpha<-  points1[j3,3]
      alpha2<- points1[j3,4]
      param<-list(sigma=sigma,theta=theta,psi=psi,alpha=alpha,alpha2=alpha2)
      
      #simulate the stochastic model 10 times with the sampled particle
      
      bt<-0
      
      for(k in 1:10){
        res<-evol(param,coln,br)
        tr.p<-res$trait
        tr.p<-(tr.p-min(tr.p))/(max(tr.p)+abs(min(tr.p)))
        geo.p<-geomat(res$geo)
        euc<-c(euc,sqrt(sum((tr.o-tr.p)^2)))
        f<-c(f,norm(abs(geo.o-geo.p),type="F"))
        f.dist<-norm(abs(geo.o-geo.p),type="F")
        
        if(euc<tol.tr[2] & f.dist<tol.geo[2]){
          bt<-bt+1
        }
      }
      
      if(bt>0){
        points2<-rbind(points2,c(sigma,psi,alpha,alpha2,bt))
      }
      
    }
  
  }
  
  return(list(first=points1,second=points2))
}  



stopCluster(clust)
saveRDS(dat,"simABC_simp.Rds")






################################################################################
#Analyze the results
################################################################################
sigmas<-c(0.1)
thetas<-0
psis<-c(0,0.05,0.1)
alphas<-c(0,0.05,0.1)
alpha2s<-c(0,1,2)
coln<-rep(0.5,5)

paramdat<-as.matrix(expand.grid(sigmas,thetas,psis,alphas,alpha2s))
colnames(paramdat)<-c("sigma","theta","psi","alpha","alpha2")


params<-paramdat[rep(1:27,each=2500),3:5]
colnames(params)<-c("psi.m","alpha.m","alpha2.m")

dat<-readRDS("simABC1.Rds")
dat2<-readRDS("simABC2.Rds")
datf<-readRDS("simABC1f.Rds")
dat2f<-readRDS("simABC2f.Rds")

eucdat<-c()
eucdat.f<-c()




for(i in 1:length(dat)){
  samp1<-dat[[i]][[1]]
  samp2<-dat[[i]][[2]]
  
  samp3<-dat2[[i]][[1]]
  samp4<-dat2[[i]][[2]]
  

  eucdat<-rbind(eucdat,cbind(samp1,samp2[,3:12]))
  eucdat.f<-rbind(eucdat.f,cbind(samp3,samp4[,4:13]))
  
}

eucdat<-as.data.frame(eucdat)
eucdat.f<-as.data.frame(eucdat.f)

f1<-function(x){
  sum(x[1:10]<1.1 & x[11:20]<3)
}

filt<-apply(eucdat[,-c(1:2)],1,f1)


eucdat1<-cbind(params,eucdat[,1:2],filt)
eucdat1<-eucdat1[eucdat1$filt>0,]
eucdat1<-eucdat1[rep((1:nrow(eucdat1)),eucdat1$filt),-6]

eucdat1$psi.m<-as.factor(eucdat1$psi.m)
eucdat1$alpha.m<-as.factor(eucdat1$alpha.m)
eucdat1$alpha2.m<-as.factor(eucdat1$alpha2.m)

colnames(eucdat1)[4:5]<-c("psi","alpha")

ggplot(eucdat1,aes(x=psi.m,y=psi))+geom_boxplot() + 
  labs(x="True \u03C8",y="Estimated \u03C8")
ggplot(eucdat1,aes(x=alpha.m,y=alpha))+geom_boxplot()+
  labs(x=bquote("True \u03B1"),y="Estimated \u03B1")

f1<-function(x){
  sum(x[1:10]<0.9 & x[11:20]<2.85)
}

filt.f<-apply(eucdat.f[,-c(1:3)],1,f1)

eucdat1f<-cbind(params,eucdat.f[,1:3],filt.f)
eucdat1f<-eucdat1f[eucdat1f$filt.f>0,]
eucdat1f<-eucdat1f[rep((1:nrow(eucdat1f)),eucdat1f$filt),-7]

eucdat1f$psi.m<-as.factor(eucdat1f$psi.m)
eucdat1f$alpha.m<-as.factor(eucdat1f$alpha.m)
eucdat1f$alpha2.m<-as.factor(eucdat1f$alpha2.m)

colnames(eucdat1f)[4:6]<-c("psi","alpha","alpha2")

getmode <- function(v) {
  his<-hist(v,seq(0,max(v)+1,length.out=10),plot=F)
  return(his$mids[which.max(his$counts)])
}

ggplot(eucdat1f,aes(x=psi.m,y=psi))+geom_boxplot() + 
  labs(x="True \u03C8",y="Estimated \u03C8")
  #stat_summary(fun=getmode,geom="point", size=2, color="red")
ggplot(eucdat1f,aes(x=alpha.m,y=alpha))+geom_boxplot()+
  labs(x=bquote("True \u03B1"),y="Estimated \u03B1")
  #stat_summary(fun=getmode,geom="point", size=2, color="red")
    
ggplot(eucdat1f,aes(x=alpha2.m,y=alpha2))+geom_boxplot()+
    labs(x=expression(paste("True ",alpha[c])),y=expression(paste("Estimated ",alpha[c])))
  #+  stat_summary(fun=getmode,geom="point", size=2, color="red")


#Bootstrap to get likihood ratios under different simulated datasets
f1<-function(x){
  sum(x[1:10]<1 & x[11:20]<4.5)
}

filt<-apply(eucdat[,-c(1:2)],1,f1)

filt.f<-apply(eucdat.f[,-c(1:3)],1,f1)
a<-b<-c()
for(i in 1:nrow(paramdat)){
  
  inds<-(2500*(i-1))+1
  samp1<-filt[inds:(inds+2499)]
  samp2<-filt.f[inds:(inds+2499)]
 
   a<-c(a,sum(samp1>0))
   b<-c(b,sum(samp2>0))
}

a
b

acc.dat<-c()
for(i in 1:nrow(paramdat)){
  
  inds<-(2500*(i-1))+1
  samp1<-filt[inds:(inds+2499)]
  samp2<-filt.f[inds:(inds+2499)]
  
  
  acc<-c()
  
  for(i1 in 1:1000){
    bootsamp<-sample(1:2500,2500,replace=TRUE)
    
    samp11<-(sum(samp1[bootsamp]>0))
    samp21<-(sum(samp2[bootsamp]>0))
    
  
    if(samp11>0){acc<-c(acc,samp21/samp11)
    }else{acc<-c(acc,2500)}
    
  }
  acc.dat<-rbind(acc.dat,acc)

}



#quantiles of null expectations:
q95<-apply(acc.dat[1:9,],1,quantile,probs=0.95)

acc.dat1<-cbind(acc.dat,rep(q95,3))

powfun<-function(x){
  return(sum(x[1:100]>x[101]))
}

power1<-apply(acc.dat1,1,powfun)


#Test the expected geographic overlap with randomized matrices

randset<-list()
for(i in 1:100){
  mat<-matrix(0,15,5)
  for(j in 1:15){
    x<-sample(c(0,1),5,replace=T)
    if(sum(x)>1){mat[j,]<-x
    }else{mat[sample(1:5,1),]<-1}
  }
  randset<-c(randset,list(mat))
}

fmat<-c()
for(i1 in 1:length(geodat)){
  fd<-c()
  geod<-geomat(geodat[[i1]])
  for(i2 in 1:100){
    geoc<-geomat(randset[[i2]])
    f<-norm(abs(geod-geoc),type="F")
   fd<-c(fd,f) 
  }
  fmat<-rbind(fmat,fd)
}

plot(density(fmat[1,]))
for(i in 2:27){
  lines(density(fmat[i,]),col=i)
}
