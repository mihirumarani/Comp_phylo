library(ape)
library(phytools)
library(geiger)
library(ggplot2)
library(cowplot)
library(RPANDA)


#A simple analysis where a trait evolution model from Nuismer and Harmon (2015) is simulated on a pure-birth phylogenetic tree with 20 tips. All the species at a given time
#are assumed to be in sympatry. Trait values do not influence speciation rates and the speciation rates do not alter the trait values of new species by a lot. 
#Trait optimum for OU process is assumed to be constant. Free parameters in the model: a. Ancestral trait value, b. OU trait optimum, c. Coefficient for OU process, d. drift variance
# e. coefficient of competition




#total time to evolve from the common ancestor (To limit the interations in the simulation)
time<-250
brt<-branching.times(tree)

#speciation events
br<-floor(time*(1-(brt/max(brt))))
br<-sort(br,decreasing=FALSE)


#Create a simulation for trait evolution
#Inputs:S, a coefficient that represents strength of selection due to interactions
#       psi, a coefficient that represents strength of stabilizing abiotic selection
#       theta, an environmental optima  
#       delta, brownian motion change (variance)
#       trait values of species
#       ancestral trait value


evol<-function(ent,br){
  ancestor<-ent[1]
  S<-ent[2]
  psi<-ent[3]
  delta<-ent[5]
  theta<-ent[4]
  tr<-c(ancestor)
  for(i in 0:time){
    if(i%in%br){
      l<-sum(i==br)
      tr<-c(tr,tr[length(tr)]+rnorm(l,0,0.01))}
    tr<-tr-S*(mean(tr)-tr)+psi*(theta-tr)+rnorm(length(tr),0,delta)
  }
  return(tr)
}

#Design a ABC model to estimate the model paramter fitting the observed data

#Input phylogenetic tree and trait data for extant species
tree<-pbtree(b=1,d=0,n=20)
Data<-rnorm(20,0,3)


#Draw values for the parameters to know from their prior distributions
#Set up a rejection threshold for euclidean distance
#Model: tr.change=tr+S*(mean(tr)-tr)+psi*(theta-tr)+rnorm(length(tr),0,delta)
#Parameters and their priors
# Ancestral trait  ancestor: dunif(1,8)
# Coefficient for interaction term  S:dunif(0,5)
# Coefficient for stabilizing selection  psi: dlnorm(0,4)
# Selection optima 1  theta: dunif(1,8)
# drift sd per time  delta: dunif(0,5)

final<-vector(length=6)
pb<-txtProgressBar(min =1, max =10000, style = 3)
for(i in 1:10000){
  Sys.sleep(0.05)
  # update progress bar
  setTxtProgressBar(pb, i)
  pars<-c(c(runif(1,1,8),runif(1,0,5),runif(1,0,4),rnorm(1,0,8),runif(1,0,10)))
  res<-evol(pars,br)
  
  #Calculate summary statistics
  stat1<-sqrt(sum((Data-res)^2))  # Simple euclidean distance between trait vectors
  #stat.sp1<-spectR_t(tree,Data)$splitter-spectR_t(tree,res)$splitter # characterstic (splitter)of spectral density plot of a tree and trait data
  #stat.sp2<-spectR_t(tree,Data)$tracer-spectR_t(tree,res)$tracer   #characterstic (tracer)of spectral density plot of a tree and trait data
  #stat.sp3<-spectR_t(tree,Data)$fragmenter-spectR_t(tree,res)$fragmenter #characterstic (fragmenter)of spectral density plot of a tree and trait data
  if(stat1<25) {final<-rbind(final,c(pars,stat1))}
  #final<-rbind(final,c(pars,stat1,stat.sp1,statsp2,statsp3))
}

final1<-final[-1,]
final1<-final1[,1:5]

#Second filtering
#create sets of parameter values based on first filtering
#Deviate param values using rnorm function
final2<-vector(length=5)
for(i in 1:length(final1)){
  for(j in 1:5){
   final2<-rbind(final2,rnorm(5,final1[i,],0.05*final1[i,])) 
  }
}

#Sequential filtering
final.second<-final2


final.third<-vector(length=6)
pb<-txtProgressBar(min = 1, nrow(final.second), style = 3)
for(i in 1:nrow(final.second)){
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i)
  pars<-c(final.second[i,])
  res<-evol(pars,br)
  stat<-sqrt(sum((Data-res)^2))
  if(stat<10){final.third<-rbind(final.third,c(pars,stat))}
}


#Plot overlap of priors and posteriors

par(mfrow=c(3,2))

prior.A<-data.frame(points=runif(1000,1,8))
post.A<-data.frame(points=final.third[,1])  
prior.A$type<-"prior for Ancestral trait"
post.A$type<-"Posterior for Ancestral trait"
A.dat<-rbind(prior.A,post.A)
ggplot(A.dat,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,10)

priorS<-data.frame(points=runif(100000,0,5))
postS<-data.frame(points=final.third[,2])  
priorS$type<-"prior for S"
postS$type<-"Posterior for S"
Sdat<-rbind(priorS,postS)
tiff('AnolS.tiff', units="in", width=5, height=5, res=300)

ggplot(Sdat,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,20)
dev.off()

prior.psi<-data.frame(points=runif(1000,0,4))
post.psi<-data.frame(points=final.third[,3])
prior.psi$type<-"Prior for delta"
post.psi$type<-"Posterior for delta"
psidat<-rbind(prior.psi,post.psi)
tiff('Anoldelta.tiff', units="in", width=5, height=5, res=300)

ggplot(psidat,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,20)
dev.off()

prior.delta<-data.frame(points=runif(1000,0,10))
post.delta<-data.frame(points=final.third[,5])  
prior.delta$type<-"Prior for sigma"
post.delta$type<-"Posterior for sigma"
delta.dat<-rbind(prior.delta,post.delta)
tiff('Anolsigma.tiff', units="in", width=5, height=5, res=300)

ggplot(delta.dat,aes(points,fill=type))+geom_density(alpha=0.2)
dev.off()
prior.theta<-data.frame(points=runif(1000,0,8))
post.theta<-data.frame(points=final.third[,4])  
prior.theta$type<-"prior for theta"
post.theta$type<-"Posterior for theta"
theta.dat<-rbind(prior.theta,post.theta)
ggplot(theta.dat,aes(points,fill=type))+geom_density(alpha=0.2)+ xlim(1,8)

ggplot(postS,aes(points,fill=type))+geom_density()+xlim(0,2)
ggplot(post.psi,aes(points,fill=type))+geom_density()+xlim(0,2)
compare<-rbind(postS,post.psi)
ggplot(compare,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,2.5)

#Summary of posterior values
means<-apply(final.third,2,mean)

mostlike<-function(x){
  a<-density(x)
  return(a$x[which.max(a$y)])
}
modes<-apply(final.third,2,mostlike)
sds<-apply(final.third,2,sd)

Summary<-cbind(c("A","S","Psi","Theta","Delta"),means,modes,sds)
Summary<-data.frame(Parameters=c("A","S","Psi","Theta","Delta"),Mean=means,Mode=modes,SD=sds)


#Model fit

simdat<-evol(means,br)

simdat1<-sort(simdat,decreasing = F)
Data1<-sort(Data,decreasing = F)
result<-data.frame(Data=Data1,Simulation=simdat1)
ggplot(result,aes(x=Data,y=Simulation))+geom_point(size=2)+geom_smooth(method="lm")



final.third<-final.third[-1,]

final1<-final.third

final1<-final[-1,]
H1_post<-order(final1[,6])[1:500]
H1Post<-matrix(ncol=3,nrow=500)
H1Post[,1]<-final[,2][H1_post]  #S
H1Post[,2]<-final[,3][H1_post]  #Stabilizing selection coeff
H1Post[,3]<-final[,5][H1_post]  #Drift variance
MLE<-H1Post[1,]

#Draw values for the parameters to know from their prior distributions
#Set up a rejection threshold for euclidean distance
#Model: tr.change=tr+S*(mean(tr)-tr)+psi*(theta-tr)+rnorm(length(tr),0,delta)
#Parameters and their priors
# Ancestral trait  ancestor: dunif(1,8)
# Coefficient for interaction term  S:dunif(0,5)
# Coefficient for stabilizing selection  psi: dlnorm(0,4)
# Selection optima 1  theta: dunif(1,8)
# drift sd per time  delta: dunif(0,5)


k   = kde(H1Post, xmin=c(0,min(H1Post[,2]),0), xmax=c(max(H1Post[,1]),max(H1Post[,2]),max(H1Post[,3])))
k0  = kde(H1Post, xmin=c(0,min(H1Post[,2]),0), xmax=c(0,max(H1Post[,2]),max(H1Post[,3])))       


# Use kernel smoothing to estimate likelihood maxima with and without competition.
k_max_index     = which(k$estimate == max(k$estimate), arr.ind = TRUE)
H1_lik          = k$estimate[k_max_index[1], k_max_index[2],k_max_index[3]]

k0_max_index    = which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
H0_lik          = k0$estimate[k0_max_index[1], k0_max_index[2],k0_max_index[3]]

LRT             = -2 * log( H0_lik / H1_lik )


#bootstrap for estimated parameter values


#Simulate datasets
Dataset<-matrix(ncol=2,nrow=100)
MLE<-matrix(ncol=3,nrow=500)
pb<-txtProgressBar(min = 1, 250, style = 3)
for(m in 1:250){
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, m)
  #Simulate data
 
  
 
  dat<-evol_non(c(0.05,0.05,a),br)
  #Run a likelihood test
  final<-vector(length=6)
  for(n in 1:5000){
    pars<-c(0,runif(1,0,5),rnorm(1,0,0.5),0,runif(1,0,1))
    res<-evol(pars,br)
    stat<-sqrt(sum((dat-res)^2))
    final<-rbind(final,c(pars,stat))
  }
  final<-final[-1,]
  H1_post<-order(final[,6])[1:500]
  H1Post<-matrix(ncol=3,nrow=500)
  H1Post[,1]<-final[,2][H1_post]  #S
  H1Post[,2]<-final[,3][H1_post]  #Stabilizing selection coeff
  H1Post[,3]<-final[,5][H1_post]  #Drift variance
  MLE[m,]<-H1Post[1,]
  
  k   = kde(H1Post, xmin=c(0,min(H1Post[,2]),0), xmax=c(max(H1Post[,1]),max(H1Post[,2]),max(H1Post[,3])))
  k0  = kde(H1Post, xmin=c(0,min(H1Post[,2]),0), xmax=c(0,max(H1Post[,2]),max(H1Post[,3])))       
  
  # Use kernel smoothing to estimate likelihood maxima with and without competition.
  k_max_index     = which(k$estimate == max(k$estimate), arr.ind = TRUE)
  H1_lik          = k$estimate[k_max_index[1], k_max_index[2],k_max_index[3]]
  
  k0_max_index    = which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
  H0_lik          = k0$estimate[k0_max_index[1], k0_max_index[2],k0_max_index[3]]
  
  LRT             = -2 * log( H0_lik / H1_lik )
  
  Dataset20[m,]<-c(a,dat,LRT)
}






