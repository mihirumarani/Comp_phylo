library(ape)
library(phytools)
library(geiger)
library(ggplot2)
library(cowplot)


#Input phylogenetic tree from Mahler et.al.(2013)
tree<-read.tree("GA_Anolis_MCC.tre")

#Input trait data
Data<-read.csv("GA_Anolis_traits.csv")

#Trait to analyze is SVL, snout to vent length
Data1<-Data[,2]

#number of tips
n<-length(tree$tip.label)

##Calculate speciation events


#total time to evolve from the common ancestor (in terms of generations~5 years)

t<-10^3
brt<-branching.times(tree)

#speciation events
br<-floor(t*(1-(brt/max(brt))))
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
  for(i in 0:250){
    
    if(i%in%br){
      l<-sum(i==br)
      tr<-c(tr,tr[length(tr)]+rnorm(l,0,0.01))}
    tr<-tr+S*(mean(tr)-tr)+psi*(theta-tr)+rnorm(length(tr),0,delta)
    tr[which(1>tr)]<-1
    tr[which(tr>8)]<-8
  }
  return(tr)
}

#Run a bootstrap with simulated data under BM+OU model to calculate the null expectation of likelihood ratios of two models
# Null: BM+OU Alternative: BM+OU+Interactions
dataset1<-fastBM(tree,a=mean(Data1),mu=mean(Data1),sig2=runif(1,0,5),bounds=c(1,8),nsim=500)
A<-mean(Data1)
Like<-c()
pb<-txtProgressBar(min =1, max =500, style = 3)


  
  Sys.sleep(0.01)
  # update progress bar
  setTxtProgressBar(pb, k)


# ABC model

#Draw values for the parameters to know from their prior distributions
#Set up a rejection threshold for euclidean distance
#Model: tr.change=tr+S*(mean(tr)-tr)+psi*(theta-tr)+rnorm(length(tr),0,delta)
#Parameters and their priors
# Ancestral trait  ancestor: fixed
# Coefficient for interaction term  S:dunif(0,5)
# Coefficient for stabilizing selection  psi: dlnorm(0,4)
# Selection optima 1  theta: fixed
# drift sd per time  delta: dunif(0,5)


Data<-dat1[2,]

final<-vector(length=6)
for(i in 1:10000){
  pars<-c(c(0,runif(1,0,5),rnorm(1,0,0.5),runif(1,0,1)))
  res<-evol(pars,br)
  stat<-sqrt(sum((Data-res)^2))
  final<-rbind(final,c(pars,stat))
}

final<-final[-1,]
H1_post<-order(final[,6])[1:500]
H1Post<-matrix(ncol=2,nrow=500)
H1Post[,1]<-final[,2][H1_post]  #S
H1Post[,2]<-final[,3][H1_post]  #Stabilizing selection coeff
H1Post[,3]<-final[,5][H1_post]  #Drift variance


k   = kde(H1Post, xmin=c(0,min(H1Post[,2]),0), xmax=c(max(H1Post[,1]),max(H1Post[,2]),max(H1Post[,3])))
k0  = kde(H1Post, xmin=c(0,min(H1Post[,2]),0), xmax=c(0,max(H1Post[,2]),max(H1Post[,3])))       

# Use kernel smoothing to estimate likelihood maxima with and without competition.
k_max_index     = which(k$estimate == max(k$estimate), arr.ind = TRUE)
H1_lik          = k$estimate[k_max_index[1], k_max_index[2],k_max_index[3]]

k0_max_index    = which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
H0_lik          = k0$estimate[k0_max_index[1], k0_max_index[2],k0_max_index[3]]

LRT             = -2 * log( H0_lik / H1_lik )

Like<-c(Like,LRT)









#Plot overlap of priors and posteriors

par(mfrow=c(3,2))

prior.A<-data.frame(points=runif(1000,1,8))
post.A<-data.frame(points=final.third[,1])  
prior.A$type<-"prior for Ancestral trait"
post.A$type<-"Posterior for Ancestral trait"
A.dat<-rbind(prior.A,post.A)
ggplot(A.dat,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,10)

priorS<-data.frame(points=rlnorm(1000,0,2))
postS<-data.frame(points=final.third[,2])  
priorS$type<-"prior for S"
postS$type<-"Posterior for S"
Sdat<-rbind(priorS,postS)
ggplot(Sdat,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,20)

prior.psi<-data.frame(points=rlnorm(1000,0,2))
post.psi<-data.frame(points=final.third[,3])
prior.psi$type<-"prior for psi"
post.psi$type<-"Posterior for psi"
psidat<-rbind(prior.psi,post.psi)
ggplot(psidat,aes(points,fill=type))+geom_density(alpha=0.2)+xlim(0,20)

prior.delta<-data.frame(points=runif(1000,0,10))
post.delta<-data.frame(points=final.third[,5])  
prior.delta$type<-"prior for delta"
post.delta$type<-"Posterior for delta"
delta.dat<-rbind(prior.delta,post.delta)
ggplot(delta.dat,aes(points,fill=type))+geom_density(alpha=0.2)

prior.theta<-data.frame(points=runif(1000,0,8))
post.theta<-data.frame(points=final.third[,4])  
prior.theta$type<-"prior for theta"
post.theta$type<-"Posterior for theta"
theta.dat<-rbind(prior.theta,post.theta)
ggplot(theta.dat,aes(points,fill=type))+geom_density(alpha=0.2)

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
ggplot(result,aes(x=Data,y=Simulation))+geom_point(size=2)+geom_abline(slope=1,intercept = 0,color="red")+geom_smooth(method="lm")




###Alternative model of trait evolution
#For the interaction component, trait change is linearly proportional to  trait matching function with thresholds
#Inputs: z= vector of all trait values, S=coefficient of rate change
traitevol<-function(z,omega,S){
  zdash<-c()
  for(i in 1:length(z)){
    x<-z-z[i]
    x[which(abs(x)>=5)]<-0
    y<-S*sum(-sign(x)*exp(-abs(x)/omega))
    zdash<-c(zdash,z[i]+y)
  }
  return(zdash)
}


traitevol2<-function(z,S){
  zdash<-z-S*(mean(z)-z)
  return(zdash)
}


#Simulate the process over multiple time steps
n<-20
traits<-rnorm(n,0,4)
primer<-traits
newdat<-traits
for(i in 1:10000){
  newdat<-traitevol(newdat,1,0.05)
  traits<-rbind(traits,newdat)
}

plot(traits[,1],type="l",ylim=c(min(traits),max(traits)),xlab="time",ylab="Traits")
for(i in 2:n){
  lines(traits[,i],col=i)
}

means<-apply(traits,1,mean)
plot(means,ylab="grand mean")

variance<-apply(traits,1,var)
plot(variance,ylab="grand variance")
