

ed<-tree$edge
bt1<-nodeHeights(tree)
bt1<-floor(bt1*times/max(bt1))
res<-as.matrix(cbind(ed,bt1))

phylo$tip.label <- 1:length(tree$tip.label)
phylomat <- dist.nodes(tree)
tipmat <- phylomat[1:length(phylo$tip.label), 1:length(phylo$tip.label)]

geo<-matrix(sample(c(0,1),25,replace=T),5,5)


results<-c()
for(j in 1:1000){
  geo<-matrix(sample(c(0,1),200,replace=T),20,10)
  
  if(sum(rowSums(geo)==0)>0){
    geo[which(rowSums(geo)==0),sample(1:ncol(geo),1)]<-1
  }
  mat1<-geomat(geo)
  mat<-tipmat*(exp(-mat1))

  x <- eigen(graph.laplacian(graph.adjacency(mat, weighted = T, 
                                             mode = "undirected"), normalized = T), only.values = F)
  x <- subset(x$values, x$values > 0.1)
  l = length(x)
  integr <- function(x, f) {
    int = 0.5 * sum((x[2:l] - x[1:(l - 1)]) * (f[2:l] + f[1:(l - 
                                                               1)]))
    return(int)
  }
  d <- density(x)
  dsc <- d$y/integr(d$x, d$y)
  splitter = max(x)
  tracer = max(dsc)
  fragmenter = (sum((x - mean(x))^3)/l)/(sum((x - mean(x))^2)/l)^(3/2)
  results<-rbind(results,c(splitter,tracer,fragmenter))
  
  }


results<-as.data.frame(results)
colnames(results)<-c("splitter","tracer","fragmenter")
results$range<-as.factor(results$range)
ggplot(results,aes(x=range,y=splitter))+geom_boxplot()
ggplot(results,aes(x=range,y=fragmenter))+geom_boxplot()
ggplot(results,aes(x=range,y=tracer))+geom_boxplot()



#Test simspat
tr<-rnorm(20,0,2)
alpha2<-1
coln<-c(rep(0.7,3),rep(0.3,7))
coln<-rep(0.5,10)

geo<-matrix(0,nrow=20,ncol=10)
for(i in 1:nrow(geo)){
  geo[i,sample(1:ncol(geo),1)]<-1
}

for(j in 1:10000){
  geo<-simspat(tr,geo,coln,alpha2)
}

geolap=geomat(geo)
plot(density((geolap[upper.tri(geolap)])))
plot(density((geolap[lower.tri(geolap)])))

mean(geolap[upper.tri(geolap)])
mean(geolap[lower.tri(geolap)])

sd(geolap[upper.tri(geolap)])
sd(geolap[lower.tri(geolap)])
#Perform a ABC-SMC model selection exercise on simulated and real finches phylogeny.
#Models to compare: 1. OU model with no explicit spatial dynamic; 2. OU + competition with 
#trait independent spatial dynamic and 3. OU + competition model with explicit trait-dependent
#spatial dynamic.

library(TESS)
library(RPANDA)
library(ape)
library(igraph)
library(phytools)
library(picante)

#Necessary functions

geomat.r<-function(geo){
  mat<-matrix(nrow=nrow(geo),ncol=nrow(geo))
  if(nrow(geo)<2){
    mat<-matrix(2,nrow=1,ncol=1)
  }else{
    for(i in 1:nrow(geo)){
      for(j in 1:nrow(geo)){
        mat[i,j]<-sum(geo[i,]*geo[j,])
      }
    }
  }
  return(as.matrix(mat))
}


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
    mat<-mat/rowSums(geo)
  }
  return(as.matrix(mat))
}

spectR_g<-function (phylo, geolap, draw = F) 
{

  if(length(phylo$tip.label) != nrow(geolap))
    stop("geo does not match phylo")
  phylo$tip.label <- 1:length(phylo$tip.label)
  phylomat <- dist.nodes(phylo)
  tipmat <- phylomat[1:length(phylo$tip.label), 1:length(phylo$tip.label)]
  difs <- geolap
  left <- phylomat[length(phylo$tip.label):dim(phylomat)[1], 
                   length(phylo$tip.label):dim(phylomat)[2]]
  difmat <- tipmat/(1+difs)
  x <- eigen(graph.laplacian(graph.adjacency(difmat, weighted = T, 
                                             mode = "undirected"), normalized = T), only.values = F)
  x <- subset(x$values, x$values > 0.1)
  l = length(x)
  integr <- function(x, f) {
    int = 0.5 * sum((x[2:l] - x[1:(l - 1)]) * (f[2:l] + f[1:(l - 
                                                               1)]))
    return(int)
  }
  d <- density(x)
  dsc <- d$y/integr(d$x, d$y)
  splitter = max(x)
  tracer = max(dsc)
  fragmenter = (sum((x - mean(x))^3)/l)/(sum((x - mean(x))^2)/l)^(3/2)
  res <- list(eigenvalues = x, splitter = splitter, tracer = tracer, 
              fragmenter = fragmenter)
  if (draw == T) {
    par(mar = c(4, 5, 1, 1))
    plot(d$x, dsc, type = "l", xlab = expression(""[n] * 
                                                   lambda), ylab = expression(paste("f(x)/", integral(f(y) * 
                                                                                                        dy)), sep = ""))
    polygon(d$x, dsc, col = colors(1)[runif(1, 1, 500)])
    return(res)
  }
  else {
    return(res)
  }
}


spectR_tg<-function (phylo, dat, geolap,draw = F) 
{
  if (length(phylo$tip.label) != length(dat)) 
    stop("dat do not match phylo")
  if (length(phylo$tip.label) != nrow(geolap)) 
    stop("geo do not match phylo")
  phylo$tip.label <- 1:length(phylo$tip.label)
  phylomat <- dist.nodes(phylo)
  tipmat <- phylomat[1:length(phylo$tip.label), 1:length(phylo$tip.label)]
  geolap <-2/(1+geolap)
  difs <- as.matrix(dist(dat))*geolap
  
  left <- phylomat[length(phylo$tip.label):dim(phylomat)[1], 
                   length(phylo$tip.label):dim(phylomat)[2]]
  difmat <- difs * tipmat
  x <- eigen(graph.laplacian(graph.adjacency(difmat, weighted = T, 
                                             mode = "undirected"), normalized = T), only.values = F)
  x <- subset(x$values, x$values > 0.1)
  l = length(x)
  integr <- function(x, f) {
    int = 0.5 * sum((x[2:l] - x[1:(l - 1)]) * (f[2:l] + f[1:(l - 
                                                               1)]))
    return(int)
  }
  d <- density(x)
  dsc <- d$y/integr(d$x, d$y)
  splitter = max(x)
  tracer = max(dsc)
  fragmenter = (sum((x - mean(x))^3)/l)/(sum((x - mean(x))^2)/l)^(3/2)
  res <- list(eigenvalues = x, splitter = splitter, tracer = tracer, 
              fragmenter = fragmenter)
  if (draw == T) {
    par(mar = c(4, 5, 1, 1))
    plot(d$x, dsc, type = "l", xlab = expression(""[n] * 
                                                   lambda), ylab = expression(paste("f(x)/", integral(f(y) * 
                                                                                                        dy)), sep = ""))
    polygon(d$x, dsc, col = colors(1)[runif(1, 1, 500)])
    return(res)
  }
  else {
    return(res)
  }
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



#function to get all summary statistics
#spectral density profile with and without spatial information

get.stat<-function(tree,tr,geo){

  geolap<-geomat(geo)
  means<-mean(geolap[lower.tri(geolap)])
  sds<-sd(geolap[lower.tri(geolap)])
  
  split.t<-spectR_t(tree,tr)$splitter
  frag.t    <-spectR_t(tree,tr)$fragmenter
  trace.t  <-spectR_t(tree,tr)$tracer
  
  
  return(list(split.t=split.t,frag.t=frag.t,trace.t=trace.t,
              geo_mean=means,geo_sd=sds))
  
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



evol<-function(param,coln,br){

  sigma  <-param$sigma
  theta  <-param$theta
  psi    <-param$psi
  alpha  <-param$alpha
  alpha2 <-param$alpha2
  
  
  tr<-as.matrix(c(0,rnorm(1,0,sigma)))
  tr<-cbind(which(br[,3]==0),tr)
  
  #randomize the initial geographic distribution for MRCA
  geo<-matrix(rep(0,2*length(coln)),2,length(coln))
  geo[1,sample(1:length(coln),1)]<-1
  geo[2,sample(1:length(coln),1)]<-1
  geo<-cbind(which(br[,3]==0),geo)
  
  tips<-1:nrow(tr)
  
  for(i in 1:200){
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
  ind<-which(br[,2]%in%(1:20))
  tr<-tr[which(tr[,1]%in%ind),]
  geo<-geo[which(tr[,1]%in%ind),]
  br1<-br[tr[,1],1:2]
  trait<-cbind(br1[,2],tr[,2])
  trait<-trait[order(trait[,1]),2]
  geo<-cbind(br1[,2],geo[,-1])
  geo<-geo[order(geo[,1]),-1]
  
  return(list(trait=trait,geo=geo))
}




evol2<-function(param,br){
  
  sigma  <-param$sigma
  theta  <-param$theta
  psi    <-param$psi
  alpha  <-param$alpha

  
  tr<-as.matrix(c(0,rnorm(1,0,sigma)))
  tr<-cbind(which(br[,3]==0),tr)
  
  
  tips<-1:nrow(tr)
  
  for(i in 1:200){
    #If there is a branching event, add trait and geo characteristics  
    if(i%in%br[,3]){  
      
      parent<-tr[which(tr[,1]%in%which(br[,4]==i)),1]
      
      for(i1 in parent){
        tips<-c(tips,nrow(tr)+1,nrow(tr)+2)
        
        new.ind<-which(br[,1]==br[i1,2])
        p.tr<-tr[which(tr[,1]==i1),2]
        tr<-rbind(tr,cbind(new.ind,c(p.tr,p.tr+rnorm(1,0,sigma))))
        
        tips<-tips[!tips%in%which(tr[,1]==i1)]
      }
    }
    newtr<-tr[tips,-1]
    newtr<-newtr + psi*(theta-newtr) + rnorm(length(newtr),0,sigma)
    for(i2 in 1:length(newtr)){
      newtr[i2]<-newtr[i2]+alpha*sum(((newtr[i2]-newtr)*exp(-abs(newtr[i2]-newtr))))
    }
    tr[tips,2]<-newtr
    
  }
  #Outputs: traits and geo data
  ind<-which(br[,2]%in%(1:20))
  tr<-tr[which(tr[,1]%in%ind),]
  br1<-br[tr[,1],1:2]
  trait<-cbind(br1[,2],tr[,2])
  trait<-trait[order(trait[,1]),2]
  
  return(trait)
}

#Sample dataset and params

tree<-tess.sim.taxa(n=1,nTaxa=20, lambda=1, mu=0, max=5)[[1]]
ed<-tree$edge
bt1<-nodeHeights(tree)
bt1<-floor(bt1*times/max(bt1))
res<-as.matrix(cbind(ed,bt1))

param<-list(sigma=0.1,theta=0,psi=0,alpha=0,alpha2=1)

coln<-rep(0.5,10)



big.g<-c()

for(j in 1:50){
  for(i in c(0,0.05,0.1)){
  param$psi<-i
  res1<-evol(param,coln,res)$geo
  geolap<-geomat.r(res1)

  st<-spectR_g(tree,geolap)

  big.g<-rbind(big.g,c(i,st$splitter,st$fragmenter,st$tracer))
  }
}

 geodat<-cbind(rep(c(0,0.1,0.2),50),geodat)
geodat<-as.data.frame(geodat)
colnames(geodat)<-c("al","mean1","mean2","sd1","sd2")
geodat$al<-as.factor(geodat$al)
ggplot(geodat,aes(x=al,y=mean1))+geom_boxplot()



big.g<-as.data.frame(big.g)
colnames(big.g)<-c("sig","splitter","fragmenter","tracer")
big.g$sig<-as.factor(big.g$sig)



ggplot(big.g,aes(x=sig,y=splitter))+geom_boxplot()


ggplot(big.g,aes(x=sig,y=fragmenter))+geom_boxplot()


ggplot(big.g,aes(x=sig,y=tracer))+geom_boxplot()

summary(lm(big1$splitter~big1$sig))

#########################################################################
#Set up a batch of simulations with different combinations of parameters

sigmas<-c(0.1)
thetas<-c(0)
psis<-c(0,0.05,0.1)
alphas<-c(0,0.05,0.1)
alpha2s<-c(0,1,2)
coln<-c(rep(0.7,3),rep(0.3,7))
coln<-rep(0.5,10)

paramdat<-as.matrix(expand.grid(sigmas,thetas,psis,alphas,alpha2s))
colnames(paramdat)<-c("sigma","theta","psi","alpha","alpha2")

#Start the cluster and split the loop on #cores
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

times<-200
trees <- tess.sim.taxa(n=50,nTaxa=20, lambda=1, mu=0, max=5)
brdat<-list()


for(i in 1:50){
  tree<-trees[[i]]
  ed<-tree$edge
  bt1<-nodeHeights(tree)
  bt1<-floor(bt1*times/max(bt1))
  res<-as.matrix(cbind(ed,bt1))
  res[which(res[,2]%in%(1:20)),4]<-200
  while(sum(res[,3]==res[,4])!=0){
    ind<-which(res[,3]==res[,4])
    for(j in ind){
      res[j,4]<-res[j,4]+1
      res[which(res[,1]==res[j,2]),3]<-res[which(res[,1]==res[j,2]),3]+1
    }
  }
  brdat<-c(brdat,list(res))
}


numcores<-detectCores()
clust<-makeCluster(numcores)
clusterExport(clust,c("geomat","simspat","evol","trees","spectR_t"
                      ,"get.stat","paramdat","coln","brdat"))
registerDoParallel(clust)




sim.res<-foreach(j1=1:length(brdat),.combine=rbind,
                  .packages=c("RPANDA","TESS","ape","igraph","phytools"))%dopar%
{
  br<-brdat[[j1]]
  res1<-c()
  for(j2 in 1:nrow(paramdat)){

      param<-as.list(paramdat[j2,])

      res<-evol(param,coln,br)
      tr<-res$trait
      geo1<-res$geo
  
    res1<-rbind(res1,c(unlist(param),unlist(get.stat(trees[[j1]],tr,geo1))))
  }
  return(res1)
}
stopCluster(clust)

write.csv(sim.res,"phylosim1_het.csv")
  


##############################################################################
#Plot results
library(ggplot2)
library(stringr)
library(plotly)


dat<-read.csv("phylosim1_hom.csv")

dat<-dat[,c(4:11)]

dat$psi<-as.factor(dat$psi)
dat$alpha<-as.factor(dat$alpha)
dat$alpha2<-as.factor(dat$alpha2)


#Perform ANOVA

model.split.t<-lm(split.t~psi*alpha*alpha2,data=dat)
summary(aov(model.split.t))
capture.output(summary(aov(model.split.t)),file="split_hom.doc")

model.frag.t<-lm(frag.t~psi*alpha*alpha2,data=dat)
summary(aov(model.frag.t))
capture.output(summary(aov(model.frag.t)),file="frag_hom.doc")

model.trace.t<-lm(trace.t~psi*alpha*alpha2,data=dat)
summary(aov(model.trace.t))
capture.output(summary(aov(model.trace.t)),file="trace_hom.doc")


model.geomean<-lm(geo_mean~psi*alpha*alpha2,data=dat)
summary(aov(model.geomean))
capture.output(summary(aov(model.geomean)),file="geomean_hom.doc")

model.geosd<-lm(geo_sd~psi*alpha*alpha2,data=dat)
summary(aov(model.geosd))
capture.output(summary(aov(model.geosd)),file="geosd_hom.doc")


colnames(dat)[3]<-"alpha2"

bxp1.t<-ggplot(dat,aes(x=psi,y=split.t,color=alpha))+geom_boxplot()+
  facet_grid(~paste0("\u03b1","c","=",alpha2))+ylab("splitter")+ylim(range(dat[,c(4)]))+
  labs(color=paste0("\u03b1"),x=paste0("\u03c8"))+
theme(
  strip.text.x = element_text(
    size = 11, face = "bold"
  ),
  strip.text.y = element_text(
    size = 10, face = "bold"
  ),
  legend.title = element_text(color = "black", size = 13,face="bold"),
  legend.key.size = unit(1, "cm"),
  legend.text=element_text(size=12),
  axis.title=element_text(size=13,face="bold")
)
bxp1.t



bxp2.t<-ggplot(dat,aes(x=psi,y=trace.t,color=alpha))+geom_boxplot()+
  facet_grid(~paste0("\u03b1","c","=",alpha2))+ylab("tracer")+ylim(range(dat[,c(6)]))+
  labs(color=paste0("\u03b1"),x=paste0("\u03c8"))+theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 13,face="bold"),
    legend.key.size = unit(1, "cm"),
    legend.text=element_text(size=12),
    axis.title=element_text(size=13,face="bold")
  )
bxp2.t


bxp3.t<-ggplot(dat,aes(x=psi,y=frag.t,color=alpha))+geom_boxplot()+
  facet_grid(~paste0("\u03b1","c","=",alpha2))+ylab("Fragmenter")+ylim(range(dat[,c(5)]))+
  labs(color=paste0("\u03b1"),x=paste0("\u03c8"))+theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 13,face="bold"),
    legend.key.size = unit(1, "cm"),
    legend.text=element_text(size=12),
    axis.title=element_text(size=13,face="bold")
  )
bxp3.t


bxp4<-ggplot(dat,aes(x=psi,y=geo_mean,color=alpha))+geom_boxplot()+
  facet_grid(~paste0("\u03b1","c","=",alpha2))+ylab("Mean geographic overlap")+ylim(range(dat[,c(7)]))+
  labs(color=paste0("\u03b1"),x=paste0("\u03c8"))+theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 13,face="bold"),
    legend.key.size = unit(1, "cm"),
    legend.text=element_text(size=12),
    axis.title=element_text(size=13,face="bold")
  )
bxp4

bxp5<-ggplot(dat,aes(x=psi,y=geo_sd,color=alpha))+geom_boxplot()+
  facet_grid(~paste0("\u03b1","c","=",alpha2))+ylab("SD in geographic overlap")+ylim(range(dat[,c(8)]))+
  labs(color=paste0("\u03b1"),x=paste0("\u03c8"))+theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 13,face="bold"),
    legend.key.size = unit(1, "cm"),
    legend.text=element_text(size=12),
    axis.title=element_text(size=13,face="bold")
  )
bxp5


dat1<-dat
dat1$model<-NA
dat1[(dat1$psi==0 & dat1$alpha==0 & dat1$alpha2==0),9]<-"BM"
dat1[(dat1$psi!=0 & dat1$alpha==0 & dat1$alpha2==0),9]<-"OU"
#dat1[(dat1$psi!=0 & dat1$alpha!=0 & dat1$alpha2==0),9]<-"OU+Competition"
dat1[(dat1$psi!=0 & dat1$alpha!=0),9]<-"OU+Competition"

dat1<-na.omit(dat1)

dat1$trace.t<-log(dat1$trace.t)

dat1$model<-as.factor(dat1$model)


box.model1<-ggplot(dat1,aes(x=model,y=split.t))+geom_boxplot()+
  ylab("Splitter")+
  theme(
    panel.background = element_blank()
    axis.title=element_text(size=13,face="bold"),
    axis.text = element_text(size = 11,face="bold")
  )
box.model1

box.model2<-ggplot(dat1,aes(x=model,y=frag.t))+geom_boxplot()+
  ylab("Fragmenter")+
  theme(
    axis.title=element_text(size=13,face="bold"),
    axis.text = element_text(size = 11,face="bold")
  )
box.model2

box.model3<-ggplot(dat1,aes(x=model,y=trace.t))+geom_boxplot()+
  ylab("Tracer")+
  theme(
    axis.title=element_text(size=13,face="bold"),
    axis.text = element_text(size = 11,face="bold")
  )
box.model3

box.model4<-ggplot(dat1,aes(x=model,y=geo_mean))+geom_boxplot()+
  ylab("Mean Geographic Overlap")+
  theme(
    axis.title=element_text(size=13,face="bold"),
    axis.text = element_text(size = 11,face="bold")
  )
box.model4

box.model5<-ggplot(dat1,aes(x=model,y=geo_sd))+geom_boxplot()+
  ylab("S.D. in Geographic Overlap")+
  theme(
    axis.title=element_text(size=13,face="bold"),
    axis.text = element_text(size = 11,face="bold")
  )
box.model5


fig1<-plot_ly(dat1,x=~trace.t,y=~frag.t,z=~split.t,
             color=~model,colors=c("red","blue","green"))
fig1 <- fig1 %>% add_markers()
fig1

fig2<-plot_ly(dat1,x=~geo_mean,y=~geo_sd,
             color=~model,colors=c("red","blue","green"))
fig2 <- fig2 %>% add_markers()
fig2

fig3<-plot_ly(dat1,x=~trace.t,y=~geo_sd,z=~geo_mean,
              color=~model,colors=c("red","blue","green"))
fig3 <- fig3 %>% add_markers()
fig3


