library("DiceDesign")
library("DiceKriging")
library("DiceView")


test<-function(model){
  resLOO <- leaveOneOut.km(model, type="UK", trend.reestim=FALSE)
  Q <- 1 - sum((resLOO$mean - model@y)^2) / sum( (model@y - mean(model@y))^2)
  M <- mean((model@y-resLOO$mean)^2)
  standardised_residuals <- (model@y-resLOO$mean)/resLOO$sd
  return(list(Q2=Q, mse= M, sd_error=sd(standardised_residuals)))
}

dist2<-function(X){
  distances<-matrix(1, nrow = nrow(X), ncol = nrow(X))
  for(i in 1:nrow(X)){
    for(j in 1:nrow(X)){
      if(i==j){distances[i,j]<-0}
      else if(distances[j,i] == 1){distances[i,j]<-norm(X[i,]-X[j,], type="2")}
      else{ distances[i,j]<-distances[j,i] }
    }
  }
  return(distances)}




generateCVT <- function(npts, dim, nite){
 
  X0<-replicate(dim, runif(npts,0,1))
  J<-rep(1,npts)
  for(i in 1:nite){
    w<-t(runif(dim,0,1))
    distances<-as.matrix(dist(rbind(w, X0), method='euclidean'))[1,-1] 
    #on récupère la distance entre x1 et x2...xn et w1
    crit<-which.min(distances)
    pt<-X0[crit,]
    nouv<-(J[crit]*pt + w)/(J[crit]+1)
    X0[crit,]<-nouv
    
    J[crit]<-J[crit]+1
  }  
  return(X0)
}



generateLHS <- function(npts, dim){
  temp<-replicate(dim, runif(npts,0,1))
  return(temp)}

evalMinDist <- function(X){ 
  D<-as.matrix(dist(X, method='euclidean'))
  temp<-matrix(NA, nrow = nrow(X), ncol = 1)
  for(i in 1:nrow(X)) { 
    s<-as.vector(D[i,-i]) #retire le pt
    ind<-which.min(s)
    temp[i,]<-s[ind] }
  d<-temp
  return( list(minDist = d, allDist = D) )
}
#minDist correspond à la plus petite distance pour chaque point Xi


MinDistAvcAutresPts<-function(X,ind.pt){
  temp<-matrix(NA, nrow = nrow(X), ncol = 1)
  for(i in 1:nrow(X)){
    if(i != ind.pt){temp[i,]<-norm(X[i,]-X[ind.pt,], type="2") }
    else{temp[ind.pt,]<-2 } #on lui attribue 2, car on sait que dans tous les cas la distance ne dépasse
    #pas sqrt(d) avec d=2 (ici distance de la diagonale)
  }
  mindist<-temp[which.min(temp)] #la plus petite distance avec les autres pts
  return(mindist)
}

#on doit remplit le plan
optimiseLHS<-function(lhs,nite,T0){
  for(i in 1:nite){
    Temp<-T0/(log(i))
     rnd<-runif(1, min=0, max=1)
     lhs.temp<-lhs
     tmp<-evalMinDist(lhs)$minDist
     ind.ptCritique<-which.min(tmp)
     distActuelle<-tmp[ind.ptCritique]
     rnd.col<-round(runif(1, min=1, max=ncol(lhs))) #numero de colonne du pt critique aleatoire
     repeat{#on repete jusqua ce que l'indice du pt aleatoire soit different à l'indice du pt critique
       rnd.pt<-round(runif(1, min=1, max=nrow(lhs)))
       if(rnd.pt!=ind.ptCritique){break}
     }
     
     #on cree des variables temporaires (ici les points apres echange)
     temp1<-lhs[ind.ptCritique,]
     temp1[rnd.col]<-lhs[rnd.pt,rnd.col]
     
     temp2<-lhs[rnd.pt,]
     temp2[rnd.col]<-lhs[ind.ptCritique,rnd.col]
     
     #on cree une matrice temporaire (ici apres substitution des coordonnees de deux points)
     lhs.temp[ind.ptCritique,]<-temp1
     lhs.temp[rnd.pt,]<-temp2
     
     #on recupere la nouvelle distance minimum de l'ancien point critique 
     distNouvelle<-evalMinDist(lhs.temp)$minDist[ind.ptCritique]
     if(distNouvelle>distActuelle) {lhs<-lhs.temp
     distActuelle<-distNouvelle}
     
     else if(rnd<=exp(-(distActuelle-distNouvelle)/Temp)){
     lhs<-lhs.temp
     distActuelle<-distNouvelle}
  } 
  return(lhs)}
 

#MinDistAvcAutresPts(lhs.temp, ind.ptCritique)

discrepancy <- function (design){
  X <- as.matrix(design)
  dimension <- dim(X)[2]
  n <- dim(X)[1]
  
  dL2 <- 0
  for (j in 1:n) {
    for (i in 1:n) {
      if (i != j) {
        t <- c()
        for (l in 1:dimension) t <- c(t, 1 - max(X[i, 
                                                   l], X[j, l]))
        t <- (prod(t))/(n^2)
      }
      else {
        t1 <- 1 - X[i, ]
        t1 <- prod(t1)
        t2 <- 1 - X[i, ]^2
        t2 <- prod(t2)
        t <- t1/(n^2) - ((2^(1 - dimension))/n) * t2
      }
      dL2 <- dL2 + t
    }
  }
  DisL2star <- sqrt(3^(-dimension) + dL2)
  return(DisL2star)
}



#Dimension 2, 10 pts
cvt2<-generateCVT(10,2,5000)
lhs2<-generateLHS(10,2)
lhs_optim2<-optimiseLHS(lhs2,5000,10) #plus T est grand plus les pts vont se deplacer

par(mfrow=c(2,2))
plot(cvt2)
plot(lhs2)
plot(lhs_optim2)
par(op)


par(mfrow=c(2,2))
hist(cvt2[,1])
hist(cvt2[,2])

hist(lhs2[,1])
hist(lhs2[,2])

hist(lhs_optim2[,1])
hist(lhs_optim2[,2])

par(op)

par(mfrow=c(2,2))
pairs(cvt2)
pairs(lhs2)
pairs(lhs_optim2)
par(op)



#critere maximin
mindistCVT<-evalMinDist(cvt2)$minDist
minValCvt2<-mindistCVT[which.min(mindistCVT)]

mindistLHS<-evalMinDist(lhs2)$minDist
minValLHS2<-mindistLHS[which.min(mindistLHS)]

mindistLHSoptim<-evalMinDist(lhs_optim2)$minDist
minValLHS2optim<-mindistLHSoptim[which.min(mindistLHSoptim)]

comparaison1<-cbind(minValCvt2,minValLHS2,minValLHS2optim)
comparaison1 #on cherche celui qui a la plus grde maximin

#discrepancy
disc_cvt2<-discrepancy(cvt2)
disc_lhs2<-discrepancy(lhs2)
disc_lhsOptim2<-discrepancy(lhs_optim2)
comparaison2<-cbind(disc_cvt2,disc_lhs2,disc_lhsOptim2)
comparaison2#on cherche la plus petite discrepancy



#Dimension 5, 70 pts
cvt5<-generateCVT(70,5,5000)
lhs5<-generateLHS(70,5)
lhs_optim5<-optimiseLHS(lhs5,5000,10) #plus T est grand plus les pts vont se deplacer
#lhs_optim50<-maximinSA_LHS(lhs5,T0=10,c=0.95,it=2000,p=50,profile="GEOM",Imax=100)
lhs_optim5<-maximinSA_LHS(lhs5,T0=10,c=0.95,it=2000,p=50,profile="GEOM",Imax=100)$design


op<-par(mfrow=c(3,2))
for(i in 1:5){hist(cvt5[,i], main="Tessellation centroïdale de Voronoï")}

op<-par(mfrow=c(3,2))

for(i in 1:5){hist(lhs5[,i], main="Hypercube latin")}
op<-par(mfrow=c(3,2))

for(i in 1:5){hist(lhs_optim5[,i], main="Hypercube latin optimisé")}
par(op)

op<-par(mfrow=c(2,2))
pairs(cvt5, main="Tessellation centroïdale de Voronoï")
pairs(lhs5, main="Hypercube latin")
pairs(lhs_optim5, main="Hypercube latin optimisé")
par(op)

  #critere maximin
mindistCVT<-evalMinDist(cvt5)$minDist
minValCvt5<-mindistCVT[which.min(mindistCVT)]

mindistLHS<-evalMinDist(lhs5)$minDist
minValLHS5<-mindistLHS[which.min(mindistLHS)]

mindistLHSoptim<-evalMinDist(lhs_optim5)$minDist
minValLHS5optim<-mindistLHSoptim[which.min(mindistLHSoptim)]

comparaison1<-cbind(minValCvt5,minValLHS5,minValLHS5optim)
rownames(comparaison1)<-'maximin' #on cherche celui qui a la plus grde maximin
comparaison1
barplot.default(comparaison1, col=4)

#discrepancy
disc_cvt5<-discrepancy(cvt5)
disc_lhs5<-discrepancy(lhs5)
disc_lhsOptim5<-discrepancy(lhs_optim5)
comparaison2<-cbind(disc_cvt5,disc_lhs5,disc_lhsOptim5)
rownames(comparaison2)<-'discrepancy'
comparaison2#on cherche la plus petite discrepancy
barplot.default(comparaison2, col=4)

#phiP criterion
phip_cvt5<-phiP(cvt5)
phip_lhs5<-phiP(lhs5)
phip_lhsOptim5<-phiP(lhs_optim5)
comparaison3<-cbind(phip_cvt5,phip_lhs5,phip_lhsOptim5)
rownames(comparaison3)<-'phiP'
comparaison3#on cherche la plus petite phiP (ca revient à la distance...)
barplot.default(comparaison3, col=4)

#mesh ration
mesh_cvt5<-meshRatio(cvt5)
mesh_lhs5<-meshRatio(lhs5)
mesh_lhsOptim5<-meshRatio(lhs_optim5)
comparaison4<-cbind(mesh_cvt5,mesh_lhs5,mesh_lhsOptim5)
rownames(comparaison4)<-'mesh ratio'

comparaison4#on chercheune mesh ratio qui s'approche de 1
barplot.default(comparaison4, col=4)

#coverage
cov_cvt5<-coverage(cvt5)
cov_lhs5<-coverage(lhs5)
cov_lhsOptim5<-coverage(lhs_optim5)
comparaison5<-cbind(cov_cvt5,cov_lhs5,cov_lhsOptim5)
rownames(comparaison5)<-'coverage'
comparaison5#on chercheune coverage qui s'approche de 0
barplot.default(comparaison5, col=4)



#Dimension 10, 150 pts
cvt10<-generateCVT(150,10,5000)
lhs10<-generateLHS(150,10)
lhs_optim10<-optimiseLHS(lhs10,5000,10) #plus T est grand plus les pts vont se deplacer




par(mfrow=c(2,2))
for(i in 1:10){hist(cvt10[,i])}
for(i in 1:10){hist(lhs10[,i])}
for(i in 1:10){hist(lhs_optim10[,i])}
par(op)

par(mfrow=c(2,2))
pairs(cvt5)
pairs(lhs5)
pairs(lhs_optim5)
par(op)

#critere maximin
mindistCVT<-evalMinDist(cvt10)$minDist
minValCvt10<-mindistCVT[which.min(mindistCVT)]

mindistLHS<-evalMinDist(lhs10)$minDist
minValLHS10<-mindistLHS[which.min(mindistLHS)]

mindistLHSoptim<-evalMinDist(lhs_optim10)$minDist
minValLHS10optim<-mindistLHSoptim[which.min(mindistLHSoptim)]

comparaison1<-cbind(minValCvt10,minValLHS10,minValLHS10optim)
comparaison1 #on cherche celui qui a la plus grde maximin

#discrepancy
disc_cvt10<-discrepancy(cvt10)
disc_lhs10<-discrepancy(lhs10)
disc_lhsOptim10<-discrepancy(lhs_optim10)
comparaison2<-cbind(disc_cvt10,disc_lhs10,disc_lhsOptim10)
comparaison2#on cherche la plus petite discrepancy





#krigeage
#plan.saved<-lhs_optim5 #meilleur repartition, minVal critere aussi
#write.table(plan.saved, "plan.csv", sep = ";",quote = FALSE, row.names = FALSE)

plan.saved<-read.csv("plan.csv", sep=';')
colnames(plan.saved) <- c("longitude", "latitude", "elevation", "radius", 'overpressure')
Jx<-compute_wls(plan.saved)
plan.saved<-as.data.frame(plan.saved)

##plot Y against the inputs
op<-par(mfrow=c(2,2)) # split the plotting area in 5
plot(plan.saved$longitude,Jx)
plot(plan.saved$latitude,Jx)
plot(plan.saved$elevation,Jx)
plot(plan.saved$radius,Jx)#ca a l'air d'etre le plus corrélé
plot(plan.saved$overpressure,Jx)
par(op)
pairs(plan.saved, main="Plan d'experience")

sobol.ind<-fast99(model = compute_wls, factors = 5, n = 1000, q = "qunif",     q.arg = list(min = 0, max = 1))
plot(sobol.ind)


n<-400
X1 <- runif(n,min = 0, max = 1 )
X2 <- runif(n,min = 0, max = 1 )
X3 <- runif(n,min = 0, max = 1 )
X4<-seq(0,1,length = 20)
X5<-seq(0,1,length = 20)
X<-cbind(X1,X2,X3,expand.grid(X4,X5))
colnames(X) <- c("longitude", "latitude", "elevation", "radius", 'overpressure')
resp<-compute_wls(X)


#Changer les kernels
m1<-km(~-1,design=plan.saved, Jx, estim.method = "MLE", covtype="gauss")
plot(m1)

m2<-km(~-1,design=plan.saved, Jx, estim.method = "MLE", covtype="powexp")
plot(m2)

m3<-km(~-1,design=plan.saved, Jx, estim.method = "MLE", covtype="exp")
plot(m3)

m4<-km(~-1,design=plan.saved, Jx, estim.method = "MLE", covtype="matern5_2")
plot(m4)

m5<-km(~-1,design=plan.saved, Jx, estim.method = "MLE", covtype="matern3_2")
plot(m5)

Comparison<-matrix(NA,nrow=3,ncol=5)
colnames(Comparison) <- c("gauss", "powexp", "exp", "matern5_2", 'matern3_2')
rownames(Comparison) <- c("Q2", "mse", "sd_error")

Comparison[1,]<-cbind(test(m1)$Q2,test(m2)$Q2,test(m3)$Q2,test(m4)$Q2,test(m5)$Q2)
Comparison[2,]<-cbind(test(m1)$mse,test(m2)$mse,test(m3)$mse,test(m4)$mse,test(m5)$mse)
Comparison[3,]<-cbind(test(m1)$sd_error,test(m2)$sd_error,test(m3)$sd_error,test(m4)$sd_error,test(m5)$sd_error)
Comparison


#Après avoir choisir le kernel (gauss), on va essayer de changer le trend
m12<-km(~-1+I((radius-elevation)^2)+(radius*elevation*overpressure),design=plan.saved, Jx, estim.method = "MLE", covtype="gauss")
plot(m12)

m13<-km(~1+(radius+elevation+overpressure)^3,design=plan.saved, Jx, estim.method = "MLE", covtype="gauss")
plot(m13)

m14<-km(~-1+radius/elevation+(radius-longitude)*(radius-latitude)*(radius-elevation)+longitude/latitude,design=plan.saved, Jx, estim.method = "MLE", covtype="gauss")
plot(m14)

m15<-km(~1+radius+overpressure+elevation/radius,design=plan.saved, Jx, estim.method = "MLE", covtype="gauss")
plot(m15)

Comparison2<-matrix(NA,nrow=3,ncol=5)
colnames(Comparison2) <- c("initial", "test1", "test2", "test3",'test4')
rownames(Comparison2) <- c("Q2", "mse", "sd_error")

Comparison2[1,]<-cbind(test(m4)$Q2,test(m12)$Q2,test(m13)$Q2,test(m14)$Q2,test(m15)$Q2)
Comparison2[2,]<-cbind(test(m4)$mse,test(m12)$mse,test(m13)$mse,test(m14)$mse,test(m15)$mse)
Comparison2[3,]<-cbind(test(m4)$sd_error,test(m12)$sd_error,test(m13)$sd_error,test(m14)$sd_error,test(m15)$sd_error)
Comparison2


#test m1, les résidus
res <- leaveOneOut.km(m13, type="UK", trend.reestim=FALSE)
std_res<-   (m13@y-res$mean)/res$sd

qqnorm(std_res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE)
qqline(std_res,col="red")

smir<-ks.test(std_res,"pnorm",mean(std_res),sd(std_res))
D<-smir$statistic
0.895-D*(sqrt(70)+0.85)/(sqrt(70)-0.01)#si positive ca passe

Box.test(std_res, type='Ljung-Box') #non correlation entre les residus
#au final suit une loi normale et non correles entre eux, donc independance
#donc modele valide



sectionview.km(m13, type='UK',center = as.matrix(plan.saved[which.min(Jx),]))


# pred 

predicted.values.model1 <- predict(m13, X, "UK",checkNames=FALSE)
plot(resp,predicted.values.model1$mean)
abline(0,1, col="red")

tst<-predicted.values.model1$sd
X[which.max(tst),]
plan.saved<-rbind(plan.saved,X[which.max(tst),])
Jx<-compute_wls(plan.saved)









#cas test

pl_exp<- plan.saved
rep_wls <- compute_wls(pl_exp)
plot(rep_wls, ylab = "Observation", main = "Obsevation des valeurs calculées par compute_wls")

sol <- unnorm_var(pl_exp)
pairs(cbind(sol,rep_wls))

fitted.model1 <- km(~1, design=pl_exp, response=rep_wls,
                    covtype="matern5_2", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model1)

library(rgenoud)
nsteps <- 20
lower <- rep(0,dim)
upper <- rep(1,dim)
oEGO <- EGO.nsteps(model=fitted.model1, fun=compute_wls, nsteps=nsteps,
                   lower=lower, upper=upper, control=list(pop.size=20, BFGSburnin=2))
print(oEGO$par)
val_ego <- oEGO$value
print(val_ego)

unn_ego<- unnorm_var(oEGO$par)
pairs(cbind(val_ego,unn_ego))

plot(rep_wls, ylab = "Observation",
     main = "Convergence de la fonction objectif \n Nouveaux points en bleu")
points(oEGO$value, col = "red")
pairs(oEGO$lastmodel@X, main = "Distribution des points dans le plan")
hist(oEGO$value)

min(oEGO$value)
ind <- which.min(oEGO$value)
argMin <- unn_ego[ind,]
argMin
