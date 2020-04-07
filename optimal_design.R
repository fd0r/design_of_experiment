library(AlgDesign)
library(ggplot2)
library(fields)

levels<-seq(-1,1,by=.1)

dat<-expand.grid(list(A=levels,B=levels))
desL<-optFederov(~quad(.), dat, nTrials=16, criterion="A", eval=TRUE)

ggplot(desL$design, aes(x=A, y=B))+geom_point(size=6, color="#69b3a2")


prediction_variance <- function(X, x0, sigma=1){
  return (1+t(x0)%*%solve(t(X)%*%X)%*%x0)*sigma*sigma
}

prediction_variance(matrix(c(desL$design$A, desL$design$B), ncol=2), c(0.6,0.7), sigma)


levels<-seq(-1,1,by=.05)
data <- expand.grid(list(X1=levels,X2=levels))

D_plan = optFederov(~quad(.), data, nTrials=16, criterion="D", eval=TRUE)
I_plan = optFederov(~quad(.), data, nTrials=16, criterion="I", eval=TRUE)

X1 <- levels
X2 <- levels

ggplot(I_plan$design, aes(x=X1, y=X2))+geom_point(size=6, color="blue")+ggtitle("I optimal plan")

D_matrix = matrix(c(rep(1, length(D_plan$design$X1)), D_plan$design$X1, D_plan$design$X2,
                    D_plan$design$X1**2, D_plan$design$X2**2, D_plan$design$X1*D_plan$design$X2), ncol=6)

var_D = matrix(rep(0, length(levels)**2), ncol=length(levels))

for (i in 1:length(levels)){
  for (j in 1:length(levels))
  var_D[i,j] <- prediction_variance(D_matrix, c(1, levels[i], levels[j], levels[i]**2, levels[j]**2, levels[i]*levels[j]))
}

I_matrix = matrix(c(rep(1, length(I_plan$design$X1)), I_plan$design$X1, I_plan$design$X2,
                    I_plan$design$X1**2, I_plan$design$X2**2, I_plan$design$X1*I_plan$design$X2), ncol=6)

var_I = matrix(rep(0, length(levels)**2), ncol=length(levels))

for (i in 1:length(levels)){
  for (j in 1:length(levels))
    var_I[i,j] <- prediction_variance(I_matrix, c(1, levels[i], levels[j], levels[i]**2, levels[j]**2, levels[i]*levels[j]))
}

image.plot(X1, X2, var_D, zlim=c(1.1, 1.7))
image.plot(X1, X2, var_I, zlim=c(1.1, 1.7))

