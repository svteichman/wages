library(rblm)
require(knitr)
require(kableExtra)

set.seed(324313)
model = m2.mini.new(10,serial = F,fixb=T)

# we set the parameters to something simple
model$A1 = seq(0,2,l=model$nf) # setting increasing intercepts
model$B1 = seq(1,2,l=model$nf) # adding complementarity (increasing interactions)
model$Em = seq(0,1,l=model$nf) # adding sorting (mean type of workers is increasing in k)

# we make the model stationary (same in both periods)
model$A2 = model$A1
model$B2 = model$B1

# setting the number of movers and stayers 
model$Ns   = array(300000/model$nf,model$nf)
model$Nm   = 10*toeplitz(ceiling(seq(100,10,l=model$nf)))

# creating a simulated data set
ad =  m2.mini.simulate(model)

ms    = grouping.getMeasures(ad,"ecdf",Nw=20,y_var = "y1")
grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=250)
ad   = grouping.append(ad,grps$best_cluster,drop=T)

