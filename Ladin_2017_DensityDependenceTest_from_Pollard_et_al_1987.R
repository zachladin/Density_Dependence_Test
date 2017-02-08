#R script to run randomization test for density dependence from Pollard et al. (1987),
#These functions were used to produce Fig. 7 from Ladin, Z. S., and C. K. Williams. 2017. Detecting and Analyzing Density Dependence in 

plotmath()

#clear environment
rm(list=ls())

##############################################################################################################
#load packages for producing figures
require(ggplot2)
require(gridExtra)

##############################################################################################################
#Simulate population data that is density-independent "random walk with drift" (Eq. 2 from Pollard et al. 1987)
#Eq. 2. x[i+1]=r + x[i]+ e[i]
#where e is a random variate from a normal distribution with mean 0 and sd=0.1
#init.pop = initial population size
#n = number of time steps (e.g., years)
#r = drift parameter

#getNdi
getNdi<-function(init.pop, n, r){
  
  #error term (mean 0, var=0.01, equal to sd=0.0001)
  e=rnorm(n,mean=0, sd=0.1)
  
  #natural log of initial pop. size
  new.x=log(init.pop)
  
  #save as pop. size at initial time step
  Nout1=new.x
  
  Nout<-list()
  for(i in 1:(n-1)){
    #error term (mean 0, var=0.01, equal to sd=0.1)
    new.x<-r+new.x+e[i]
    Nout<-rbind(Nout, new.x)
  }
  Nout<-data.frame(logN=unlist(Nout))
  
  row.names(Nout)<-NULL
  dataSet<-data.frame(rbind(Nout1, Nout))
  return(dataSet)
}

###################################################################################################################################
#Simulate population data that is density-dependent (Eq. 3 from Pollard et al. 1987)
#Eq. 3. x[i+1]=r + Bx[i]+ e[i]
#where e is a random variate from a normal distribution with mean 0 and sd=0.1
#init.pop = initial population size
#n = number of time steps (e.g., years)
#r = drift parameter
#B = slope parameter, where B < 1

#get N at t+n function
getNdd<-function(init.pop, n, r,B){
  
  #error term (mean 0, var=0.01, equal to sd=0.0001)
  e=rnorm(n,mean=0, sd=0.1)
  
  #initial pop size
  
  new.x=log(init.pop)
  #save as N at initial time step
  Nout1=new.x
  
  Nout<-list()
  for(i in 1:(n-1)){
    #error term (mean 0, var=0.01, equal to sd=0.0001)
    new.x<-r+(B*new.x)+e[i]
    Nout<-rbind(Nout, new.x)
  }
  Nout<-data.frame(logN=unlist(Nout))
  
  row.names(Nout)<-NULL
  dataSet<-data.frame(rbind(Nout1, Nout))
  return(dataSet)
}

###################################################
#m1
m1func<-function(x){
  n=length(x)
  x.out1<-list()
  for(i in 1:(n-1)){
          x.temp1<-x[i]/(n-1)
          x.out1<-rbind(x.out1,x.temp1)
  }
  df.out1<-data.frame(x.out1=unlist(x.out1))
  m1=sum(df.out1$x.out1)
  return(m1)
}

m1=m1func(x)
###################################################
#m2
m2func<-function(x){
  n=length(x)
  x.out2<-list()
  for(i in 1:(n-1)){
    x.temp2<-x[i+1]/(n-1)
    x.out2<-rbind(x.out2,x.temp2)
  }
  df.out2<-data.frame(x.out2=unlist(x.out2))
  m2=sum(df.out2$x.out2)
  return(m2)
}

m2=m2func(x)
##############################################################################################################
#b
bfunc<-function(x){
  n=length(x)
  m1=m1func(x)
  m2=m2func(x)
  
  b.out1<-list()
  b.out2<-list()
  for(i in 1:(n-1)){
    temp1<-(x[i] -m1)*(x[i+1]-m2)
    temp2<-(x[i]-m1)^2
  
    b.out1<-rbind(b.out1, temp1)
    b.out2<-rbind(b.out2, temp2)
  
    }
  df.b.out1<-data.frame(b.out1=unlist(b.out1))
  df.b.out2<-data.frame(b.out2=unlist(b.out2))

  b=sum(df.b.out1$b.out1)/sum(df.b.out2$b.out2)
  return(b)
}

b=bfunc(x)
##############################################################################################################
#Test statistic T(2,3)

Tfunc=function(x){
  n=length(x)
  m1=m1func(x)
  m2=m2func(x)
  b=bfunc(x)
  
  temp1.out<-list()
  temp2.out<-list()
  temp3.out<-list()
  for(i in 1:(n-1)){
    temp1<-(x[i+1]-m2)^2 - b
    temp2<-(x[i]-m1)*(x[i+1]-m2)
    temp3<-((x[i+1]-x[i])^2) - (x[n]-x[1])^2/(n-1)

    temp1.out<-rbind(temp1.out, temp1)
    temp2.out<-rbind(temp2.out, temp2)
    temp3.out<-rbind(temp3.out, temp3)
  }

  Tout<-(sum(unlist(temp1.out))*sum(unlist(temp2.out)))/sum(unlist(temp3.out))
return(Tout)
}

#However, T(2,3), the test statistic equal to the likelihood ratio test between model 2 (density-independent) and 
#model 3 (density-dependent) = rdx (see Appendix in Pollard et al. 1987, where rdx is the correlation coef between di and xi).

#function to estimate rdx, the correlation coefficient between di and xi of observered data
rdxfunc<-function(x){
  n=length(x)
  di.out<-list()
  for(i in 1:(n-1)){
    d.temp<-x[i+1]-x[i]
    di.out<-rbind(di.out, d.temp)
  }
  df.di.out<-data.frame(di.out=unlist(di.out))
  new.df<-data.frame(x[2:n],df.di.out)
  colnames(new.df)<-c("xi","di")
  new.df$xi<-as.numeric(new.df$xi)
  new.df$di<-as.numeric(new.df$di)
  
  #get correlation coef
  rdx=cor(new.df$di,new.df$xi)
  return(rdx)
}


#Function to estimate rdx using randomization process
rdxfuncRand<-function(x){
  n=length(x)
  di.out<-list()
  for(i in 1:(n-1)){
    d.temp<-x[i+1]-x[i]
    di.out<-rbind(di.out, d.temp)
    }
  df.di.out<-data.frame(di.out=unlist(di.out))
  
  #randomize di order
  rand.order<-sample(seq(1:(n-1)),(n-1),replace=FALSE)
  new.di.out<-data.frame(rand.order=rand.order, di.out=df.di.out)
  new.di.out<-new.di.out[order(new.di.out$rand.order, decreasing=FALSE),]
  
  di.rand<-as.numeric(new.di.out$di.out)
  x.1=x[1]
  x.temp=x.1
  xi.out<-list()
  for(i in 1:(n-1)){
    x.temp<-x.temp+di.rand[i]
    xi.out<-rbind(xi.out, x.temp)
  }
  
  new.df<-data.frame(unlist(xi.out),new.di.out[,2])
  colnames(new.df)<-c("xi","di")
  new.df$xi<-as.numeric(new.df$xi)
  new.df$di<-as.numeric(new.df$di)

  #get correlation coef
  rdx=cor(new.df$di,new.df$xi)
  return(rdx)
}

##################################################################################################################################
#Function to simulate density-independent data loop to generate (n) datasets.
sim.data.di.func<-function(init.pop, nPops, n, r){
  sim.data.di<-list()
    for(j in 1:nPops){
      new.data.set<-getNdi(init.pop=init.pop, n=n, r=r)
      row.names(new.data.set)<-NULL
      df<-data.frame(dataSet=j, timeStep=seq(1:n), logN=new.data.set)
      sim.data.di<-rbind(sim.data.di, df)
  }
return(sim.data.di)
}

######################################################################################################################################
#Function to simulate density-dependent data loop to generate (n) datasets.
sim.data.dd.func<-function(init.pop, nPops, n, r, B){
  sim.data.dd<-list()
  for(j in 1:nPops){
    new.data.set<-getNdd(init.pop=init.pop, n=n, r=r, B=B)
    row.names(new.data.set)<-NULL
    df<-data.frame(dataSet=j, timeStep=seq(1:n), logN=new.data.set)
    sim.data.dd<-rbind(sim.data.dd, df)
    
  }
  return(sim.data.dd)
}
######################################################################################################################################
#Simulate "observed" data, and then perform randomization test, and plot for density-independent data.

#Step 1) generate population data sets to test (nPops=200)
new.di.data.1<-sim.data.di.func(init.pop=3, nPops=1, n=30, r=0.4)

#Step 2) compute rdx test statistic T=rdx
rdx.1<-rdxfunc(new.di.data.1$logN)

#Step 3) generate 200 simulated population data
sim.di.data<-sim.data.di.func(init.pop=3, nPops=200, n=30, r=0.4)

#Step 4) get rdx for each simulated population
rdx.out<-list()
for(i in 1:200){
  
  #subset data
  temp.data<-subset(sim.di.data, dataSet==i)
  
  #get rdx (correlation coef)
  rdx=rdxfuncRand(temp.data$logN)
  
  #save rdx output
  rdx.out<-rbind(rdx.out, rdx)
}
#Step 5) compare each to rdx.1 and tallyt frequency to test null hypothesis

#make dataframe with all rdx values (including first one)
rdx.all<-data.frame(rdx=c(rdx.1,unlist(rdx.out)))

#loop to compare for test
test.out<-list()
for(i in 1:200){
  
  result<-ifelse(rdx.all$rdx[i+1]<=rdx.all$rdx[1],1,0)
  test.out<-rbind(test.out, result)
  
}
#compute mean of test results.
final.result.di<-ifelse(mean(unlist(test.out))<=0.05,"Density-Dependent","Density-Independent")
final.result.di


#plot with ggplot
plot.di<-ggplot(data=sim.di.data, aes(x=timeStep, y=logN,group=dataSet))+
  #geom_point(color="black",size=1)+
  geom_path(color="gray60",alpha=0.3)+
  geom_path(data=new.di.data.1, aes(x=timeStep, y=logN),color="black",size=1)+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(color="black",fill=NA))+
  theme(panel.border = element_rect(color="black",fill=NA))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))+
  theme(plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=(c(1, (length(unique(sim.di.data$timeStep))/2), length(unique(sim.di.data$timeStep)) )))+
  labs(x=expression(paste("Time step (",italic("t"),")")), y=expression(ln~italic(N)))+
  ggtitle(paste("a)",final.result.di))
plot.di
plot.di.out<-plot.di+annotate("text", x=19, y=2.5, label="italic(x[i+1] == r + x[i] + e[i])",parse=TRUE,size=8)
plot.di.out
######################################################################################################################################
#Simulate "observed" data, and then perform randomization test, and plot for density-dependent data.

#Step 1) get population data to test
new.dd.data<-sim.data.dd.func(init.pop=3, nPops=1, n=30, r=0.4,B=0.8)
#new.dd.data.1<-sim.data.di.func(init.pop=11, nPops=1, n=10, r=0.4)

#Step 2) compute rdx test statistic T=rdx
rdx.1<-rdxfunc(new.dd.data$logN)

#Step 3) generate 200 simulated population data
sim.dd.data<-sim.data.dd.func(init.pop=3, nPops=200, n=30, r=0.4, B=0.8)

#Step 4) get rdx for each simulated population
rdx.out<-list()
for(i in 1:200){
  
  #subset data
  temp.data<-subset(sim.dd.data, dataSet==i)
  
  #get rdx (correlation coef)
  rdx=rdxfuncRand(temp.data$logN)
  
  #save rdx output
  rdx.out<-rbind(rdx.out, rdx)
}
#Step 5) compare each to rdx.1 and tallyt frequency to test null hypothesis

#make dataframe with all rdx values (including first one)
rdx.all<-data.frame(rdx=c(rdx.1,unlist(rdx.out)))

#loop to compare for test
test.out<-list()
for(i in 1:200){
  
  result<-ifelse((rdx.all$rdx[i+1]<=rdx.all$rdx[1]),1,0)
  test.out<-rbind(test.out, result)
  
}
#compute mean of test results.
final.result.dd<-ifelse(mean(unlist(test.out))<=0.05,"Density-Dependent","Density-Independent")
final.result.dd


#plot with ggplot
plot.dd<-ggplot(data=sim.dd.data, aes(x=timeStep, y=logN,group=dataSet))+
  #geom_point(color="black",size=1)+
  geom_path(color="gray60",alpha=0.3)+
  geom_path(data=new.dd.data, aes(x=timeStep, y=logN),color="black",size=1)+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(color="black",fill=NA))+
  theme(panel.border = element_rect(color="black",fill=NA))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))+
  theme(plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=(c(1, (length(unique(sim.dd.data$timeStep))/2), length(unique(sim.dd.data$timeStep)) )))+
  labs(x=expression(paste("Time step (",italic("t"),")")), y=expression(ln~italic(N)))+
  ggtitle(paste("b)",final.result.dd))
plot.dd
plot.dd.out<-plot.dd+annotate("text", x=19, y=1.2, label="italic(x[i+1] == r + Beta*x[i] + e[i])",parse=TRUE,size=8)
plot.dd.out

#combine plots using grid extra
plot.out<-grid.arrange(plot.di.out, plot.dd.out, ncol = 1, nrow = 2)

#save plot to Desktop
tiff(filename=paste(getwd(),"Desktop",'PollardTest_figure.tiff',sep="/"), width=600, height=600)
plot.out<-grid.arrange(plot.di.out, plot.dd.out, ncol = 1, nrow = 2)
dev.off()
