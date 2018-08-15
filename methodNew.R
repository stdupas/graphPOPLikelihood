library(markovchain)
library(MASS)


setMethod(
  f = "ncellA",
  signature = "RasterLayer",
  definition = function(object){
    length(na.omit(values(object)))
  }
)
setMethod(
  f = "ncellA",
  signature = "RasterStack",
  definition = function(object){
    ncellA(object[[1]])
  }
)
setMethod(
  f = "ncellA",
  signature = "RasterBrick",
  definition = function(object){
    ncellA(object[[1]])
  }
)

setMethod(
  f = "valuesA",
  signature = "RasterLayer",
  definition = function(object){
    #x=data.frame(variable=na.omit(values(object)))
    select <- !is.na(values(object))
    x=values(object)[select]
    names(x) <- which(select)
    #colnames(x)=names(object)
  x
  }
  )

setMethod(
  f = "valuesA",
  signature = "RasterBrick",  
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(object)
    rownames(x)=cellNumA(object)
    x
    }
)

setMethod(
  f = "valuesA",
  signature = "RasterStack",
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(x)
    rownames(x) <- cellNumA(object)
    x
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterLayer",
  definition = function(object){
    df=xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
    df
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterStack",
  definition = function(object){
    df =xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
    df
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterBrick",
  definition = function(object){
    df=xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
    df
  }
)



setMethod(
  f = "cellNumA",
  signature = "RasterLayer",
  definition = function(object){
    which(!is.na(values(object)))
  }
)

setMethod(
  f = "cellNumA",
  signature = "RasterStack",
  definition = function(object){
    cellNumA(object[[1]])
  }
)

setMethod(
  f = "cellNumA",
  signature = "RasterBrick",
  definition = function(object){
    cellNumA(object[[1]])
  }
)

    
    
# copier la fonction distanceMatrix  dans la method en  l'adaptant --> idem pour transitionMatrix        


setMethod(
  f="distanceMatrixA",
  signature="RasterLayer",
  definition = function(object) {
      coords = xyFromCellA(object)
      distance = as.matrix(dist(coords)) 
      dimnames(distance) <- list(which(!is.na(values(object))),which(!is.na(values(object))))
      distance
  }
)

setMethod(
  f="transitionMatrixA",
  signature=c("RasterLayer","prior"),
  definition = function(object1,object2){
    object2 <- sampleP(object2)
    transitionMatrixA(object1,object2)
  })

setMethod(
  f="transitionMatrixA",
  signature=c("character","prior"),
  definition = function(object1,object2){
    object2 <- sampleP(object2)
    transitionMatrixA(object1,object2)
  })

setMethod(
  f="transitionMatrixA",
  signature=c("character","parameters"),
  definition = function(object1,object2){
    switch(object2[[1]]$model,
           stepwise={
             
           })
    object2[["mutation_rate"]]
    transitionMatrixA(object1,object2)
  })


setMethod(
  f="transitionMatrixA",
  signature=c("RasterLayer","parameters"),
  definition = function(object1,object2,Option="demes"){
    if(Option=="demes"){
      X=valuesA(object1)
      K = ReactNorm(X,object2$K[[names(object1)]]$p,object2$K[[names(object1)]]$model)[,"Y"]
      r = ReactNorm(X,object2$R[[names(object1)]]$p,object2$R[[names(object1)]]$model)[,"Y"] 
      migration <- migrationMatrixA(object1,object2$dispersion[[names(object1)]]$model, object2$dispersion[[names(object1)]]$p[[1]])
      if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
      if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
      if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
      if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
      if ((length(r)>1)&(length(K)>1)) {transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
      transition = transition / rowSums(transition)  # Standardisation
      transition = transition * as.numeric(transition>1e-6) # removal of values below 1e-6
      transition = transition / rowSums(transition)  # Standardisation again
      transition
    } else {
      transition = migrationMatrixA(object1,object2$mutation_rate[[names(object1)]]$model[Option], object2$mutation_rate[[names(object1)]]$p[[Option]])
      names(transition)
      transition
    }
  } 
  )

setMethod(
  f="varnames",
  signature="parameters",
  definition=function(object){
    names(object$K)
  })

setMethod(
  f="varnames",
  signature="prior",
  definition=function(object){
    names(object$K)
  })

setMethod(
  f="sampleP",
  signature="prior",
  definition=function(Prior) {
  Result=list()
  for(parametreBio in names(Prior)) {
    for (variableEnvironnemental in names(Prior[[parametreBio]])) {
      Result[[parametreBio]] <- list() 
      Result[[parametreBio]][[variableEnvironnemental]]$p <- switch(Prior[[parametreBio]][[variableEnvironnemental]]$a$distribution,
               uniform=data.frame(variableEnvironnemental=runif(1,min=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
               fixed =data.frame(variableEnvironnemental=Prior[[parametreBio]][[variableEnvironnemental]]$a$p),
               normal=data.frame(variableEnvironnemental=rnorm(1,mean=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],sd=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
               loguniform={
                 out=list()
                 for (locus in names(Prior[[parametreBio]][[variableEnvironnemental]]$a$p)){
                   out[[locus]] <- exp(runif(1,
                                             min=log(Prior[[parametreBio]][[variableEnvironnemental]]$a$p[[locus]]["min"]),
                                             max=log(Prior[[parametreBio]][[variableEnvironnemental]]$a$p[[locus]]["max"])))
                   }  
                 out},
               uniform=data.frame(variableEnvironnemental= runif(1,min=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
               fixed =data.frame(variableEnvironnemental= Prior[[parametreBio]][[variableEnvironnemental]]$a$p),
               normal=data.frame(variableEnvironnemental =rnorm(1,mean=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],sd=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
               loguniform=data.frame(variableEnvironnemental= log(runif(1,min=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=Prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])))
               )
      Result[[parametreBio]][[variableEnvironnemental]]$model <-  Prior[[parametreBio]][[variableEnvironnemental]]$model
      names(Result[[parametreBio]][[variableEnvironnemental]]$model)=names(Result[[parametreBio]][[variableEnvironnemental]]$p )
      if(class(Result[[parametreBio]][[variableEnvironnemental]]$p)=="data.frame") 
      {colnames(Result[[parametreBio]][[variableEnvironnemental]]$p)=variableEnvironnemental
       rownames(Result[[parametreBio]][[variableEnvironnemental]]$p)=c("a")
       names(Result[[parametreBio]][[variableEnvironnemental]]$model)=variableEnvironnemental
      }
    }
  }
  #  names(Result) <- names(prior)
  return(new("parameters",Result))
}
)

setMethod(
  f="nodes",
  signature="listOfNodes",
  definition=function(object){
    unlist(lapply(object,function(sub) sub@nodeNo))
  }
)

setMethod(
  f="currentNodes",
  signature=c("listOfNodes","numeric"),
  definition=function(object,age){
    which((sapply(object,function(object) object@ancestorAge)>=(age+1E-15))&(sapply(object,function(object) object@age)<=age))
  }
)
setMethod(
  f="currentNodes",
  signature=c("graphicalListOfNodes","numeric"),
  definition=function(object,age){
    which((sapply(object,function(object) object@ancestorAge)>=(age+1E-15))&(sapply(object,function(object) object@age)<=age))
  }
)
setMethod(
  f="nodesByStates",
  signature=c("listOfNodes","numeric","character"),
  definition=function(object,age,Which){
    switch(Which,
           allByDemeAndAllele = {
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             lapply(NBS,function(x) as.numeric(row.names(x)))},
           allByDeme = {
             states=state(object,age,"Demes")
             lapply(split(states,states),function(x) as.integer(names(x)))},
           allAllele = {
             states=state(object,age,"Alleles")
             lapply(split(states,states),function(x) as.integer(names(x)))},
           notAloneByDemeAndAllele = {
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             whichNBS<- which(lapply(NBS,nrow)>1)
             lapply(whichNBS,function(x) as.integer(row.names(NBS[[x]])))},
           notAloneAndDemeAndAllele = {
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             whichNBS<- which(lapply(NBS,nrow)>1)
             lapply(whichNBS,function(x) NBS[[x]])},
           notAloneByDeme = {
             states=state(object,age,"Demes")
             NBS=split(states,states)
             whichNBS<- which(lapply(NBS,length)>1)
             lapply(whichNBS,function(x) as.integer(names(NBS[[x]])))},
           notAloneByAllele = {
             states=state(object,age,"Alleles")
             NBS=split(states,states)
             whichNBS<- which(lapply(NBS,length)>1)
             lapply(whichNBS,function(x) as.integer(names(NBS[[x]])))},
           statesWhereCoalesce ={
             states=state(object,age,"DemesAndAlleles")
             NBS=split(states,paste(states$statusDemes,states$statusAlleles,sep="."))
             whichNBS<- which(lapply(NBS,nrow)>1)
             lapply(whichNBS,function(x) NBS[[x]][1,])})
  }
  )    


setMethod(
  f="state",
  signature=c("listOfNodes","numeric","character"),
  definition=function(object,age=numeric(),type){
    if (length(age)==0) Nodes = names(coalescent) else Nodes = currentNodes(object,age)
    switch(type,
           DemesAndAlleles = {df=data.frame(
             statusDemes = unlist(lapply(Nodes,function(i) {
               demes <- object[[i]]@statusDemes[(object[[i]]@agesDemes<=age)]
               demes[length(demes)]})),
             statusAlleles = unlist(lapply(Nodes,function(i) {
               alleles <- object[[i]]@statusAlleles[(object[[i]]@agesAlleles<=age)]
               alleles[length(alleles)]})))
             df$statusDemes=as.character(df$statusDemes)
             df$statusAlleles=as.character(df$statusAlleles)
             df},
           Demes = as.character(unlist(lapply(Nodes,function(i) {
             demes <- object[[i]]@statusDemes[(object[[i]]@agesDemes<=age)]
             demes[length(demes)]}))),
           Alleles = as.character(unlist(lapply(Nodes,function(i) {
             alleles <- object[[i]]@statusAlleles[(object[[i]]@agesAlleles<=age)]
             alleles[length(alleles)]}))),
           aparitionDemes = unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@agesDemes[(object[[i]]@agesDemes<=age)]
             apparition[length(apparition)]})),
           aparitionAlleles = unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@agesAlleles[(object[[i]]@agesAlleles<=age)]
             apparition[length(apparition)]})),
           age =  {tmp=unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@age}))
             names(tmp)<-Nodes
             tmp},
           ancestorAge =  {tmp= unlist(lapply(Nodes,function(i){
             tmp <- object[[i]]@ancestorAge}))
             names(tmp)=Nodes
             tmp},
           nodeNo =  {tmp=unlist(lapply(Nodes,function(i){
             object[[i]]@nodeNo}))
             names(tmp)=Nodes
             tmp})
  }
)

setMethod(
  f="state",
  signature=c("listofnodEs","numeric","character"),
  definition=function(object,age=numeric(),type){
    if (length(age)==0) Nodes = names(coalescent) else Nodes = currentNodes(object,age)
    switch(type,
           age =  {tmp=unlist(lapply(Nodes,function(i){
             apparition <- object[[i]]@age}))
             names(tmp)<-Nodes
             tmp},
           ancestor =  {tmp= unlist(lapply(Nodes,function(i){
             tmp <- object[[i]]@ancestor}))
             names(tmp)=Nodes
             tmp},
           descendants =  {tmp= unlist(lapply(Nodes,function(i){
             tmp <- object[[i]]@descendants}))
             names(tmp)=Nodes
             tmp},
           nodeNo =  {tmp=unlist(lapply(Nodes,function(i){
             object[[i]]@nodeNo}))
             names(tmp)=Nodes
             tmp})
  }
)


setMethod("setState",
          signature=c("listOfNodes","integer","character","list"),
          definition=function(object,Nodes,attribut,newValues){
            switch(attribut,
                   ancestorAge={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@ancestorAge = newValues[[i]]}},
                   age={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@age = newValues[[i]]}},
                   statusDemes={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@statusDemes = append(object[[Nodes[i]]]@statusDemes,newValues[[i]])}},
                   ageDemes={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@ageDemes = append(object[[Nodes[i]]]@ageDemes,newValues[[i]])}},
                   statusAlleles={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@statusAlleles = append(object[[Nodes[i]]]@statusAlleles,newValues[[i]])}},
                   ageAlleles={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@ageAlleles = append(object[[Nodes[i]]]@ageAlleles,newValues[[i]])}},
                   nodeNo={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@nodeNo = newValues[[i]]
                     names(object)[Nodes[i]] <- newValues[[i]]}},
                   descendants={for (i in 1:length(Nodes)){
                     object[[Nodes[i]]]@descendants = newValues[[i]]}})
          object
          }
)
            


setMethod(
  f="currentState",
  signature=c("listOfNodes","character","integer"),
  definition=function(object,type,nodes){
    switch(type,
           allStatus = data.frame(statusDemes = unlist(lapply(nodes,function(i) coalescent[[i]]@statusDemes[length(coalescent[[i]]@statusDemes)])),
                                  statusAlleles = unlist(lapply(nodes,function(i) coalescent[[i]]@statusAlleles[length(coalescent[[i]]@statusAlleles)]))),
           statusDemes = unlist(lapply(nodes,function(i) coalescent[[i]]@statusDemes[length(coalescent[[i]]@statusDemes)])),
           statusAlleles = unlist(lapply(nodes,function(i) coalescent[[i]]@statusAlleles[length(coalescent[[i]]@statusAlleles)])),
           agesDemes = unlist(lapply(nodes,function(i) coalescent[[i]]@agesDemes[length(coalescent[[i]]@agesDemes)])),
           agesAlleles = unlist(lapply(nodes,function(i) coalescent[[i]]@agesAlleles[length(coalescent[[i]]@agesAlleles)])))
  }
)




setMethod(
  f="migrationMatrixA",
  signature="RasterLayer",
  definition = function(object,shapeDisp,pDisp){
      distMat<-distanceMatrixA(object)
      Ndim = 1+all(ncellA(object)!=dim(object)[1:2])
      migration = apply(distMat, c(1,2), 
                        function(x)(switch(shapeDisp,
                                           # 1: alphaDisp   2: betaDisp ; note: esperance = 1/alphaDisp
                                           fat_tail1 = 1/(1+x^pDisp[2]/pDisp[1]), # Molainen et al 2004
                                           # 1: sigmaDisp
                                           gaussian = (dnorm(x, mean = 0, sd = pDisp[1], log = FALSE)),
                                           # 1: sigma Disp
                                           exponential = (dexp(x, rate = 1/pDisp[1], log = FALSE)),
                                           # 1: sigmaDisp
                                           contiguous = (x==0)*(1-pDisp[1])+((x>0)-(x>1.4*res(object)[1]))*(pDisp[[1]]/(2*Ndim)),
                                           # 1: sigmaDisp
                                           stepwise = (x==0)*(1-pDisp[1])+((x>0)-(x>1.4*res(object)[1]))*(pDisp[1]/(Ndim)),
                                           # 1: sigmaDisp
                                           contiguous8 = (x==0)*(1-pDisp[[1]])+((x>0)-(x>2*res(object)[1]))*(pDisp[[1]]/(4*Ndim)),
                                           #island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]/(nCell-1)),
                                           island = (x==0)*(1-pDisp[[1]])+(x>0)*(pDisp[[1]]),
                                           #1: sigmaDisp    2: gammaDisp
                                           fat_tail2 = x^pDisp[2]*exp(-2*x/(pDisp[[1]]^0.5)),
                                           #1: sigmaDisp, 2: 
                                           contiguous_long_dist_mixt = pDisp["plongdist"]/ncellA(object)+(x==0)*(1-pDisp["pcontiguous"]-pDisp["plongdist"])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp["pcontiguous"]/2),
                                         gaussian_long_dist_mixt = pDisp[2]/ncellA(object) + (dnorm(x, mean = 0, sd = pdisp[1], log = FALSE))
                        )))
      migration
    }
    )

setMethod(
  f="absorbingTransitionA",
  signature="transitionModel",
  definition = function(object){
    #which= 'dem' or 'allel'
    absorbingTransitionA(object[paste(Which,"icTransition",sep="")],object["Ne"])
  })

    
setMethod(
  f="absorbingTransitionA",
  signature="matrix",
  definition = function(object,N){
    N[N<1]=1
    Ndeme <- dim(object)[1]
    # States of pairs of genes
    # N Homodemic not coalesced states
    # Heterodemic states (ij)!=(kl)
    Nhetero <- Ndeme*(Ndeme-1)/2
    # Ndeme Coalesced demic state, can be merged to 1 if we don't mind the deme
    #
    # First: 
    # Calculate the transient states (not coalesced) transition matrix
    # for transition among Q among N*(N+1)/2 not coalesced states 
    # It is composed of a submatrixes of transition between not coalesced heterodemic and
    # homodemic states
    #
    #    /          \
    #    |HeHe  HeHo|
    # Q= |HoHe  HoHo|
    #    \          /
    # where He are heterodemic states, and Ho are homodemic states
    # In submatrix HeHe:
    # lines ij are ordered for (i in 2:Ndeme) for (j in 1:(i-1))
    # columns kl are ordered for (k in 2:Ndeme) for (l in 1:(k-1))
    # In submatrix HoHe the lines are from 1 to Ndeme
    Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
    # contains the line number or column number 
    # in transient Q matrix for each pair of deme {i,j} 
    # presented as in the transition matrix
    # 
    HeHe <- matrix(NA,Nhetero,Nhetero)
    ij=0;kl=0
    Check=TRUE
    for (i in 2:Ndeme){
      for(j in 1:(i-1)) {
        ij=ij+1
        Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
        # is in Q matrix lines or columns
        #      i_j_Q_lines[ij,]<-c(i,j)
        kl=0
        for (k in 2:Ndeme){
          for (l in 1:(k-1)){
            kl=kl+1
            HeHe[ij,kl] <- object[i,k]*object[j,l]
            #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
          }
        }
      }
    }
    
    # then transient matrix 
    HoHe <- matrix(NA,Ndeme,Nhetero)
    kl=0
    Check=TRUE
    for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
      Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
      # in the lines of Q matrix
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){ # only heterodemic targets
          kl=kl+1
          HoHe[i,kl] <- object[i,k]*object[i,l]    # i=j (homodemic sources)
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
      kl=0
    }
    HeHo <- matrix(NA,Nhetero,Ndeme)
    HeCo <- matrix(NA,Nhetero,Ndeme)
    ij=0
    for (i in 2:Ndeme){
      for (j in 1:(i-1)){ # only heterodemic sources
        ij=ij+1
        for(k in 1:Ndeme){ # only homodemic targets
          HeHo[ij,k] <- object[i,k]*object[j,k]*(1-1/(2*N[k]))
          HeCo[ij,k] <- object[i,k]*object[j,k]*1/(2*N[k])
          # homodemic targets that have not coalesced
        }
      }
    }
    HoHo <- object*object*matrix(1-1/(2*N),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
    HoCo <- object*object*matrix(1/(2*N),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
    CoHo <- matrix(0,nrow=length(N),ncol=length(N))
    CoHe <- matrix(0,nrow=length(N),ncol=nrow(HeHe))
    CoCo <- object
    Q <- cbind(rbind(HeHe,HoHe,CoHe),rbind(HeHo,HoHo,CoHo),rbind(HeCo,HoCo,CoCo))
    dimnames(Qline) <- list(names(Ne),names(Ne))
    tmp <- rev(c(Qline))
    names(tmp) <- rev(paste(rep(rownames(Qline),nrow(Qline)), rep(colnames(Qline),each=nrow(Qline)),sep="_"))
    tmp <- names(tmp[order(tmp)][!duplicated(tmp[order(tmp)])])
    tmp <- append(tmp,tmp[(length(tmp)-length(N)+1):length(tmp)])
    dimnames(Q) <- list(tmp,tmp)
    result=list(Q=Q,Qline=Qline)
    result
    #
    #
    #    /                \
    #    |HeHe  HeHo  HeCo |
    # Q= |HoHe  HoHo  HoCo |
    #    |CoHe  CoHo  CoCo |  
    #    \                /
    #
    #
    
  } 
  )

setAs("demeTransition", "matrix",
       function(from , to){
       })


setMethod(f="plot", 
          signature="demeTransition",
          definition = function(x, y , ...){
            
})

setMethod(
  f="absorbingTransitionA",
  signature="demeTransition",
  definition = function(object){
    N=valuesA(object)
    N[N<1]=1
    Ndeme <- dim(object)[1]
    # States of pairs of genes
    # N Homodemic not coalesced states
    # Heterodemic states (ij)!=(kl)
    Nhetero <- Ndeme*(Ndeme-1)/2
    # Ndeme Coalesced demic state, can be merged to 1 if we don't mind the deme
    #
    # First: 
    # Calculate the transient states (not coalesced) transition matrix
    # for transition among Q among N*(N+1)/2 not coalesced states 
    # It is composed of a submatrixes of transition between not coalesced heterodemic and
    # homodemic states
    #
    #    /          \
    #    |HeHe  HeHo|
    # Q= |HoHe  HoHo|
    #    \          /
    # where He are heterodemic states, and Ho are homodemic states
    # In submatrix HeHe:
    # lines ij are ordered for (i in 2:Ndeme) for (j in 1:(i-1))
    # columns kl are ordered for (k in 2:Ndeme) for (l in 1:(k-1))
    # In submatrix HoHe the lines are from 1 to Ndeme
    Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
    # contains the line number or column number 
    # in transient Q matrix for each pair of deme {i,j} 
    # presented as in the transition matrix
    # 
    HeHe <- matrix(NA,Nhetero,Nhetero)
    ij=0;kl=0
    Check=TRUE
    for (i in 2:Ndeme){
      for(j in 1:(i-1)) {
        ij=ij+1
        Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
        # is in Q matrix lines or columns
        #      i_j_Q_lines[ij,]<-c(i,j)
        kl=0
        for (k in 2:Ndeme){
          for (l in 1:(k-1)){
            kl=kl+1
            HeHe[ij,kl] <- object[i,k]*object[j,l]
            #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
          }
        }
      }
    }
    
    # then transient matrix 
    HoHe <- matrix(NA,Ndeme,Nhetero)
    kl=0
    Check=TRUE
    for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
      Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
      # in the lines of Q matrix
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){ # only heterodemic targets
          kl=kl+1
          HoHe[i,kl] <- object[i,k]*object[i,l]    # i=j (homodemic sources)
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
      kl=0
    }
    HeHo <- matrix(NA,Nhetero,Ndeme)
    ij=0
    for (i in 2:Ndeme){
      for (j in 1:(i-1)){ # only heterodemic sources
        ij=ij+1
        for(k in 1:Ndeme){ # only homodemic targets
          HeHo[ij,k] <- object[i,k]*object[j,k]*(1-1/(2*N[k,]))
          # homodemic targets that have not coalesced
        }
      }
    }
    HoHo <- object*object*matrix(1-1/(2*N[,1]),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
    Q <- cbind(rbind(HeHe,HoHe),rbind(HeHo,HoHo))
    result=list(Q=Q,Qline=Qline)
    result
  } 
  
)

setMethod(
  f="plotgenealogy",
  signature="listOfNodes",
  definition = function(object,tipcols=NA)    {
    if(is.na(tipcols)) { tipcols=1}
    plot(coalescent_2_phylog(object),direction="downward",tip.color=tipcols)
    }
)

setMethod(
  f="plotLandG",
  signature="LandGenetGenealogy",
  definition=function(object,rasK=NULL) {
    par(mfrow=c(1,2))
    plot(object)
    tipcells <- genotypes[,2][as.numeric(coalescent_2_phylog(object,direction="downward",tip.color=tip.colors)$tip.label)]
    tipcols = rainbow(ncell(rasK))[tipcells]
    plotgenealogy(object@genealogy,tip.color=tipcols)
    legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
    
    
  }
)


setMethod(
  f="simul_coalescent",
  signature="transitionModel",
  definition = function(transitionMod,Option="allSimul")#transitionList,geneticData)
  {
    if (Option=="demeSimul"){
      Ne <- round(transitionMod@Ne);Ne[Ne==0]<-1 
      coalescent <- list()
      for (i in 1:length(transitionMod@demicStatusOfStartingIndividuals)){
        coalescent[[i]] <- new("Node",nodeNo=i,descendant=integer(),ancestor=integer(),
                               new("branchTransition",age=0,ancestorAge=Inf,
                                   statusDemes=transitionMod@demicStatusOfStartingIndividuals[i],
                                   agesDemes=0,
                                   statusAlleles=transitionMod@allelicStatusOfStartingIndividuals[i],
                                   agesAlleles=0))
        names(coalescent)[i]=i
      }
      numberOfNodes <- length(coalescent)
      coalescent= new("listOfNodes",coalescent)
      Age=0
      single_coalescence_events=0
      single_and_multiple_coalescence_events=0
      notCoalesced <- names(which(state(object = coalescent,age = Age,type = "ancestorAge")==Inf))
      #
      # Coalescing
      #
      while(length(notCoalesced)>1){
        status_of_nodes <- state(object = coalescent,age = Age,type = "Demes")
        parent_status_of_nodes <- status_of_nodes
        for (node in notCoalesced)#node = nodes[1];node = nodes[2];node = nodes[3]# parent_status_of_nodes
        {
          # migrations and mutations
          parent_status_of_nodes[node] = sample(transitionMod@demesNames,size=1,prob=c(transitionMod@demicTransition[as.character(status_of_nodes[node,"statusDemes"]),]))
        }
        changes <- status_of_nodes!=parent_status_of_nodes;changes
        # we go to the previous generation: parents become current generation
        Age=Age+1
        for (node in as.character(names(which(changes,arr.ind=TRUE)))){
          coalescent[[as.character(node)]]@agesDemes=
            append(coalescent[[as.character(node)]]@agesDemes,Age)
          coalescent[[as.character(node)]]@statusDemes=
            append(coalescent[[as.character(node)]]@statusDemes,as.vector(parent_status_of_nodes[node,"statusDemes"]))
        }
        #        nodesByStates(object = coalescent,age = Age, Which = "allByDemeAndAllele")
        nodesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "notAloneByDeme")
        statesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "statesWhereCoalesce")
        nodesDfStatesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "notAloneAndDemeAndAllele")
        for (States in names(nodesDfStatesThatCanCoalesce)){
          # get get the Deme of the nodes that can coalesce in the list
          current = unlist(strsplit(States,"\\."))# current[1] = deme, current[2]=allele
          nbNodesThatCanCoalesceInTheDeme=length(nodesThatCanCoalesce[[1]])
          smp = sample(Ne[current[1]],nbNodesThatCanCoalesceInTheDeme,replace=TRUE)
          parentoffspringmatrix <- matrix(smp,nrow=nbNodesThatCanCoalesceInTheDeme,ncol=Ne[current[1]])==matrix(1:Ne[current[1]],nrow=nbNodesThatCanCoalesceInTheDeme,ncol=Ne[current[1]],byrow=TRUE)
          #        colnames(parentoffspringmatrix) <- nodes_remaining_in_the_cell
          # the columns  column of parentoffspringmatrix that have more than one TRUE
          # identifies the individuals that coalesce with the lines of the TRUEs
          if (any(colSums(parentoffspringmatrix)>1) )
          {
            rownames(parentoffspringmatrix) <- nodesThatCanCoalesce[[1]]
            #  which(colSums(parentoffspringmatrix)>1)) gives the column names 
            #  of parentoffspringmatrix that have coalescence information
            for (multiple in which(colSums(parentoffspringmatrix)>1)) # multiple<-which(colSums(parentoffspringmatrix)>1)[1]
            {
              # there is coalescence 
              single_coalescence_events = single_coalescence_events +1
              # which(parentoffspringmatrix[,multiple]) identifies which node in the column coalesce
              nodes_that_coalesce = names(which(parentoffspringmatrix[,multiple]))
              # attibutes new node number to the ancestor, adds this to the nodes vector, removes the nodes that coalesced from the node vector
              numberOfNodes=numberOfNodes+as.integer(1)
              new_node <- as.integer(numberOfNodes)
              # updating the data frame parent_status_of_nodes (adding the cell number of the new node and removing the nodes that disapeared)
              parent_status_of_nodes[new_node,] <- current
              parent_status_of_nodes <- parent_status_of_nodes[!(rownames(parent_status_of_nodes)%in%nodes_that_coalesce),]
              # adds the event to the list coalescent: time, which node coalesced, and the number of the new node
              coalescent[[new_node]] <- new("Node",nodeNo=as.integer(new_node),
                                            descendant=as.integer(nodes_that_coalesce),
                                            ancestor=integer(),
                                            age=Age,ancestorAge=Inf,
                                            statusDemes=current[1],agesDemes=Age,
                                            statusAlleles=current[2],agesAlleles=Age)
              names(coalescent)[new_node] <- new_node
              # adds information on the ancestor for coalescing nodes
              for (Nod in as.integer(nodes_that_coalesce)) {
                coalescent[[Nod]]@ancestor <- as.integer(new_node)
                coalescent[[Nod]]@ancestorAge <- Age
              }
              # updates the number of coalescent events 
              single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
              notCoalesced <- as.integer(names(which(state(object = coalescent,age = Age,type = "ancestorAge")==Inf)))
            }
          }
        }
        notCoalesced <- names(which(state(object = coalescent,age = Age,type = "ancestorAge")==Inf))
        status_of_nodes <- parent_status_of_nodes
      }
      coalescent
    }
    if (Option=="allSimul"){
      # transitionList =  list of transition matrix 
      #                   sublist demes contains list of demic transitions
      #                   sublist alleles contains list of allelic transitions for each allele
      # Ne = a data.frame with number of individuals in each deme
      # demes = all the possible deme status (attibutted cells in the raster lanscape of population sizes)
      # alleles
      # demeStatus= deme status of the nodes
      # alleleStatus = allele status of the nodes
      Ne <- round(transitionMod@Ne);Ne[Ne==0]<-1 
      coalescent <- list()
      for (i in 1:length(transitionMod@demicStatusOfStartingIndividuals)){
        coalescent[[i]] <- new("Node",nodeNo=i,descendant=integer(),ancestor=integer(),
                               new("branchTransition",age=0,ancestorAge=Inf,
                                   statusDemes=transitionMod@demicStatusOfStartingIndividuals[i],
                                   agesDemes=0,
                                   statusAlleles=transitionMod@allelicStatusOfStartingIndividuals[i],
                                   agesAlleles=0))
        names(coalescent)[i]=i
      }
      numberOfNodes <- length(coalescent)
      coalescent= new("listOfNodes",coalescent)
      Age=0
      single_coalescence_events=0
      single_and_multiple_coalescence_events=0
      notCoalesced <- names(which(state(object = coalescent,age = Age,type = "ancestorAge")==Inf))
      #
      # Coalescing
      #
      while(length(notCoalesced)>1){
        status_of_nodes <- state(object = coalescent,age = Age,type = "DemesAndAlleles")
        parent_status_of_nodes <- status_of_nodes
        for (node in notCoalesced)#node = nodes[1];node = nodes[2];node = nodes[3]# parent_status_of_nodes
        {
          # migrations and mutations
          parent_status_of_nodes[node,] = as.vector(c(sample(transitionMod@demesNames,size=1,prob=c(transitionMod@demicTransition[as.character(status_of_nodes[node,"statusDemes"]),])),
                                                      sample(transitionMod@allelesNames,size=1,prob=c(transitionMod@allelicTransition[as.character(status_of_nodes[node,"statusAlleles"]),]))))
        }
        changes <- status_of_nodes!=parent_status_of_nodes;changes
        # we go to the previous generation: parents become current generation
        Age=Age+1
        for (node in as.character(names(which(changes[,"statusDemes"],arr.ind=TRUE)))){coalescent[[as.character(node)]]@agesDemes=
                                                                                         append(coalescent[[as.character(node)]]@agesDemes,Age)
                                                                                       coalescent[[as.character(node)]]@statusDemes=
                                                                                         append(coalescent[[as.character(node)]]@statusDemes,as.vector(parent_status_of_nodes[node,"statusDemes"]))
        }
        for (node in as.character(names(which(changes[,"statusAlleles"],arr.ind=TRUE)))){coalescent[[as.character(node)]]@agesAlleles=
                                                                                           append(coalescent[[as.character(node)]]@agesAlleles,Age)
                                                                                         coalescent[[as.character(node)]]@statusAlleles=
                                                                                           append(coalescent[[as.character(node)]]@statusAlleles,as.vector(parent_status_of_nodes[node,"statusAlleles"]))
        }
        #        nodesByStates(object = coalescent,age = Age, Which = "allByDemeAndAllele")
        nodesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "notAloneByDemeAndAllele")
        statesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "statesWhereCoalesce")
        nodesDfStatesThatCanCoalesce <- nodesByStates(object = coalescent,age = Age, Which = "notAloneAndDemeAndAllele")
        for (States in names(nodesDfStatesThatCanCoalesce)){
          # get get the Deme of the nodes that can coalesce in the list
          current = unlist(strsplit(States,"\\."))# current[1] = deme, current[2]=allele
          nbNodesThatCanCoalesceInTheDeme=length(nodesThatCanCoalesce[[1]])
          smp = sample(Ne[current[1]],nbNodesThatCanCoalesceInTheDeme,replace=TRUE)
          parentoffspringmatrix <- matrix(smp,nrow=nbNodesThatCanCoalesceInTheDeme,ncol=Ne[current[1]])==matrix(1:Ne[current[1]],nrow=nbNodesThatCanCoalesceInTheDeme,ncol=Ne[current[1]],byrow=TRUE)
          #        colnames(parentoffspringmatrix) <- nodes_remaining_in_the_cell
          # the columns  column of parentoffspringmatrix that have more than one TRUE
          # identifies the individuals that coalesce with the lines of the TRUEs
          if (any(colSums(parentoffspringmatrix)>1) )
          {
            rownames(parentoffspringmatrix) <- nodesThatCanCoalesce[[1]]
            #  which(colSums(parentoffspringmatrix)>1)) gives the column names 
            #  of parentoffspringmatrix that have coalescence information
            for (multiple in which(colSums(parentoffspringmatrix)>1)) # multiple<-which(colSums(parentoffspringmatrix)>1)[1]
            {
              # there is coalescence 
              single_coalescence_events = single_coalescence_events +1
              # which(parentoffspringmatrix[,multiple]) identifies which node in the column coalesce
              nodes_that_coalesce = names(which(parentoffspringmatrix[,multiple]))
              # attibutes new node number to the ancestor, adds this to the nodes vector, removes the nodes that coalesced from the node vector
              numberOfNodes=numberOfNodes+as.integer(1)
              new_node <- as.integer(numberOfNodes)
              # updating the data frame parent_status_of_nodes (adding the cell number of the new node and removing the nodes that disapeared)
              parent_status_of_nodes[new_node,] <- current
              parent_status_of_nodes <- parent_status_of_nodes[!(rownames(parent_status_of_nodes)%in%nodes_that_coalesce),]
              # adds the event to the list coalescent: time, which node coalesced, and the number of the new node
              coalescent[[new_node]] <- new("Node",nodeNo=as.integer(new_node),
                                            descendant=as.integer(nodes_that_coalesce),
                                            ancestor=integer(),
                                            age=Age,ancestorAge=Inf,
                                            statusDemes=current[1],agesDemes=Age,
                                            statusAlleles=current[2],agesAlleles=Age)
              names(coalescent)[new_node] <- new_node
              # adds information on the ancestor for coalescing nodes
              for (Nod in as.integer(nodes_that_coalesce)) {
                coalescent[[Nod]]@ancestor <- as.integer(new_node)
                coalescent[[Nod]]@ancestorAge <- Age
              }
              # updates the number of coalescent events 
              single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
              notCoalesced <- as.integer(names(which(state(object = coalescent,age = Age,type = "ancestorAge")==Inf)))
            }
          }
        }
        notCoalesced <- names(which(state(object = coalescent,age = Age,type = "ancestorAge")==Inf))
        status_of_nodes <- parent_status_of_nodes
      }
      coalescent
    }
  }
)







setMethod ("myPlot", 
           signature=c("listOfNodes","transitionModel"),
           definition = function(coalescent,transitionMod){
             alleles <- 1:length(transitionMod@allelesNames);names(alleles)<-transitionMod@allelesNames
             demes <- 1:length(transitionMod@demesNames);names(demes)<-transitionMod@demesNames
             maxAge<-max(state(coalescent,Inf,"age"))
             tips <- as.character(rev(tipsInOrder(coalescent)))
             ageNodes <- state(coalescent,numeric(),"age")
             ageInternalNodes <- ageNodes[ageNodes>0]
             internalNodesFromYoungest <- names(ageInternalNodes[order(ageInternalNodes)])
             Nodes <- append(tips,internalNodesFromYoungest)
             Y=unlist(lapply(Nodes,function(i){coalescent[[i]]@age/maxAge}))
             names(Y)<-Nodes
             X=.8*(1:length(tips)-.5)/length(tips)
             names(X)<-tips
             for (node in internalNodesFromYoungest){
               X[node] <- mean(unlist(lapply(as.character(coalescent[[node]]@descendant),function(i){X[i]})))
             }
             names(X)<-Nodes
             dev.off()
             plot.new()
             for (Node in Nodes){#Node="3"
               numberOfDemesInTheBranch <- length(coalescent[[Node]]@statusDemes)
               ancestor <- coalescent[[Node]]@ancestor
               if (length(ancestor)==0) ancestor=length(Nodes)
               yDemBranch <- c(coalescent[[Node]]@agesDemes/maxAge,rep(Y[ancestor],2))
               xDemBranch <- c(rep(X[Node],1+length(coalescent[[Node]]@agesDemes)),X[ancestor])
               segments(x0=xDemBranch[-length(xDemBranch)], 
                        y0=yDemBranch[-length(yDemBranch)], 
                        x1 =xDemBranch[-1], 
                        y1 =yDemBranch[-1],
                        col = rainbow(length(transitionMod@demesNames))[c(demes[coalescent[[Node]]@statusDemes],
                                                                          demes[coalescent[[Node]]@statusDemes[length(coalescent[[Node]]@statusDemes)]])],
                        lwd=4)
             }
             text(x=rep(.88,length(transitionMod@demesNames)),y=.9*(1:length(transitionMod@demesNames))/length(transitionMod@demesNames),labels=paste("d",transitionMod@demesNames,sep=""),col=rainbow(length(transitionMod@demesNames)),add=TRUE,cex=.8)
             text(x=.88,y=1,labels="dem",cex=.8)
             numberOfAlleles=length(transitionMod@allelesNames)
             for (Node in as.character(1:length(coalescent))){#Node="3"
               numberOfAllelesInTheBranch <- length(coalescent[[Node]]@statusAlleles)
               yDemBranch <- c(coalescent[[Node]]@agesAlleles/maxAge,rep(Y[coalescent[[Node]]@ancestor],2))
               if(length(yDemBranch)<length(xDemBranch)) yDemBranch[length(yDemBranch)+1]<-yDemBranch[length(yDemBranch)]
               xDemBranch <- c(rep(X[Node],1+length(coalescent[[Node]]@agesAlleles)),X[coalescent[[Node]]@ancestor])
               segments(x0=xDemBranch[-length(xDemBranch)]+.02, 
                        y0=yDemBranch[-length(xDemBranch)], 
                        x1 =xDemBranch[-1]+0.02, 
                        y1 =yDemBranch[-1],
                        col = rainbow(length(transitionMod@allelesNames))[alleles[coalescent[[Node]]@statusAlleles]],
                        lwd=4)
             }
             segments(x0=.1,x1=.1,y0=.9,y1=.9-round(maxAge/10)/maxAge,lwd=4)
             text(x=.14,y=0.85,paste(round(maxAge/10)),cex=.8)
             text(x=0.1,y=1,maxAge,cex=.8)
             text(x=rep(1,length(transitionMod@allelesNames)),y=.9*(1:length(transitionMod@allelesNames))/length(transitionMod@allelesNames),labels=paste("a",transitionMod@allelesNames,sep=""),col=rainbow(length(transitionMod@allelesNames)),add=TRUE,cex=.8)
             text(x=1,y=1,labels="allel",cex=.8)
           }
           )

setMethod ("coalescent_2_newick",
           signature="listOfNodes",
           definition = function(coalescent){
             tree=paste(" ",coalescent[[length(coalescent)]]@nodeNo," ",sep="")
             ageNodes <- state(coalescent,numeric(),"age")
             ageInternalNodes <- ageNodes[ageNodes>0]
             internalNodesFromOldestToYoungest <- rev(names(ageInternalNodes[order(ageInternalNodes)]))
             for (i in internalNodesFromOldestToYoungest)
             {
               Time = coalescent[[i]]@age
               coalesc <- as.character(coalescent[[i]]@descendant)
               br_length <-  unlist(lapply(as.character(coalescent[[i]]@descendant),function(node){coalescent[[i]]@age-coalescent[[node]]@age}))
               tree <- str_replace(tree,paste(" ",as.character(coalescent[[i]]@nodeNo)," ",sep=""),paste(" ( ",paste(" ",coalesc," :",br_length,collapse=" ,",sep=""),") ",sep=""))
             }
             tree <- gsub(" ","",paste(tree,";",sep=""))
             tree
           }
)

setMethod ("coalescent_2_newick",
           signature="listofnodEs",
           definition = function(coalescent){
             tree=paste(" ",coalescent[[length(coalescent)]]@nodeNo," ",sep="")
             ageNodes <- state(coalescent,numeric(),"age")
             ageInternalNodes <- ageNodes[ageNodes>0]
             internalNodesFromOldestToYoungest <- rev(names(ageInternalNodes[order(ageInternalNodes)]))
             for (i in internalNodesFromOldestToYoungest)
             {
               Time = coalescent[[i]]@age
               coalesc <- as.character(coalescent[[i]]@descendant)
               br_length <-  unlist(lapply(as.character(coalescent[[i]]@descendant),function(node){coalescent[[i]]@age-coalescent[[node]]@age}))
               tree <- str_replace(tree,paste(" ",as.character(coalescent[[i]]@nodeNo)," ",sep=""),paste(" ( ",paste(" ",coalesc," :",br_length,collapse=" ,",sep=""),") ",sep=""))
             }
             tree <- gsub(" ","",paste(tree,";",sep=""))
             tree
           }
)

setMethod ("coalescent_2_phylog",
           signature="listOfNodes",
           definition = function(coalescent){
             read.tree(text=coalescent_2_newick(coalescent))
           })


setMethod ("tipsInOrder",
           signature="character",
           definition = function(tree){
             Nodes <- integer();i=0
             positions <- which(strsplit(tree, "")[[1]]==":")
             for (position in positions){
               options(warn=-1)
               beforeTwoPoints <- as.numeric(substr(tree,position-1,position-1))
               if (is.numeric(beforeTwoPoints)&!is.na(beforeTwoPoints)){
                 i=i+1
                 Nodes[i] <- beforeTwoPoints
                 position <- position-1
                 beforeBefore <- as.numeric(substr(tree,position-1,position-1))
                 while(is.numeric(beforeBefore)&!is.na(beforeBefore)){
                   Nodes[i] <- paste(Nodes[i],beforeBefore)
                   position <- position-1
                   beforeBefore <- as.numeric(substr(tree,position-1,position-1))}
               }
             }
             Nodes
           }
)


setMethod ("tipsInOrder",
           signature="listOfNodes",
           definition = function(tree){
             tree <- coalescent_2_newick(tree)
             tipsInOrder(tree)
             }
)


setMethod("modifyBranch", 
          signature=c("listOfNodes","transitionModel"),
          definition = function(coalescent,transitionMod){
            node=sample(length(coalescent), 1)
            coalescent[[node]]@ageDemes=coalescent[[node]]@ageDemes[1]
            coalescent[[node]]@ageAlleles=coalescent[[node]]@ageAlleles[1]
            deme= coalescent@statusDemes=coalescent@statusDemes[1]
            allele=coalescent@statusAlleles=coalescent@statusAlleles[1]
            notCoalesced=TRUE
            age=coalescent[[node]]@ageAlleles
            while (notCoalesced){
              age=age+1
              newDeme=transitionMod@demeNames[sample(transitionMod@demeTransition[deme,],)]
              if (deme != newDeme) {coalescent@statusDeme=append(coalescent@statusDeme,newDeme)
                                    deme=newDeme
              }
              newAllele=transitionMod@allelesNames[sample(transitionMod@alleleTransition[allele,],)]
              if (allele != newAllele) {coalescent@statusAllele=append(coalescent@statusAllele,newAllele)
                                        allele =newAllele
                                        coalescent@ageAllele=append(coalescent@ageAllele,newAllele)
              }
              potentialCoalescent<-state(coalescent, age, "Node")[ state(coalescent, age, "allele")==allele]
            }
actualCoalescent = sample (c(potentialCoalescent, NA), c(rep(length(potentialCoalescent, 1/transitionMod@Ne[deme]), 1-length(potentialCoalescent) / transitionMod@Ne[deme])))
if (!is.na(actualCoalescent)){
  coalescent[[node]]@ancestor=coalescent[[actualCoalescent]]@ancestor
  coalescent[[coalescent[[node]]@ancestor]]@descendants[which(coalescent[[coalescent[[node]]@ancestor]]@descendants==actualCoalescent)]=node
  
  }})

transitionModel <- function(demicTransition,allelicTransitionList,demeNames,alleleNamesList,alleleTipsDataFrame,demeTips,Ne){
  if(class(Ne)=="RasterLayer") Ne <- round(valuesA(Ne))
  Ne[Ne==0] <- 1
  if (dim(demicTransition)[1]!=dim(demicTransition)[2]) stop("demicTransition not square")
  for (locus  in names(alleleNamesList)) {
    if (dim(allelicTransitionList[[locus]]))[1]!=dim(allelicTransitionList[[locus]]))[2]) stop(paste("locus", locus,"allelicTransition not square"))
    if (is.null(colnames(allelicTransitionList[[locus]])))) colnames(allelicTransitionList[[locus]]))=alleleNamesList[[locus]]
    if (is.null(rownames(allelicTransitionList[[locus]])))) rownames(allelicTransitionList[[locus]]))=alleleNamesList[[locus]]
    if (length(alleleNamesList[[locus]])!=dim(allelicTransitionList[[locus]]))[1]) stop(paste("locus", locus,"length of alleleNames differs from dimention of allelic transition matrix"))
    if (as.integer(colnames(allelicTransitionList[[locus]])))==1:length(alleleNames)) colnames(allelicTransitionList[[locus]]))=alleleNames
    if (as.integer(rownames(allelicTransitionList[[locus]])))==1:length(alleleNames)) rownames(allelicTransitionList[[locus]]))=alleleNames
    if (class(alleleTipsList[[locus]])=="integer") alleleTipsList[[locus]]=as.character(alleleTipsList[[locus]])
    if (!all((alleleTipsList[[locus]]%in%alleleNamesList[[locus]]))) stop("alleleTips are not included in alleleNames")
  }
  if (length(demeNames)!=dim(demicTransition)[1]) stop("length of demeNames differs from dimention of demic transition matrix")
  if (is.null(colnames(demicTransition))) colnames(demicTransition)=demeNames
  if (is.null(rownames(demicTransition))) rownames(demicTransition)=demeNames
  if (as.integer(colnames(demicTransition))==1:length(demeNames)) colnames(demicTransition)=demeNames
  if (as.integer(rownames(demicTransition))==1:length(demeNames)) rownames(demicTransition)=demeNames
  if (is.null(names(Ne))) names(Ne)<-demeNames
  if (length(demeTips)!=length(alleleTips)) stop("number of tips differs between alleleTips and demetips")
  if (class(demeTips)=="integer") demeTips=as.character(demeTips)
  if (!all((demeTips%in%demeNames))) stop("demeTips are not included in demeNames")
  new("transitionModel",
      demicTransition=demicTransition,
      allelicTransition=allelicTransition,
      Ne=Ne, 
      demesNames=demeNames,
      allelesNames=alleleNames,
      demicStatusOfStartingIndividuals=demeTips,
      allelicStatusOfStartingIndividuals=alleleTips)
}


# laplaceMatrix returns Laplacian matrix from transition matrix
setMethod("laplaceMatrix", 
          signature=c("transitionModel"),
          definition= function(transitionMod){
            matrixDeme = diag(rep(1,dim(transitionMod@demicTransition)[1])) # diagonal equals to 1
            laplacianMatrixDeme = matrixDeme - transitionMod@demicTransition
            laplacianMatrixDeme[is.na(laplacianMatrixDeme)]<-0 # replace NA by 0
            matrixAllele = diag(rep(1,dim(transitionMod@allelicTransition)[1])) # diagonal equals to 1
            laplacianMatrixAllele = matrixAllele - transitionMod@allelicTransition
            laplacianMatrixAllele[is.na(laplacianMatrixAllele)]<-0 # replace NA by 0
            list(demic=laplacianMatrixDeme,allelic=laplacianMatrixAllele)
          })
# ordinary laplacian
# Boley et al 2011
#

setMethod("ordinary_laplacian",
          definition = function(transitionMod){
            for (i in c("dem","allel"))
            markovB<-new("markovchain", states=transitionMod@demesNames, transitionMatrix=transitionMod@demicTransition)
            PI<-diag(steadyStates(markovB)[1,])
            PI - PI%*%transitionMod@demicTransition
          })

# Calcul of resistance between two points of the graph
# with the Moore-Penrose generalized inverser matrix.
# ref : Bapat et all, A Simple Method for Computing Resistance Distance (2003)
# ref : Courrieu, Fast Computation of Moore-Penrose Inverse Matrices (2005)
#r <- inverseMP
#for (i in 1:dim(r)[1]){
#  for (j in 1:dim(r)[2]){
#    r[i,j] <- inverseMP[i,i]+inverseMP[j,j]-inverseMP[i,j]-inverseMP[j,i]
#  }
#}

setMethod("commute_time_undigraph",
          definition = function(transitionMod){
            laplacian = laplaceMatrix(transitionMod)
            commute_time<-list()
            for (i in names(laplacian)){
              inverseMP = ginv(laplacian[[i]]) # generalized inverse matrix  (Moore Penrose)
              diag = diag(inverseMP) # get diagonal of the inverse matrix
              mii = matrix(diag, nrow =dim(inverseMP), ncol = dim(inverseMP))
              mjj = t(mii)
              mij = inverseMP
              mji = t(mij)
              commute_time[[i]]= mii + mjj - mij - mji
            }
            commute_time
          })


# for weighted digraphs
# Boley et al 2011
# Hitting time is given by 
# H = 1 · [diag(Z )]T − Z 
# Where Z, the scaled Fundamental Matrix is defined by
# Z = (L + pi_ pi_ T)⁻¹ = (PI ( 1 - P π T))⁻¹
# where P is the transition probability matrix, 
# pi_ (or π) is the stationary probablity distribution of P
# and PI is the diagonal matrix of π

setMethod("hitting_time_digraph",
          definition = function(transitionMod){
          H <- list()
          for (i in c("dem","allel")){
            Ones <- rep(1,length(transitionMod[paste(i,"esNames",sep="")]))
            markovB<-new("markovchain", states=transitionMod[paste(i,"esNames",sep="")], transitionMatrix=transitionMod[paste(i,"icTransition",sep="")])
            pi_<-steadyStates(markovB)[1,]
            PI <- diag(pi_)
            L <- PI - PI%*%transitionMod[paste(i,"icTransition",sep="")]
            Z <- ginv(L + pi_%*%t(pi_))
            H[[paste(i,"ic",sep="")]] <- Ones%*%t(diag(Z))-Z
          }
          H})

setMethod("commute_time_digraph",
          definition = function(transitionMod){
            C <- hitting_time_digraph(transitionMod)
            for (i in names(C)){
              C[[i]]=C[[i]]+t(C[[i]])
            }
            C})

setReplaceMethod(f = "$", signature = "transitionModel",
          definition = function(x, name, value) {
            if(name=="demicTransition"){x@demicTransition<-value}
            if(name=="allelicTransition"){x@allelicTransition<-value}
            if(name=="allelicTransition"){x@allelicTransition<-value}
            if(name=="Ne"){x@Ne<-value}
            if(name=="demeNames"){x@demeNames<-value}
            if(name=="alleleNames"){x@alleleNames<-value}
            if(name=="demicStatusOfStartingIndividuals"){x@demicStatusOfStartingIndividuals<-value}
            if(name=="allelicStatusOfStartingIndividuals"){x@allelicStatusOfStartingIndividuals<-value}
            validObject(x)
            return (x)
          }
)

setMethod(f = "$", signature = "transitionModel",
          definition = function(x, name) {
            slot(x,name)
          }
)

setMethod(f = "[", signature = c("model","character"),
          definition = function(x,i) {
            slot(x,i)
          }
)

setMethod(f = "[", signature = c("transitionModel","character"),
          definition = function(x,i) {
            slot(x,i)
          }
)

setMethod(f = "[", signature = c("transitionModel","character"),
          definition = function(x,i) {
            slot(x,i)
          }
)
setReplaceMethod(f = "[", signature = "transitionModel",
                 definition = function(x, i,j, value) {
                   if(i=="demicTransition"){x@demicTransition<-value}
                   if(i=="allelicTransition"){x@allelicTransition<-value}
                   if(i=="allelicTransition"){x@allelicTransition<-value}
                   if(i=="Ne"){x@Ne<-value}
                   if(i=="demeNames"){x@demeNames<-value}
                   if(i=="alleleNames"){x@alleleNames<-value}
                   if(i=="demicStatusOfStartingIndividuals"){x@demicStatusOfStartingIndividuals<-value}
                   if(i=="allelicStatusOfStartingIndividuals"){x@allelicStatusOfStartingIndividuals<-value}
                   validObject(x)
                   return (x)
                 }
)


setMethod("as",
          signature="listOfNodes",
          function(object,Class){
            newObject <- list()
            for (i in names(object)){
              newObject[[i]] <- new("nodE",
                                    nodeNo=object[[i]]@nodeNo,
                                    ancestor=object[[i]]@ancestor,
                                    age=object[[i]]@age,
                                    descendants=object[[i]]@descendant)
            }
            new(Class,newObject)
            }
          )
