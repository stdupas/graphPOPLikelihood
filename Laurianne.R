library(markovchain)


cleanerData <-function(data) 
  # take off rows wich have missing data by replace this missing data to false data
{
  data[is.na(data)] <- as.integer(1E9)
  data=data[rowSums(data[,grep("Locus",colnames(data),value=TRUE)])<7E9,]
  data[data==1E9]=NA
  data <- TwoCols2OneCol(data)
  
}


TwoCols2OneCol <- function(tip_genotype)
{
  nonGenotypeData <- tip_genotype[,!grepl("ocus",colnames(tip_genotype))]
  
  matrix_pair = tip_genotype[,grep("\\.2", colnames(tip_genotype))]
  colnames(matrix_pair) <- grep("ocus",unlist(strsplit(colnames(matrix_pair),"\\.")),value=TRUE)
  matrix_pair <- cbind(nonGenotypeData,matrix_pair)
  rownames(matrix_pair)<-paste(rownames(matrix_pair),".1",sep="")
  matrix_impair = tip_genotype[,grep("\\.1", colnames(tip_genotype))]
  colnames(matrix_impair) <- grep("ocus",unlist(strsplit(colnames(matrix_impair),"\\.")),value=TRUE)  
  matrix_impair <- cbind(nonGenotypeData,matrix_impair)
  rownames(matrix_impair)<-paste(rownames(matrix_impair),".2",sep="")  
  tip_genotype <- rbind(matrix_pair,matrix_impair)[rep(1:nrow(matrix_pair),each=2)+rep(c(nrow(matrix_pair),0),nrow(matrix_pair)),]
  tip_genotype
}
#


proportional <- function(X,p,Log=FALSE)
{
  if (Log) {log(p[rep("a",dim(X)[1]),colnames(X)]*X)} else {
    p[rep("a",dim(X)[1]),colnames(X)]*X
  }
}


ReactNorm <- function(X,p,shapes)
  
  # X is a matrix with the different environemental variables in column and the map cells in line
  # or a numeric vector
  #p=c(p["Xmin"]=10,p["Xmax"]=20,p["Xopt"]=18,p["Ymax"]=0.1)
  ## shapeDisps in c("enveloppe","envelin","envloglin","loG","conquadratic","conquadraticskewed","conquadraticsq","conquadraticskewedsq")
{
  if (class(p)=="numeric") {p=t(as.matrix(p));rownames(p)="a"}
  if (class(X)=="numeric") {X=data.frame(X=X);colnames(X)=colnames(p)}
  Y=X
  if (!all(colnames(p)%in%names(shapes))) {stop ("variable names do not correspond between parameters 'p' ans 'shape'")}
  for (shape in as.character(levels(as.factor(shapes))))
  {
    variables = colnames(p)[which(shapes==shape)]
    Y[,variables]=switch(shape,
                         constant=p,
                         proportional = proportional(subset(X,select=variables),p), # erreur subset
                         linear = linear(subset(X,select=variables),p),
                         enveloppe=enveloppe(subset(X,select=variables),p),
                         envelin=envelinear(subset(X,select=variables),p),
                         envloglin=envelinear(subset(X,select=variables),p,log=TRUE),
                         loG = log(subset(X,select=variables)),
                         conquadratic=conquadratic(subset(X,select=variables),p),
                         conquadraticskewed=conquadraticskewed(subset(X,select=variables),p),
                         conquadraticsq=conquadraticsq(subset(X,select=variables),p),
                         conquadraticskewedsq=conquadraticskewedsq(subset(X,select=variables),p)
    )
  }
  Y=cbind(Y,Y=apply(Y, 1, prod)^(1/dim(p)[2])) # geometric mean
  Y
}


# populationSize: uses K_Function to obtain the populationsize landscape raster
#
#
populationSize <- function(donneesEnvironmentObs, p, shapes)
{
  populationSize <- donneesEnvironmentObs
  values(populationSize) <- ReactNorm(valules(donneesEnvironmentObs), p, shapes)[1,]
  populationSize
}



# (optional) distanceMatrix return distance between all cells of raster
distanceMatrix <- function(rasterStack){
  #get x and y coordinates for each cell of raster object put in parameters
  coords = xyFromCell(rasterStack, 1:length(values(rasterStack[[1]])), spatial=FALSE)
  distance = as.matrix(dist(coords)) # distance matrix of coordinates
  return(distance)
}




# migrationMatrix return matrix of migration rate for each cell
# From the raster Stack and the dispersal parameter pDisp (estimated),
# this function calculates distance between all cells of raster
# To this matrix of distance is applied a shapeDisp of migration.
# And this function allows to choose the shapeDisp of study:
# shapeDisp:  "gaussian" a simple normal density distribution, 
#           "exponential" density distribution,
#           "fat_tail1", "fat_tail2": ref :Chapman et all, Journal of Animal Ecology (2007) 76 , 36– 44
#           "island" probability 1-m to stay, else homogen dispersion,
#           "contiguous" near dispersal.
migrationMatrix <- function(rasterStack,shapeDisp, pDisp){
  distMat<-distanceMatrix(rasterStack)
  Ndim = 1+all(ncell(rasterStack)!=dim(rasterStack)[1:2])
  migration = apply(distMat, c(1,2), 
                    function(x)(switch(shapeDisp,
                                       # 1: alphaDisp   2: betaDisp ; note: esperance = 1/alphaDisp
                                       fat_tail1 = 1/(1+x^pDisp[,2]/pDisp[,1]), # Molainen et al 2004
                                       # 1: sigmaDisp
                                       gaussian = (dnorm(x, mean = 0, sd = pDisp[,1], log = FALSE)),
                                       # 1: sigma Disp
                                       exponential = (dexp(x, rate = 1/pDisp[,1], log = FALSE)),
                                       # 1: sigmaDisp
                                       contiguous = (x==0)*(1-pDisp[,1])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp[,1]/(2*Ndim)),
                                       # 1: sigmaDisp
                                       contiguous8 = (x==0)*(1-pDisp[,1])+((x>0)-(x>2*res(rasterStack)[1]))*(pDisp[,1]/(4*Ndim)),
                                       #island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]/(nCell-1)),
                                       island = (x==0)*(1-pDisp[,1])+(x>0)*(pDisp[,1]),
                                       #1: sigmaDisp    2: gammaDisp
                                       fat_tail2 = x^pDisp[,2]*exp(-2*x/(pDisp[,1]^0.5)),
                                       #1: sigmaDisp, 2: 
                                       contiguous_long_dist_mixt = pDisp["plongdist"]/ncellA(rasterStack)+(x==0)*(1-pDisp["pcontiguous"]-pDisp["plongdist"])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp["pcontiguous"]/2),
                                       gaussian_long_dist_mixt = pDisp[,2]/ncellA(rasterStack) + (dnorm(x, mean = 0, sd = pDisp[,1], log = FALSE))
                    )))
  return(migration)
}




enveloppe <- function(X,p)
{
  p[rep("Yopt",dim(X)[1]),colnames(X)]*((X>p[rep("Xmin",dim(X)[1]),])&(X<p[rep("Xmax",dim(X)[1]),colnames(X)]))
}

# envelope = linear response within an envelope
# X : matrix or data frame providing the values of independent variable to calculate reaction norm
# p : matrix parameter values for the reaction norm
# line names of p used for caclculation : c("Xmin","Xmax","Yxmin","Yxmax"), 
# column names of p : names of the independent variables of X used for the reaction norm calculation
# [Xmin, Xmax] is the enveloppe, 
# Yxmin and Yxmax are the values at Xmin and Xmax



# # transitionMatrix obtained with an isotropic migration hypothesis for a backward model

transitionMatrixBackward <- function(rasterStack=environmentalData,prior){
  listeSample=sampleP(prior)
  X=valuesA(rasterStack)
  K = ReactNorm(X,listeSample$K$busseola$p,listeSample$K$busseola$model)[,"Y"]
  r = ReactNorm(X,listeSample$R$busseola$p,listeSample$R$busseola$model)[,"Y"] 
  migration <- migrationMatrix(rasterStack,listeSample$dispersion$busseola$model, listeSample$dispersion$busseola$p)
  if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)>1)&(length(K)>1)) {transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  transition = transition / rowSums(transition)  # Standardisation
  transition = transition * as.numeric(transition>1e-6) # removal of values below 1e-6
  transition = transition / rowSums(transition)  # Standardisation again
  transition
}



# il faut que $p soit une data frame  colonne n de variable  une seule ligne a --> effectue dans Initialisation
### function to sample  
sampleP <- function(prior) {
  Result=list()
  for(parametreBio in names(prior)) {
    for (variableEnvironnemental in names(prior[[parametreBio]])) {
        Result[[parametreBio]] <- list() 


         Result[[parametreBio]][[variableEnvironnemental]]$p <- switch(prior[[parametreBio]][[variableEnvironnemental]]$a$distribution,
                                uniform=data.frame(variableEnvironnemental=runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                fixed =data.frame(variableEnvironnemental=prior[[parametreBio]][[variableEnvironnemental]]$a$p),
                                normal=data.frame(variableEnvironnemental=rnorm(1,mean=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],sd=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                loguniform=data.frame(variableEnvironnemental=log(runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1]))),
                                uniform=data.frame(variableEnvironnemental= runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                fixed =data.frame(variableEnvironnemental= prior[[parametreBio]][[variableEnvironnemental]]$a$p),
                                normal=data.frame(variableEnvironnemental =rnorm(1,mean=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],sd=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1])),
                                loguniform=data.frame(variableEnvironnemental= log(runif(1,min=prior[[parametreBio]][[variableEnvironnemental]]$a$p[1,1],max=prior[[parametreBio]][[variableEnvironnemental]]$a$p[2,1]))))
         
     Result[[parametreBio]][[variableEnvironnemental]]$model <-  prior[[parametreBio]][[variableEnvironnemental]]$model
     colnames(Result[[parametreBio]][[variableEnvironnemental]]$p)=variableEnvironnemental
     rownames(Result[[parametreBio]][[variableEnvironnemental]]$p)=c("a")
     names(Result[[parametreBio]][[variableEnvironnemental]]$model)=variableEnvironnemental
     
    }
  }
#  names(Result) <- names(prior)
  return(new("parameters",Result))
}


# Aggregate_and_adjust_raster_to_data change resolution and extent of environmental stacked layers
# according to data geographic range and extension zone outside geographic range of data
# ARGUMENTS:
# Envir_raster_stack = raster file 
# xy = genetics data optionnal
# extend_band_size =

Aggregate_and_adjust_raster_to_data <- function(Envir_raster_stack,xy=NULL,extend_band_size=NA,aggregate_index){
  resDivAbs<-dim(Envir_raster_stack)[1]%%aggregate_index
  resDivOrd<-dim(Envir_raster_stack)[2]%%aggregate_index
  # dimention new windows [dim(Envir-raster-stack)[1]/aggregate_index][dim(Envir-raster-stack)[1]/aggregate_index]
  # note = dans ?extend on met les coordonnee
    if (aggregate_index > 1) {
      Envir_raster_stack <-aggregate(crop(Envir_raster_stack,extent(as.vector(extent(xy)+extend_band_size *c(-1*res(Envir_raster_stack)[1],1*res(Envir_raster_stack)[1],
                                                                                                                     -1*res(Envir_raster_stack)[2],1*res(Envir_raster_stack)[2])))), fact=aggregate_index, fun=mean, expand=TRUE, na.rm=TRUE)
    }

  Envir_raster_stack
}



#
# Absorbing transition matrix
# Estimates probability of transition among transient and absorbant states
# for pairs of genes in demes newtork
# dimension is N(N+3)/2
# arg: 
# transition : backward transition probability matrix of genes among demes
# N : population sizes of the demes
# Mcrae 2006
# Hey 2001


AbsorbingTransition <- function(transition,N)
{
  N[N<1]=1
  Ndeme <- dim(transition)[1]
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
  QheteroHetero <- matrix(NA,Nhetero,Nhetero)
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
          QheteroHetero[ij,kl] <- transition[i,k]*transition[j,l]
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
    }
  }
  
  # then transient matrix 
  QhomoHetero <- matrix(NA,Ndeme,Nhetero)
  kl=0
  Check=TRUE
  for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
    Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
    # in the lines of Q matrix
    for (k in 2:Ndeme){
      for (l in 1:(k-1)){ # only heterodemic targets
        kl=kl+1
        QhomoHetero[i,kl] <- transition[i,k]*transition[i,l]    # i=j (homodemic sources)
        #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
      }
    }
    kl=0
  }
  QheteroHomo <- matrix(NA,Nhetero,Ndeme)
  ij=0
  for (i in 2:Ndeme){
    for (j in 1:(i-1)){ # only heterodemic sources
      ij=ij+1
      for(k in 1:Ndeme){ # only homodemic targets
        QheteroHomo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k,]))
        # homodemic targets that have not coalesced
      }
    }
  }
  QhomoHomo <- transition*transition*matrix(1-1/(2*N[,1]),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
  Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
  result=list(Q=Q,Qline=Qline)
  result
}

#to modify structure genealogy

change_genealogy <- function(genealogy) {
  for (event in 1:length(genealogy)) {
    genealogy[[event]]$colesc/home/egce/Téléchargements/coalescent.rdaing<-list(c(genealogy[[event]]$coalescing,c("vector of position","x","y"),c("changement of position","x")))
  }
genealogy
  }

# to sample node in tree genealogy


change_node <- function(genealogy){
  coalescent_event <- sample(length(genealogy),1)
  new_node <- genealogy[[coalescent_event]]$new_node
  # get the coalescing modified
  coalescing_modified <- sample(genealogy[[coalescent_event]]$coalescing,1)
  # where is new node branched ?
  coalescent_event_where_new_node_of_the_brother_will_be_branched <-  grep(new_node,lapply(genealogy,"[[","coalescing"))
  # replace new node in this coalescent event with the coalescing node that is not modified
  # get the coalescing node that is not modified
  for( i in 1:length(genealogy[[coalescent_event]]$coalescing)) {
    if (genealogy[[coalescent_event]]$coalescing[i]!=coalescing_modified){
      coalescing_brother_no_modifie=genealogy[[coalescent_event]]$coalescing[i]
    }
  }
  genealogy[[coalescent_event_where_new_node_of_the_brother_will_be_branched]]$new_node=coalescing_brother_no_modifie
  # modify the information of the coalescent event that has been removed
  
  # which coalescing the node will be associated with
  all_coalescing <- unlist(lapply(genealogy,"[[","coalescing"))
  potential_new_associated_node <- all_coalescing[which(all_coalescing!=coalescing_modified)]
  new_associated_node<-sample(potential_new_associated_node,1)
  # where the node will be
  # modify the brother of the coalescing modified
  genealogy[[coalescent_event]] <-
    # take the same event put 
  genealogy[[coalescent_event]]$new_node= 
  genealogy[[coalescent_event]]$coalescing=c(new_associated_node,coalescing_modified)
}







# Probability distribuiton of coalescence times
#

coalescenceTimeProbDistrib <- function(Qlist){
  Ndeme <- (2*dim(Qlist[[1]])[1]+1/4)^.5-.5
  cTPD <- array(NA,dim=c(Ndeme,Ndeme,1))
  condition=TRUE
  t=1
  QPower <- Q <- Qlist[["Q"]]
  cTPD[,,t] <- matrix(colSums((matrix(1,nrow(QPower),ncol(QPower))-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme)
  while (condition){
    t=t+1
    QPowerBefore <- QPower
    QPower <- Q%*%QPowerBefore
    cTPD <- abind(cTPD,matrix(colSums((QPowerBefore-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme),along=3)
    condition=!all(colSums(cTPD,3)>0.9999)
  }
  cTPD
}


#
# distribution of coalescence time among pair of genes depending on their deme of sampling at t=0
# arguments : max_time_interval : maximum time interval to jump over generations
#                                 in the probability calculation
#
#
# value: list 
#        -first : array of coalescence times for deme pairs
#        -second : matrix of expected coalescence times (mean of the distribution)
#


# coalescent_2_newick
# fuinction that converts coalescent
# 

coalescent_2_newick <- function(coalescent)
{
  #  tree=paste("(",coalescent[[length(coalescent)]]$new_node,")",sep="")
  tree=paste(" ",coalescent[[length(coalescent)]]@nodeNo," ",sep="")
  for (i in length(coalescent):(1+length(state(coalescent,0,"nodeNo"))))
  {
    Time = coalescent[[i]]@age
    coalesc <- as.character(coalescent[[i]]@descendant)
    br_length <-  unlist(lapply(as.character(coalescent[[i]]@descendant),function(node){coalescent[[i]]@age-coalescent[[node]]@age}))
    tree <- str_replace(tree,paste(" ",as.character(coalescent[[i]]@nodeNo)," ",sep=""),paste(" ( ",paste(" ",coalesc," :",br_length,collapse=" ,",sep=""),") ",sep=""))
  }
  tree <- gsub(" ","",paste(tree,";",sep=""))
  tree
}

# to an other type of genealogy

coalescent_2_newick_new <- function(coalescent)
{
  #  tree=paste("(",coalescent[[length(coalescent)]]$new_node,")",sep="")
  tree=paste(" ",coalescent[[length(coalescent)]]@nodeNo," ",sep="")
  for (i in length(coalescent):1)
  {
    Time = coalescent[[i]]@age
    coalesc <- as.character(coalescent[[i]]@descendant)
    br_length <-  lapply(as.character(coalescent[["9"]]@descendant),function(node){coalescent[["9"]]@age-coalescent[[node]]@age})
    tree <- str_replace(tree,paste(" ",as.character(coalescent[[i]]@nodeNo)," ",sep=""),paste(" ( ",paste(" ",coalescent[[i]]@nodeNo," :",br_length,collapse=" ,",sep=""),") ",sep=""))
  }
  tree <- gsub(" ","",paste(tree,";",sep=""))
  tree
}

coalescent_2_phylog <- function(coalescent)
{
  read.tree(text=coalescent_2_newick(coalescent))
}


coalescence_prob_time_distribution_matrix <- function(transition,max_time_interval=4,rasK,threshold=1E-6)
{
  # calculates the probability distribution of coalescence among pairs of
  # genes in different cells (cell 1 as line, cell 2 as column)
  # knowing "transition" = gene backward transition probability among demes
  #         "time"
  #         "max_time_interval" = maximum number of generations to group 
  #                           the computation of coalescence
  #         "rasK" = population sizes
  # uses function "matrix.power" of package matrixcalc
  # value:
  # a list of matrix of coalescence probability among genes in the deme graph 
  # at different times in the past starting in the following sequence
  # 1,2,4,16,32,...,2^floor(log2(max_time_grouping)),2*2^floor(log2(max_time_grouping))
  # 3*2^floor(log2(max_time_grouping)) ... n*2^floor(log2(max_time_grouping))
  # untill probability of coalescence is lower than a threshold and reducing everywhere in 
  # the matrix
  time_points <- NA # this vector will contain all the time points where probabilities
  # will be calculated
  occurence_prob <- array(NA, dim=c(dim(transition),1))# occurence probabilities among demes of individuals at the generation
  # given in the time_points values at the same position
  # in time_point vector as in the list
  occurence_and_coalescence <- array(NA, dim=c(dim(transition),1))# occurence and coalescence probabilities among demes of 
  # second  individual at the time considered in the 
  # time_point vector
  popSizes_receiving <- matrix(valuesA(rasK),nrow=dim(transition)[1],ncol=dim(transition)[2],byrow=TRUE)
  coal_column <- 1/(2*popSizes_receiving)
  coal_column[coal_column>1]<-1
  coal_column_tiMax <- 1-(1-coal_column)^max_time_interval
  not_coalesced_prob <- array(NA, dim=c(dim(transition),1))
  coalescence_prob <- array(NA, dim=c(dim(transition),1))#array(NA,dim=)
  time_points <- 1
  occurence_prob[,,1] <-  diag(1,dim(transition)[1])
  occurence_and_coalescence[,,1] <- occurence_prob[,,1] / (2*popSizes_receiving)
  coalescence_prob[,,1] <- occurence_prob[,,1] %*% t(occurence_and_coalescence[,,1])
  coalescence_prob[,,1][coalescence_prob[,,1]>1]=1
  not_coalesced_prob[,,1] <- 1-coalescence_prob[,,1]
  condition=TRUE
  i=1
  maxtransition <- matrix.power(transition,max_time_interval)
  while (condition){
    i=i+1
    time_interval <- time_points[i-1]
    if (time_interval>max_time_interval) {
      time_interval <- max_time_interval;transition_powered = maxtransition
      coal_column_ti <- coal_column_tiMax
    } else {
      transition_powered = matrix.power(transition,time_interval)
      coal_column_ti <- 1-(1-coal_column)^time_interval
    }
    time_points[i] <- time_points[i-1]+time_interval
    occurence_prob <- abind(occurence_prob,occurence_prob[,,i-1] %*% transition_powered,along=3)
    occurence_and_coalescence <- abind(occurence_and_coalescence,occurence_prob[,,i] *coal_column_ti,along=3)
    coalescence_prob <- abind(coalescence_prob,occurence_prob[,,i] %*% t(occurence_and_coalescence[,,i]),along=3)
    coalescence_prob[,,i][coalescence_prob[,,i]>1]=1
    not_coalesced_prob <- abind(not_coalesced_prob,not_coalesced_prob[,,i-1]*(1-coalescence_prob[,,i]),along=3)
    coalescence_prob[,,i] <- coalescence_prob[,,i]*not_coalesced_prob[,,i-1]
    #condition = any(coalescence_prob[[i]]>=coalescence_prob[[i-1]])|any(coalescence_prob[[i]]>threshold)
    condition = !all(not_coalesced_prob[,,i]<threshold)
  }
  dimnames(coalescence_prob) <- list(1:dim(transition)[1],1:dim(transition)[2],time_points)
  expected_coalescence_times <- matrix(NA,nrow=dim(transition)[1],ncol=dim(transition)[2])
  for (i in 1:dim(transition)[1]){
    for (j in 1:dim(transition)[2]){
      expected_coalescence_times[i,j] <- coalescence_prob[i,j,]%*%time_points
    }
  }
  #expected_coalescence_times<-array(unlist(coalescence_prob),dim=c(dim(transition),length(coalescence_prob)))
  #expected_coalescence_times<-unlist(coalescence_prob,dim=c(dim(matrix),length(coalescence_prob))
  list(coalescent_prob=coalescence_prob,exp_times=expected_coalescence_times)
}


combn = function(n,x){factorial(n)/(factorial(n-x)*factorial(x))}
