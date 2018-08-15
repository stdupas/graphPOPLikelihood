model <- setClass("model",
                  slots=c(model="character",type="character",variableName="character",parametersNames="character"),
                  validity = function(object){
                    if (!object@type%in%c("niche","mutation","dispersion")) stop("type must be one of c('niche', 'mutation', 'dispersion'")
                    if ((object@model=="proportional")&(length(object@parametersName)!=1)) stop("proportional model has one parameter and one only")
                    if ((object@model=="proportional")&(length(object@variableName)!=1)) stop("proportional model has one variable and one only")
                    if ((object@model=="stepwise")&(length(object@parametersName)!=1)) stop("stepwise model has one parameter and one only")
                    })

instanceOfModel <- setClass("instanceOfModel",
                            contains="model",
                            slots=c(parameters="numeric"),
                            validity = function(object){
                              if (length(object@parametersNames)!=(length(object@parameters))) stop("length of @parameterNames differs form lenth of @parameters")
                            })

hyperModel <- setClass("priorModel",
                       contains = "model",
                       slots=c(hyperModel="character",hyperParameters="numeric"),
                       validity = function(object){
                         if (length(object@models)!=length(object@parameters)) stop("length of @models vector differs from length of @parameters list")
                       })


nicheModel <- setClass("nicheModel",
                       slots=c("model","variableEnvironnemental","parameters"),
                       validity = function(object){
                         if (length(object@models)!=length(object@parameters)) stop("length of @models vector differs from length of @parameters list")
                       })

mutationModel <- setClass("mutationModel",
                          slots=c("model","parameter"))

dispersionModel <- setClass("dispersionModel",
                            slots=c("model","parameter"))

model <- setClass("model",
                  slots=c("nicheModel","mutationModel","dispersionModel")
                  )

prior <- setClass("prior",
                  contains="list",
                  validity = function(object){
                    
                  }) # add validity (contains K, R and mutation_rate)

parameters <- setClass("parameters",
                       contains="list") # add validity (contains K, R and mutation_rate)

transitionList <- setClass("transitionList",
                              contains="list",
                              validity=function(object){
                                if(!lapply(object,"class")=="matrix") stop("transitionList is expected to be a list of matrix")})

populationTypeList <- setClass("populationTypeList",
                      contains="list",
                      validity=function(object){
                        if(!lapply(object,"class")=="character") stop("pouplationTypeList is expected to be a list of matrix")})


transitionModel <- setClass("transitionModel",
                           slots = c(transitionList="transitionList",
                                     Ne="numeric", 
                                     populationTypeList="populationTypeList", 
                                     populationTypeOfStartingIndividuals = "data.frame"),
                           validity = function(object){
                             if (any(dim(object@transitionList)!=length(object@populationTypeList))) stop("one dimension of transitionList slot differs from length of populationTypeList slot")
                             if (ncol(object@populationTypeOfStartingIndividuals)!=length(object@populationTypeList)) stop("ncol of populationTypeOfStartingIndividuals data.frame slot differs from length of populationTypeList slot")
                             if (!("demic"%in%names(object@populationTypeList))) stop("'demic' is absent from names of object@populationTypeList")
                             if (!("demic"%in%colnames(object@poulationTypeOfStartingIndividuals))) stop("'demic' is absent from in colnames of object@populationTypeOfStartingIndoviduals data.frame")
                             if (length(object@populationTypeList[["demic"]])!=length(object@Ne)) stop ("length of Ne slot differs from length of populationTypeList slot")
                             if (any(object@populationTypeList[["demic"]])!=names(object@Ne)) stop ("names of object@Ne differ from object@populationTypeList[['demic']]")
                             if (!(all(names(object@transitionList)%in%names(object@populationTypeList))&
                                     all(names(object@populationTypeList)%in%names(object@transitionList)))) stop ("names of object@transitionList differs from names of object@populationTypeList")
                             if (any(names(object@transitionList)!=names(object@populationTypeList))) stop ("names of object@transitionList are not in the same order as names object@populationTypeList")
                             for (Name in names(object@transitionList)){
                               if (!all(object@populationTypeOfStartingIndividuals[,"Name"]%in%object@populationTypeList[[Name]])) stop (paste(which(!(object@populationTypeOfStartingIndividuals[,"Name"]%in%object@populationTypeList[[Name]]),"th starting individuals population type is not defined in  a object@populationType[['",Name,"']]")))
                             }}
) 

# nrow replace dim() because object is a matrix 


demeTransition <- setClass("demeTransition",
                           slots=c(transition="matrix"),
                           contains="RasterLayer")

setClassUnion("numericOrNULL", c("numeric", "NULL"))

nodE <- setClass("nodE",
                   slots=c(nodeNo="integer",ancestor="integer",age="numeric",descendants="integer"),
                   validity=function(object){
                     if (length(object@ancestor)>1) stop("a node has at most one known ancestor")
                   })

listofnodEs <- setClass("listofnodEs",
                   contains="list",
                   validity=function(object){
                     if (all(lapply(object,"class")=="nodE")) TRUE else FALSE
                     })


branchTransition <- setClass("branchTransition",
                             slots=c(age="numeric",ancestorAge="numeric",
                                     statusDemes="character",agesDemes="numeric",
                                     statusAlleles="character",agesAlleles="numeric"))


Node <- setClass("Node",
                      contains="branchTransition",
                      slots = c(nodeNo="integer",descendant="integer",ancestor="integer")
                      )

listOfNodes <- setClass("listOfNodes",
                      contains ="list",
                      validity = function(object){
                        if (all(lapply(object,"class")=="Node")) TRUE else FALSE
                      }
                      )
graphicalNode <- setClass("graphicalNode", contains="Node",
                                 slots=c(x="numeric",y="numeric")
                                 )

graphicalListOfNodes <- setClass("graphicalListOfNodes", contains="list",
                                 validity = function(object){
                                   if (all(lapply(object,"class")=="graphicalNode")) TRUE else FALSE
                                 }
)

spatialListOfNodes <- setClass("spatialListOfNodes",
                        contains ="listOfNodes",
                        slots=c(populationsSizes="RasterLayer"),
                        validity = function(object){
                          if (all(state(object,Inf,"Demes")%in%1:ncellA(object@populationsSizes))) TRUE else FALSE
                        }
)

# start to remove
LandGenealogy <- setClass("LandGenealogy",
                           slots=c(genealogy="Node"),
                           contains="RasterLayer")
#stop to remove

genetic <- setClass("genetic",
                    slots = c(ploidy="integer",ploidyByrow="logical"),
                    contains="data.frame",
                    prototype = prototype(data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=as.integer(2), ploidyByrow=FALSE),
                    validity = function(object){
                      if (all(grepl("ocus",names(object)))) TRUE else stop("col names of genetic data.frame do not contain 'ocus'")
                      if ((object@ploidy==2)&(object@ploidyByrow==FALSE)) {#tip.color=tipcols
                        if (length(grep("\\.1",names(object)))==0|length(grep("\\.2",names(object)))==0) {
                          if ((grep("\\.1",names(object))%%2!=1)|(grep("\\.2",names(object))%%2!=0)){
                            stop("Columns of diploid by row FALSE data frame have to be named as follows: 'c('.1','.2','.1','.2')'")
                          }
                        }
                      }
                    }
)


spatialGenetic <- setClass("spatialGenetic",
                           slots = c(x="numeric", y="numeric",Cell_numbers="integer"),
                           contains = "genetic",
                           prototype = prototype(genetic(),x=c(1,2),y=c(1,1),Cell_numbers=as.integer(c(1,2))),
                           validity = function(object){
                             if (length(object@x)!=length(object@y)) stop("slots x and y do not have the same length")
                             if (length(object@x)!=nrow(object)) stop("slots x and genetic do not have the same number of individuals")
                           }
)

GenealPopGenet <- setClass("GenealPopGenet",
                   contains="listOfNodes",
                   slots=c(Genet="spatialGenetic",Pop="RasterBrick"),

                   )

# to remove start
LandGenetGenealogy <- setClass( "LandGenetGenealogy",
                                contains = "LandGenealogy",
                                slot= c(Genotype="spatialGenetic"))
# stop to remove

genetic <- function(df=data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=NA, ploidyByrow=NA){
  if (is.na(ploidy)) { 
    if (any(grep("\\.2", colnames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",colnames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
    if (any(grep("\\.2", rownames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",rownames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
  }
  if (is.na(ploidyByrow)) ploidyByrow = !(any(grep(paste("\\.",ploidy,sep=""), colnames(df))))
  new("genetic",df,ploidy=ploidy,ploidyByrow=ploidyByrow)
}

spatialGenetic <- function(df=NA,x=NA,y=NA,Cell_numbers=NA)
{
  if(is.na(df)) df=cbind(data.frame(genetic()),x=c(1,2),y=c(1,1),Cell_numbers=c(1,2))
  if (!(all(c("x","y","Cell_numbers")%in%colnames(df)))){
    if (is.na(Cell_numbers)){
      if ("Cell_numbers"%in%colnames(df)) {
        Cell_numbers=df$Cell_numbers
        df <- df[,-which(colnames(df)=="Cell_numbers")]
      }
    }
    if (is.na(x)){
      if ("x"%in%colnames(df)) {
        x=df$x
        df <- df[,-which(colnames(df)=="x")]
      }
    }
    if (is.na(y)){
      if ("x"%in%colnames(df)) {
        y=df$y
        df <- df[,-which(colnames(df)=="y")]
      }
    }
  }
  new("spatialGenetic",genetic(df[,-which(colnames(df)%in%c("x","y","Cell_numbers"))]),x=df$x,y=df$y,Cell_numbers=df$Cell_numbers)
}

