#Initialisation



# fixed uniform normal loguniform lognormal



wd="/M1bi 2015 2016/stage/busseola/"

wd="~/Documents/Lauriane/busseola/"

#  wd="/home/dupas/busseola"

# wd= "~/busseola/"
setwd(wd)
library(raster)
library(rgdal) # necessary to use function raster()
library(stringr)
source("class.R")
source("generic.R");source("method.R")
source("Laurianne.R")
load("Prior.rda")

#
# Create priors
#

# Prior distributions
Prior=list()
Prior$K$busseola$a$distribution = "uniform"
Prior$K$busseola$model = data.frame(busseola="proportional")
Prior$K$busseola$a$p = data.frame(busseola=c(min=0.001,max=0.5))
Prior$R$busseola$a$distribution = "fixed"
Prior$R$busseola$model =data.frame(busseola="constant")
Prior$R$busseola$a$p = 20
Prior$mutation_rate$busseola$model =data.frame(busseola= "stepwise")
Prior$mutation_rate$busseola$a$distribution = "loguniform"
Prior$mutation_rate$busseola$a$p =data.frame(busseola=c(min=1E-6,max=1E-2))
Prior$dispersion$busseola$model=data.frame(busseola="contiguous")
Prior$R$busseola$model = c(K="constant")

Prior$K$busseola$model =c("proportional")
Prior$K$busseola$a$p = data.frame(busseola=c(min=0.001,max=0.5))
Prior$R$busseola$a$distribution = "fixed"
Prior$R$busseola$model =c("constant")
Prior$R$busseola$a$p = 20
Prior$mutation_rate$busseola$model =data.frame(busseola= "stepwise")
Prior$mutation_rate$busseola$a$distribution = "loguniform"
Prior$mutation_rate$busseola$a$p =data.frame(busseola=c(min=1E-6,max=1E-2))

Prior$dispersion$busseola$model="contiguous"
Prior$dispersion$busseola$model=data.frame(busseola="contiguous")
Prior$R$busseola$model = c("constant")
Prior$R$busseola$a$p = 20
Prior$mutation_rate$busseola$model =c("stepwise")
Prior$mutation_rate$busseola$a$distribution = "loguniform"
Prior$mutation_rate$busseola$a$p =data.frame(busseola=c(min=1E-6,max=1E-2))
Prior$dispersion$busseola$model=c("contiguous")
Prior$dispersion$busseola$a$distribution="loguniform"
Prior$dispersion$busseola$a$p=data.frame(busseola=c(min=0.001,max=0.5))
Prior<-new("prior",Prior)
save(Prior,file="Prior.rda")

populations <- raster(matrix(c(20,30,40,10,NA,NA,20,30,40),nrow=3),xmn=0,xmx=3,ymn=0,ymx=3)
populations <- raster(matrix(c(8,4,2,4,8),ncol=5),xmn=0,xmx=5,ymn=0,ymx=1)
populations <- raster(matrix(c(4,8,3,NA,4,8),ncol=3,nrow=2),xmn=0,xmx=3,ymn=0,ymx=2)
plot(populations,main=names(populations))
names(populations)<-"busseola"
save(populations,file= "populations.rda")

#
# Load data
#


#genealogie 
load("genealogy.rda")

# donnes genetique pour le mini jeu de donnees genotypes
load("genotypes.rda")
#genotypes<-data.frame(Cell_numbers=as.integer(c(4,7,7,9,6)),x=c(.5,.5,.5,2.5,2.5),y=c(1.5,.5,.5,.5,1.5),locus1=as.integer(c(122,122,122,120,120)))
#genotypes<-data.frame(Cell_numbers=as.integer(c(1,4,4,3,6)),x=c(.5,.5,.5,2.5,2.5),y=c(1.5,.5,.5,1.5,.5),locus1=as.integer(c(122,122,122,120,120)))
#save(genotypes,file="genotypes.rda")
#genotypes<-data.frame(Cell_numbers=as.integer(c(1,2,3,3,3)),x=c(.5,1.5,.5,.5,.5),y=c(1.5,1.5,.5,.5,.5),locus1=as.integer(c(122,122,120,120,120)))

# population data
load("populations.rda")
#populations <- raster(matrix(c(4,8,3,NA,4,8),ncol=3,nrow=2),xmn=0,xmx=3,ymn=0,ymx=2)
#populations <- raster(matrix(c(4,8,12,2,NA,NA,4,6,8),ncol=3,nrow=3),xmn=0,xmx=3,ymn=0,ymx=3)
#names(populations)<-"busseola"
#save(populations,file="populations.rda")
#populations <- raster(matrix(c(4,8,3,NA),ncol=2,nrow=2),xmn=0,xmx=2,ymn=0,ymx=2)
#populations <- raster(matrix(c(4,8,3,NA),ncol=2,nrow=2),xmn=0,xmx=2,ymn=0,ymx=2)
plot(populations)
alleles <- raster(matrix(1,ncol=2,nrow=1),xmn=0,xmx=2,ymn=0,ymx=1)
#save(alleles,file="alleles.rda")
#load("alleles.rda")
# prior data
load("Prior.rda")

# spatial genetic classes
gen <- new("genetic",data.frame(Locus1=genotypes[,"locus1"]),ploidy=as.integer(1), ploidyByrow=FALSE)
spgen <- new("spatialGenetic",gen,x=genotypes[,"x"],y=genotypes[,"y"],Cell_numbers=as.integer(genotypes[,"Cell_numbers"]))
#popgg <- new("GenealPopGenet",genealogy,Genet=spgen,Pop=populations)


#
# Simulate genealogy
#
Parameters <- sampleP(Prior)
transition=transitionMatrixA(object1 = populations,object2 = Parameters)
#demeTrans <- new("demeTransition",populations,transition=transition)
allelicTransition <- matrix(c(1-mutRate,mutRate,mutRate,1-mutRate),nrow=2)
dimnames(allelicTransition) <- list(c(120,122),c(120,122))

transitionMod=transitionModel(demicTransition =transitionMatrixA(object1 = populations,object2 = Parameters),
                              allelicTransition = transitionMatrixA(object1 = alleles,object2 = Parameters,option="allele"),
                              demeNames = names(valuesA(populations)),
                              alleleNames = c("120","122"),
                              alleleTips = genotypes$locus1,
                              demeTips = genotypes$Cell_numbers,
                              Ne = populations)

coalescent=simul_coalescent(transitionMod)
myPlot(coalescent,transitionMod)
statesdf <- data.frame(nodeNo=1:5,alleles=c(122,122,120,120,122),demes=spgen@Cell_numbers,time=0)
states <- statesdf[,c("alleles","demes")]
alleles=transitionMod@allelesNames
names(alleles)<-1:length(alleles)
alleleStatus = as.integer(statesdf$alleles)
demes=rownames(transitionMod@demesNames)
demeStatus=as.integer(statesdf$demes)
Ne=valuesA(populations)


a=list()
a[[1]]=  new("graphicalNode",
      age=0,
      ancestorAge=5,
      statusDemes=as.character(1,4),
      agesDemes=c(0,4),
      statusAlleles=as.character(122),
      agesAlleles=c(0),
      nodeNo=as.integer(1),
      descendant=integer(0),
      ancestor=as.integer(4),
      x=0,
      y=0
      )
a[[2]]=  new("graphicalNode",
      age=0,
      ancestorAge=5,
      statusDemes=as.character(4),
      agesDemes=c(0),
      statusAlleles=as.character(122),
      agesAlleles=c(0),
      nodeNo=as.integer(2),
      descendant=integer(0),
      ancestor=as.integer(4),
      x=.5,
      y=0)
a[[3]]=  new("graphicalNode",
      nodeNo=as.integer(3),
      age=0,
      ancestorAge=140,
      statusDemes=as.character(as.integer(c(6,3,2,1,4))),
      agesDemes=c(0,10,135,136,140),
      statusAlleles=as.character(as.integer(c(120,122))),
      agesAlleles=c(0,112),
      descendant=integer(0),
      ancestor=as.integer(5),
      x=1,
      y=0)
a[[4]]=  new("graphicalNode",
      nodeNo=as.integer(4),
      age=5,
      ancestorAge=140,
      statusDemes=as.character(4),
      agesDemes=5,
      statusAlleles=as.character(122),
      agesAlleles=c(5),
      descendant=as.integer(c(1,2)),
      ancestor=as.integer(5),
      x=.25,
      y=5/140)
a[[5]]=  new("graphicalNode",
      age=140,
      ancestorAge=Inf,
      statusDemes=as.character(4),
      agesDemes=5,
      statusAlleles=as.character(122),
      agesAlleles=c(5),
      nodeNo=as.integer(4),
      descendant=as.integer(c(3,4)),
      ancestor=as.integer(),
      x=(.25+1)/2,
      y=1
      )
names(a)<- 1:5
Coalescent <- new("graphicalListOfNodes",a)


# 
environmentalData <- raster("busseola.tif")
genetData <- read.table("WBf16genelandcoord.txt")
genetData <- cbind(genetData,read.table("WBf16genelandgeno.txt"))
genetSP <- SpatialPoints(genetData[,c("x","y")])
bbox(genetSP)
environmentalData <- crop(environmentalData, bbox(genetSP)+res(environmentalData)*c(-2,-2,2,2))
genetData$Cell_numbers <- cellFromXY(environmentalData,genetData)
environmentalData[environmentalData==0] <- NA
ncellA(environmentalData)
any(is.na(extract(environmentalData,genetData[,c("x","y")])))
plot(environmentalData)
plot(genetSP,add=TRUE)
min(extract(environmentalData,genetData[,c("x","y")]))
genetData[which(extract(environmentalData,genetData[,c("x","y")])<0.01),]
genetData=cleanerData(genetData)
listePrior=sampleP(Prior) 
matriceMigration=migrationMatrixA(environmentalData,listePrior$dispersion$busseola$model,listePrior$dispersion$busseola$p)
matriceTrans=transitionMatrixA(environmentalData,Prior)
listePrior=sampleP(Prior) 
plot(raster(matriceMigration))
plot(raster(matriceTrans))
#diminution des pixels
environmentalDataSimple=Aggregate_and_adjust_raster_to_data(environmentalData,xy=genetData[,c("x","y")],extend_band_size=1,aggregate_index=4)
plot(environmentalDataSimple)
matriceTrans=transitionMatrixA(environmentalDataSimple,Prior)
plot(raster(matriceTrans))
matriceAbsorbante=absorbingTransitionA(matriceTrans,valuesA(environmentalDataSimple))
plot(raster(matriceAbsorbante[[1]]))

# avec l utilisation de la class

Matrice_de_transition=new("Transition_Matrix",matriceTrans,populationsize=environmentalDataSimple)
absorbante=absorbingTransitionA(Matrice_de_transition)
# pour simplifier 

ma=raster(matrix(c(20,30,40,10,NA,NA,20,30,40),nrow=3,ncol=3))
extent(ma)=c(0,3,0,3)
names(ma)="busseola"
plot(ma)
dim(ma)
matriceAbsMa=migrationMatrixA(ma,listePrior$dispersion$busseola$model,listePrior$dispersion$busseola$p)

plot(raster(matriceAbsMa))
matriceTransMa=transitionMatrixA(ma,Prior)
plot(raster(matriceTransMa))
matriceAbsorbanteMa=absorbingTransitionA(matriceTransMa,valuesA(ma))
plot(raster(matriceAbsorbanteMa[[1]]))
pdf("exemplesimple.pdf")
par(mfrow=c(2,2))
plot(ma,main="raster de depart")
plot(raster(matriceAbsMa),main="matrice de migration")
plot(raster(matriceTransMa),main="matrice de transition")
dev.off()

pdf("comparaison.pdf")
par(mfrow=c(4,2))
plot(ma,main="raster de depart simple")
plot(environmentalDataSimple,main="raster de depart busseola")
plot(raster(matriceAbsMa),main="matrice de migration simple")
plot(raster(matriceMigration),main="matrice de migration busseola")
plot(raster(matriceTransMa),main="matrice de transition simple")
plot(raster(matriceTrans),main="matrice de transition busseola")
plot(raster(matriceAbsorbanteMa[[1]]),main="matrice absorbante simple")
plot(raster(matriceAbsorbante[[1]]),main="matrice absorbante busseola")
dev.off()


#genealogie 
load("genealogy.rda")

# donnes genetique pour le mini jeu de donnees genotypes
load("genotypes.rda")

# population data
load("populations.rda")
gen <- new("genetic",data.frame(Locus1=genotypes[,"Locus1"]),ploidy=as.integer(1), ploidyByrow=FALSE)
spgen <- new("spatialGenetic",gen,x=genotypes[,"x"],y=genotypes[,"y"],Cell_numbers=genotypes[,"Cell_numbers"])
popgg <- new("GenealPopGenet",genealogy,Genet=spgen,Pop=populations)

genealogy[[2]]
library(ape)
library(stringr)
plot(coalescent_2_phylog(genealogy))
coalescent=new("Genealogy",genealogy,genotypes[,2])
plotgenealogy(coalescent,tipcols = NA)
geneasimple=new("LandGenealogy",ma,genealogy=coalescent)
plotLandG(geneasimple,rasK = NULL)
geneaandgenet=new("LandGenetGenalogy",geneasimple,Genotype=spgen)


########## modify gealogy to colorate this

genealogymod=change_genealogy(genealogy)


# try to use phtools for colorate but no my data
# our data = ma (raster simplify dim (3,3,1))

## set seed for reproducibility
set.seed(100)
## load phytools
library(phytools)
## simulate a tree & some data
tree <- pbtree(n = 26, scale = 1)
tree$tip.label <- LETTERS[26:1]
x <- fastBM(tree)

Q <- matrix(c(-2, 1, 1, 1, -2, 1, 1, 1, -2), 3, 3)
colnames(Q) <- rownames(Q) <- letters[1:3]
x <- sim.history(tree, Q, anc = "a")$states
print(x)
# estimate ancestral states under a ER model
fitER <- rerootingMethod(tree, x, model = "ER")
print(fitER)
plotTree(tree, setEnv = TRUE, offset = 0.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow"), cex = 0.6)
tiplabels(pie = to.matrix(x, sort(unique(x))), piecol = c("blue", "red", "yellow"), 
          cex = 0.3)


# Mutation landscape
#

mutland <- raster(matrix(1:4,nrow=2))
extent(mutland) <- c(0,2,0,2)
plot( mutland,col=c("blue","red","yellow","green"))
plot(SpatialPointsDataFrame(SpatialPoints(xyFromCell(mutland,1:4)),data.frame(Genotype=c("A","T","G","C"))),legend=c("A","T","G","C"),add=TRUE)


load("coalescent.rda")
coalescent
