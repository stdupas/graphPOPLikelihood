setGeneric(
  name = "ncellA",
  def = function(object) { return(standardGeneric("ncellA"))})

setGeneric(
  name = "valuesA",
  def = function(object) { return(standardGeneric("valuesA"))})

setGeneric(
  name="nodes",
  def= function(object) { return(standardGeneric("nodes"))})

setGeneric(
  name="state",
  def=function(object,age,type) { return(standardGeneric("state"))})

setGeneric(
  name="currentState",
  def= function(object,type,nodes) { return(standardGeneric("currentState"))})

setGeneric(
  name="currentNodes",
  def= function(object,age) { return(standardGeneric("currentNodes"))})

setGeneric(
  name="varnames",
  def=function(object){return(standardgeneric("varnames"))})

setGeneric(
  name = "cellNumA",
  def = function(object) { return(standardGeneric("cellNumA"))})


setGeneric(
  name = "xyFromCellA",
  def = function(object,cellNum) { return(standardGeneric("xyFromCellA"))})

setGeneric(
  name="distA",
  def=function(object)
  { return(standardGeneric("distA"))
  }
)

setGeneric( 
  name="distanceMatrixA",
  def = function(object) { return(standardGeneric("distanceMatrixA"))}
)

setGeneric(
  name="migrationMatrixA",
  def= function(object,shapeDisp,pDisp) { return(standardGeneric("migrationMatrixA"))}
)

setGeneric(
  name="transitionMatrixA",
  def= function(object1,object2,...)  { return(standardGeneric("transitionMatrixA"))}
)

setGeneric(
  name="absorbingTransitionA",
  def= function(object,...)  { return(standardGeneric("absorbingTransitionA"))}
)

setGeneric(
  name="plotgenealogy",
  def=function(object,tipcols) {return(standardGeneric("plotgenealogy"))}
)

setGeneric(
  name = "plotLandG",
  def=function(object,rasK) {return(standardGeneric("plotLandG"))}
)

setGeneric(
   name="nodesByStates",
   def=function(object,age,Which) {return(standardGeneric("nodesByStates"))}
 )

setGeneric(name="setState",
           def=function(object,Nodes,attribut,newValues){return(standardGeneric("setState"))}
)

setGeneric(name="sampleP",
           def=function(prior){return(standardGeneric("sampleP"))})
          
setGeneric(name="simul_coalescent",
           def=function(transitionMod,...)
             {return(standardGeneric("simul_coalescent"))})

setGeneric(name="myPlot",
           def=function(coalescent,transitionMod){return(standardGeneric("myPlot"))})

setGeneric(name="tipsInOrder",
           def=function(tree){return(standardGeneric("tipsInOrder"))})

setGeneric(name="modifyBranch",
          def=function(coalescent,transitionMod){return("modifyBranch")})

setGeneric(name="coalescent_2_newick",
           def=function(coalescent){return("coalescent_2_newick")})

setGeneric(name="coalescent_2_phylog",
           def=function(coalescent){return("coalescent_2_phylog")})

setGeneric(name="laplaceMatrix",
           def=function(transitionMod){return("laplaceMatrix")})

setGeneric(name="ordinary_laplacian",
           def=function(transitionMod){return("ordinary_laplacian")})

setGeneric(name="commute_time_undigraph",
           def=function(transitionMod){return("commute_time_undigraph")})

setGeneric(name="commute_time_digraph",
           def=function(transitionMod){return("commute_time_digraph")})

setGeneric(name="hitting_time_digraph",
           def=function(transitionMod){return("hitting_time_digraph")})

commute_time_undigraph
           