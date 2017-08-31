# Copyright (C) 2016 VRVis.
# All rights reserved.
# Contact: VRVis Forschungs-GmbH (office@vrvis.at)
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. All advertising materials mentioning features or use of this software
#    must display the following acknowledgement:
#    This product includes software developed by the VRVis Forschungs-GmbH.
# 4. Neither the name of the VRVis Forschungs-GmbH nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.

# This script is for maximizing contrast of gene expression synergy for gene sets (sets need to be computed by get_gene_expression_for_genesets.m)

# The code for the rbga.int function is from the R package genalg (general code) and gramEvol (genalg extended for integer chromosome) adapted for parallel computing
# genalg: https://github.com/egonw/genalg by Egon Willighagen and Michel Ballings
# gramEvol: https://github.com/fnoorian/gramEvol/ by Farzad Noorian and Anthony Mihirana

useCores<-20                                  #cores that will be used for the genetic algorithm
#If one chooses 1, no parallel computing will be used

setwd("/home/ubuntu/") 

testSetName<-"storage//test_genesets"                               #Name of the testset (without file ending) 
#(what will be stored by get_gene_expression_for_genesets)

resultsFolder<-"storage//results//test_genesets_optimized"   #Folder where the results will be stored

#create results folder
dir.create("storage//results")
dir.create(resultsFolder)


library(R.matlab)
library(rjson)
library(circlize)
library(foreach)
library(doParallel)
library(gplots)
library(png)
library(iterators)
library(fields)
library(igraph)
library(fastcluster)
library(clue)
library(GA)
library(R.utils)

`%op%` <- if (useCores>1) `%dopar%` else `%do%`

#ontology file (gives information for every region in the Allen brain atlas)
document <- fromJSON(file="storage//ontology.json", method='C')

#gets acronym of a region given by a region id (from ontology file)
getAcronymByID <- function(child,ID){
  if(!is.null(child)){
    if(child$id==ID){
      return(child)
    }
  }
  
  for(actchild in child$children){
    res<-getAcronymByID(actchild,ID)
    if(!is.null(res)){
      if(res$id==ID){
        return(res)
      }
    }
  }
  return(NULL)
}

#gives a vector with the size of the atlas (3D volume rescaled to 1D vector) where every position
#is zero, execpt for positions where the region (given by the ID) can be found
getAtlasRegionsOfID <- function(atlasRegions,ID){
  newAtlasRegions<-atlasRegions==ID
  
  childrens<-getAcronymByID(document$msg[[1]],ID)$children
  while(length(childrens)>0){
    newChildren<-c()
    for(i in 1:length(childrens)){
      newAtlasRegions<-newAtlasRegions|(atlasRegions==childrens[[i]]$id)
      newChildren<-c(newChildren,getAcronymByID(document$msg[[1]],childrens[[i]]$id)$children)
    }
    
    childrens<-newChildren
  }
  return(newAtlasRegions)
}

#Genetic algorithm for optimizing sets according to a criteria 
#The code is from the R package genalg (general code) and gramEvol (genalg extended for integer chromosome) adapted for parallel computing
#genalg: https://github.com/egonw/genalg by Egon Willighagen and Michel Ballings
#gramEvol: https://github.com/fnoorian/gramEvol/ by Farzad Noorian and Anthony Mihirana
rbga.int <- function (ints=c(0,1), size = 10, suggestions = NULL, popSize = 200, iters = 100,
                      mutationChance = NA, elitism = NA, monitorFunc = NULL, ints.probs=rep(1/length(ints),length(ints)),
                      evalFunc = NULL, showSettings = FALSE, verbose = FALSE,useCores=1,dataMatrix=NULL,equalsIterationsForConvergence=10) {
  `%op%` <- if (useCores>1) `%dopar%` else `%do%`
  newIters<-0
  if (is.null(evalFunc)) {
    stop("A evaluation function must be provided. See the evalFunc parameter.")
  }
  vars = size
  if (is.na(mutationChance)) {
    mutationChance = 1/(vars + 1)
  }
  if (is.na(elitism)) {
    elitism = floor(popSize/5)
  }
  if (verbose)
    cat("Testing the sanity of parameters...\n")
  if (popSize < 5) {
    stop("The population size must be at least 5.")
  }
  if (iters < 1) {
    stop("The number of iterations must be at least 1.")
  }
  if (!(elitism < popSize)) {
    stop("The population size must be greater than the elitism.")
  }
  if (length(ints)<2) {
    stop("The number of integers must be at least 2.")
  }
  if (showSettings) {
    if (verbose)
      cat("The start conditions:\n")
    result = list(size = size, suggestions = suggestions,
                  popSize = popSize, iters = iters, elitism = elitism,
                  mutationChance = mutationChance)
    class(result) = "rbga"
    cat(summary(result))
  } else {
    if (verbose)
      cat("Not showing GA settings...\n")
  }
  if (vars > 0) {
    if (!is.null(suggestions)) {
      if (verbose)
        cat("Adding suggestions to first population...\n")
      population = matrix(nrow = popSize, ncol = vars)
      suggestionCount = dim(suggestions)[1]
      for (i in 1:suggestionCount) {
        population[i, ] = suggestions[i, ]
      }
      if (verbose)
        cat("Filling others with random values in the given domains...\n")
      for (child in (suggestionCount + 1):popSize) {
        population[child, ] = sample(ints, vars, rep = TRUE, prob=ints.probs)
      }
    } else {
      if (verbose)
        cat("Starting with random values in the given domains...\n")
      population = matrix(nrow = popSize, ncol = vars)
      for (child in 1:popSize) {
        population[child, ] = sample(ints, vars, rep = TRUE, prob=ints.probs)
      }
    }
    bestEvals = rep(NA, iters)
    meanEvals = rep(NA, iters)
    evalVals = rep(NA, popSize)
    for (iter in 1:iters) {
      if (verbose)
        cat(paste("Starting iteration", iter, "\n"))
      if (verbose)
        cat("Calucating evaluation values... ")
      
      clComputing <- c()
      if(useCores>1){
        clComputing <- makeCluster(20,type="PSOCK")
        registerDoParallel(clComputing)
      }
      evalVals<-foreach(allPop=1:popSize, i=icount(),.combine="c", .inorder=TRUE,.multicombine=TRUE,.export=c("geneExpressionSynergyOrg","expressionIndex","expressionOfGenes","trimExpression")) %op%{
        return(evalFunc(population[i,]))
      }
      if(useCores>1){
        stopCluster(clComputing)
      }
      
      bestEvals[iter] = min(evalVals)
      meanEvals[iter] = mean(evalVals)
      if (verbose)
        cat(" done.\n")
      if (!is.null(monitorFunc)) {
        if (verbose)
          cat("Sending current state to rgba.monitor()...\n")
        result = list(type = "integer chromosome", size = size,
                      popSize = popSize, iter = iter, iters = iters,
                      population = population, elitism = elitism,
                      mutationChance = mutationChance, evaluations = evalVals,
                      best = bestEvals, mean = meanEvals)
        class(result) = "rbga"
        monitorFunc(result)
      }
      
      newIters<-newIters+1
      
      if(iter>equalsIterationsForConvergence){
        if(length(unique(bestEvals[(iter-equalsIterationsForConvergence+1):iter]))==1){
          print("Converged!")
          break;
        }
      }
      
      if (iter < iters) {
        if (verbose)
          cat("Creating next generation...\n")
        newPopulation = matrix(nrow = popSize, ncol = vars)
        newEvalVals = rep(NA, popSize)
        if (verbose)
          cat("  sorting results...\n")
        sortedEvaluations = sort(evalVals, index = TRUE)
        sortedPopulation = matrix(population[sortedEvaluations$ix,], ncol = vars)
        if (elitism > 0) {
          if (verbose)
            cat("  applying elitism...\n")
          newPopulation[1:elitism, ] = sortedPopulation[1:elitism,]
          newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
        }
        if (vars > 1) {
          if (verbose)
            cat("  applying crossover...\n")
          for (child in (elitism + 1):popSize) {
            parentProb = dnorm(1:popSize, mean = 0, sd = (popSize/3))
            parentIDs = sample(1:popSize, 2, prob = parentProb)
            parents = sortedPopulation[parentIDs, ]
            # Crossover probability???
            crossOverPoint = sample(0:vars, 1)
            if (crossOverPoint == 0) {
              newPopulation[child, ] = parents[2, ]
              newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
            } else if (crossOverPoint == vars) {
              newPopulation[child, ] = parents[1, ]
              newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
            } else {
              newPopulation[child, ] = c(parents[1, ][1:crossOverPoint],
                                         parents[2, ][(crossOverPoint + 1):vars])
            }
          }
        }
        else {
          if (verbose)
            cat("  cannot crossover (#vars=1), using new randoms...\n")
          newPopulation[(elitism + 1):popSize, ] = sortedPopulation[sample(1:popSize,
                                                                           popSize - elitism), ]
        }
        population = newPopulation
        evalVals = newEvalVals
        
        if (mutationChance > 0) {
          if (verbose)
            cat("  applying mutations... ")
          mutatedCells=which(runif((popSize-elitism)*vars)<mutationChance)
          m_rows=elitism+ceiling(mutatedCells/vars)
          m_cols=mutatedCells%%vars
          m_cols=ifelse(m_cols==0,vars,m_cols)
          mutationCount=length(mutatedCells)
          for (muts in 1:mutationCount) {
            if (length(ints)==2) {
              population[m_rows[muts], m_cols[muts]]=ints[ints!=population[m_rows[muts], m_cols[muts]]]
            } else {
              population[m_rows[muts], m_cols[muts]] = sample(ints[ints!=population[m_rows[muts], m_cols[muts]]], 1, prob=ints.probs[ints!=population[m_rows[muts], m_cols[muts]]])
            }
          }

          if (verbose)
            cat(paste(mutationCount, "mutations applied\n"))
        }
      }
      
    }
  }
  result = list(type = "integer chromosome", size = size, popSize = popSize,
                iters = newIters, suggestions = suggestions, population = population,
                elitism = elitism, mutationChance = mutationChance, evaluations = evalVals,
                best = bestEvals, mean = meanEvals)
  class(result) = "rbga"
  return(result)
}



#############################START OF THE METHOD#########################################


print("Load atlas")
atlasRegions<-readMat("storage//atlasRegions.mat")

#atlasRegions is a 3D volume reformated to a 1D vector
#this is also done for gene expression volumes. In this script
#we do not work with 3D volumes, since it's easier in R to work 
#with vectors. The variable index maps [x,y,z] coordinates to the
#1D vector 
index<-atlasRegions$index
atlasRegions<-atlasRegions$atlasRegions

#Remove fiber tracts (and subregions) from brain
atlasRegions[getAtlasRegionsOfID(atlasRegions,8)==0]<-0 
atlasRegions[atlasRegions==8]<-0

#Basically all atlas regions (everything that is within the brain)
atlasRegionsBiggerZero<-atlasRegions[atlasRegions>0]

indexX=c()
indexY=c()
indexZ=c()

for(x in 1:(dim(index)[1])){
  for(y in 1:(dim(index)[2])){
    for(z in 1:(dim(index)[3])){
      indexX[index[x,y,z]]=x
      indexY[index[x,y,z]]=y
      indexZ[index[x,y,z]]=z
    }
  }
}

mY<-max(indexY)
mZ<-max(indexZ)



amountOfSets<-length(setnames)   #Amount of gene-sets in testSetName.

#load all relevant data from all gene-expression sets and computes gene-expression synergies
print("Load gene-expression sets and computes gene-expression synergies")
start <- Sys.time ()
for(actSet in 1:amountOfSets){
  fname = setnames[[actSet]][[1]][1]
  print(paste0("",Sys.time()," ","Load: ",fname))
  
  conData<-readMat(paste0(testSetName,"/",fname,".mat"))
  
  expressionIndex<-conData$expressionIndex
  
  expressionOfGenes<-as.matrix(conData$expressionMatrix)
  expressionOfGenes<-as.matrix(expressionOfGenes[,expressionIndex])
  
  expressionOfGenes[expressionOfGenes==Inf]<-NA
  expressionOfGenes[expressionOfGenes==-Inf]<-NA
  
  #Stanardize every gene by its mean and standard deviation
  expressionOfGenes<-apply(expressionOfGenes,2,function(x){
    y<-x[atlasRegions>0]
    y[is.na(y)]<-0
    (x-mean(y))/sd(y)
  })
  
  #gene expression synergy of the original set
  repCount<-min((1:ceiling(1/trimExpression))[(length(expressionIndex)*(1:ceiling(1/trimExpression)))%%ceiling(1/trimExpression)==0])
  geneExpressionSynergyOrg<-apply(expressionOfGenes[,expressionIndex],1,function(x){mean(rep(x, repCount),trim=trimExpression,na.rm=TRUE)})
  
  GAmodel = rbga.int(size=length(expressionIndex), ints=expressionIndex, popSize=500, iters=500, evalFunc = function(y) {
    if(length(unique(y))<=1){ #y is the subset
      return(Inf)
    }
    repCount<-min((1:ceiling(1/trimExpression))[(length(unique(y))*(1:ceiling(1/trimExpression)))%%ceiling(1/trimExpression)==0])
    synergy<-apply(expressionOfGenes[,unique(y)],1,function(x){mean(rep(x, repCount),trim=trimExpression,na.rm=TRUE)}) #gene expression synergy of the subset
    
    if(cor(geneExpressionSynergyOrg,synergy,method="pearson",use="pairwise.complete.obs")<0.95){
      return(Inf)
    }
    
    return(-mad(synergy,na.rm=TRUE)*sqrt(length(unique(y)))) #GA will optimize towards the smallest value, therefore we take the negative corrected MAD
  },monitorFunc=function(x){ 
    if(x$iter==1 || (x$iter%%10)==0){
      lengthOfBestSet<-0
      if(sum(x$evaluations==min(x$evaluations,na.rm=TRUE))==1){
        lengthOfBestSet<-length(unique(x$population[x$evaluations==min(x$evaluations),]))
        print(paste0(x$iter,': ',-x$best[x$iter]," length: ",lengthOfBestSet))
      }else{
        lengthOfBestSet<-min(apply(x$population[x$evaluations==min(x$evaluations),],1,function(y){length(unique(y))}))
        print(paste0(x$iter,': ',-x$best[x$iter]," length: ",lengthOfBestSet))
      }
      
    }
  },useCores=useCores,equalsIterationsForConvergence=5)
  
  lengthOfBestSet<-0 # Get the best set from the output of GA
  newExpressionIndex<-c()
  if(sum(GAmodel$evaluations[!is.na(GAmodel$evaluations)]==min(GAmodel$evaluations,na.rm=TRUE))==1){
    newExpressionIndex<-sort(unique(GAmodel$population[!is.na(GAmodel$evaluations)][GAmodel$evaluations[!is.na(GAmodel$evaluations)]==min(GAmodel$evaluations,na.rm=TRUE)]))
    lengthOfBestSet<-length(newExpressionIndex)
  }else{
    bestSets<-GAmodel$population[!is.na(GAmodel$evaluations),][GAmodel$evaluations[!is.na(GAmodel$evaluations)]==min(GAmodel$evaluations,na.rm=TRUE),]
    
    lengthOfBestSet<-min(apply(bestSets,1,function(y){length(unique(y))}))
    newExpressionIndex<-sort(unique(as.matrix(t(bestSets)[,lengthOfBestSet==apply(bestSets,1,function(y){length(unique(y))})])[,1]))
  }
  
  
  if(min(GAmodel$evaluations,na.rm=TRUE)<Inf){ #if a better subset is found, take it
    print(paste0("Finished after ",GAmodel$iters,' iterations: ',-GAmodel$best[GAmodel$iters]," length: ",lengthOfBestSet))
    expressionIndex<-newExpressionIndex
  }else{
    #if there is no better subset, take the original one
    print(paste0("Finished after ",GAmodel$iters,' iterations: ',-GAmodel$best[GAmodel$iters]," length: ",length(expressionIndex)))
  }
  
  write.table(data.frame(expressionindex=expressionIndex,genenames=unlist(conData$genenames[expressionIndex]),entrez=unlist(conData$entrez[expressionIndex])),file=paste0(resultsFolder,'//',"GA_sets//",fname,".csv"),sep = ",",row.names = FALSE)
}
print(paste0("Finished after: ",difftime(Sys.time(),start,units="mins"),"min"))



