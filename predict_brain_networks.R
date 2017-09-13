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

useCores<-20                                  #cores that will be used for p-value computation and clustering. 
#If one chooses 1, no parallel computing will be used

precision<-1000                                  #Amount of random gene-sets (synergies) that will be drawn for p-value computation

splitInBatches<-10                                #Random gene-sets (synergies) will be stored in one matrix with "precision" cols.
#This can be a problem regarding the size of the memory, therefore
#you can split the computation into batches. Batchsize = precision/splitInBatches
#The batchsize shouldn't be bigger than 1000

setwd("/home/ubuntu/") 

testSetName<-"storage//test_genesets"                               #Name of the testset (without file ending) 
#(what will be stored by get_gene_expression_for_genesets)

resultsFolder<-"storage//results//test_genesets"   #Folder where the results will be stored

#create results folder
dir.create("storage//results")
dir.create(resultsFolder)

minimumPValueToVisualize<-10^(-6)  #minimumPval is the negative Power of the minimum p-value that the color bar (for mri-style plots and cluster plots) provides

randomdata<-"storage//random_genes"                                 #Folder with random genes
#(what will be stored by get_gene_expression_for_random_genes)

amountOfRandomFiles<-10                                             #Amount of random-genes files in randomdata (not all of the files in random data
#need to be used!)

trimExpression<-0.05                                                #Trim-factor for trimmed-mean gene-expression calculation of a gene set
#0.05 --> The upper 5% and lower 5% percent gene-expression for every grid point
#will be filtered. 


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

#Resizes a 2D image
resizeImage = function(im, w.out, h.out) {
  
  # initial width/height
  w.in = nrow(im)
  h.in = ncol(im)
  
  # Create empty matrix
  im.out = matrix(rep(0,w.out*h.out), nrow =w.out, ncol=h.out )
  
  #Compute ratios 
  w_ratio = w.in/w.out
  h_ratio = h.in/h.out
  
  #resizing 
  im.out <- im[ floor(w_ratio* 1:w.out), floor(h_ratio* 1:h.out)]
  
  return(im.out)
}

#plots (stores a png to path_to_file) of spatial p-values (plotPvals) into several slices using the jet-color palette
#For a better overview, instead of viewing al slices in z-Direction, 5 slices each are interpolated
#to one slice. A color bar will be added, that shows significance at fdr005 and fdr01.
plot_mri_style_image<-function(path_to_file,plotPvals,fdr005,fdr01){
  
  #volume that will be plotted 
  m<-array(0, dim(index))
  plotPvals[is.na(plotPvals)]<-1
  
  
  plotPvals[plotPvals<minimumPValueToVisualize]<-minimumPValueToVisualize
  m[atlasRegions>0]<--log10(plotPvals)
  
  if(length(unique(m))>1){
    png(paste0(path_to_file),width=3000,height=2000)  
    par(mfrow=c(4,3),mar = rep(0.1, 4)) 
    m[is.nan(m)]<-0
    m[is.na(m)]<-0
    max(m)
    m[m<0]<-0
    
    
    m<-m*100
    m[atlasRegions<=0]<--10
    #builds multiple slice view
    for (z in 1:11){
      im<-apply(m[,,(5*z-2):(5*z+2)],c(1,2),max)
      im<-im[,dim(im)[2]:1]
      
      
      rb<-jet.colors(-log10(minimumPValueToVisualize)*100*100)
      rb<-rb[1:round(max(im)/(-log10(minimumPValueToVisualize)*100)*length(rb))]
      
      im<-resizeImage(im,dim(im)[1]*2,dim(im)[2]*2)
      im<-image.smooth(im, theta=1)$z  #smoothes image
      image(im, xaxt= "n", yaxt= "n",col=c("#000000",rb) )
      
    }
    image.scale <- function(z, zlim, col = heat.colors(12),
                            breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
      if(!missing(breaks)){
        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
      }
      if(missing(breaks) & !missing(zlim)){
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
      }
      if(missing(breaks) & missing(zlim)){
        zlim <- range(z, na.rm=TRUE)
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
      }
      poly <- vector(mode="list", length(col))
      for(i in seq(poly)){
        poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
      }
      xaxt <- ifelse(horiz, "s", "n")
      yaxt <- ifelse(horiz, "n", "s")
      if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
      if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
      if(missing(xlim)) xlim=XLIM
      if(missing(ylim)) ylim=YLIM
      plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i",cex.axis=8, ...)  
      for(i in seq(poly)){
        if(horiz){
          polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
        }
        if(!horiz){
          polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
        }
      }
      if(!(fdr005==0)){
        lines(c(-log10(fdr005),-log10(fdr005)),c(0,1),lty=3,lwd=5)
      }
      if(!(fdr01==0)){
        lines(c(-log10(fdr01),-log10(fdr01)),c(0,1),lwd=5)
      }
    }
    
    par(mgp=c(5,6,0))
    par(mar=c(15,5,15,5))
    image.scale(m/100,col=rb)
    dev.off()
    
  }
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


#indizes of the injection volumes and its IDs to load the 
#from the hard-drive (generated with get_connectivity_data())
conData<-readMat("storage//injectedIndizes.mat")


injectedIDs<-unique(conData$injectedIDs) #IDs (IDs correspond to the filenames of the connectivity data)  of the injection sites
injectedIDsIndex<-c()
injectedIDsIndex[injectedIDs]<-1:length(injectedIDs) #Position of an ID in injectedIDs

injectedIndizesToID<-as.matrix(conData$injectedIndizesToID) #maps indizes (for the 1D vector space) to ID
injectedIndizes<-unique(injectedIndizesToID[,1])     #all indizes of injection sites (all grid points that are parts of injection sites)

remove(conData)

#injectedIndizes are the indizes within the complete 1D vector. injectedIndizesWithinBrain are
#the inizies of the injection sites within a 1D vector that only contains voxels that are within the brain (so without background)
injectedIndizesWithinBrain<-rep(0,length(atlasRegions))  
injectedIndizesWithinBrain[atlasRegions>0]<-1:sum(atlasRegions>0)
injectedIndizesWithinBrain<-injectedIndizesWithinBrain[injectedIndizes[atlasRegions[injectedIndizes]>0]]

print(paste0("Amount of injections: ",length(unique(injectedIDs))))
print(paste0("Grid points that are parts of the injection sites: ",length(injectedIndizes)))


#Load the injection data into the memory. Every file represents the targets of an injection site. Injection data is binarized with a threshold of 10^-4.5 
print("Load the injection data into the memory")
pb <- txtProgressBar(min=0, max=100, initial=0,style=3)
act=1;
injectionData=matrix(as.raw(0), nrow = length(injectedIDs), ncol = sum(atlasRegions>0))
for(injIndbG in injectedIDs){
  
  conInd<-readMat(paste0("storage//id_connectivity/",injIndbG,".mat"),sep="")
  
  conRow<-(conInd$con)[atlasRegions>0]
  conRow[conRow<0]<-0
  
  conRow[is.nan(conRow)]<-0
  conRow[conRow==Inf]<-0
  injectionData[act,]=as.raw(conRow>(10^-4.5))
  
  act=act+1;
  setTxtProgressBar(pb, act/(length(injectedIDs))*100)
  
}
close(pb)


#Build connectivity matrix. Every injection data represents the targets of an injection site, but since an injection site consists
#of multiple grid points, some lines of injection data must be present multiple times in connectivityMatrix           
print("Build connectivity matrix")
pb <- txtProgressBar(min=0, max=100, initial=0,style=3)
act=1;
connectivityMatrix=matrix(as.raw(0), nrow = length(injectedIndizes[atlasRegions[injectedIndizes]>0]), ncol = sum(atlasRegions>0))
for(injIndbG in injectedIndizes[atlasRegions[injectedIndizes]>0]){
  idSofIndex<-unique(injectedIndizesToID[injectedIndizesToID[,1]==injIndbG,2]) #idSofIndex is the file id number of injection site given by injIndbG
  if(length(idSofIndex)>1){
    conRow<-colSums(injectionData[injectedIDsIndex[idSofIndex],]>0)>0
  }else{
    conRow<-injectionData[injectedIDsIndex[idSofIndex],]>0
  }
  connectivityMatrix[act,]=as.raw(conRow>0)
  
  act=act+1;
  setTxtProgressBar(pb, act/(length(injectedIndizes[atlasRegions[injectedIndizes]>0]))*100)
  
}
close(pb)

print(paste0("Dimensions of connectivity matrix: ",dim(connectivityMatrix)[1],"x",dim(connectivityMatrix)[2]))
remove(injectionData)
gc()


#gives a subset of connectivityMatrix
getSubsetOfConnectivity<-function(from,to){
  tryCatch({
    return(connectivityMatrix[from,to])
  }, error = function(e) {
    if(length(from)<=1){
      return(connectivityMatrix[from,][to])
    }else{
      return(connectivityMatrix[from,][,to])
    }
  })
}

#Load random gene expression data (will be used fro random gene-sets computation)
expressionOfRandomGenes<-c()
for(actRand in 1:amountOfRandomFiles){
  print(paste0("Load randomdata ",actRand,"/",amountOfRandomFiles))
  randData<-readMat(paste0(randomdata,"/random_genes_",actRand,".mat"))
  hexpressionOfRandomGenes<-randData$expressionMatrix
  hexpressionOfRandomGenes<-as.matrix(hexpressionOfRandomGenes)
  hexpressionOfRandomGenes[atlasRegions<=0,]<-0
  hexpressionOfRandomGenes[hexpressionOfRandomGenes==-Inf]<-NA
  hexpressionOfRandomGenes[hexpressionOfRandomGenes==Inf]<-NA
  remove(randData)
  expressionOfRandomGenes<-cbind(expressionOfRandomGenes,hexpressionOfRandomGenes[atlasRegions>0,])
}
remove(hexpressionOfRandomGenes)

print(paste0("Check distribution mean (should be close to 0): ",mean(apply(expressionOfRandomGenes,1,function(x){mean(x,na.rm=TRUE)}))))
print(paste0("Check distribution std (should be close to 1): ",mean(apply(expressionOfRandomGenes,1,function(x){sd(x,na.rm=TRUE)}))))


#Stanardize every random gene by its mean and standard deviation
expressionOfRandomGenes<-apply(expressionOfRandomGenes,2,function(x){
  y<-x
  y[is.na(y)]<-0
  (x-mean(y))/sd(y)
})


meanExprs<-apply(expressionOfRandomGenes,1,function(x){mean(x,na.rm=TRUE)}) #mean gene expression for every grid point
sdExprs<-apply(expressionOfRandomGenes,1,function(x){sd(x,na.rm=TRUE)})     #standard deviation of gene expression for every grid point



amountOfSets<-length(setnames)   #Amount of gene-sets in testSetName.

geneExpressionSynergyMatrix<-matrix(0,nrow=sum(atlasRegions>0),ncol = amountOfSets) #stores gene-expression synergy (grid-point-wise trimmed-mean gene-expression)
ratiosMatrix<-list() #every gene in the gene-expression set can have a weighting ratio that gives them more weight in the gene-expression synergy , usually 1 for all
amountOfGenesVector<-c()
genenamesList<-list()
geneExressionMatrix<-list()

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
  
  ratios<-conData$ratios
  if(length(ratios)==0){
    ratios<-rep(1,length(expressionIndex))
  }else{
    ratios<-ratios[expressionIndex]
  }
  
  genenames<-conData$genenames[expressionIndex][ratios>0]
  genenamesList[[actSet]]<-genenames
  ratios<-ratios[ratios>0]
  
  ratiosMatrix[[actSet]]<-ratios
  
  expressionOfGenes[atlasRegions<=0,]<-0
  
  
  expressionOfGenes<-t(t(expressionOfGenes)*ratios)
  amountOfGenes<-dim(expressionOfGenes)[2]
  
  
  geneExressionMatrix[[actSet]]<-expressionOfGenes
  amountOfGenesVector[actSet]<-amountOfGenes
  
  #calculate gene expression synergy (trimmed mean)
  if(amountOfGenes>1){
    geneExpressionSynergy<-apply(expressionOfGenes,1,function(x){mean(rep(x, ceiling(1/(trimExpression))),trim=trimExpression,na.rm=TRUE)})
    #the repetition (rep) of genes is necessary, since a 5% lower and upper cutoff would only work for at least 20 genes. 
    #example: one could not compute the 10% trimmed mean of 5 numbers [0,1,2,3,10], but taking those 5 numbers 2 times [0,1,2,3,10,0,1,2,3,10] 
    #its possible (it doesn't trim the "10" perfectly, but still better than without)
  }else{
    geneExpressionSynergy<-expressionOfGenes
  }
  
  
  geneExpressionSynergy<-geneExpressionSynergy[atlasRegions>0]
  
  geneExpressionSynergy[is.na(geneExpressionSynergy)]<-0
  
  geneExpressionSynergyMatrix[,actSet]<-geneExpressionSynergy
  
}
print(paste0("Finished after: ",difftime(Sys.time(),start,units="mins"),"min"))


#Calculate incoming and outgoing node-strength

inBetweenStart<-Sys.time() 

chunkSize<-min(1000,ceiling(dim(connectivityMatrix)[2]/useCores)) #splits the connectivityMatrix into chunks column wise, to fit the computation into the memory
print("Calculate incoming node-strength for all gene-sets")

pb <- txtProgressBar(min=0, max=100, initial=0,style=3)
IncomingStrengthOfWeightedExpressionMatrices<-foreach(conMatSubset=iter(connectivityMatrix, by = "col",chunksize =  chunkSize), i=icount(),.combine="cbind",.multicombine=TRUE, .inorder=TRUE) %do%{
  setTxtProgressBar(pb, ceiling(i/(dim(connectivityMatrix)[2]/chunkSize)*100))
  
  return((crossprod(geneExpressionSynergyMatrix[injectedIndizesWithinBrain,],conMatSubset>0)))
}
close(pb)

print(paste0("Finished after: ",difftime(Sys.time(),inBetweenStart,units="mins"),"min"))
inBetweenStart<-Sys.time() 
print("Calculate outgoing node-strength for all gene-sets")

pb <- txtProgressBar(min=0, max=100, initial=0,style=3)
OutgoingStrengthOfWeightedExpressionMatrices<-foreach(conMatSubset=iter(connectivityMatrix, by = "row",chunksize =  ceiling(dim(connectivityMatrix)[1]/dim(connectivityMatrix)[2]*chunkSize)), i=icount(),.combine="cbind",.multicombine=TRUE, .inorder=TRUE) %do%{
  setTxtProgressBar(pb, ceiling(i/(dim(connectivityMatrix)[2]/chunkSize)*100))
  return((crossprod(geneExpressionSynergyMatrix,t(conMatSubset)>0)))
}
close(pb)

print(paste0(difftime(Sys.time(),inBetweenStart,units="mins"),"min"))



for(actSet in 1:amountOfSets){
  fname = setnames[[actSet]][[1]][1]
  print(paste0("",Sys.time()," ","Start pvalue calculation for: ",fname," (",actSet,"/",amountOfSets,")"))
  
  
  dir.create(paste0(resultsFolder,'//',paste0(fname)))
  dir.create(paste0(resultsFolder,'//',"pvalues"))
  
  start <- Sys.time () 
  
  #get current gene set (actSet) specific data 
  geneExpressionSynergy<-geneExpressionSynergyMatrix[,actSet]
  IncomingStrengthOfWeightedExpressionMatrix<-IncomingStrengthOfWeightedExpressionMatrices[actSet,]
  OutgoingStrengthOfWeightedExpressionMatrix<-OutgoingStrengthOfWeightedExpressionMatrices[actSet,]
  expressionOfGenes<-geneExressionMatrix[[actSet]]
  amountOfGenes<-amountOfGenesVector[actSet]
  genenames<-genenamesList[[actSet]]
  ratios<-ratiosMatrix[[actSet]]
  
  #if p-values havent been computed yet, compute them now
  if(!file.exists(paste0(resultsFolder,'//','pvalues/',fname,'_i_ttest.rds')) | !file.exists(paste0(resultsFolder,'//','pvalues/',fname,'_ttest_geneExpr.rds'))){
    
    
    #those are the mean/standard deviation of incoming/outgoing node strength and gene-expression, computed
    #from random-gene set synergies
    IncomingStrengthMeans<-rep(0,sum(atlasRegions>0)) 
    IncomingStrengthSDs<-rep(0,sum(atlasRegions>0))
    OutgoingStrengthMeans<-rep(0,length(injectedIndizesWithinBrain))
    OutgoingStrengthSDs<-rep(0,length(injectedIndizesWithinBrain))
    geneExpressionMeans<-rep(0,sum(atlasRegions>0))
    geneExpressionSDs<-rep(0,sum(atlasRegions>0))
    
    
    computationStart<-Sys.time()
    for(actBatch in 1:splitInBatches){
      
      inBetweenStart<-Sys.time() 
      
      
      chunkSize<-min(1000,ceiling(dim(connectivityMatrix)[2]/useCores)) #splits the connectivityMatrix into chunks column wise, to fit the computation into the memory
      
      if(useCores>1){
        print("Setup parallel computing...")
        clComputing <- makeCluster(useCores,outfile="parallel_log.txt",type="PSOCK")
        registerDoParallel(clComputing)
        print("Parallel computing setup done")
      }
      
      
      print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5)),": Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Generate ",ceiling(precision/splitInBatches)," random gene-sets..."))
      
      #Generates an amount (precision/splitInBatches) of random gene-sets and computes the gene expression synergy (trimmed-mean)
      expressionOfRandomGenesRowSumsMatrix<-foreach(o=1:ceiling(precision/splitInBatches),.combine="cbind",.multicombine=TRUE, .inorder=FALSE) %op%{
        randnum<-sample(1:dim(expressionOfRandomGenes)[2],amountOfGenes)
        expressionOfRandomGenesRowSums<-0
        if(amountOfGenes>1){
          expressionOfRandomGenesRowSums<-apply(t(expressionOfRandomGenes[,randnum])*ratios,2,function(x){mean(rep(x,ceiling(1/(trimExpression))),trim=trimExpression,na.rm=TRUE)})
          #the repetition (rep) of genes is necessary, since a 5% lower and upper cutoff would only work for at least 20 genes. 
          #example: one could not compute the 10% trimmed mean of 5 numbers [0,1,2,3,10], but taking those 5 numbers 2 times [0,1,2,3,10,0,1,2,3,10] 
          #its possible (it doesn't trim the "10" perfectly, but still better than without)
        }else{
          expressionOfRandomGenesRowSums<-expressionOfRandomGenes[,randnum]
        }
        return(expressionOfRandomGenesRowSums)
      }
      print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5))," Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Random gene-sets generated after ",difftime(Sys.time(),inBetweenStart,units="mins"),"min"))
      
      
      if(useCores>1){
        print(paste0("",Sys.time()," ","Restart parallel computing (that's normal)"))
        stopCluster(clComputing) 
        Sys.sleep(10)
        gc()
        clComputing <- makeCluster(useCores,outfile="parallel_log.txt",type="PSOCK")
        registerDoParallel(clComputing)
        print("Parallel computing setup done")
      }
      
      inBetweenStart<-Sys.time() 
      print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5))," Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Node-strength calculation..."))
      
      geneExpressionMeans<-geneExpressionMeans+apply(expressionOfRandomGenesRowSumsMatrix,1,function(x){mean(x,na.rm=TRUE)})
      geneExpressionSDs<-geneExpressionSDs+apply(expressionOfRandomGenesRowSumsMatrix,1,function(x){sd(x,na.rm=TRUE)})
      
      expressionOfRandomGenesRowSumsMatrix[is.na(expressionOfRandomGenesRowSumsMatrix)]<-0
      
      print("-----calculate incoming node-strength for all random gene-sets")
      strength_MEANs_and_SDs<-foreach(conMatSubset=iter(connectivityMatrix, by = "col",chunksize =  chunkSize), i=icount(),.combine="c", .inorder=TRUE,.multicombine=TRUE) %op%{
        randomStrengthSample<- crossprod(expressionOfRandomGenesRowSumsMatrix[injectedIndizesWithinBrain,],conMatSubset>0)
        return(c(apply(randomStrengthSample,2,function(x){mean(x,na.rm=TRUE)}),apply(randomStrengthSample,2,function(x){sd(x,na.rm=TRUE)})))
      }
      
      mean_sd_indizes<-foreach(conMatSubset=split(1:dim(connectivityMatrix)[2],ceiling((1:dim(connectivityMatrix)[2])/chunkSize)),.combine="c", .inorder=TRUE)%do%{
        return(c(rep(1,length(conMatSubset)),rep(2,length(conMatSubset))))
      }
      
      IncomingStrengthMeans<-IncomingStrengthMeans+strength_MEANs_and_SDs[mean_sd_indizes==1]
      IncomingStrengthSDs<-IncomingStrengthSDs+strength_MEANs_and_SDs[mean_sd_indizes==2]
      
      
      print(paste0("",Sys.time()," ","...finished"))
      if(useCores>1){
        print(paste0("",Sys.time()," ","Restart parallel computing (that's normal)"))
        stopCluster(clComputing) 
        Sys.sleep(10)
        gc()
        clComputing <- makeCluster(useCores,outfile="parallel_log.txt",type="PSOCK")
        registerDoParallel(clComputing)
        print("Parallel computing setup done")
      }
      
      
      print("-----calculate outgoing node-strength for all random gene-sets")
      strength_MEANs_and_SDs<-foreach(conMatSubset=iter(connectivityMatrix, by = "row",chunksize =  ceiling(dim(connectivityMatrix)[1]/dim(connectivityMatrix)[2]*chunkSize)), i=icount(),.combine="c", .inorder=TRUE,.multicombine=TRUE) %op%{
        randomStrengthSample<- crossprod(expressionOfRandomGenesRowSumsMatrix,t(conMatSubset)>0)
        
        return(c(apply(randomStrengthSample,2,function(x){mean(x,na.rm=TRUE)}),apply(randomStrengthSample,2,function(x){sd(x,na.rm=TRUE)})))
      }
      
      mean_sd_indizes<-foreach(conMatSubset=split(1:dim(connectivityMatrix)[1],ceiling((1:dim(connectivityMatrix)[1])/ceiling(dim(connectivityMatrix)[1]/dim(connectivityMatrix)[2]*chunkSize))),.combine="c", .inorder=TRUE)%do%{
        return(c(rep(1,length(conMatSubset)),rep(2,length(conMatSubset))))
      }
      
      OutgoingStrengthMeans<-OutgoingStrengthMeans+strength_MEANs_and_SDs[mean_sd_indizes==1]
      OutgoingStrengthSDs<-OutgoingStrengthSDs+strength_MEANs_and_SDs[mean_sd_indizes==2]
      
      
      print(paste0("",Sys.time()," ","...finished"))
      
      if(useCores>1){
        print("Parallel computing finished")
        stopCluster(clComputing) 
      }
      
      rm(expressionOfRandomGenesRowSumsMatrix)
      rm(strength_MEANs_and_SDs)
      rm(mean_sd_indizes)
      gc()
      
      print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5))," Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Node-strengths computed after ",difftime(Sys.time(),inBetweenStart,units="mins"),"min"))
      
      
    }
    
    print("Calculate p-values...")
    pvals_Incoming_ttest<-(1-pnorm(IncomingStrengthOfWeightedExpressionMatrix,mean=IncomingStrengthMeans/splitInBatches,sd=IncomingStrengthSDs/splitInBatches))
    pvals_Outgoing_ttest<-rep(NA,sum(atlasRegions>0))
    pvals_Outgoing_ttest[injectedIndizesWithinBrain]<-(1-pnorm(OutgoingStrengthOfWeightedExpressionMatrix,mean=OutgoingStrengthMeans/splitInBatches,sd=OutgoingStrengthSDs/splitInBatches))
    
    pvalsExpr_ttest<-(1-pnorm((geneExpressionSynergy),mean=geneExpressionMeans/splitInBatches,sd=geneExpressionSDs/splitInBatches))
    
    zscores_Incoming<-(IncomingStrengthOfWeightedExpressionMatrix-(IncomingStrengthMeans/splitInBatches))/(IncomingStrengthSDs/splitInBatches)
    zscores_Outgoing<-rep(NA,sum(atlasRegions>0))
    zscores_Outgoing[injectedIndizesWithinBrain]<-(OutgoingStrengthOfWeightedExpressionMatrix-(OutgoingStrengthMeans/splitInBatches))/(OutgoingStrengthSDs/splitInBatches)
    zscores_Expr<-(geneExpressionSynergy-(geneExpressionMeans/splitInBatches))/(geneExpressionSDs/splitInBatches)
    
    
    pvalsExpr_ttest[geneExpressionMeans==0 | geneExpressionSDs==0]<-NA
    
    print("p-values calculated!")
    saveRDS(pvals_Incoming_ttest, file = paste0(resultsFolder,'//','pvalues/',fname,'_i_ttest.rds'))
    saveRDS(pvals_Outgoing_ttest, file = paste0(resultsFolder,'//','pvalues/',fname,'_o_ttest.rds'))
    saveRDS(pvalsExpr_ttest, file = paste0(resultsFolder,'//','pvalues/',fname,'_ttest_geneExpr.rds'))
    
    saveRDS(zscores_Incoming, file = paste0(resultsFolder,'//','pvalues/',fname,'_i_zscores.rds'))
    saveRDS(zscores_Outgoing, file = paste0(resultsFolder,'//','pvalues/',fname,'_o_zscores.rds'))
    saveRDS(zscores_Expr, file = paste0(resultsFolder,'//','pvalues/',fname,'_zscores_geneExpr.rds'))
    
  }
  
  #Load calculated p-values from hard-drive
  subfolder<-""
  print("TTest: ")
  pvalsIncoming<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_i_ttest.rds'))
  pvalsOutgoing<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_o_ttest.rds'))
  pvalsExpr<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_ttest_geneExpr.rds'))
  subfolder<-"//TTest"
  
  #Calculate FDR (q-values)
  qvalsIncoming<-rep(1,length(pvalsIncoming))
  qvalsIncoming[!is.na(pvalsIncoming)]<-p.adjust(pvalsIncoming[!is.na(pvalsIncoming)], method="BH")
  qvalsOutgoing<-rep(1,length(pvalsOutgoing))
  qvalsOutgoing[!is.na(pvalsOutgoing)]<-p.adjust(pvalsOutgoing[!is.na(pvalsOutgoing)], method="BH")
  qvalsExpr<-rep(1,length(pvalsExpr))
  qvalsExpr[!is.na(pvalsExpr)]<-p.adjust(pvalsExpr[!is.na(pvalsExpr)], method="BH")
  
  dir.create(paste0(resultsFolder,'//',fname,subfolder))
  
  pvalsIncoming[is.na(pvalsIncoming)]<-1
  qvalsIncoming[is.na(qvalsIncoming)]<-1
  pvalsOutgoing[is.na(pvalsOutgoing)]<-1
  qvalsOutgoing[is.na(qvalsOutgoing)]<-1
  pvalsExpr[is.na(pvalsExpr)]<-1
  qvalsExpr[is.na(qvalsExpr)]<-1
  
  print("Percentage of significant regions in the brain: ")
  print(paste0("First order FDR(0.05): ",sum(qvalsExpr<=0.05)/length(qvalsIncoming)))
  print(paste0("First order FDR(0.1): ",sum(qvalsExpr<=0.1)/length(qvalsIncoming)))
  
  print(paste0("Second orderIncoming FDR(0.05): ",sum(qvalsIncoming<=0.05)/length(qvalsIncoming)))
  print(paste0("Second orderIncoming FDR(0.1): ",sum(qvalsIncoming<=0.1)/length(qvalsIncoming)))
  
  print(paste0("Second orderOutgoing FDR(0.05): ",sum(qvalsOutgoing<=0.05)/length(qvalsOutgoing)))
  print(paste0("Second orderOutgoing FDR(0.1): ",sum(qvalsOutgoing<=0.1)/length(qvalsOutgoing)))
  
  print("Plotting results (MRI style plot)...")
  plot_mri_style_image(paste0(resultsFolder,'//',fname,subfolder,'/mri_style_pvals_combined_minimum_pvalue.png'),pmin(pvalsExpr,pvalsOutgoing,pvalsIncoming),0,0)
  plot_mri_style_image(paste0(resultsFolder,'//',fname,subfolder,'/mri_style_pvals_max10_incoming_secondOrder.png'),pvalsIncoming,max(pvalsIncoming[qvalsIncoming<=0.05]),max(pvalsIncoming[qvalsIncoming<=0.1]))
  plot_mri_style_image(paste0(resultsFolder,'//',fname,subfolder,'/mri_style_pvals_max10_firstOrder.png'),pvalsExpr,max(pvalsExpr[qvalsExpr<=0.05]),max(pvalsExpr[qvalsExpr<=0.1]))
  plot_mri_style_image(paste0(resultsFolder,'//',fname,subfolder,'/mri_style_pvals_max10_outgoing_secondOrder.png'),pvalsOutgoing,max(pvalsOutgoing[qvalsOutgoing<=0.05]),max(pvalsOutgoing[qvalsOutgoing<=0.1]))
  
  
  writeMat(paste0(resultsFolder,'//',fname,'//','pvalues.mat'),pvals_firstOrder=pvalsExpr,pvals_secondOrder_Incoming=pvalsIncoming,pvals_secondOrder_Outgoing=pvalsOutgoing,qvals_firstOrder=qvalsExpr,qvals_secondOrder_Incoming=qvalsIncoming,qvals_secondOrder_Outgoing=qvalsOutgoing)
  
  
  
  
  print("Plotting heatmap...")
  colorsOfRegions<-c()
  qvalThresh<-0.1
  
  expressionOfGenes<-cbind(geneExressionMatrix[[actSet]])
  amountOfGenes<-amountOfGenesVector[actSet]
  
  
  
  colsAndRows<-unique(atlasRegionsBiggerZero[qvalsIncoming<=qvalThresh | qvalsOutgoing<= qvalThresh | qvalsExpr<=qvalThresh])
  atlasRegionExpression<-matrix(0,nrow=dim(expressionOfGenes)[2]+3,ncol=length(colsAndRows))
  
  #Calculate the trimmed gene expression for every gene in the set
  #-->the gene expression for every gene that is left after trimming the set
  if(trimExpression>0 && amountOfGenes>1){
    expressionMatrixTrimmed<-t(apply(expressionOfGenes[atlasRegions>0,],1,function(x){
      nam<-1:length(x)
      repX<-rep(x,1/trimExpression)
      repNam<-rep(nam,1/trimExpression)
      ordX<-order(repX)
      repX<-repX[ordX[(length(ordX)*trimExpression+1):(length(ordX)*(1-trimExpression))]]
      repNam<-repNam[ordX[(length(ordX)*trimExpression+1):(length(ordX)*(1-trimExpression))]]
      
      newX<-c()
      for(i in nam){
        newX<-c(newX,sum(repX[repNam==i])/((1-(2*trimExpression))/trimExpression))
      }
      
      return(newX)
    }))
  }else{
    expressionMatrixTrimmed<-cbind(expressionOfGenes[atlasRegions>0,])
  }
  
  #P-values for trimmed gene-expression compared to all random genes
  for(k in 1:dim(expressionMatrixTrimmed)[1]){
    expressionMatrixTrimmed[k,]<-(((1-pnorm(expressionMatrixTrimmed[k,],mean=meanExprs[k],sd=sdExprs[k]))))
  }
  for(j in 1:dim(expressionMatrixTrimmed)[2]){
    expressionMatrixTrimmed[,j]<-(-log10(expressionMatrixTrimmed[,j]))
  }
  
  #Add additional rows to heatmap (First order, second order node strength/gene-expression)
  expressionMatrixTrimmed<-cbind(-log10(pvalsOutgoing),expressionMatrixTrimmed)
  expressionMatrixTrimmed<-cbind(-log10(pvalsIncoming),expressionMatrixTrimmed)
  expressionMatrixTrimmed<-cbind(-log10(pvalsExpr),expressionMatrixTrimmed)
  
  expressionMatrixTrimmed[is.na(expressionMatrixTrimmed)]<-0
  expressionMatrixTrimmed[is.nan(expressionMatrixTrimmed)]<-0
  expressionMatrixTrimmed[expressionMatrixTrimmed>min(2,max(expressionMatrixTrimmed))]<-min(2,max(expressionMatrixTrimmed))
  
  for(idindex in 1:length(colsAndRows)){
    tryCatch({
      
      atlasRegionExpression[,idindex]<-apply(expressionMatrixTrimmed[atlasRegionsBiggerZero==colsAndRows[idindex],],2,mean)
    }, error = function(e) {
      atlasRegionExpression[idindex]<-mean(expressionOfGenes[atlasRegions==colsAndRows[idindex]])
    })
    
    colorsOfRegions<-c(colorsOfRegions,paste0("#",getAcronymByID(document$msg[[1]],colsAndRows[idindex])$color_hex_triplet))
    
    colsAndRows[idindex]<-getAcronymByID(document$msg[[1]],colsAndRows[idindex])$acronym
  }
  
  colsAndRows<-colsAndRows[order(colorsOfRegions)]
  atlasRegionExpression<-atlasRegionExpression[,order(colorsOfRegions)]
  
  colnames(atlasRegionExpression)<-c(colsAndRows)
  rownames(atlasRegionExpression)<-capitalize(c("First Order","Second Ord. In","Second Ord. Out",unlist(genenames)))
  
  colorsOfRegions<-sort(colorsOfRegions)
  
  
  png(paste0(resultsFolder,'//',fname,subfolder,'//Heatmap.png'),width=5000,height=8000)  
  
  par(mfrow=c(1,1),mar = c(1,1,1,1))
  heatmap.2(t(atlasRegionExpression),dendrogram="none",key.xlab="0  -log10(pvalue) 2",key.xtickfun=function() {
    return(list(labels=FALSE, tick=FALSE))
  },key.title=NA,key.par = list(mar=c(10,40,120,1),cex.lab=7,cex.axis=0.01),scale="none",keysize = 1,col=colorpanel(100, "white", "orange", "red"),margins=c(70,20),cexRow=2,cexCol=10,Rowv=FALSE,Colv=FALSE,srtCol=45,RowSideColors=colorsOfRegions,density.info="none", trace="none")
  
  dev.off() 
  
  print("Write brain-region wise csv...")
  signRegions<-unique(atlasRegionsBiggerZero[qvalsIncoming<=qvalThresh | qvalsOutgoing<=qvalThresh | qvalsExpr<=qvalThresh])
  signRegions<-signRegions[!is.na(signRegions)]
  if(length(signRegions)>1){
    regions_via_genes<-data.frame(region=c(),regionShort=c(),regionName=c(),regionCol=c(),regionPvalFirstOrder=c(),regionPvalIncomingSecondOrder=c(),regionPvalOutgoingSecondOrder=c(),regionPvalCombined=c(), stringsAsFactors=FALSE)
    for(sR in signRegions){
      if(sR>0 ){
        brainRegion<-getAcronymByID(document$msg[[1]],sR)    
        combPVal<-min(c(min(pvalsExpr[atlasRegionsBiggerZero==sR]),min(pvalsIncoming[atlasRegionsBiggerZero==sR]),min(pvalsOutgoing[atlasRegionsBiggerZero==sR])))
        regions_via_genes<-rbind(regions_via_genes,data.frame(sR,brainRegion$acronym,brainRegion$name,brainRegion$color_hex_triplet,min(pvalsExpr[atlasRegionsBiggerZero==sR]),min(pvalsIncoming[atlasRegionsBiggerZero==sR]),min(pvalsOutgoing[atlasRegionsBiggerZero==sR]),combPVal, stringsAsFactors=FALSE))
        
      }
      
    }
    
    suppressWarnings(colnames(regions_via_genes)<-c("region","regionShort","regionName","regionCol","regionPvalFirstOrder","regionPvalIncomingSecondOrder","regionPvalOutgoingSecondOrder","regionPvalCombined"))
    
    regions_via_genes<-unique(regions_via_genes[order(regions_via_genes$regionPvalCombined),])
    write.csv2(regions_via_genes, paste0(resultsFolder,'//',fname,subfolder,'//Regions.csv'))
    
    
  }
  
  
}


amountOfSets<-length(setnames)

print("")
print("Calculate clusters...")
print("")

#Generate a brain-mask for the right-hemisphere
#Plots of the brain will only contain p-values 
#from the left-hemisphere, since they would overlap
#in a z-projection (they are equal anyway, due to mirroring
#of the data)
atlasRegionsRight<-atlasRegions
for(x in 1:(dim(index)[1])){
  for(y in 1:(dim(index)[2])){
    for(z in 1:ceiling((dim(index)[3])/2)){
      atlasRegionsRight[index[x,y,z]]=0
    }
  }
}

getIndexOtherSideWithinBrain<-function(indicesPos){
  atlasRegionHelp<-atlasRegions<0
  atlasRegionHelp[indicesPos]<-1
  
  atlasRegionHelpMir<-atlasRegions<0
  
  atlasRegionHelpMir[index[1:(dim(index)[1]),1:(dim(index)[2]),(ceiling((dim(index)[3])/2)+1):(dim(index)[3])]]=atlasRegionHelp[index[1:(dim(index)[1]),1:(dim(index)[2]),((dim(index)[3])+1)-(ceiling((dim(index)[3])/2)+1):(dim(index)[3])]]
  return((1:sum(atlasRegions>0))[atlasRegionHelpMir[atlasRegions>0]>0])
}


for(actSet in 1:amountOfSets){ 
  
  
  fname<-setnames[[actSet]][[1]][1]
  print(paste0("Cluster data for set ",actSet,"/",amountOfSets,": ",fname))
  
  #Load calculated p-values from hard-drive
  pvalsIncoming<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_i_ttest.rds'))
  pvalsOutgoing<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_o_ttest.rds'))
  pvalsExpr<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_ttest_geneExpr.rds'))
  
  #Calculate FDR (q-values)
  qvalsIncoming<-rep(1,length(pvalsIncoming))
  qvalsIncoming[!is.na(pvalsIncoming)]<-p.adjust(pvalsIncoming[!is.na(pvalsIncoming)], method="BH")
  qvalsOutgoing<-rep(1,length(pvalsOutgoing))
  qvalsOutgoing[!is.na(pvalsOutgoing)]<-p.adjust(pvalsOutgoing[!is.na(pvalsOutgoing)], method="BH")
  qvalsExpr<-rep(1,length(pvalsExpr))
  qvalsExpr[!is.na(pvalsExpr)]<-p.adjust(pvalsExpr[!is.na(pvalsExpr)], method="BH")
  
  
  pvalsIncoming[is.na(pvalsIncoming)]<-1
  qvalsIncoming[is.na(qvalsIncoming)]<-1
  pvalsOutgoing[is.na(pvalsOutgoing)]<-1
  qvalsOutgoing[is.na(qvalsOutgoing)]<-1
  pvalsExpr[is.na(pvalsExpr)]<-1
  qvalsExpr[is.na(qvalsExpr)]<-1
  
  
  pvalsCombined<-pmin(pvalsIncoming,pvalsOutgoing,pvalsExpr)
  qvalsCombined<-pmin(qvalsIncoming,qvalsOutgoing,qvalsExpr)
  
  #remove p-values from the right-hemisphere, since we
  #want to plot only p-values from the left-hemisphere
  #(they would overlap in a z-projection and they are 
  #equal anyway, due to mirroring of the data))
  pvalsCombined[atlasRegionsRight[atlasRegions>0]>0]<-1
  qvalsCombined[atlasRegionsRight[atlasRegions>0]>0]<-1
  
  #Cluster regions with an FDR of 0.05. If that is
  #more than 25% of the brain, just plot the
  #most significant values (but maximal 25% of 
  #the brain). Otherwise it would be a bit to complex/chaotic
  #to visualize
  
  qvalThresh<-0.1
  if(sum(qvalsCombined<=qvalThresh)>(length(qvalsCombined)*0.25/2)){
    qvalThresh<-quantile(qvalsCombined,0.25/2)
    print(paste0("Bigger than 25%, Thresh=",qvalThresh))
  }
  #Cluster regions only if there are at least a few
  #significant grid point
  if(sum(qvalsCombined<=qvalThresh)>(length(qvalsCombined)*0.001)){   
    
    
    #indizes of injection sites (corresponds to rows of connectivityMatrix) 
    #that are on grid points with significant p-value
    fromIndex<-1:length(atlasRegions)
    fromIndex<-fromIndex[atlasRegions>0]
    fromIndex<-fromIndex[qvalsCombined<=qvalThresh]
    fromIndex<-intersect(fromIndex,injectedIndizes)    
    
    #indizes of targets (corresponds to cols of connectivityMatrix)   
    #that are on grid points with significant p-value
    toIndex<-1:length(atlasRegions)
    toIndex<-toIndex[atlasRegions>0]
    toIndex<-toIndex[qvalsCombined<=qvalThresh] 
    
    #indizes of injection sites of grid points that are within the brain
    fromIndexWithinBrain<-rep(NA,length(atlasRegions))
    fromIndexWithinBrain[injectedIndizes[atlasRegions[injectedIndizes]>0]]<-1:length(injectedIndizes[atlasRegions[injectedIndizes]>0])
    fromIndexWithinBrain<-fromIndexWithinBrain[fromIndex]
    
    #indizes of targets of grid points that are within the brain
    toIndexWithinBrain<-(1:sum(atlasRegions>0))
    toIndexWithinBrain<-toIndexWithinBrain[qvalsCombined<=qvalThresh]
    
    pvalsCombined<-pvalsCombined[qvalsCombined<=qvalThresh] 
    
    pvalsCombined[pvalsCombined<(10^-10)]<-(10^-10)
    
    atlasIndizes<-1:length(atlasRegions)
    atlasIndizes[TRUE]<-0
    atlasIndizes[atlasRegions>0]<-1:sum(atlasRegions>0)
    
    result = tryCatch({
      #get a subset of the connectivity matrix with only connectivity between 
      #grid points that are significant
      connectivityClus<-getSubsetOfConnectivity(fromIndexWithinBrain,toIndexWithinBrain)
      
      
      start <- Sys.time ()
      
      print(paste0("Matrixsize: ",dim(connectivityClus)[2]," x ",dim(connectivityClus)[2]))
      
      #take only grid points with at least 10 incoming connections (=true in validData)
      #since we are clustering by incoming connections, grid points
      #with only a few incoming connections would be clustered together, although they
      #don't have necessarily much in common. Therefore we remove them from clustering
      validData<-(colSums(connectivityClus>0)>10) 
      connectivityClus<-connectivityClus[,validData]
      
      
      if(file.exists(paste0(resultsFolder,"//",fname,"//hclustCorrelation_25p.rds"))==FALSE){
        
        clComputing<-c()
        if(useCores>1){
          print("Setup parallel computing...")
          clComputing <- makeCluster(useCores,outfile="parallel_log.txt",type="PSOCK")
          registerDoParallel(clComputing)
          print("Parallel computing setup done")
        }
        
        getDistIndex<-function(x,y,dimension){
          ret<-(((dimension-1)^2+(dimension-1))/2-((dimension-y)^2+(dimension-y))/2+(x-y))
          ret[y==1]<-x-1
          ret[y>=x]<-NA
          return(ret)
        }
        
        start <- Sys.time ()
        print(paste0("Start clustering with computation of the distance map: ",start))
        
        sumValid<-sum(validData)
        dd<-dist(rep(0,sum(validData)))
        
        chunkSize<-1000
        
        #Parallelized calculation of a distance map based on the correlation coefficient 
        #The distance map dd is based on the correlation of incoming connections
        #Why do we use incoming connections? Not every grid point is an injection site and
        #has outgoing connections. But nearly every grid point has incoming connections
        pb <- txtProgressBar(min=0, max=100, initial=0,style=3)
        a<-foreach(conMatSubset2=iter(connectivityClus, by = "col",chunksize =  chunkSize),i=icount(), .inorder=TRUE) %do%{
          setTxtProgressBar(pb, i/ceiling(dim(connectivityClus)[2]/chunkSize)*100)
          actCol<-(foreach(conMatSubset=iter(connectivityClus, by = "col",chunksize =  chunkSize),.combine="rbind", .inorder=TRUE) %op%{
            return(1-cor(conMatSubset>0,conMatSubset2>0))
          })
          
          distIndex<-rep(0,sumValid*dim(actCol)[2])
          act<-1
          for(rowI in 1:dim(actCol)[2]+(((i-1)*chunkSize))){  
            
            distIndex[act:(act+sumValid-1)]<-getDistIndex(rowI,1:dim(actCol)[1],sumValid)
            act<-act+sumValid
            
          }
          dd[distIndex[!is.na(distIndex)]]<-actCol[!is.na(distIndex)]
          return(NULL)
        }
        
        close(pb)
        
        if(useCores>1){
          print("Parallel computing finished")
          stopCluster(clComputing) 
        }
        
        print(paste0("Clustering took: ",difftime(Sys.time(),start,units="hours")))
        
        
        print("Clustering distance map: ")
        
        
        hc<-hclust(dd,"ward.D2")
        print("Save clustering...")
        
        saveRDS(hc,paste0(resultsFolder,"//",fname,"//hclustCorrelation_25p.rds"))
        saveRDS(validData,paste0(resultsFolder,"//",fname,"//validData_25p.rds"))
        print(paste0("Clustering ended!"," ",Sys.time()))
        
      }
      hc<-readRDS(paste0(resultsFolder,"//",fname,"//hclustCorrelation_25p.rds"))
      validData<-readRDS(paste0(resultsFolder,"//",fname,"//validData_25p.rds"))
      
      
      toIndex<-toIndex[validData]
      pvalsCombined<-pvalsCombined[validData]
      toIndexReverse<-c()
      toIndexReverse[toIndex]<-1:length(toIndex)
      fromIndex<-intersect(toIndex,fromIndex)
      
      print(paste0("Start kClustering"," ",Sys.time()))
      #Generate a connectivity graph for different k
      for(kOfClustering in 3:10){
        
        print(paste0("Clustering k=",kOfClustering))
        
        dir.create(paste0(resultsFolder,"//",fname,'//',"clustering_k",kOfClustering))
        
        clustering<-cutree(hc,k=kOfClustering)
        allClusters<-unique(clustering)
        
        clusteringPosFrom<-clustering[toIndexReverse[fromIndex]]
        
        #Generate a brain with z-projection for every cluster (visCluster)
        for(visCluster in allClusters){
          
          png(paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",visCluster,".png"),1000,1000,bg = "transparent")
          laststartx<-0
          laststarty<-0.75
          
          startXY<-c()
          for(i in 1:10){
            startXY<-rbind(startXY,c(laststartx,laststarty))
            startx<-cos(2*pi/10)*laststartx-sin(2*pi/10)*laststarty
            starty<-sin(2*pi/10)*laststartx+cos(2*pi/10)*laststarty
            laststartx<-startx
            laststarty<-starty
          }
          startXY[,1]<-startXY[,1]*-1
          startXY[3,1]<-startXY[3,1]-0.075
          startXY[4,1]<-startXY[4,1]-0.075
          startXY[8,1]<-startXY[8,1]+0.05
          startXY[9,1]<-startXY[9,1]+0.07
          startXY[1,2]<-startXY[1,2]+0.05
          startXY[6,2]<-startXY[6,2]-0.05
          
          
          allRegionsOfCluster<-(atlasRegions[toIndex][clustering==visCluster])
          factorsName<-c()
          act<-1
          
          
          #Mean Position is necessary to know where the label line should point to
          meanPvalsOfRegions<-c()
          regionsOfCluster<-c()  #significant regions of the cluster 
          meanPosition<-c()
          
          #gets the most significant regions of a cluster, their p-values and their mean position (center)
          #Region must be at least 0.1% of the cluster
          for(rid in unique(allRegionsOfCluster)){
            if(sum(allRegionsOfCluster==rid)>=(length(allRegionsOfCluster)/1000)){
              regionsOfCluster<-c(regionsOfCluster,rid)
              meanPvalsOfRegions<-c(meanPvalsOfRegions,mean(pvalsCombined[clustering==visCluster][atlasRegions[toIndex][clustering==visCluster]==rid]))
              meanPosition<-rbind(meanPosition,c((mean(indexX[toIndex][clustering==visCluster][atlasRegions[toIndex][clustering==visCluster]==rid])/dim(atlasRegions)[1]-0.5)*1.29-0.02,mean(mY-indexY[toIndex][clustering==visCluster][atlasRegions[toIndex][clustering==visCluster]==rid])/dim(atlasRegions)[2]*0.865-0.365))
            }
          }
          
          #The next few lines are basically about ordering
          #and arange the labels and label lines arround the
          #z-projection of the brain to point at the 10 most
          #significant regions
          orderOfPvals<-order(meanPvalsOfRegions)
          if(length(orderOfPvals)>10){
            orderOfPvals<-orderOfPvals[1:10]
          }
          regionsOfCluster<-(regionsOfCluster)[orderOfPvals]
          meanPosition<-meanPosition[orderOfPvals,]
          
          
          if(length(orderOfPvals)<10){
            regionsOfCluster<-c(regionsOfCluster,rep(0,10-length(orderOfPvals)))
            for(i in (length(orderOfPvals)+1):10){
              meanPosition<-rbind(meanPosition,c(0,0))
            }
          }
          
          labelingCosts<-matrix(0,10,10)
          for(start in 1:10){
            for(end in 1:10){
              labelingCosts[start,end]<-sqrt((startXY[start,1]-meanPosition[end,1])^2+(startXY[start,2]-meanPosition[end,2])^2)
            }
          }
          y <- solve_LSAP(labelingCosts)
          regionsOfCluster<-regionsOfCluster[y]
          meanPosition<-meanPosition[y,]
          
          
          for(rid in (regionsOfCluster)){
            region<-getAcronymByID(document$msg[[1]],rid)
            if(length(region)==0){
              factorsName<-c(factorsName,"")
            }
            else{
              factorsName<-c(factorsName,region$acronym)
              
            }
            act<-act+1
          }
          
          
          par(mar = c(0, 0, 0, 0))
          
          #Plotting region labels for the z-projection of the brain
          circos.initialize(1:length(factorsName), xlim=c(0,1))
          
          circos.trackPlotRegion(c(4:length(factorsName),1:3),ylim=c(0,1),bg.border="transparent",panel.fun = function(x, y){
            sector.index = get.cell.meta.data("sector.index")
            circos.text(x=0.5,y=0.5,labels=factorsName[c(4:length(factorsName),1:3)][as.numeric(sector.index)],facing="downward",cex=3.5,adj=c(0.47,0.5),font=2)
          })
          
          circos.clear()
          
          
          
          zProj<-array(0, c(dim(index)[1],dim(index)[2]))
          
          for (x in 1:(dim(index)[1]-1)){
            for(y in 1:(dim(index)[2]-1)){
              for (z in 1:(dim(index)[3]-(57/(0+1)))){ 
                
                zProj[x,mY-y]<-(zProj[x,mY-y]+(atlasRegions[index[x,y,z]]>0)*3.3)
              }  
            }
          }
          zProj<-zProj/max(zProj)*99
          
          
          contoursZProjH<-zProj<0
          for(x in 1:(dim(zProj)[1])){
            for(y in 1:(dim(zProj)[2])){
              
              if((zProj[min(x+1,dim(zProj)[1]),y]>0|
                  zProj[min(x+1,dim(zProj)[1]),min(y+1,dim(zProj)[2])]>0|
                  zProj[x,min(y+1,dim(zProj)[2])]>0|
                  zProj[x,max(y-1,1)]>0|
                  zProj[max(x-1,1),max(y-1,1)]>0|
                  zProj[max(x-1,1),min(y+1,dim(zProj)[2])]>0|
                  zProj[min(x+1,dim(zProj)[1]),max(y-1,1)]>0|
                  zProj[max(x-1,1),y]>0)&zProj[x,y]<=0){
                contoursZProjH[x,y]<-1
              }
            }
            
          }
          contoursZProjH[130,7:59]<-1
          contoursZProjH[132,]<-0
          contoursZProj<-image.smooth( contoursZProjH, theta=1)$z
          
          pvalsCombinedCols<--log10(pvalsCombined)
          pvalsCombinedCols[is.na(pvalsCombinedCols)]<-0
          pvalsCombinedCols[pvalsCombinedCols>(-log10(minimumPValueToVisualize))]<-(-log10(minimumPValueToVisualize))
          
          colorsClus<-(jet.colors(round((-log10(minimumPValueToVisualize))/max(pvalsCombinedCols[clustering==visCluster])*800)))[1:800]
          colors=c(rep(rgb(1,1,1,alpha=0),10),colorpanel(190,rgb(1,1,1,alpha=0), rgb(0,0,0,alpha=1)),colorsClus)
          
          
          leftPosTo<-toIndex[clustering==visCluster]
          leftPvalsCols<-pvalsCombinedCols[clustering==visCluster]
          leftPosToIndex<-c()
          leftPosToIndex[leftPosTo]<-1:length(leftPosTo)
          
          for(x in 1:(dim(contoursZProj)[1])){
            for(y in 1:(dim(contoursZProj)[2])){
              intersectionPos<-intersect(index[x,y,],leftPosTo)
              if(length(intersectionPos)>0){
                contoursZProj[x,mY-y]<-max(leftPvalsCols[leftPosToIndex[intersectionPos]])
              }
            }
          }
          
          
          #Plotting z-projection of the brain
          contoursZProj<-resizeImage(contoursZProj,dim(contoursZProj)[1]*5,dim(contoursZProj)[2]*5)
          add.image(-0.02,0.07,image.smooth( contoursZProj, theta=2)$z,image.width=0.40,image.height = 0.6,col=colors)
          
          #Plotting lines between labels and region center
          for(i in 1:10){
            if(nchar(factorsName[i])>0){
              lines(c(startXY[i,1],meanPosition[i,1]),c(startXY[i,2],meanPosition[i,2]),lwd=3)
            }
          }
          
          dev.off()
        }
        
        #Calculating the weights (normalized by the injection volume) between the clusters
        edgesRaw<-c()
        weightsRaw<-c()
        for(i in 1:length(allClusters)){
          for(j in 1:length(allClusters)){
            edgesRaw<-rbind(edgesRaw,c(i,j))
            weightsRaw<-c(weightsRaw,(sum(connectivityClus[clusteringPosFrom==i,clustering==j]>0)+sum(connectivityMatrix[fromIndexWithinBrain[clusteringPosFrom==i],getIndexOtherSideWithinBrain(toIndex[clustering==j])]>0))/sum(clusteringPosFrom==i))
            
          }
        }
        
        weights<-(weightsRaw)
        
        #scaling the weights, so they look nicer in the plots
        #Maximum weight is the 90% percentile of all weights
        #(to be robust versus outliers)
        weights<-weights/quantile(weights,0.9, na.rm = TRUE)*100  
        weights[is.na(weights)]<-0
        weights[weights>100]<-100
        edges<-edgesRaw[weights>10,]
        weights<-weights[weights>10]
        curves=c()
        
        for(i in 1:dim(edges)[1]){
          found<-0
          for(j in 1:dim(edges)[1]){
            if(edges[i,1]==edges[j,2] & edges[i,2]==edges[j,1]){
              found<-found+1
            }
          }
          if(found>=1){
            curves=c(curves,0.2)
          }else{
            curves=c(curves,0)
          }
          
        }
        
        
        #Build gml file and plot the png
        nodes<-allClusters
        if(length(edges)>0){
          G <- graph( t(edges), directed = TRUE )
          E(G)$weight<-(weights)
          E(G)$color<-gray.colors(100,start = 0.9, end = 0.2)[round(weights)]
          
          for(visCluster in allClusters){
            img <- readPNG(paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",visCluster,".png"),native=TRUE)
            
            V(G)$raster[visCluster]<-replicate(1,img, simplify=FALSE)
            
          }
          
          graphName="graph"
          
          png(paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".png"),4000,4000)
          
          G$layout <- layout_in_circle
          
          plot(G,vertex.label=NA,edge.curved=curves,vertex.shape="raster",vertex.size=80, vertex.size2=80,margin=0,edge.loop.angle=1.25*pi,edge.width=15,edge.arrow.size=5,edge.arrow.width=2)
          
          dev.off()
          
          
          write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=FALSE,x="graph")
          write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="[")
          write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="  directed 1")
          for(i in 1:length(nodes)){
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="node")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="[")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("id ",nodes[i]))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("label \"",i,"\""))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="graphics")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="[")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("fill \"#ffffff\""))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="type  \"ellipse\"")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("w ",1000))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("h ",1000))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="]")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="]")
          }
          for(i in 1:dim(edges)[1]){
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="edge")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="[")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("source ",edges[i,1]))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("target ",edges[i,2]))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="graphics")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="[")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("width ",25))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x=paste0("fill \"",E(G)$color[i],"\""))
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="targetArrow \"standard\"")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="]")
            write(file=paste0(resultsFolder,"//",fname,"//clustering_k",kOfClustering,"//",graphName,".gml"),append=TRUE,x="]")
          }
        }
      }
      
      print(paste0("kClustering finished!"," ",Sys.time()))
      
      
    }, error = function(e) {
      print(paste0("Error@ ",fname,": ",e))
      
    })
  }
}




