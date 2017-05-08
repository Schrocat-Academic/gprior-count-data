# noncontdata.R
#
# Created by Michael Walker 13/10/2014
# Emails: walkerm1@student.unimelb.edu.au, m.walker@aip.org.au
#
# Contains code for
# a) Converting non-continuous code into suitable latent variables
# b) Updating latent variables at every iteration
#
# Written for use with ARTiVa but ideally could be used with any g-prior
# package with a single insertion
#

library(truncnorm)

#########GLOBAL VARIABLES#########

noncontinit = FALSE		# Have the latent variables been initiated yet?
binary = c()	# Vector of binary variables
isbinarydata = FALSE	#Is the target data binary?
isintegerdata = FALSE	#Is the target data integer-valued?
originalData = array()	# Store original values of integer data
oldZdata = array()  # Store Z-values from previous iteration for moving cutoffs
Mtimepos = c()      # Vector of changing timepoints in binary data
cutoff = 0.5        # Cutoff point between enhanced binary values
cutoffs = c(-Inf,.5,Inf)

######### FUNCTIONS #########

print("Initiating noncontdata.R")

noncontARTIVA <- function(s, E, Sall, Ball, parentData, targetData, Mphase){
	# Functions in this module called via this function
	
	# Sample of extracting required data within ARTIVA
	
	if (!noncontinit){
		# Identify binary data fields automatically
        targetData = init_noncont(targetData)
	}
    else if (isintegerdata) {
    #else if (isbinarydata) {
        ## For each current phase
		for (phase in E[1:(s+1)]){
	    	posPhase = which(E==phase)
			parents = which(Sall[posPhase,]==1)
			nb_parents = length(parents)
 		}
        # Ensure vanishing coeffs for non-edges
        regresscoeffs = as.matrix(Ball * Sall)
        targetData = update_latentvar(parentData,regresscoeffs,E)
	}
    oldZdata <<- targetData
    #print("exiting noncont")
	targetData
}


#########Initiate latent variables#################

init_noncont <- function(targetData,interval=0.4){
	# Initiate necessary variables to replace binary variables
	# with latent variables
	
	# Identify integer data fields automatically
	if (is.vector(targetData)) {
        #if (is_binary(targetData)){
        if (is_integer(targetData)){
		    print("Initiating integer data")
            #print("targetData")
            originalData <<- targetData
            #isintegerdata <<- TRUE
            noncontinit <<- TRUE
            
            cutoffs <<- c(-Inf,seq((min(originalData)+.5),(max(originalData)-.5)),Inf)
            #lwrcutoffs = c(seq((min(originalData)-interval),(max(originalData)-interval)))
            #upprcutoffs = c(seq((min(originalData)+interval),(max(originalData)+interval)))
            #cutoffs <<- unique(sort(c(lwrcutoffs,upprcutoffs)))
            
            cat("cutoffs",cutoffs,'\n')
        
            originalnames = names(originalData)
            timepoints = as.integer(substr(originalnames,nchar(originalnames),nchar(originalnames)))
            timechanges = timepoints[1:length(timepoints)-1] != timepoints[2:length(timepoints)]
            Mtimepos <<- c(1,grep(TRUE,timechanges)+1,length(timepoints)+1)
		}
	}
    else {  # For packages where targets  are not processed one-at-a-time.
        binary <<- find_binary(targetData)
	}	
    if (isintegerdata){
        #if (isbinarydata){
            # Store original data values. GLOBAL. 
		originalData <<- targetData
	
        # Replace binary variables with latent ones
        #targetData = subst_binary()
        # Replace integer variables with latent ones
        targetData = subst_integer()
	}
	targetData
}

is_binary <- function(targetData){
	#Checks that vector targetData is binary valued
    length(unique(targetData)) == 2
}

is_integer <- function(targetData){
    # Checks that vector targetData is integer valued
    notinteger = FALSE %in% (targetData == round(targetData))
    !notinteger
}

find_binary <- function(targetData){
	# Identify which predictors in targetData are binary-valued
	if (is.null(rownames(targetData))) rownames(targetData) = "target"
	if (is.vector(targetData)){
		binaryvals = targetData %in% c(0,1)
		print(rownames(binaryvals))
	}
	else {
		binaryvals = array(targetData %in% c(0,1),dim= c(nrow(targetData),ncol(targetData)))
		rownames(binaryvals) = rownames(targetData)
	}
	# Identify non-binary observables
	nonbinary = vector('logical',ncol(binaryvals))
	for (observe in colnames(binaryvals)) nonbinary[observe] = FALSE %in% binaryvals[,observe]
	
	# Return binary observables
	colnames(binaryvals)[!nonbinary]
}

######### Substitute latent variables in place of non-continuous data ########

subst_integer <- function(){
    #subst_binary <- function(){
	#Replaces binary/integer data fields in target data with latent variables
	# originalData is the original data for target variables
	
	# Call truncnorm to generate latent values for binaries
    integerdata = call_truncnorm()
	integerdata
}


#########Update latent variables according to sample and regresscoeffs#######

update_latentvar <- function(parentData,regresscoeffs,E){
	#Generates new values for latent variables according to changes in 
	# DBN and the regression coefficients.
    targetData = call_truncnorm(parentData,regresscoeffs,E) # binaryData
    #targetData = binaryData
	targetData
}

call_truncnorm <- function(parentData=NA, regresscoeffs,E,trunsd=1,interval=0.4){
	# Generates latent variable values, choosing values between appropriate cutoffs
	# Uses R package truncnorm.R
	
	#Default values for regression coefficients

    newZdata = array()
    
    originalnames = names(originalData)
    sampletimes = as.integer(unique(substring(originalnames,nchar(originalnames),nchar(originalnames))))

    #datarow = originalData
    
	# Choose limits to put 0,1 variables below/above cutoff respectively
    #interval = 0.1
    #lwrlimits = originalData - interval
    #upprlimits = originalData + interval

    lwrlimits = cutoffs[originalData+1] #- interval
    upprlimits = cutoffs[originalData+2] #+ interval
    
	# if initiating, no coefficients available
    if (length(parentData)==1 & (TRUE %in% is.na(parentData))) {
        newZdata = as.vector(t(rtruncnorm(1, lwrlimits, upprlimits,mean=originalData,sd=trunsd)))
    }
    else {	# if updating after iteration
        newZdata = c()
        phase=1
		for (sampletimepos in sampletimes){
            if (sampletimepos == E[phase+1]) phase = phase + 1
            phasecoeffs = as.matrix(regresscoeffs[phase,])
            
            parentRow = as.matrix(parentData[seq(Mtimepos[sampletimepos-1],Mtimepos[sampletimepos]-1),])

            means = as.vector(t(parentRow %*% phasecoeffs))
            meannames = substr(rownames(parentRow),1,nchar(rownames(parentRow))-1)
            newmeannames = paste0(meannames,sampletimepos)
            names(means) = newmeannames

            timelwr = lwrlimits[Mtimepos[sampletimepos-1]:(Mtimepos[sampletimepos]-1)]
            timeuppr = upprlimits[Mtimepos[sampletimepos-1]:(Mtimepos[sampletimepos]-1)]

            newZdata[Mtimepos[sampletimepos-1]:(Mtimepos[sampletimepos]-1)] = rtruncnorm(1, timelwr, timeuppr, mean = means,sd=trunsd)
            
            if (TRUE %in% is.na(newZdata)){
                cat("means",means,'\n')
                cat("newZdata",newZdata,'\n')
                cat("timelwr",timelwr,'\n')
                cat("timeuppr",timeuppr,'\n')
                stop("NA appeared in newZdata")
            }
            #else print("FALSE")
		}
	}
    names(newZdata) = originalData
	newZdata
}

sample_cutoffs <- function(Zdata,E){
    # Choose next set of cutoffs
    
    integervals = seq(min(originaldata),max(originaldata))
    #minZ = c()
    #maxZ = c()
    newcutoffs = c()
    for (val in integervals){
        #minZ = c(minZ,min(Zdata[originaldata==val]))
        #maxZ = c(maxZ,max(Zdata[originaldata==val]))
        Zval = Zdata[originaldata==val]
        lwrrange = c(cutoffs[2*val],min(Zval))
        # Generate random number from uniform distribution in lwrrange
        newcutoffs = c(newcutoffs) 
        upprrange = c(max(Zval),cutoffs[2*val+3])
        # Generate random number from uniform distribution in upprrange
        newcutoffs = c(newcutoffs) 
    }
}

#########Call module with source() from main loop########

