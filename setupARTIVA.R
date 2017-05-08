source("work/code/functions.R")
	
call_ARTIVA <- function(xlsdata = NULL,niter = 50000, correctGaussianPrior=TRUE, meanImpute = FALSE, rerun = FALSE, name_only = FALSE){	
	# xlsdata - data.frame, CAS data for inferring networks
	# niter - integer, ARTIVA parameter, number of iterations in MCMC
	# correctGaussianPrior - boolean, use Albert and Chib variables for binary/count
	#targets
	# meanImpute - boolean, meanImpute missing data if TRUE, discard sample if FALSE
	# rerun - boolean, if TRUE ignore whether network has been found already
	# name_only, boolean, if TRUE find name of network but then stop

	timesteps = 5
	firsttime = 0
	alphacp = 1		# alpha for change point prior distribution
	betacp = 0.5		# beta for change point prior distribution
	alphaedges = 1		# alpha for change point prior distribution
	betaedges = 0.5		# beta for change point prior distribution
	parentfolder = ""    # list of parents in foldername
	targetfolder = ""    # list of targets in foldername
	bothfolder = ""    # list of both parents and targets in foldername
	conditionfolder = ""    # list of conditions in foldername
	simfolder = "/Users/walkerm1/work/data/simulations/AlbertChibb/"  # folder with sim files
	suboutfolder = ""
	simfile <- file('/Users/walkerm1/work/code/testnodes.txt','r') 
	simname <- readtest()
	simname = paste0(simfolder,simname,'.txt')
	#print(c("simname",simname))
	fileline = ""
	simfile <- file(simname,'r') 
	fileline <- strsplit(readLines(simfile,1),'\t')
	#print(c('fileline',fileline))
	datafiles = c()
	while (length(fileline)){    # read lines until EOF
		if (fileline[[1]][1] == 'datafiles') datafiles = c(datafiles,fileline[[1]][2])
		else if (fileline[[1]][1] == 'outfolder') {
			outfolder = fileline[[1]][2] 
			#if (substring(outfolder,nchar(outfolder)) != "/") outfolder = paste0(outfolder,"/")
		}
		else if (fileline[[1]][1] == 'simnodes') simnodes = fileline[[1]][2:length(fileline[[1]])] 
		else if (fileline[[1]][1] == 'parents') parents <- fileline[[1]][2:length(fileline[[1]])]
		else if (fileline[[1]][1] == 'targets') targets = fileline[[1]][2:length(fileline[[1]])] 
		else if (fileline[[1]][1] == 'timesteps') timesteps = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'suboutfolder') suboutfolder = fileline[[1]][2]
		else if (fileline[[1]][1] == 'firsttime') firsttime = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'alphaCP') alphacp = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'betaCP') betacp = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'alphaEdges') alphaedges = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'betaEdges') betaedges = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'conditions') conditions = fileline[[1]][2] 
		fileline <- strsplit(readLines(simfile,1),'\t')
	}
	close(simfile)

	outname = paste0(outfolder,'/ARTIVA_FinalNetwork.txt')
	
	# Put data in suitable form and run ARTIVA
	# Check first that this particular network has not been done already
	if ((!file.exists(outname) | rerun) & !name_only){	#Do not regenerate preexisting network 
		library(ARTIVA)	#,lib.loc="/Users/walkerm1/work/ARTIVA/code/ARTIVA/") 
		#source("/Users/walkerm1/work/ARTIVA/code/ARTIVA/R/ARTIVAnet.R")
		#stop('DEBUG only')

		print("Reading .xls file")
		if (length(xlsdata)== 0) xlsdata <- readframe(datafiles)    # Call function readframe from functions.R
		xlsnames = names(xlsdata)
		#print(length(rownames(xlsdata)))
		nb_entries = length(rownames(xlsdata))
		#print(nb_entries)
		boundaries <- c()
		dynaroots <- c()
		newboundaries <- c()
		newdynaroots <- c()
		for (simnode in simnodes){
			#cat("DEBUG: simnode",simnode,'\n')
			if (simnode %in% xlsnames) boundaries <- c(boundaries,simnode)
			else if (TRUE %in% grepl(simnode,substr(xlsnames,1,nchar(simnode)))){	#substr(xlsnames,1,nchar(simnode)) 
				dynaroots <- c(dynaroots,simnode)
				appendsimnode = TRUE
				grepnames = grep(simnode,xlsnames,value = TRUE)
				for (grepname in grepnames){
					if (nchar(grepname) <= nchar(simnode)+1 | simnode == strsplit(grepname,'6m')[[1]][1] | simnode == strsplit(grepname," ")[[1]][1]){
						appendsimnode = FALSE 
						break
					}
				}
				if (appendsimnode){
					xlsdata = appendnewdata(xlsdata,simnode)
					xlsnames = names(xlsdata)
				}
			}
			else{
				print(c(simnode,"not currently in xlsdata"))
				compound = strsplit(simnode,"_")[[1]]
				if (length(compound) >1){	# compound nodes
					if (compound[1] %in% xlsnames & compound[3] %in% xlsnames) newboundaries <- c(newboundaries,simnode)  #compound of boundaries is boundary  #TODO: fix to take boundaries from newnodes.txt
					else newdynaroots <- c(newdynaroots,paste0(simnode,' ')) #otherwise dynamic
					print("Updating xlsdata")
					xlsdata <- appendnewdata(xlsdata,simnode)
					xlsnames = names(xlsdata)
				}
				else{
					cat("Reading in datafield from newnodes",simnode,'\n')
					xlsdata <- appendnewdata(xlsdata,simnode)
					xlsnames = names(xlsdata)
					grepnames = grep(simnode,xlsnames,value = TRUE)
					if (length(grepnames)==0) stop("ERROR: node not found in any dataset.")
					for (grepname in grepnames){
						if (grepname == simnode) newboundaries <- c(newboundaries,simnode)
						else if (nchar(grepname) == nchar(simnode)+1 | simnode == strsplit(grepname," ")[[1]][1]) newdynaroots <- c(newdynaroots,simnode)
					}
				}
			}
		}
	#print("dynamic")
    #print(dynamic)
	#print("parents")
	#print(parents)
	dynamic = c()
	adddynaroots = c(dynaroots,newdynaroots)
	print(c("firsttime",firsttime,"timesteps",timesteps))
	if (length(adddynaroots)>0){
		calltimesteps = timesteps(adddynaroots,xlsnames)
		dynamic = as.vector(unlist(calltimesteps$dynamic))
		steps = as.vector(unlist(calltimesteps$timesteps))
		timesteps = length(steps)
		#cat("DEBUG: steps",steps,"timesteps",timesteps,'\n')
		if (FALSE){
		firsttimeold = firsttime
		timestepsold = timesteps
		for (loop in firsttimeold:timestepsold){	# Construct vector of timepoint name roots
			dyna <- paste(adddynaroots,as.character(loop),sep="")
			for (dynaloop in dyna){
				#cat("DEBUG: dynaloop",dynaloop,'\n')
				if (!dynaloop %in% names(xlsdata)){
					#print(c("Timestep",loop,"missing in datafield",dynaloop),quote=FALSE)
					dynaroot = substr(dynaloop,1,nchar(dynaloop)-1)
					if (firsttime == loop){
						#if (dynaroot %in% parents){
							firsttime = loop + 1	# Leave out this timestep
						#}
					}
					else if (firsttime < loop & loop <= timesteps){
						#if (dynaroot %in% targets){
							timesteps = loop - 1
						#}
					}
				}
				else dynamic <- c(dynamic,dynaloop)
			}
		}
		}
	}
	#else steps = seq(5)
	boundaries = c(boundaries,newboundaries)
	datanodes = unique(c(boundaries,dynamic))
	#print("targets")
	#print(targets)
	print("Finding reduced data set") 
	reduceddata <- subset(xlsdata, eval(parse(text=conditions)), select = datanodes) 
	#print(reduceddata)
	
	#datamissing = vector(mode="integer",length=length(datamean))
	#names(datamissing) <- names(datamean)
	
	# Identify fields with missing entries
	missingnames = unique(dynamic)	# seed names(datamissing) #TODO: why is dynamic not unique?
	for (step in steps) missingnames = c(missingnames,paste0(boundaries,step))	#names of missing data vector
	datamissing = vector(mode="integer",length=length(missingnames))
	names(datamissing) = missingnames
	rootnodes = c(boundaries,dynaroots,newdynaroots)
	
	# Either impute missing data or remove defective records
	if (meanImpute){
		print("Finding means")
		datamean = colMeans(reduceddata,na.rm = TRUE)
		if (correctGaussianPrior){
			# Round-off count data
			integers = find_integers(reduceddata)
			datamean[integers] = as.integer(round(datamean[integers]))
		}
	}
	else {
		## Identify records (now columns) with missing data
		##removerows = apply(reduceddata,1,anyNA) #!grepl(0,rowSums(isnamatrix))
		##keeprows = !removerows
		
		# Remove records with missing data
		##reduceddata = reduceddata[keeprows,]
		reduceddata = remove_na(reduceddata)
	}
	
	# Find number of included records
	nb_records = nrow(reduceddata)
	
	#TODO: implement shortening with timesteps
	#timemean = array(dim = c(length(rootnodes),length(steps)))	#timesteps-(firsttime-1)))	# 1:5*length(datanodes), 
	#print(timemean)
	datadyna = as.data.frame(matrix(ncol=length(rootnodes)))
	boundstep <- reduceddata[boundaries]
	names(datadyna) = rootnodes
	for (step in steps){	# Adding each timestep
		stepind = grep(step,steps)
		stepentries = grep(step,substr(dynamic,(1+nchar(dynamic)-nchar(step)),nchar(dynamic)))
		stepnames = dynamic[stepentries]	# timestep entries from last character
		#cat("DEBUG: step",step,'stepnames',stepnames,'\n')
		dynastep <- reduceddata[stepnames]
		dynastep <- cbind(boundstep,dynastep)
		#cat("DEBUG: dynamic",dynamic,'\n')
		names(dynastep) = names(datadyna) # rootnodes[rootnodes %in% c(boundaries,substr(stepnames,1,(nchar(stepnames)-1)))]	 #TODO: recognise parents/targets
		rownames(dynastep) = paste(rownames(dynastep),stepind,sep="-t")
		datadyna <- rbind(datadyna,dynastep)
	}
	cleandata <- datadyna[2:length(rownames(datadyna)),]
	
	# Either impute missing data or remove incomplete entries

	# Replace missing data with mean of data
	isnamatrix <- is.na(cleandata)	# locate NAs in cleandata
	if (meanImpute){
		for (col in boundaries){   # Substitute means into boundary NAs
			narows <- grep(TRUE,isnamatrix[,col])
			cleandata[narows,col] = datamean[[col]]
			#for (step in firsttime:timesteps) missingnames = c(missingnames,paste0(col,step))
			#datamissing[[col]] = length(narows)/(timesteps+1-firsttime) 
		}
		#datamissing = vector(mode="integer",length=length(missingnames))
		for (col in c(dynaroots,newdynaroots)){	# Substitute means into dynamic NAs
			narows <- grep(TRUE,isnamatrix[,col])
			if (length(narows)) {
				timestep = strsplit(rownames(isnamatrix)[narows],"-t")
				for (ind in 1:length(narows)){
						stepindex = as.integer(timestep[[ind]][[2]])
						stepname = paste(col,steps[stepindex],sep="")
						cleandata[narows[[ind]],col] = datamean[[stepname]]
				}
				#datamissing[[paste0(col,timestep[[ind]][[2]])]] = datamissing[[paste0(col,timestep[[ind]][[2]])]] + 1
			}
		}
		for (col in rootnodes){	# Generate datamissing array
			narows <- grep(TRUE,isnamatrix[,col])
			if (length(narows)) {
				timestep = strsplit(rownames(isnamatrix)[narows],"-t")
				for (ind in 1:length(narows)){
					stepindex = as.integer(timestep[[ind]][[2]])
					stepname = paste(col,steps[stepindex],sep="")
					datamissing[[stepname]] = datamissing[[stepname]] + 1
				}
			}
		}
	}
	
	# Construct final data frame
	# Row labels
	#print(c('datamissing',datamissing))
	datadescr = c()
	timestep = strsplit(rownames(cleandata),"-t")
	times = matrix(0,ncol=timesteps)
	for (step in timestep) times[[as.numeric(step[[2]])]] = times[[as.numeric(step[[2]])]] + 1	
	for (step in firsttime:timesteps) datadescr = c(datadescr,matrix(step,ncol=times[[step]]))	# datadescr: how many samples in each timestep?, gives value of Datadescription in ARTIVAnet
	# Frame construction in form required by ARTIVA
	print("Construct final data frame")
	finaldata <- t(cleandata)
	compndparents = c()
	compndtargets = c()

	#Find compndparents
	for (parent in parents) if (!parent %in% rownames(finaldata)) compndparents = c(compndparents,grep(parent,parents)) 
	# Avoid multiple corrections from similar compound names
	compndparents = unique(compndparents)

	for (ind in compndparents) parents[[ind]] = paste0(parents[[ind]]," ")
	#Find compndtargets

	for (target in targets) if (!target %in% rownames(finaldata)) compndtargets = c(compndtargets,grep(target,targets)) #Find compndtargets
	# Avoid multiple corrections from similar compound names
	compndtargets = unique(compndtargets)
	for (ind in compndtargets) targets[[ind]] = paste0(targets[[ind]]," ")

	#print(targets)
	#print(parents)
	#cat("DEBUG: rownames(finaldata)")
	finaltargets <- finaldata[targets,]
	finalparents <- finaldata[parents,]
	#print(head(finalparents))
	#print(head(finaltargets))
	noCP = length(steps) - 1
	steps[steps == 'cd'] = 'birth'
	print(c("timesteps",steps))
	ARTIVAnet(targetData = finaltargets, parentData = finalparents, segMinLength = 1, edgesThreshold = 0.4, dataDescription = datadescr, niter = niter, dyn = 1, maxCP=noCP, targetNames= targets, parentNames = parents, outputPath = outfolder,alphaCP=alphacp,betaCP=betacp,alphaEdges=alphaedges,betaEdges=betaedges,correctGaussianprior=correctGaussianPrior)

	dictmissing = ""  # String to generate dictionary of missing values
	for (node in names(datamissing)) dictmissing = paste0(dictmissing,node,':',datamissing[node],',')
	dictmissing = substr(dictmissing,1,nchar(dictmissing)-1) # Remove final ','

	#firsttime = steps[1]
	#timesteps = length(steps)
	outfile <- file(outname,'a') 
	writeLines('ADDENDA',outfile) # Adding extra data to ARTIVA output.
	writeLines(paste0(c('time steps',steps),collapse='\t'),outfile)
	writeLines(paste('start time',firsttime,sep='\t'),outfile)
	writeLines(paste('end time',timesteps,sep='\t'),outfile)
	writeLines(paste('nb records',nb_records,sep='\t'),outfile)
	writeLines(paste('nb entries',nb_entries,sep='\t'),outfile)
	writeLines(paste('data missing',dictmissing,sep='\t'),outfile)
	writeLines(paste('suboutfolder',suboutfolder,sep='\t'),outfile)
	close(outfile)
	}
#	else print(c("Network",outfolder,"already exists"))
	outfolder
}