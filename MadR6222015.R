# This is a modified version of a talk given at the MadR meeting on 06/17/2015
# The original datafiles for this talk are not included, but a step is provided
# to download comparable data from the Paleobiology Database into R using
# the new API.

################################# Load Necessary Packages #############################
# non-function objects are PascalCase and custom functions are camelCase
# Dependency functions are not embedded in master functions

# Required Packages for today's examples
library(RCurl) # Package for downloading https
library(eHOF) # Model selection program for various non-linear logistic regressions
library(nnet) # Package for multinomial logistic regression
library(doParallel) # Extension of packages parallel and iterations
library(rbenchmark) # Benchmarking package

# Shout-outs to the following packages/alternatives
library(pbapply) # Package for adding progress bars to apply() family functions
library(MASS) # Package capable of ordinal logistic regression
library(plyr) # Package with optimized apply() family functions
library(microbenchmark) # Alternative benchmarking package
library(pscl) # Package for computing pseudo-R2 for logistic regressions

# Alternative packages that I do not especially recommend
library(foreach) # Alternative Parallelization package
library(mlogit) # Alternative package for multinomial logistic regression
library(mnlogit) # Alternative package for multinomial logistic regression

####################### Download Data from the Paleobiology database###################

# Download from PBDB by Class and Period
# Object Taxa must be a vector of taxa
downloadPBDB<-function(Taxa="Bivalvia",Period="Pleistocene") {
	# Only one taxon
 	if (length(Taxa)==1) {
		URL<-paste("http://paleobiodb.org/data1.1/occs/list.csv?base_name=",Taxa,"&interval=",Period,"&show=paleoloc,phylo&limit=all",sep="")
		File<-read.csv(URL,header=T)
		return(File)
		}
 	# Vector of multiple taxa
	else {
		Taxa<-paste(Taxa,collapse=",")
		URL<-paste("http://paleobiodb.org/data1.1/occs/list.csv?base_name=",Taxa,"&interval=",Period,"&show=paleoloc,phylo&limit=all",sep="")
		File<-read.csv(URL,header=T)
		return(File)
		}
	}

Epoch<-downloadPBDB("Bivalvia","Pleistocene")

####################### Calculate Maximum Great Circle Distance #######################
# Examine a sample dataset of Pleistocene Bivalves (clams)
# from the Paleobiology Database
head(Epoch)
    	genus paleolat paleolng
	1  Barnea    37.06   -76.05
	2 Anadara    37.06   -76.05
	3 Mulinia    37.06   -76.05
	4   Ensis    37.06   -76.05
	5 Angulus    37.06   -76.05
	6    Abra    37.06   -76.05

# Check the dimensions
dim(Epoch)
	[1] 11059     3

# Check the number of unique genera in the dataset
length(unique(Epoch$genus))
	[1] 533 # For comparison, the full OBIS dataset has ~5500 genera

# A function for converting degrees to radians
convertRadians<-function(deg) {
	return(deg*pi/180)
	}

# Convert latitude and longitude from degrees to radians
Epoch$paleolat<-convertRadians(Epoch$paleolat)
Epoch$paleolng<-convertRadians(Epoch$paleolng)

# Calculates the geodesic distance between two points specified by 
# radian Latitude/Longitude using the Haversine formula - r-bloggers.com
# The Haversine formula is a special case of the spherical law
# of cosines, which is a special case of the law of cosines.
calcHaversine<-function(Long1,Lat1,Long2,Lat2) {
	Radius<-6371 # radius of the Earth (km)
	DeltaLong<-(Long2-Long1)
	DeltaLat<-(Lat2-Lat1)
	A<-sin(DeltaLat/2)^2+cos(Lat1)*cos(Lat2)*sin(DeltaLong/2)^2
	C<-2*asin(min(1,sqrt(A)))
	Distance<-Radius*C
	return(Distance) # Distance in km
	}

# Make a matrix of all unique latitude and longitude combinations
# Based around use of the combn() function in base R
# combn() and its cousin expand.grid() are brute force functions
# They create all possible combinations elements in two vectors
parseCombine<-function(GenusSubset) {
	# Collapse the latitude and longitude of a point into a single string
	CollapsedCoords<-apply(GenusSubset[,c("paleolat","paleolng")],1,paste,collapse=",")
	# Make every unique combinaton of coordinate pairs using combn()
	CombinedCoords<-combn(CollapsedCoords,2)
	# Rebind the output list into a dataframe
	CombinedCoords<-do.call(rbind,CombinedCoords)
	# Split the string of coordinates back into unique lat and long values
	FirstCoordsSplit<-strsplit(CombinedCoords[,1],",")
	FirstCoords<-do.call(rbind,FirstCoordsSplit)
	SecondCoordsSplit<-strsplit(CombinedCoords[,2],",")
	SecondCoords<-do.call(rbind,SecondCoordsSplit)
	# Rebind in a final dataset
	Coordinates<-cbind(FirstCoords,SecondCoords)
	Coordinates<-apply(Coordinates,2,as.numeric)
	colnames(Coordinates)<-c("Lat1","Long1","Lat2","Long2")
	return(Coordinates)
	}

# Subset the Pleistocene PBDB dataset to look at only one genus, "Anadara"
Anadara<-subset(Epoch,Epoch$genus=="Anadara")

# View the subset data structure
head(Anadara)
	     genus    paleolat  paleolng
	2  Anadara 0.011289122 -1.327323
	17 Anadara 0.011249521 -1.326101
	36 Anadara 0.009741665 -1.414938
	37 Anadara 0.009741665 -1.414938
	38 Anadara 0.009741665 -1.414938
	47 Anadara 0.009723388 -1.414938

# Eliminate duplicated rows; 
# Equivalent to using unique(), but more explicit
# There is no statistical difference in speed
Anadara<-Anadara[!duplicated(Anadara),]

# Make  a matrix of all unique combinations of Anadara points
# This is a brute force approach
AnadaraCoordinates<-parseCombine(Anadara)
head(AnadaraCoordinates)
	           Lat1     Long1        Lat2     Long2
	[1,] 0.01128912 -1.327323 0.011249521 -1.326101
	[2,] 0.01128912 -1.327323 0.009741665 -1.414938
	[3,] 0.01128912 -1.327323 0.009723388 -1.414938
	[4,] 0.01128912 -1.327323 0.008419625 -1.403245
	[5,] 0.01128912 -1.327323 0.010125483 -2.085145
	[6,] 0.01128912 -1.327323 0.008005346 -1.407259

# Use a for() loop to calculate the calcHaversine formula
AnadaraDistances<-vector("numeric",length=nrow(AnadaraCoordinates))
for (j in 1:nrow(AnadaraCoordinates)) {
	AnadaraDistances[j]<-calcHaversine(
		Long1=AnadaraCoordinates[j,2],
		Lat1=AnadaraCoordinates[j,1],
		Long2=AnadaraCoordinates[j,4],
		Lat2=AnadaraCoordinates[j,3]
		)
	}

# Find the maximum great circle distance of Anadara points
max(AnadaraDistances)
	[1] 19959.55 # in km
# For perspective the greatest circle distance possible on the earth
# is just over 20,000 km

################# Calculate Max Great Circle Distance of all Genera ###################
####################### Discuss for(),apply(),parApply() Merits #######################
# Vectorization is also an option, but we are going to ignore that for the moment
# We want to compare for(), apply() family functions, and parallelized apply()

# This is a very common in problem in R, and there are numerous packages with functions
# to handle it quickly and efficiently - e.g., plyr, reshape, data.table, sqldf, and doBy
# But, I want to talk a little a bit about principles, so we are going to use
# native R functions for this demo.

# "I feel the need for speed." - Maverick

# First function
# Calculate maximum great circle distance of all genera
# Uses for() loops without assigned memory
forRangeGenus<-function(Epoch=Pleistocene) {
	Genera<-unique(Epoch$genus)
	FinalVector<-vector()
	# Subset the dataset by each genus
	for (i in 1:length(Genera)) {
		GenusSubset<-subset(Epoch,Epoch$genus==Genera[i])
		Coordinates<-parseCombine(GenusSubset)
		GenusRanges<-vector()
		# Calculate great circle distance of all genus points
		for (j in 1:nrow(Coordinates)) {
			GenusRanges[j]<-calcHaversine(
				Long1=Coordinates[j,2],
				Lat1=Coordinates[j,1],
				Long2=Coordinates[j,4],
				Lat2=Coordinates[j,3]
				)
			}
		# Find the maximum distance of each set of distances
		FinalVector[i]<-max(GenusRanges,na.rm=T)
		}
	names(FinalVector)<-Genera
	return(FinalVector)
	}

# Second Function
# Calculate maximum great circle distance of all genera
# Uses split() and sapply()
applyRangeGenus<-function(Epoch=Pleistocene) {
	# Subset dataset by each genus
	Genera<-split(Epoch,Epoch$genus)
	# Calculate great circle distance of all genus points and find max
	RangeList<-sapply(Genera,maxRange)
	names(RangeList)<-names(RangeList)
	return(RangeList)
	}

# Third function
# Calculate maximum great circle distance of all genera
# Uses ONLY for() Loops, but with pre-defined arrays
assignRangeGenus<-function(Epoch=Pleistocene) {
	Genera<-unique(Epoch$genus)
	FinalVector<-vector("numeric",length=length(Genera)) # Assign memory
	# Subset the dataset by each genus
	for (i in 1:length(Genera)) {
		GenusSubset<-subset(Epoch,Epoch$genus==Genera[i])
		Coordinates<-parseCombine(GenusSubset)
		GenusRanges<-vector("numeric",length=nrow(Coordinates)) # Assign memory
		# Calculate great circle distance of all genus points
		for (j in 1:nrow(Coordinates)) {
			GenusRanges[j]<-calcHaversine(
				Long1=Coordinates[j,2],
				Lat1=Coordinates[j,1],
				Long2=Coordinates[j,4],
				Lat2=Coordinates[j,3]
				)
			}
		# Find the maximum distance of each set of distances
		FinalVector[i]<-max(GenusRanges,na.rm=T)
		}
	names(FinalVector)<-Genera
	return(FinalVector)
	}

# Fourth function
# Calculate maximum great circle distance of all genera
# Uses by(), the native forebear of ddply() in the plyr package
byRangeGenus<-function(Epoch=Pleistocene) {
	# Calculate great circle distance of all genus points
	# Find the maximum distance of each set of distances
	RangeList<-by(Epoch,Epoch$genus,maxRange)
	# Convert by() output to a vector with names
	RangeList<-setNames(as.vector(RangeList),names(RangeList)) # setNames saves memory!
	return(RangeList)
	}

# Perform a benchmark test of different Haversine functions functions
benchmark(
	forRangeGenus(Epoch),
	applyRangeGenus(Epoch),
	assignRangeGenus(Epoch),
	byRangeGenus(Epoch),
	)
	                      test replications elapsed relative
	3   applyRangeGenus(Epoch)          100 547.528    1.037 # 3% slower
	4      byRangeGenus(Epoch)          100 527.816    1.000 # By is the fastest
	2  assignRangeGenus(Epoch)          100 583.822    1.106 # 10% slower
	1     forRangeGenus(Epoch)          100 750.532    1.422 # 42% slower

# So you're probably asking, since by() is the fastest, why not just use
# an optimized version of by(), like ddply(), and be done with it?

# Fifth Function
# Calculate maximum great circle distance of all genera
# Using a parallelized version of sapply(), parSapply()
# There is no parallelized version of by() in the doParallel package
parallelRangeGenus<-function(Epoch=Pleistocene) {
	# Make a cluster for parallelization
	cl<-makeCluster(4)
	# Export functions to be used in parallel
	clusterExport(cl=cl, varlist=c("calcHaversine","maxRange","parseCombine"))
	# Subset dataset by each genus
	Genera<-split(Epoch,Epoch$genus)
	# Calculate great circle distance of all genus points
	# Find the maximum distance of each set of distances
	RangeList<-parSapply(cl,Genera,maxRange)
	stopCluster(cl) # Don't forget to do this!
	names(RangeList)<-names(Genera)
	return(RangeList)
	}

# Perform a benchmark test of different great circle functions
benchmark(
	forRangeGenus(Epoch),
	applyRangeGenus(Epoch),
	assignRangeGenus2(Epoch),
	byRangeGenus(Epoch),
	parallelRangeGenus(Epoch)
	)
	                       test replications elapsed relative
	3    applyRangeGenus(Epoch)          100 537.487    1.974
	4       byRangeGenus(Epoch)          100 533.517    1.960 # by is twice as long!
	2  definedRangeGenus(Epoch)          100 580.627    2.133
	1      forRangeGenus(Epoch)          100 773.786    2.842
	5 parallelRangeGenus(Epoch)          100 272.262    1.000

# So now you're probably asking, since parApply() is the fastest, why not just use
# parallelization and be done with it?

# "Beware of arguments related to programming speed. 
# All things being equal, faster is better. 
# But all things are never equal." - Josh Tyrangiel

############################## Parallelizaton Caveats ##################################
# There are lots of deprecated or otherwise out-of-date packages for parallelization
# Formerly popular parallelization packages (i.e., snow, multicore) are now part of 
# base R. In fact, the packages we're using here, parallel and doParallel, are based 
# on snow and multicore. This can make it very confusing to learn parallelization 
# because many online tutorials are out of date.

# I recommend using the doParallel package and parApply() family functions if you can 
# get away with it. If you absolutely cannot, then try using the foreach package, though
# foreach operates in largely the same maner as parLapply().

# All iterations of a parallelized process must be independent of each other. 
# Meaning that iteration i cannot affect or depend upon iteraton i++ or i--
# Parallelization where the different processes talk to each other is a whole
# different ballgame. You don't want to deal with that unless if you're a 
# high end computer scientist

# If the parallelized operations are too few in number, it can actually run slower than
# not using it. This is because it takes effort for the cluster to assign
# jobs among each processor, and if there aren't enough jobs to justify that effort then
# it's actually going to be slower.

################################### Memory Considerations ##############################
# The talk used a dataset downloaded from iobis.org, that dataset is not included with this
# file in GitHub. I recommend downloading a comparable dataset from the paleobiology database

# options(timeout=300)
# CanonicalTaxa<-c("Bivalvia","Gastropoda","Brachiopoda","Anthozoa","Bryozoa","Crinoidea")
# PhanerozoicPBDB<-downloadPBDB(CanonicalTaxa,"Phanerozoic")
# OBISData<-PhanerozoicPBDB

# Let's look at the size of an OBIS dataset
OBISData<-read.csv("OBISData.csv")
dim(OBISData)
	[1] 1373885      35 # 1.37e6 rows
# Check the memory usage of OBISData
object.size(OBISData)
	297177040 bytes # 0.29 GB

# Case 1: Subset the matrix by genus using split()
OBISGenera<-split(OBISData,OBISData$genus)

# Check memory
object.size(OBISData)+object.size(OBISGenera)
	29470143240 bytes # ~29.5 GB!!

# Case 2: Iteratively subset the matrix with for()
for (i in 1:length(unique(OBISData$genus))) {
	# Overwrites the last iteration of GenusSubset
	GenusSubset<-subset(OBISData,OBISData$genus==unique(OBISData$genus)[i])
	}

# Maximum possible memory usage during any one iteration is OBIS plus the largest
# Genus subset - in this case the genus Mytilus with 35,591 rows
Mytilus<-subset(OBISData,OBISData$genus=="Mytilus")
object.size(Mytilus)+object.size(OBISData)
	310778112 bytes # 0.3 GB 

# Effectively two orders of magnitude less RAM than using split()

# The RAM usage will be even worse when you use parallelization
# You are effectively multiplying the RAM burden of each loop by the number
# of parallel processes

# For example, even though the parallelRangeGenus() function is significantly
# faster than the other versions, I am incapable of running it even on a
# machine with 64 GB of RAM.

# Bottom line, there are trade-offs between speed, memory, and readability.
# Don't be dogmatic about how you approach the problem.

####################### Data Analysis: Age and Area #######################
# Look at our merged data
head(AgeAreaMatrix)
  		            LatRange	GridCount   Circle EarlyAge   AgeBin
	Acanthotrophon      5.12            2   737.13        8  Neogene
	Acorylus           12.21            3  1558.64       12  Neogene
	Adamussium         34.38           17  7631.04       12  Neogene
	Adeonella          99.33           10 13780.22        6  Neogene
	Adipicola         104.00            8 19875.02       16  Neogene
	Aesopus            78.82            7  8932.28       12  Neogene

# The most obvious test is to correlate area with age because both are 
# continuous variables
cor(AgeAreaMatrix$Circle,AgeAreaMatrix$EarlyAge,use="pairwise.complete.obs",method="spearman")
	[1] 0.3470801

# It's a positivie correlation, but nothing to write home about

# If we look at it in a plot we see that there are numerous ties for EarlyAge
plot(x=log(AgeAreaMatrix$Circle),y=AgeAreaMatrix$EarlyAge,las=1,pch=16,ylab="Age",xlab="Range")

# This is because most fossil ages are actually derived from time bins 
# This means that age is a de facto discrete variable, not continuous

####################### Demonstration of Binary Logistic Regression ###############
# We are going to use the glm() function to test whether the Circle Distance of a taxon
# Affects the probability that a taxon did or did not originate in the Neogene
# The Neogene is the past 23 million years.

# Neogene vs. Non-Neogene origination is our binary variable, which we represent as [0,1]
# Our choice of which state is 1 and which is 0 will affect our interpretation

# Logistic regression is just like traditional regression, except that the dependent
# variable is expressed as the logit.
	y = ax + b
	log(p/(1-p)) = ax + b
		
# Let's take a look at our data table again, but with a new column added named binary
# Where 1 represents a neogene age and 0 represents a non-Neogene age
AgeAreaMatrix[448:452,c("Circle","Binary","AgeBin")]

	               Circle   Binary     AgeBin
	Zeidora      13735.63        1    Neogene
	Zemitrella   10458.26        1    Neogene
	Zethalia      1392.69        1    Neogene
	Abra         18639.52        0  Paleogene
	Acanthastrea 19789.72        0  Paleogene

# Let's peform a logistic regression using glm(), which operates similar to lm
AgeAreaReg<-glm(AgeAreaMatrix$Binary~AgeAreaMatrix$Circle,family="binomial")

# Let's look at the output
summary(AgeAreaReg)

	Call:
	glm(formula = AgeAreaMatrix$Binary ~ AgeAreaMatrix$Circle, family = "binomial")

	Deviance Residuals: 
 	   Min       1Q   Median       3Q      Max  
	-1.1773  -0.7913  -0.6230   1.2178   1.9199  

	Coefficients:
	                          Estimate Std. Error z value Pr(>|z|)    
	(Intercept)              7.339e-04  1.097e-01   0.007    0.995    
	AgeArea$circle 			-8.354e-05  8.546e-06  -9.776   <2e-16 ***
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

	(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 1929.2  on 1642  degrees of freedom
	Residual deviance: 1830.4  on 1641  degrees of freedom
	AIC: 1834.4

	Number of Fisher Scoring iterations: 4

# For every kilometer of distance the log-odds of a Neogene
# origination decreases by -8.354e-05

# Some workers prefer to interpret log-odds as odds, in which case
# you just need to take the exponent of log-odds. For every kilometer
# of distance the odds of a Neogene origination changes by a
# factor of 0.9999165
exp(-8.354e-05)
	[1] 0.9999165

# I think that's still hard to conceptualize, so I like to convert the odds 
# into probability as a function of the independent variable(s).
	p = exp(ax+b)/(1+exp(ax+b))

# Plot probability for linear logistic regression
plotProbability<-function(Regression,Data) {
	Coefficients<-coefficients(Regression)
	A<-Coefficients[1]
	B<-Coefficients[2]
	X<-1:max(Data)
	Y<-exp(A+B*X)/(1+exp(A+B*X))
	quartz()
	plot(y=Y,x=X,las=1,type="l",lwd=3,ylim=c(0,1),ylab="Probability",xlab="Range")
	}

######################## Variations of Binary Logistic Regression ################
# Logistic regression is very similar to regular regression in that we can add
# additional or non-linear terms to the right hand side of the equation
NonLinear <- glm(y ~ x + I(x^2), family="binomial")
Multiple <- glm(y ~ x + z + q), family="binomial")

# If you want to try a complex model selection process under a variety of linear
# and non-linear models, you can use the eHOF package, which evaluates the following
# seven equations using bootstrapped and/or AIC model selection
y<-switch(Model,
	"I"=rep(1/(1+exp(a)),length(x)),
	"II"=1/(1+exp(a+b*x)),
	"III"=(1/(1+exp(a+b*x)))*(1/(1+exp(c))),
	"IV"=(1/(1+exp(a+b*x)))*(1/(1+exp(c-b*x))),
	"V"=(1/(1+exp(a+b*x)))*(1/(1+exp(c-d*x))),
	"VI"=1/(1+exp(a+b*x))*1/(1+exp(c-b*x))+1/(1+exp(a+b*(x-d)))*1/(1+exp(c-b*(x-d))),
	"VII"=1/(1+exp(a+b*x))*1/(1+exp(c-b*x))+1/(1+exp(a+b*(x-d)))*1/(1+exp(c-f*(x-d)))
	)

# Incidentally, the above switch() function of different logistic models
# is the equivalent of a series of if/elseif() statements.
# Switch is faster and easier to read. Excellent for making Macros.
# But it's not as flexible as if/else
	if (Model=="I") {
		y<-1/(1+exp(a))
		return(rep(y,length(TimeVector)))
		}
	else if (Model=="II") {
		y<-1/(1+exp(a+b*x))
		return(y)
		}
	else if (Model=="III") {
		y<-(1/(1+exp(a+b*x)))*(1/(1+exp(c)))
		return(y)
		}
	else if (Model=="IV") {
		y<-(1/(1+exp(a+b*x)))*(1/(1+exp(c-b*x)))
		return(y)
		}
	else if (Model=="V") {
		y<-(1/(1+exp(a+b*x)))*(1/(1+exp(c-d*x)))
		return(y)
		}
	else if (Model=="VI") {
		y<-1/(1+exp(a+b*x))*1/(1+exp(c-b*x))+1/(1+exp(a+b*(x-d)))*1/(1+exp(c-b*(x-d)))
		return(y)
		}
	else {
		y<-1/(1+exp(a+b*x))*1/(1+exp(c-b*x))+1/(1+exp(a+b*(x-d)))*1/(1+exp(c-f*(x-d)))
		return(y)
		}
	}

# Example of HOF in action
# HOF doesn't take the inputs as a formula
RegModels<-HOF(AgeAreMatrix$Binary,AgeAreaMatrix$Circle)

	Response of:  
	Deviances and information criteria:
    	  Deviance     logLik       AICc  AICc.Diff AICc.W BIC.Diff
	I   1929.17411 -964.58706 1931.17655  110.65732         99.8610
	II  1877.36230 -938.68115 1881.36961   60.85038         55.4534
	III 1814.50459 -907.25229 1820.51923    0.00000          0.0000
	IV  1834.92532 -917.46266 1840.93997   20.42073         20.4207
	V   1813.38469 -906.69235 1821.40911    0.88988          6.2844
	VI  1834.92532 -917.46266 1842.94974   22.43051         27.8250
	VII 1825.07992 -912.53996 1835.11657   14.59734         25.3839
	Percentage of model types after bootstrapping:
	III  IV   V VII 
 	64   2  32   2 
	Sum of bootstrapped model weights:
  	I  II III  IV   V  VI VII 
  	0   0   0   0   0   0   0 

	Suggested best model (AICc, bootselect.lower): III # A sigmoidal relationship

# I like the eHOF package a lot because it evaluates logistic functions that I 
# personally use a lot. That said, you can achieve similar results just by writing 
# your own functions in glm()

####################### Demonstration of Multinomial Logistic Regression ###############
# What if we want our dependent variable to be more than binary, i.e., multinomial?
# We can use the multinom() function from the nnet package.
# Great tutorial on this at http://www.ats.ucla.edu/stat/r/dae/mlogit.htm

# Three points I want to make before we proceed.

# As I mentioned earlier, I do not recommend the mlogit and mnlogit packages.
# There is another alternative, which I didn't mention - the VGAM package. 
# I've never used it though, so up to you.

# Multinomial logistic regression is sometimes also called
# discrete choice model, polychotomous, or polytomous logistic regression.
# Also, on the point of naming conventions, nobody calls it binary logistic regression.

# Technically for what we are doing we'd want ordered logistic regression,
# which is a special case of the multinomial, but that requires
# an entirely different function polr() in an entirely different package, i.e., MASS

# Remind ourselves of what the data looks like
head(AgeAreaMatrix)
  		            LatRange	GridCount   Circle EarlyAge   AgeBin
	Acanthotrophon      5.12            2   737.13        8  Neogene
	Acorylus           12.21            3  1558.64       12  Neogene
	Adamussium         34.38           17  7631.04       12  Neogene
	Adeonella          99.33           10 13780.22        6  Neogene
	Adipicola         104.00            8 19875.02       16  Neogene
	Aesopus            78.82            7  8932.28       12  Neogene

# Relevel the data so that our base variable is on the top
# Using order() will not work, must be relevel
AgeAreaMatrix$AgeBin<-factor(AgeAreaMatrix$AgeBin)
AgeAreaMatrix$AgeBin<-relevel(AgeAreaMatrix$AgeBin,ref="Neogene")

# Multinomial logistic regression in the nnet package works very similarly to
# all other types of regression in R
AgeAreaReg<-multinom(AgeBin~GridCount,data=AgeAreaMatrix)

# Examine a summary of the multinomial regression
summary(AgeAreaReg)
	
	Call: multinom(formula = AgeBin ~ GridCount, data = AgeAreaMatrix)

	Coefficients:
	          (Intercept) 	  GridCount
	Mesozoic   -0.8314961    0.07863459 # Notice that Neogene is not listed!
	Paleogene  -0.1687549    0.04330495
	Paleozoic  -3.5125508    0.10157700

	Std. Errors:
	          (Intercept)     GridCount
	Mesozoic   0.10407993   0.007405377
	Paleogene  0.09587166   0.007475564
	Paleozoic  0.21651278   0.008869269

	Residual Deviance: 3767.166 
	AIC: 3779.166 

# Just like with regular regression we can convert the log-odds into
# odds by taking the exponential
exp(coefficients(AgeAreaReg))
   		       (Intercept) 	  GridCount
	Mesozoic   0.43539741      1.081809
	Paleogene  0.84471593      1.044256 # Not very impressive!
	Paleozoic  0.02982075      1.106915

# We can also get associated confidence intervals using confint()
exp(confint(AgeAreaReg))

	, , Mesozoic
    	              2.5 %    97.5 %
	(Intercept)   0.3550527 0.5339233
	GridCount	  1.0662207 1.0976251 # between 6% and 9%

	, , Paleogene
    	              2.5 %   97.5 %
	(Intercept)   0.7000105 1.019335
	GridCount	  1.0290676 1.059669 # between 3% and 6%

	, , Paleozoic
    	               2.5 %     97.5 %
	(Intercept)   0.01950844 0.04558423
	GridCount	  1.08783943 1.12632537 # between 9% and 13%


# Unlike in the glm() function, multinom() does not give p-values
# So we have to calculate them using a Z-test
ZScore<-summary(AgeAreaReg)$coefficients/summary(AgeAreaReg)$standard.errors
PValues<-(1-pnorm(abs(ZScore),0,1))*2

round(PValues,4)
	          (Intercept) 	  GridCount
	Mesozoic       0.0000             0 # p is effectively zero
	Paleogene      0.0784             0
	Paleozoic      0.0000             0

# Convert to probabilities using predict()
AgeAreaProbs<-predict(AgeAreaReg,,"probs")

# Examine what the probabilities look like
head(AgeAreaProbs)
	                 Neogene  Mesozoic Paleogene  Paleozoic
	Acanthotrophon 0.4035596 0.2070969 0.3744939 0.01484957
	Acorylus       0.3899807 0.2164066 0.3777354 0.01587726
	Adamussium     0.2178281 0.3612456 0.3843825 0.03654386
	Adeonella      0.2985682 0.2864188 0.3903378 0.02467516
	Adipicola      0.3238416 0.2656854 0.3886107 0.02186234
	Aesopus        0.3367690 0.2555086 0.3871744 0.02054801

# Check that the probabilities always sum to 1 along each row
ProbSums<-apply(AgeAreaProbs,1,sum)
length(ProbSums[ProbSums!=1])
	[1] 0