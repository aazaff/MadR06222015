# Custom functions are camelCase. Arrays and Arguments are PascalCase
# Dependency functions are not embedded in master functions
# []-notation is used wherever possible, and $-notation is avoided.

######################################### Load Required Libraries ###########################################
# Load velociraptr
if (suppressWarnings(require("velociraptr"))==FALSE) {
        install.packages("velociraptr",repos="http://cran.cnr.berkeley.edu/");
        library("velociraptr");
        }
        
# Load vegan
if (suppressWarnings(require("vegan"))==FALSE) {
        install.packages("vegan",repos="http://cran.cnr.berkeley.edu/");
        library("vegan");
        }

#############################################################################################################
######################################### GENERATE DATA FUNCTIONS, SQS ######################################
#############################################################################################################
# Use the Jones, M. C. and Pewsey A. (2009). Sinh-arcsinh distributions. Biometrika 96: 761–780.
# To generate a half-normal distribution with increasingly changing kurtoses
sinArc<-function(Points,Epsilon=0,Delta=1) {
        return(dnorm(sinh(Delta*asinh(Points)-Epsilon))*Delta*cosh(Delta*asinh(Points)-Epsilon)/sqrt(1+Points^2))
        }

# Generate a matrix of different "samples" with increasing un-evenness by row
populateMatrix<-function(Species,Start=0.1,End=2,Increments=0.1) {
        Deltas<-seq(Start,End,Increments)
        FinalMatrix<-matrix(NA,nrow=length(Deltas),ncol=length(Species))
        for (i in 1:length(Deltas)) {
                FinalMatrix[i,]<-sinArc(Species,0,Deltas[i])
                }
        return(ceiling(FinalMatrix*100)) # Use ceiling to ensure integer values
        }

theilEntropy<-function(Abundances) {
        Richness<-length(Abundances)
        Top<-Abundances/Richness
        Bottom<-1/Richness
        Quotient<-Top/Bottom
        Exponent<-Quotient^Abundances
        Solution<-log(Exponent^(1/Richness))
        return(Solution)
        }

########################################## GENERATE DATA SCRIPTS, SQS #######################################
# Determine how many species you'd like to simulate
Species<-seq(0,2,0.01)
Richness<-length(Species)

# Generate samples according to increasingly kurtotic half-normal
# This does not guarantee a consistent sample size, but you could manually approximate
# Relatively constant sample size across samples by multiplying with a constant
# This should not change the results, however.
PopulatedMatrix<-populateMatrix(Species,0.1,3,0.1)

# Calculate the standardized richness of each sample
BiodiversitySQS<-apply(PopulatedMatrix,1,velociraptr::subsampleEvenness,0.75,ExcludeDominant=TRUE)

# Calculate pielou's evenness for each sample. You could also look at the deltas as a more direct measure.
Pielou<-vegan::diversity(PopulatedMatrix)/log(vegan::specnumber(PopulatedMatrix))
# Calculate the Theil Entropy (evenness) for each sample.
Theil<-apply(PopulatedMatrix,1,theilEntropy)

# Plot the relationship
plot(y=BiodiversitySQS,x=Pielou,xlab="pielou's evennes",ylab="standardized biodiversity",pch=16,cex=1.5,las=1)
plot(y=BiodiversitySQS,x=Theil,xlab="pielou's evennes",ylab="standardized biodiversity",pch=16,cex=1.5,las=1)
