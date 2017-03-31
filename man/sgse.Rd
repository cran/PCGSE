\name{sgse}
\alias{sgse}
\title{Spectral gene set enrichment (SGSE) algorithm}
\description{
Implementation of the SGSE algorithm. 
Computes the statistical association between gene sets and the spectra of the specified data set. 
The association between each gene set and each PC is first computed using the \code{pcgse} function. 
The PC-specific p-values are then combined using the weighted Z-method with weights set to either the 
PC variance or the PC variance scaled by the lower-tailed p-value calculated for the variance according 
to the Tracey-Widom distribution.
}
\usage{
    sgse(data, prcomp.output, gene.sets, 
         gene.statistic="z", transformation="none", 
         gene.set.statistic="mean.diff", gene.set.test="cor.adj.parametric", 
         nperm=999, pc.selection.method="all", pc.indexes=NA, rmt.alpha=.05, 
         pcgse.weight="rmt.scaled.var")    
}
\arguments{
\item{data}{Empirical data matrix, observations-by-variables. Must be specified. Cannot contain missing values.}
\item{prcomp.output}{Output of prcomp(data,center=T,scale=T). If not specified, it will be computed.}
\item{gene.sets}{See documentation for \code{gene.sets} argument for \code{pcgse} function.}
\item{gene.statistic}{See documentation for \code{gene.statistic} argument for \code{pcgse} function.}
\item{transformation}{See documentation for \code{transformation} argument for \code{pcgse} function.}
\item{gene.set.statistic}{See documentation for \code{gene.set.statistic} argument for \code{pcgse} function.}
\item{gene.set.test}{See documentation for \code{gene.set.test} argument for \code{pcgse} function.}
\item{nperm}{See documentation for \code{nperm} argument for \code{pcgse} function.}
\item{pc.selection.method}{Method used to determine the PCs for which enrichment will be computed. 
Must be one of the following:
         \itemize{
         \item "all": All PCs with non-zero variance will be used.
         \item "specific": The set of PCs specified by pc.indexes will be used.
         \item "rmt": The set of PCs with significant eigenvalues according to the Tracy-Widom distribution for a white Wishart at the alpha specified by the
                      "rmt.alpha" parameter.
         }
}
\item{pc.indexes}{Indices of the PCs for which enrichment should be computed. Must be specified if pc.selection.method is "specific".}
\item{rmt.alpha}{Significance level for selection of PCs according to the Tracy-Widom distribution. Must be specified if pc.selection.method is "rmt".}
\item{pcgse.weight}{Type of weight to use with the weighted Z-method to combine the p-values from the PCGSE tests on all PCs selected according to the 
pc.selection.method parameter value. Must be one of the following:
        \itemize{
        \item "variance": The PC variance is used as the weight. NOTE: this should only be used for evaluation and testing.
        \item "rmt.scaled.var": The product of the PC variance and the Tracey-Widom lower-tailed p-value for the eigenvalue associated with the PC is used as the weight. 
        }
}
}
\value{
List with the following elements:
    \itemize{
      \item "pc.indexes": Indices of the PCs on which enrichment was performed. 
      \item "pcgse": Output from \code{pcgse} function on the PCs identified by pc.indexes. 
      \item "sgse": Vector of combined p-values for all PCs identified by pc.indexes.
      \item "weights": Vector of PC-specific weights for the PCs identified by pc.indexes.      
     }
}
\examples{

   library(MASS)

   p=200 ## number of genomic variables
   n=50 ## number of observations   
   f=20 ## number of gene sets

   ## Create annotation matrix with disjoint gene sets
   gene.sets = matrix(0, nrow=f, ncol=p)
   for (i in 1:f) {
      gene.sets[i, ((i-1)*p/f + 1):(i*p/f)] = 1 
   }

   ## Simulate MVN data where the
   ## first population PC loadings are
   ## associated with the first gene set.
   var=2 ## variance of first population PC
   default.var=.1 ## background variance of population PCs
   load = sqrt(.1) ## value of population loading vector for gene set 1 on PC 1
   
   ## Creates a first PC with loadings for just the first 20 genes and a 
   loadings = c(rep(load,p/f), rep(0,p-p/f))

   ## Create the population covariance matrix
   sigma = var * loadings \%*\% t(loadings) + diag(rep(default.var, p))

   ## Simulate MVN  data
   data = mvrnorm(n=n, mu=rep(0, p), Sigma=sigma)  
 
   ## Perform PCA on the standardized data
   prcomp.output = prcomp(data, center=TRUE, scale=TRUE)
 
   ## Execute SGSE using Fisher-transformed correlation coefficients as 
   ## the gene-level statistics, the standardized mean difference as the 
   ## gene set statistic and a correlation adjusted two-sided, 
   ## two-sample t-test for the determination of statistical significance,
   ## all PCs with non-zero eigenvalues for spectral enrichment and 
   ## variance weights
   sgse.results = sgse(data=data, 
                       prcomp.output=prcomp.output, 
                       gene.sets=gene.sets,
                       gene.statistic="z", 
                       transformation="none",
                       gene.set.statistic="mean.diff",
                       gene.set.test="cor.adj.parametric",
                       pc.selection.method="all",
                       pcgse.weight="variance")
   
   ## Display the PCGSE p-values for the first 5 gene sets for PC 1 
   sgse.results$pcgse$p.values[1:5,1]
   
   ## Display the SGSE weights for the first 5 PCs 
   sgse.results$weights[1:5]   
   
   ## Display the SGSE p-values for the first 5 gene sets 
   sgse.results$sgse[1:5]   
   
   ## Execute SGSE again but using RMT scaled variance weights
   sgse.results = sgse(data=data, 
                       prcomp.output=prcomp.output, 
                       gene.sets=gene.sets,
                       gene.statistic="z", 
                       transformation="none",
                       gene.set.statistic="mean.diff",
                       gene.set.test="cor.adj.parametric",
                       pc.selection.method="all",
                       pcgse.weight="rmt.scaled.var")

   ## Display the SGSE weights for the first 5 PCs 
   sgse.results$weights[1:5]   
                       
   ## Display the SGSE p-values for the first 5 gene sets 
   sgse.results$sgse[1:5]                          
   
   ## Execute SGSE again using RMT scaled variance weights and  
   ## all RMT-significant PCs at alpha=.05
   sgse.results = sgse(data=data, 
                       prcomp.output=prcomp.output, 
                       gene.sets=gene.sets,
                       gene.statistic="z", 
                       transformation="none",
                       gene.set.statistic="mean.diff",
                       gene.set.test="cor.adj.parametric",
                       pc.selection.method="rmt",
                       rmt.alpha=.05,
                       pcgse.weight="rmt.scaled.var")
                       
   ## Display the indexes of the RMT-significant PCs
   sgse.results$pc.indexes
                       
   ## Display the SGSE p-values for the first 5 gene sets 
   sgse.results$sgse[1:5]                             

}       
\seealso{\code{\link{pcgse}}}                  
\keyword{file}