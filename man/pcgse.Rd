\name{pcgse}
\alias{pcgse}
\title{Principal component gene set enrichment (PCGSE) algorithm}
\description{
Implementation of the PCGSE algorithm. Computes the statistical association between gene sets and the principal components of 
experimental data using a two-stage competitive test. Supported gene-level test statistics include the
PC loadings for each genomic variable, the Pearson correlation coefficients between each genomic variable and each PC
and the Fisher-transformed correlation coefficients. 
The input data is centered and scaled so that eigendecomposition is computed on the sample
correlation matrix rather than the sample covariance matrix. 
Because the PC loadings for PCA on a correlation matrix are proportional to the
Pearson correlation coefficients between each PC and each variable, all supported gene-level statistics provide a measure of correlation
between genomic variables and PCs.
Each gene set is quantified using either a standardized mean difference statistic or a standardized rank sum statistic. 
The statistical significance of each gene set test statistic is computed according to a competitive null hypothesis using 
either a parametric test, a correlation-adjusted parametric test or a permutation test.
}
\usage{
    pcgse(data, prcomp.output, pc.indexes=1, gene.sets, gene.statistic="z",
        transformation="none", gene.set.statistic="mean.diff", 
        gene.set.test="cor.adj.parametric", nperm=9999)    
}
\arguments{
\item{data}{Empirical data matrix, observations-by-variables. Must be specified. Cannot contain missing values.}
\item{prcomp.output}{Output of prcomp(data,center=T,scale=T). If not specified, it will be computed.}
\item{pc.indexes}{Indices of the PCs for which enrichment should be computed. Defaults to 1.}
\item{gene.sets}{Data structure that holds gene set membership information. 
Must be either a binary membership matrix or a list of gene set member indexes. 
For the member matrix, rows are gene sets, columns are genes, elements are binary membership values.
For the membership index list, each element of the list represents a gene set and holds a vector of indexes of genes that are members. 
Must be a matrix if gene.set.test is set to "permutation".}
\item{gene.statistic}{The gene-level statistic used to quantify the association between each genomic variable and each PC.
Must be one of the following (default is "z"):
    \itemize{
    \item "loading": PC loading associated with the genomic variable.
    \item "cor": Pearson correlation coefficient between the PC and the genomic variable.
    \item "z": Fisher-transformed Pearson correlation coefficient.
    }
}
\item{transformation}{Optional transformation to apply to the gene-level statistics. 
Must be one of the following (default is "none"):
    \itemize{
    \item "none": No transformations are applied to the gene-level statistics.
    \item "abs.value": The absolute value of the gene-level statistics is used.
    }
}
\item{gene.set.statistic}{The gene set statisic computed from the gene-level statistics. 
Must be one of the following (default is "mean.diff"):
    \itemize{
    \item "mean.diff": The standardized difference between the mean of the gene-level statistics for members of the gene set and the mean
             of the gene-level statistics for genomic variables not in the gene set. Equivalent to the U_D statistic from Barry et al.
    \item "rank.sum": The standardized Wilcoxon rank sum statistic computed from the gene-level statistics for members of the gene set. 
            Equivalent to the U_W statistic from Barry et al.
            } 
}
\item{gene.set.test}{The statistical test used to compute the significance of the gene set statistics under a competitive null hypothesis.
The "parametric" test is in the "class 1" test category according to Barry et al., 
the "cor.adj.parametric" and "permutation" tests are in the "class 2" test category according to Barry et al.
Must be one of the following (default is "cor.adj.parametric"):
    \itemize{             
    \item "parametric": If the mean difference is being used as the gene set statistic, corresponds to a two-sided, two-sample t-test with equal variances.
             If the rank sum is being used as the gene set statistic, this corresponds to a two-sided,
             two-sample z-test based on the standardized rank sum statistic. NOTE: both of these tests incorrectly assume the gene-level statistics are i.i.d. 
             and should therefore only be used for comparative purposes.
    \item "cor.adj.parametric": Tests statistical significance of the standardized and correlation-adjusted gene set statistic using a two-sided t-test or z-test.  
             Similar to the CAMERA method by Wu et al., standardization of either the mean different statistic or rank sum statistic is performed   
             using a variation inflation factor based on the average pair-wise correlation between the gene-level statistics for members of the gene set. 
             Per Barry et al., this is approximated by the average correlation between the genomic variables. 
             Per Wu et al., the significance of the correlation-adjusted t-statistic is tested using a two-sided t-test with n-2 df and 
             the significance of the correlation-adjusted rank sum z-statistic is tested using a two-sided z-test. 
    \item "permutation": Tests gene set enrichment via the permutation distribution of the gene set statistic.
             The permutation distribution is computed via permutation of the sample labels, which, in this case, is equivalent to permutation 
             of the elements of the target PC. This test is realized using the safe() function from the R safe package.
             The number of permutations is controlled by the "nperm" parameter. The gene.statistic cannot be set to "loadings" with this option.
             Per Barry et al., this correlation is approximated by the average correlation between the genomic variables.
             This option can be extremely computationally expensive so should not be used for most applications.
             }
}
\item{nperm}{Number of permutations to perform. Only relevant if gene.set.test is set to "permutation".}
}
\value{
List with the following elements:
    \itemize{
      \item p.values: Matrix with one row per gene set and one column for each tested PC. 
             Elements are the two-sided competitive enrichment p-values. Multiple hypothesis correction is NOT applied to these p-values.
      \item statistics: Matrix with one row per gene set and one column for each tested PC. 
             Elements are the gene set test statistics for each gene set.
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

   ## Simulate MVN data with two population PCs whose loadings are
   ## associated with the first and second gene sets, respectively.
   var1=2 ## variance of first population PC
   var2=1 ## variance of second population PC 
   default.var=.1 ## background variance of population PCs
   load = sqrt(.1) ## value of population loading vector for gene set 1 on PC 1 and set 2 on PC 2
   
   ## Creates a first PC with loadings for just the first 20 genes and a 
   ## second PC with loadings for just the second 20 genes 
   loadings1 = c(rep(load,p/f), rep(0,p-p/f))
   loadings2 = c(rep(0,p/f), rep(load, p/f), rep(0, p-2*p/f))

   ## Create the population covariance matrix
   sigma = var1 * loadings1 \%*\% t(loadings1) + var2 * loadings2 \%*\% t(loadings2) + 
    diag(rep(default.var, p))

   ## Simulate MVN  data
   data = mvrnorm(n=n, mu=rep(0, p), Sigma=sigma)  
 
   ## Perform PCA on the standardized data
   prcomp.output = prcomp(data, center=TRUE, scale=TRUE)
 
   ## Execute PCGSE using Fisher-transformed correlation coefficients as the gene-level statistics,
   ## the standardized mean difference as the gene set statistic and an unadjusted two-sided, 
   ## two-sample t-test for the determination of statistical significance.
   pcgse.results = pcgse(data=data, 
                         prcomp.output=prcomp.output, 
                         pc.indexes=1:2,
                         gene.sets=gene.sets,
                         gene.statistic="z", 
                         transformation="none",
                         gene.set.statistic="mean.diff",
                         gene.set.test="parametric")
   
   ## Apply Bonferroni correction to p-values
   for (i in 1:2) {
    pcgse.results$p.values[,i] = p.adjust(pcgse.results$p.values[,i], method="bonferroni")
   }
   
   ## Display the p-values for the first 5 gene sets for PCs 1 and 2 
   pcgse.results$p.values[1:5,]
   
   ## Execute PCGSE again but using a correlation-adjusted t-test
   pcgse.results = pcgse(data=data, 
                         prcomp.output=prcomp.output,
                         pc.indexes=1:2,
                         gene.sets=gene.sets,
                         gene.statistic="z", 
                         transformation="none",
                         gene.set.statistic="mean.diff",
                         gene.set.test="cor.adj.parametric")                         
      
   ## Apply Bonferroni correction to p-values
   for (i in 1:2) {
    pcgse.results$p.values[,i] = p.adjust(pcgse.results$p.values[,i], method="bonferroni")
   }
   
   ## Display the p-values for the first 5 gene sets for PCs 1 and 2 
   pcgse.results$p.values[1:5,]   
}
\keyword{file}