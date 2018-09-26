\name{apm}
\alias{apm}
\title{Additive Partitioning of Net Biodiversity Effects}
\description{This function uses the two-way method of \cite{Loreau and Hector (2001)} or the tripartite method of \cite{Fox (2005)} to partition net biodiversity effects into complementarity and dominance/selection effects.}
\usage{
apm(mix, mono, ry=NULL, method="loreau")
}
\arguments{
  \item{mix}{A matrix or a data frame containing the biomass of each plant species (columns) in each mixture plot (lines). Column (species) names must be provided.}
  \item{mono}{A matrix or a data frame containing the biomass of each plant species (columns) in each monoculture plot (lines). Column (species) names must be provided.}
  \item{ry}{A matrix or a data frame giving the expected relative yield of each species (columns) in each mixture plot (lines). Column (species) names must be provided. See details.}
  \item{method}{The additive partitioning method to be used. Could be either \dQuote{loreau} (two-way partitioning method of \cite{Loreau and Hector (2001)}) or \dQuote{fox} (tripartite partitioning method of \cite{Fox (2005)}).}
}
\details{
The expected relative yield of each species growing in mixtures should be provided with the \code{ry} argument. By default (\code{ry=NULL}), all species in a given mixture are assumed to have the same expected relative yield (i.e., all species are assumed to have been sown/planted in equal proportions). If a matrix or a data frame is provided, it should exclusively contain values comprised between 0 and 1.

The two-way method of \cite{Loreau and Hector (2001)} partitions net biodiversity effects into a complementarity effect and a selection effect. The tripartite method of \cite{Fox (2005)}, however, partitions net biodiversity effects into a dominance effect, a trait-dependent complementarity effect, and a trait-independent complementarity effect. The selection effect of \cite{Loreau and Hector (2001)} is equal to the sum of the dominance and trait-dependent effects of \cite{Fox (2005)}, while the trait-independent complementarity effect of \cite{Fox (2005)} is equal to the complementarity effect calculated by the method of \cite{Loreau and Hector (2001)}. With the Fox method, a total complementarity effect can be calculated by summing trait-independent and trait-dependent complementarity effects. I strongly recommend the reading of \cite{Loreau and Hector (2001)} and \cite{Fox (2005)} for more information about how interpreting results given by these two partitioning methods.

For each species, the monoculture biomass is calculated as the mean value of all monoculture plots.
}
\value{Returns a matrix with three (Loreau method) or four (Fox method) columns.

If \code{method="loreau"}, the matrix has the following columns:
\itemize{
\item \emph{NBE}: net biodiversity effect
\item \emph{CE}: complementarity effect
\item \emph{SE}: selection effect}

If \code{method="fox"}, the matrix has the following columns:
\itemize{
\item \emph{NBE}: net biodiversity effect
\item \emph{TICE}: trait-independent complementarity effect
\item \emph{TDCE}: trait-dependent complementarity effect
\item \emph{DE}: dominance effect}
}
\references{
Loreau M., Hector A. (2001) Partitioning selection and complementarity in biodiversity experiments, \emph{Nature}, DOI: 10.1038/35083573.

Fox J. (2005) Interpreting the selection effect of biodiversity on ecosystem function, \emph{Ecology Letters}, DOI: 10.1111/j.1461-0248.2005.00795.x.
}
\author{Benjamin M. Delory
}
\examples{
path<-system.file("extdata", package = "bef")
datamix<-read.table(paste(path, "datamix.txt", sep="/"), header=TRUE, sep=",")
datamono<-read.table(paste(path, "datamono.txt", sep="/"), header=TRUE, sep=",")
addpart<-apm(mix=datamix, mono=datamono, method="fox")
}