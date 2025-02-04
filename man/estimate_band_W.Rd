% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nhpp_functions.R
\name{estimate_band_W}
\alias{estimate_band_W}
\title{Estimate bandwidth}
\usage{
estimate_band_W(
  tree.discr,
  data.path,
  band.width = c("bw.nrd0", "bw.nrd0", "bw.ucv", "bw.bcv", "bw.SJ")
)
}
\arguments{
\item{tree.discr}{simmap or phylo object. A discretized tree using the 'discr_Simmap' function.}

\item{data.path}{list. A list of path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.}

\item{band.width}{character. Bandwidth selectors for the KDE, as in density.}
}
\value{
A numeric vector.
}
\description{
Estimate the bandwidth for the Markov KDE.
}
\examples{
data("hym_hm", "hym_tree")
# Get reference tree.
tree_discr <- discr_Simmap(hym_tree, res = 4000)
# Get hamming data from the head characters.
hm <- hym_hm$head
# Make NHPP path data.
nhpp <- make_data_NHPP_KDE_Markov_kernel(hm)
# Add pseudo data to path data.
psd <- lapply(nhpp, function(x) -x[x < 100] )
edge_groups <- as.list(1:length(hym_tree$edge.length))
nhpp_psd <- add_pseudodata(Edge.groups = edge_groups, Pseudo.data = psd, Path.data = nhpp)
# Calculate bandwidth.
bdw <- estimate_band_W(tree_discr, nhpp_psd, band.width = "bw.nrd0")
mean(bdw)

}
\author{
Sergei Tarasov
}
