% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/spader.R
\name{Diversity}
\alias{Diversity}
\title{Etimation of species diversity}
\usage{
Diversity(X, datatype = c("abundance", "incidence"))
}
\arguments{
\item{X}{a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).}

\item{datatype}{a character of data type,"abundance" or "incidence".}
}
\value{
a list of species diversity estimator with order q and its confidence interval.
}
\description{
\code{Diversity} This part features various diversity indices including the Shannon??s index and its effective number of species (diversity of order 1, or Shannon diversity), the Simpson??s index and its effective number of species (diversity order 2, or Simpson diversity), species richness (diversity of order 0).
}
\examples{
\dontrun{
data(DivDemoAbu)
Diversity(DivDemoAbu,datatype="abundance")
}
}
\author{
Anne Chao, K. H. Ma and T. C. Hsieh
}

