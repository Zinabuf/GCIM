% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/b_gweis.R
\name{b_gweis}
\alias{b_gweis}
\title{Perform GWEIS for binary outcome variables.}
\usage{
b_gweis(plink_path, dis_snp, bp_dis_phen, bp_dis_cov)
}
\arguments{
\item{plink_path}{Path to the PLINK executable application.}

\item{dis_snp}{Prefix for binary files for the discovery dataset.}

\item{bp_dis_phen}{File path for the phenotype data in the discovery dataset.}

\item{bp_dis_cov}{File path for covariate data in the discovery dataset.}

\item{temp_dir}{Directory to save output files.}
}
\value{
None. Results are saved to files.
}
\description{
This function performs genome-wide-by-environemnt interaction studies (GWEIS) and processes the results
to generate files for downstream analysis.
}
