% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbr_gweis.R, R/bqr_gweis.R
\name{br_gwas}
\alias{br_gwas}
\title{Perform GWAS for covariates.}
\usage{
br_gwas(plink_path, dis_snp, br_dis_phen, output_dir)

br_gwas(plink_path, dis_snp, br_dis_phen, output_dir)
}
\arguments{
\item{plink_path}{Path to the PLINK executable application.}

\item{dis_snp}{Prefix for binary files for the discovery dataset.}

\item{br_dis_phen}{File path for covariate data in the discovery dataset.}

\item{output_dir}{Directory to save output files.}

\item{br_dis_cov}{File path for covariate data in the discovery dataset.}
}
\value{
None. Results are saved to files.

None. Results are saved to files.
}
\description{
Perform GWAS for covariates.

Perform GWAS for covariates.
}
