% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GCIM_quantitative.R
\name{extract_prs_data_q}
\alias{extract_prs_data_q}
\title{Extract PRS data from GCIM quantitative results for further analysis}
\usage{
extract_prs_data_q(gcim_result)

extract_prs_data_q(gcim_result)
}
\arguments{
\item{gcim_result}{Result object from gcim_q function}
}
\value{
List containing the three PRS datasets (add_prs, int_prs, cov_prs)

List containing the three PRS datasets (add_prs, int_prs, cov_prs)
}
\description{
This helper function extracts the individual PRS datasets from GCIM quantitative results
for further analysis or inspection.

This helper function extracts the individual PRS datasets from GCIM quantitative results
for further analysis or inspection.
}
\examples{
\dontrun{
result <- gcim_q("Qphe_target.txt", "Qexp_target.txt", 
                 Add_PRS = q, Int_PRS = r, Cov_PRS = p)
prs_data <- extract_prs_data_q(result)
add_prs <- prs_data$add_prs
int_prs <- prs_data$int_prs  
cov_prs <- prs_data$cov_prs
}
\dontrun{
result <- gcim_q("Qphe_target.txt", "Qexp_target.txt", 
                 Add_PRS = q, Int_PRS = r, Cov_PRS = p)
prs_data <- extract_prs_data_q(result)
add_prs <- prs_data$add_prs
int_prs <- prs_data$int_prs  
cov_prs <- prs_data$cov_prs
}
}
