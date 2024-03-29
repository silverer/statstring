% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_strings.R
\name{format_confint}
\alias{format_confint}
\title{Wrapper to format raw point estimate and lower and upper CI}
\usage{
format_confint(num, lower.ci, upper.ci)
}
\arguments{
\item{num}{A point estimate}

\item{lower.ci}{Lower bound of a confidence interval}

\item{upper.ci}{Upper bound of a confidence interval}
}
\value{
An APA-formatted string of CI results
}
\description{
This function takes a point estimate, lower CI, and upper CI and returns
a formatted string with the results
}
