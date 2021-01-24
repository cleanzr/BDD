
#' Test data for Record Linkage
#'
#' This table contains artificial personal data for the 
#' evaluation of Record Linkage procedures. 50 records are duplicates
#' with randomly generated errors. 
#'
#' @format `RLdata500` is a  data frame with 500 records and 7 variables:
#' \describe{
#'  \item{fname_c1}{First name, first component}
#'  \item{fname_c2}{First name, second component}
#'  \item{lname_c1}{Last name, first component}
#'  \item{lname_c2}{Last name, second component}
#'  \item{by}{Year of birth}
#'  \item{bm}{Month of birth}
#'  \item{bd}{Day of birth}
#' }
#' `identity.RLdata500` is a vector representing the true record ids of the 
#' data set. A pair of records are duplicates, if and only if their 
#' corresponding values in the identity vector agree.
#' 
#' @source 
#' Generated with the data generation component of Febrl (Freely Extensible 
#' Biomedical Record Linkage), version 0.3 
#' <https://sourceforge.net/projects/febrl/>.
#' 
#' The following data sources were used (all relate to Germany):
#' 
#' * <https://blog.beliebte-vornamen.de/2009/02/prozentuale-anteile-2008/>, a 
#'   list of the frequencies of the 20 most popular female names in 2008.
#' * <https://www.beliebte-vornamen.de/760-alle_jahre.htm>, a list of the 
#'   100 most popular first names since 1890. The frequencies found in
#'   the source above were extrapolated to fit this list.
#' * <http://www.ahnenforschung-in-stormarn.de/geneal/nachnamen_100.htm>, 
#'   a list of the 100 most frequent family names with frequencies.
#'   
#' * Age distribution as of Dec 31st, 2008, statistics of Statistisches 
#'   Bundesamt Deutschland, taken from the GENESIS database 
#'   <https://www-genesis.destatis.de/genesis/online/logon>.
#' 
#' Web links as of October 2009.
#' 
#' @author Andreas Borg
#' @name RLdata500
NULL

#' @rdname RLdata500
"RLdata500"

#' @rdname RLdata500
"identity.RLdata500"
