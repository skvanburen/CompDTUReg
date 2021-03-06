% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{maketx2gene}
\alias{maketx2gene}
\title{Create a tx2gene file from a .gff or .gtf file that maps each transcript to its respective gene}
\usage{
maketx2gene(txdb_loc, save_loc = NULL)
}
\arguments{
\item{txdb_loc}{is the file location of the .gtf or .gff annotation}

\item{save_loc}{An optional directory location to save the tx2gene data.frame (will save to current working directory otherwise)}
}
\value{
a data.frame with columns gene_id, tx_id (transcript id) and NTrans (the number of transcripts in the current transcript's gene)
}
\description{
\code{maketx2gene} creates a tx2gene file from a .gff or .gtf file that maps each transcript to its respective gene
}
