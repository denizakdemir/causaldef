#' Hematopoietic Cell Transplantation Outcomes
#'
#' Simulated survival data comparing conditioning intensities in HCT.
#' Illustrates competing risks (Relapse vs Death).
#'
#' @format A data frame with 800 rows and 7 variables:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{age}{Patient age in years}
#'   \item{disease_status}{Disease stage (Early, Intermediate, Advanced)}
#'   \item{kps}{Karnofsky Performance Score (0-100)}
#'   \item{donor_type}{HLA matching status}
#'   \item{conditioning_intensity}{Treatment group (Myeloablative vs Reduced)}
#'   \item{time_to_event}{Time to first event or censoring}
#'   \item{event_status}{Type of event (Death, Relapse, Censored)}
#' }
#' 
#' @details 
#' This dataset mimics a retrospective registry study where treatment assignment
#' is confounded by patient health status (e.g., sicker patients get Reduced intensity).
#' 
#' @usage data(hct_outcomes)
"hct_outcomes"

#' Gene Expression Perturbation Data
#'
#' Simulated high-throughput screening data with potential batch effects.
#' Contains a negative control outcome ("housekeeping_gene").
#'
#' @format A data frame with 500 rows and 6 variables:
#' \describe{
#'   \item{sample_id}{Sample identifier}
#'   \item{batch}{Experimental batch (1-4)}
#'   \item{library_size}{Total read count}
#'   \item{knockout_status}{Intervention (Control vs Knockout)}
#'   \item{target_expression}{Expression level of target gene (outcome)}
#'   \item{housekeeping_gene}{Expression of unrelated housekeeping gene (negative control)}
#' }
#' 
#' @details 
#' Used to demonstrate negative control diagnostics. The housekeeping gene shares
#' batch effects and library size confounders with the target gene but is not
#' biologically affected by the knockout.
#' 
#' @usage data(gene_perturbation)
"gene_perturbation"
