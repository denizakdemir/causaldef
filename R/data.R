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

#' Lalonde National Supported Work (NSW) Benchmark
#'
#' The classical causal inference benchmark dataset constructed by Dehejia and Wahba (1999)
#' from the original Lalonde (1986) study. This version combines the experimental
#' NSW groups with the observational CPS and PSID comparison groups, enabling precise
#' calculation of the "deficiency" between observational and experimental inference.
#'
#' @format A data frame with 2915 rows and 11 variables:
#' \describe{
#'   \item{treat}{Treatment indicator (1 = Job Training, 0 = Control)}
#'   \item{age}{Age in years}
#'   \item{education}{Years of schooling}
#'   \item{black}{Indicator for Black race}
#'   \item{hispanic}{Indicator for Hispanic race}
#'   \item{married}{Indicator for marital status}
#'   \item{nodegree}{Indicator for no high school degree}
#'   \item{re74}{Real earnings in 1974 (pre-treatment)}
#'   \item{re75}{Real earnings in 1975 (pre-treatment)}
#'   \item{re78}{Real earnings in 1978 (outcome)}
#'   \item{sample_id}{Source of the observation: "nsw_treated", "nsw_control", 
#'     "cps_control", or "psid_control".}
#' }
#' 
#' @details 
#' This dataset allows you to verify causal methods by using the experimental subset 
#' ("nsw_treated" vs "nsw_control") as the ground truth, and comparing the results 
#' obtained by adjusting the observational subsets ("nsw_treated" vs "cps_control").
#'
#' @references
#' LaLonde, R. J. (1986). Evaluating the econometric evaluations of training programs with experimental data. 
#' The American Economic Review, 604-620.
#'
#' Dehejia, R. H., & Wahba, S. (1999). Causal effects in nonexperimental studies: 
#' Reevaluating the evaluation of training programs. Journal of the American statistical Association, 94(448), 1053-1062.
#' 
#' @usage data(nsw_benchmark)
"nsw_benchmark"

#' Right Heart Catheterization (RHC) Dataset
#'
#' Data from the SUPPORT study (Connors et al., 1996) examining the effectiveness of 
#' Right Heart Catheterization (RHC) in the management of critically ill patients.
#' This dataset is a classic example of high-dimensional confounding in medical observational studies.
#'
#' @format A data frame with 5735 rows and 63 variables.
#' 
#' @details 
#' The original study found that RHC was associated with higher mortality, contradicting 
#' potential benefits. Careful adjustment for the rich set of covariates (indicating sickness severity) 
#' is required. This dataset serves as a testbed for sensitivity analysis and policy bounds.
#'
#' Key variables include \code{swang1} (treatment), \code{dth30} (30-day mortality),
#' \code{t3d30} (survival time up to 30 days), and \code{aps1} (APACHE III score).
#'
#' @references
#' Connors, A. F., Speroff, T., Dawson, N. V., Thomas, C., Harrell, F. E., Wagner, D., ... & Goldman, L. (1996). 
#' The effectiveness of right heart catheterization in the initial care of critically ill patients. 
#' JAMA, 276(11), 889-897.
#'
#' @usage data(rhc)
"rhc"
