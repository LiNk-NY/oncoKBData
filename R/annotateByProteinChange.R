#' Pull annotations for mutations by protein change
#'
#' @examples
#' library(cBioPortalData)
#' cbio <- cBioPortal()
#' mutations <- getDataByGenes(
#'     cbio, studyId = "acc_tcga", genes = 1:3, by = "entrezGeneId",
#'     molecularProfileIds = "acc_tcga_mutations",
#'     sampleListId = "acc_tcga_all"
#' )[["acc_tcga_mutations"]]
#'
#' ## Annotate the protein change
#' ok <- oncoKB(token = Sys.getenv("ONCOKB_TOKEN"))
#' annotateByProteinChange(
#'     ok, mutations, by = "entrezGeneId", getOncoTreeCODES = TRUE
#' )
#'
#' @export
annotateByProteinChange <-
    function(
        api,
        mutations,
        by = c("entrezGeneId", "hugoSymbol"),
        getOncoTreeCODES = TRUE
    )
{
    studyid <- unique(mutations[["studyId"]])
    if (!identical(length(studyid), 1L))
        stop("Only one 'studyId' is allowed")

    if (getOncoTreeCODES) {
        clin <- clinicalData(cBioPortal(), studyId = studyid)
        mutations <- merge(mutations, clin, by = "patientId")
    }

    annotations <- apply(mutations, 1L, function(row) {
       httr::content(
            api$annotateMutationsByProteinChangeGetUsingGET_1(
                entrezGeneId = row[by],
                referenceGenome = row["ncbiBuild"],
                proteinStart = row["proteinPosStart"],
                proteinEnd = row["proteinPosEnd"],
                tumorType = row["ONCOTREE_CODE"]
            )
        )
    })
    annotations <- lapply(annotations, function(x) {
        x <- x[names(x) != "query"]
        x <- Filter(length, x)
        x[nzchar(x)]
    })
    annotations <- dplyr::bind_rows(annotations)
    dplyr::bind_cols(mutations, annotations)
}
