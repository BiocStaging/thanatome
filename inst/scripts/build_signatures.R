## Assisted by Claude Code (Anthropic); see README.md, section Provenance.
## Builds data/thanatome_signatures.rda from the consensus tables of the
## analysis repository (results/firmas/*_consenso.tsv). See the README.md of
## this folder.
source("inst/scripts/record_session.R")
repo <- Sys.getenv("GTEX_CONTROLS")
if (!nzchar(repo)) {
    stop("set GTEX_CONTROLS to a checkout of the analysis repository ",
         "(or to the folder with the deposited consensus tables)")
}
min_tissues <- 5
## |median lfc|: 0.5 for ventilator and agony (validation threshold); 0.25 for
## the ischemia component, whose contrast is per standard deviation of ischemic
## time (with 0.5 only 21 almost unexpressed genes remain). Ischemia is not a
## third axis: 61 of its 86 genes are in the ventilator signature with opposite
## sign; it is reported as the ischemic-time-sensitive component of the
## ventilator axis.
min_lfc <- c(ventilator = 0.5, ischemia = 0.25, agony = 0.5)
axes <- c(ventilator = "ventilador", ischemia = "isquemia", agony = "agonia")
read_axis <- function(axis) {
    d <- read.delim(
        file.path(
            repo, "results", "firmas",
            paste0(axes[[axis]], "_consenso.tsv")
        ),
        stringsAsFactors = FALSE
    )
    d <- d[d$n_deg_concordantes >= min_tissues &
        abs(d$lfc_mediano) > min_lfc[[axis]] &
        !grepl("_PAR_Y$", d$ensembl), ]
    data.frame(
        ensembl = sub("\\..*$", "", d$ensembl), symbol = d$Symbol, axis = axis,
        lfc_consensus = d$lfc_mediano, z_stouffer = d$z_stouffer,
        n_tissues = d$n_deg_concordantes,
        direction = ifelse(d$lfc_mediano > 0, "up", "down"),
        stringsAsFactors = FALSE
    )
}
signatures <- do.call(rbind, lapply(names(axes), read_axis))
signatures$t_median <- NA_real_
## Shared ischemia: protein-coding genes whose moderated t for ischemic time
## (limma on log2-CPM, adjusted for Hardy, age and sex) exceeds 2 in absolute
## value with the same sign in at least 60 percent of the non-brain GTEx
## tissues (21 of 34). Defined by consistency of sign, not by fold change:
## the effect per gene is small and a fold-change threshold removes it.
shared <- read.delim(
    file.path(repo, "results", "firmas", "isquemia_atlas_consenso.tsv"),
    stringsAsFactors = FALSE
)
shared <- shared[shared$signo != 0, ]
shared_summary <- read.delim(
    file.path(
        repo, "results", "final", "tables", "isquemia_atlas_resumen.tsv"
    )
)
signatures <- rbind(signatures, data.frame(
    ensembl = shared$ensembl, symbol = shared$symbol,
    axis = "shared_ischemia", lfc_consensus = NA_real_, z_stouffer = NA_real_,
    n_tissues = pmax(shared$n_sube, shared$n_baja),
    direction = ifelse(shared$signo > 0, "up", "down"),
    t_median = shared$t_mediano, stringsAsFactors = FALSE
))
rownames(signatures) <- NULL
vent <- signatures$ensembl[signatures$axis == "ventilator"]
signatures$in_ventilator <- signatures$axis %in%
    c("ischemia", "shared_ischemia") & signatures$ensembl %in% vent
## Technical flag sets (neighbouring-tissue contamination; Nieuwenhuis et al.
## 2020)
genes <- read.delim(
    file.path(repo, "data", "processed", "ovary", "genes.tsv"),
    stringsAsFactors = FALSE
)
genes$ensembl <- sub("\\..*$", "", genes$ensembl)
flag_sets <- list(
    pancreas = c("PNLIP", "CELA3A", "PRSS1", "CPA1", "CELA2A", "CTRB1"),
    surfactant = c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2")
)
flags <- do.call(rbind, lapply(names(flag_sets), function(f) {
    s <- flag_sets[[f]]
    e <- genes$ensembl[match(s, genes$Symbol)]
    data.frame(
        ensembl = e, symbol = s, flag = f,
        stringsAsFactors = FALSE
    )[!is.na(e), ]
}))
rownames(flags) <- NULL
thanatome_signatures <- list(
    signatures = signatures, flags = flags,
    thresholds = list(
        min_tissues = min_tissues, min_lfc = min_lfc,
        shared_ischemia = c(
            min_abs_t = 2, min_tissues = shared_summary$k,
            tissues = shared_summary$tejidos_no_cerebro
        )
    ),
    source = paste0(
        "GTEx v8 (portal gene read counts), ten tissues; DESeq2 per tissue ",
        "and Stouffer consensus; shared ischemia: limma per tissue over the ",
        "non-brain tissues and consensus by sign; ",
        "hachepunto/thanatome-manuscript results/firmas"
    )
)
print(table(signatures$axis, signatures$direction))
print(table(signatures$axis, signatures$in_ventilator))
save(thanatome_signatures,
    file = "data/thanatome_signatures.rda", compress = "xz"
)
record_session("build_signatures", repo)
