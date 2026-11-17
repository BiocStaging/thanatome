## Assisted by Claude Code (Anthropic); see README.md, section Provenance.

test_that("death_gene_sets splits each axis by direction", {
    gs <- death_gene_sets()
    expect_type(gs, "list")
    expect_named(gs, c("ventilator_up", "ventilator_down",
                        "agony_up", "agony_down"))
    sg <- death_signatures(c("ventilator", "agony"))
    up <- sg$ensembl[sg$axis == "ventilator" & sg$direction == "up"]
    expect_setequal(gs$ventilator_up, up)
    expect_length(death_gene_sets(split = FALSE), 2)
    sym <- death_gene_sets("ischemia", id = "symbol")
    expect_false(any(grepl("^ENSG", unlist(sym))))
})

test_that("death_gene_sets exports to GSEABase", {
    skip_if_not_installed("GSEABase")
    gsc <- death_gene_sets(format = "GeneSetCollection")
    expect_s4_class(gsc, "GeneSetCollection")
    expect_identical(names(gsc), names(death_gene_sets()))
})

test_that("shared ischemia is a signed signature without fold change", {
    s <- death_signatures("shared_ischemia")
    expect_true(all(s$direction %in% c("up", "down")))
    expect_true(all(is.na(s$lfc_consensus)))
    expect_true(all(sign(s$t_median) == ifelse(s$direction == "up", 1, -1)))
    gs <- death_gene_sets("shared_ischemia")
    expect_named(gs, c("shared_ischemia_up", "shared_ischemia_down"))
    expect_error(annotate_de(data.frame(ensembl = "ENSG00000000003",
                                        log2FoldChange = 1),
                            axes = "shared_ischemia"))
})

test_that("shared ischemia score follows ischemic time in the example", {
    data(thanatome_example)
    ctl <- thanatome_example[, thanatome_example$group == "control"]
    sc <- thanatome_score(ctl, "ovary")
    expect_gt(cor(sc$shared_ischemia, ctl$ischemic_time,
                    method = "spearman"), 0.5)
})
