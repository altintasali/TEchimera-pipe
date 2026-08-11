#!/usr/bin/env Rscript
# Chimera sample-QC: normalize the gene-TE junction counts matrix and produce
# the PCA + sample-distance views shipped by sample_qc.smk.
#
# Two modes, selected by the script's argument vector:
#   --transform counts.tsv samples.csv transform min_samples_present \
#                 min_total_counts out_matrix.tsv
#       Apply the chosen transformation (vst / rlog / log2) to the counts
#       matrix (event_id x sample, as written by chimera_counts.py) and write
#       the transformed matrix for the plot rule to read. vst/rlog use
#       DESeq2's blind normalization (independent of sample labels, so the QC
#       view can't be overfit); log2 is log2(x + 1) without DESeq2. Events
#       seen in fewer than min_samples_present samples, or with fewer than
#       min_total_counts supporting reads overall, are dropped for this view
#       only.
#   --plots transformed.tsv samples.csv min_events transform \
#                 out_pca_mqc.json out_heatmap_mqc.json
#       PCA scatter + sample-to-sample Euclidean-distance heatmap of the
#       transformed matrix, written as MultiQC custom-content JSON documents
#       (the `_mqc.json` suffix makes multiqc_report.html render them
#       interactively). PCA points are colored by the sample sheet's
#       "condition" column; samples without a condition form one group.
#
# QC filters (chimera.qc in the config) apply ONLY to this view: the all_events
# catalog and the counts matrix are never reduced.
# Libraries are loaded per-mode (DESeq2 for the transform); the plots mode
# needs no extra packages -- the JSON payloads are small and built by hand.
suppressMessages(library(DESeq2))

read_counts <- function(path) {
    m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
    storage.mode(m) <- "integer"
    m
}

transform_matrix <- function(counts, method) {
    if (method == "log2") {
        return(log2(counts + 1))
    }
    if (ncol(counts) < 2) {
        stop("need >= 2 samples for vst/rlog")
    }
    if (method == "vst") {
        tryCatch(
            vst(counts, blind = TRUE),
            error = function(e) transform_fallback(counts, method, e)
        )
    } else if (method == "rlog") {
        tryCatch(
            rlog(counts, blind = TRUE),
            error = function(e) transform_fallback(counts, method, e)
        )
    } else {
        stop(paste("unknown transform:", method))
    }
}

transform_fallback <- function(counts, method, orig_error) {
    # vst()/rlog() fail on small matrices (fewer rows than DESeq2's nsub);
    # the direct *Transformation() calls handle those but can still choke when
    # every event's gene-wise dispersion sits at the minimum (curve fitting
    # gives up), which happens easily with a handful of chimeric junctions.
    direct <- if (method == "vst") {
        tryCatch(
            varianceStabilizingTransformation(counts, blind = TRUE),
            error = function(e) NULL
        )
    } else {
        tryCatch(
            rlogTransformation(counts, blind = TRUE),
            error = function(e) NULL
        )
    }
    if (!is.null(direct)) {
        return(as.matrix(direct))
    }
    message(sprintf(
        "%s unavailable for %d events (both vst/rlog and the direct varianceStabilizingTransformation/rlogTransformation failed: %s); falling back to log2(counts + 1) for the QC view",
        method, nrow(counts), conditionMessage(orig_error)
    ))
    log2(counts + 1)
}

load_conditions <- function(samples_path) {
    sheet <- read.csv(samples_path, comment.char = "#", check.names = FALSE)
    if ("condition" %in% colnames(sheet)) {
        cond <- sheet$condition
    } else {
        cond <- rep("all", nrow(sheet))
    }
    cond[is.na(cond) | cond == ""] <- "all"
    setNames(as.character(cond), as.character(sheet$sample))
}

# --- MultiQC custom-content JSON writers -------------------------------------
# The plots mode emits two documents MultiQC picks up by the `_mqc.json`
# suffix (results/chimera/qc/*_mqc.json) and renders as an interactive
# scatter (PCA) and heatmap (sample distances). Hand-built JSON: the payloads
# are small and fixed, and the QC env needs no extra JSON dependency.

json_escape <- function(x) gsub('([\\\\"])', '\\\\\\1', x)
json_num <- function(x) sprintf("%.6g", as.numeric(x))
json_str_arr <- function(x) paste(sprintf('"%s"', json_escape(x)), collapse = ", ")

write_pca_mqc <- function(path, samples, x, y, colors, pc1, pc2, transform, note = NULL) {
    pts <- vapply(seq_along(samples), function(i) {
        sprintf('"%s": {"x": %s, "y": %s, "color": "%s"}',
                json_escape(samples[i]), json_num(x[i]), json_num(y[i]), colors[i])
    }, character(1))
    desc <- if (is.null(note)) {
        sprintf(paste0("Principal-component analysis of the %s-transformed ",
                       "chimera counts matrix, colored by sample condition."),
                transform)
    } else {
        note
    }
    body <- paste0(
        '{\n',
        '  "id": "chimera_sample_qc_pca",\n',
        '  "section_name": "Chimera sample-QC: PCA",\n',
        sprintf('  "description": "%s",\n', json_escape(desc)),
        '  "plot_type": "scatter",\n',
        '  "pconfig": {\n',
        '    "id": "chimera_pca_plot",\n',
        '    "title": "Chimera counts: PCA",\n',
        sprintf('    "xlab": "PC1 (%s%%)",\n', json_num(pc1)),
        sprintf('    "ylab": "PC2 (%s%%)"\n', json_num(pc2)),
        '  },\n',
        sprintf('  "data": {%s}\n', paste(pts, collapse = ", ")),
        '}\n'
    )
    writeLines(body, path)
}

write_heatmap_mqc <- function(path, samples, d, transform, note = NULL) {
    desc <- if (is.null(note)) {
        sprintf(paste0("Pairwise Euclidean distances between samples on the ",
                       "%s-transformed chimera counts matrix."),
                transform)
    } else {
        note
    }
    rows <- apply(d, 1, function(r) paste0('[', paste(json_num(r), collapse = ", "), ']'))
    body <- paste0(
        '{\n',
        '  "id": "chimera_sample_qc_heatmap",\n',
        '  "section_name": "Chimera sample-QC: sample distances",\n',
        sprintf('  "description": "%s",\n', json_escape(desc)),
        '  "plot_type": "heatmap",\n',
        '  "pconfig": {\n',
        '    "id": "chimera_heatmap_plot",\n',
        '    "title": "Euclidean distance between samples"\n',
        '  },\n',
        sprintf('  "ycats": [%s],\n', json_str_arr(samples)),
        sprintf('  "xcats": [%s],\n', json_str_arr(samples)),
        sprintf('  "data": [%s]\n', paste(rows, collapse = ",\n    ")),
        '}\n'
    )
    writeLines(body, path)
}

write_empty_mqc <- function(path, kind) {
    # Valid custom-content documents for the no-data case, so the multiqc rule
    # always sees its inputs and the report still documents why nothing is
    # plotted.
    note <- "No chimeric events passed the QC-view filters; plot skipped."
    body <- if (kind == "scatter") {
        paste0(
            '{\n',
            '  "id": "chimera_sample_qc_pca",\n',
            '  "section_name": "Chimera sample-QC: PCA",\n',
            sprintf('  "description": "%s",\n', json_escape(note)),
            '  "plot_type": "scatter",\n',
            '  "pconfig": {"id": "chimera_pca_plot", "title": "Chimera counts: PCA"},\n',
            '  "data": {}\n',
            '}\n'
        )
    } else {
        paste0(
            '{\n',
            '  "id": "chimera_sample_qc_heatmap",\n',
            '  "section_name": "Chimera sample-QC: sample distances",\n',
            sprintf('  "description": "%s",\n', json_escape(note)),
            '  "plot_type": "heatmap",\n',
            '  "pconfig": {"id": "chimera_heatmap_plot", "title": "Euclidean distance between samples"},\n',
            '  "ycats": [],\n',
            '  "xcats": [],\n',
            '  "data": []\n',
            '}\n'
        )
    }
    writeLines(body, path)
}

do_transform <- function(argv) {
    counts_path <- argv[2]
    samples_path <- argv[3]
    method <- argv[4]
    min_samples_present <- as.integer(argv[5])
    min_total_counts <- as.integer(argv[6])
    out_matrix <- argv[7]

    counts <- read_counts(counts_path)
    if (nrow(counts) == 0) {
        file.create(out_matrix)
        message("no chimeric events; wrote empty transformed matrix")
        quit(save = "no", status = 0)
    }
    # QC-view filter: keep events seen in >= min_samples_present samples with
    # >= min_total_counts supporting reads overall. Only affects this view.
    present <- rowSums(counts > 0)
    keep <- present >= min_samples_present &
        rowSums(counts) >= min_total_counts
    counts <- counts[keep, , drop = FALSE]
    if (nrow(counts) == 0) {
        file.create(out_matrix)
        message("no events passed the QC-view filters; empty transformed matrix")
        quit(save = "no", status = 0)
    }
    message(sprintf("transform=%s on %d events x %d samples", method,
                    nrow(counts), ncol(counts)))
    transformed <- as.matrix(transform_matrix(counts, method))
    write.table(transformed, out_matrix, sep = "\t", quote = FALSE)
}

do_plots <- function(argv) {
    matrix_path <- argv[2]
    samples_path <- argv[3]
    min_events <- as.integer(argv[4])
    transform <- argv[5]
    out_pca <- argv[6]
    out_heatmap <- argv[7]

    sz <- file.info(matrix_path)$size
    tmat <- if (!is.na(sz) && sz > 0) {
        as.matrix(read.delim(matrix_path, row.names = 1, check.names = FALSE))
    } else {
        # the transform writes a 0-byte file when nothing passed its QC filter
        matrix(nrow = 0, ncol = 0)
    }
    empty <- function(msg) {
        message(msg)
        write_empty_mqc(out_pca, "scatter")
        write_empty_mqc(out_heatmap, "heatmap")
        quit(save = "no", status = 0)
    }
    if (ncol(tmat) < 2) {
        # a single sample can't support a PCA or a sample-distance view
        empty("only one sample; skipping PCA/heatmap")
    }
    # Events that transform to a constant (or NaN/Inf) column carry no signal
    # and make prcomp(scale.=TRUE) fail outright, so drop them for the QC view.
    finite_rows <- apply(tmat, 1, function(row) all(is.finite(row)))
    var_rows <- apply(tmat, 1, var, na.rm = TRUE) > 0
    dropped <- sum(!(finite_rows & var_rows))
    if (dropped > 0) {
        message(sprintf("dropping %d zero-variance/non-finite event(s) from the QC view", dropped))
    }
    tmat <- tmat[finite_rows & var_rows, , drop = FALSE]
    if (nrow(tmat) < max(min_events, 2)) {
        # too few events for a meaningful PCA/distance view; ship empty custom
        # content rather than a cryptic prcomp error.
        empty(sprintf("fewer than %d events; skipping PCA/heatmap", min_events))
    }

    conditions <- load_conditions(samples_path)
    samples <- colnames(tmat)
    groups <- conditions[samples]
    groups[is.na(groups)] <- "all"

    pca <- prcomp(t(tmat), center = TRUE, scale. = TRUE)
    pc1 <- round(100 * summary(pca)$importance[2, 1], 1)
    pc2 <- round(100 * summary(pca)$importance[2, 2], 1)

    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7")
    ugroups <- unique(groups)
    group_colors <- setNames(palette[seq_along(ugroups)], ugroups)
    point_colors <- unname(group_colors[groups])

    write_pca_mqc(out_pca, samples, pca$x[, 1], pca$x[, 2], point_colors,
                  pc1, pc2, transform)

    write_heatmap_mqc(out_heatmap, samples, as.matrix(dist(t(tmat))), transform)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("usage: chimera_sample_qc.R --transform|--plots ...")
}
mode <- args[1]
if (mode == "--transform") {
    do_transform(args)
} else if (mode == "--plots") {
    do_plots(args)
} else {
    stop(paste("unknown mode:", mode))
}
