#!/usr/bin/env Rscript
# Chimera sample-QC: normalize the gene-TE junction counts matrix and produce
# the PCA + clustering views shipped by sample_qc.smk.
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
#   --plots transformed.tsv samples.csv min_events out_pca.svg out_heatmap.svg
#       PCA + pheatmap of the transformed matrix, colored by the sample sheet's
#       "condition" column (all samples colored by their own group; samples
#       without a condition are treated as one group).
#
# QC filters (chimera.qc in the config) apply ONLY to this view: the all_events
# catalog and the counts matrix are never reduced.
# Libraries are loaded per-mode (DESeq2 for the transform, pheatmap for the
# plots) so a transform-only run doesn't need pheatmap and vice versa.
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
        # vst() on a matrix with fewer rows than DESeq2's internal nsub
        # errors out; the direct varianceStabilizingTransformation() handles
        # small event sets (common for gene-TE chimeras) fine.
        tryCatch(
            vst(counts, blind = TRUE),
            error = function(e)
                varianceStabilizingTransformation(counts, blind = TRUE)
        )
    } else if (method == "rlog") {
        tryCatch(
            rlog(counts, blind = TRUE),
            error = function(e) rlogTransformation(counts, blind = TRUE)
        )
    } else {
        stop(paste("unknown transform:", method))
    }
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

write_svg <- function(path, width, height, expr) {
    svg(path, width = width, height = height)
    expr
    dev.off()
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
    out_pca <- argv[5]
    out_heatmap <- argv[6]

    sz <- file.info(matrix_path)$size
    tmat <- if (!is.na(sz) && sz > 0) {
        as.matrix(read.delim(matrix_path, row.names = 1, check.names = FALSE))
    } else {
        # the transform writes a 0-byte file when nothing passed its QC filter
        matrix(nrow = 0, ncol = 0)
    }
    placeholder <- function(msg) {
        message(msg)
        for (p in c(out_pca, out_heatmap)) {
            svg(p, width = 6, height = 6)
            plot.new()
            text(0.5, 0.5, sprintf("< %d events, QC view skipped", min_events))
            dev.off()
        }
        quit(save = "no", status = 0)
    }
    if (ncol(tmat) < 2) {
        # a single sample can't support a PCA or a clustering view
        placeholder("only one sample; skipping PCA/heatmap")
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
        # too few events for a meaningful PCA/clustering view; ship empty
        # outputs rather than a cryptic DESeq2/pheatmap error.
        placeholder(sprintf("fewer than %d events; skipping PCA/heatmap", min_events))
    }

    conditions <- load_conditions(samples_path)
    samples <- colnames(tmat)
    groups <- conditions[samples]
    groups[is.na(groups)] <- "all"

    pca <- prcomp(t(tmat), center = TRUE, scale. = TRUE)
    pc1 <- round(100 * summary(pca)$importance[2, 1], 1)
    pc2 <- round(100 * summary(pca)$importance[2, 2], 1)
    scores <- as.data.frame(pca$x)
    scores$sample <- rownames(scores)
    scores$condition <- groups

    write_svg(out_pca, 7, 6, {
        palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                     "#D55E00", "#CC79A7")
        cols <- setNames(palette[seq_along(unique(groups))], unique(groups))
        plot(scores[, 1], scores[, 2], col = cols[groups], pch = 19, cex = 1.4,
             xlab = sprintf("PC1 (%.1f%%)", pc1),
             ylab = sprintf("PC2 (%.1f%%)", pc2),
             main = "Chimera sample-QC: PCA")
        text(scores[, 1], scores[, 2], labels = scores$sample, pos = 3, cex = 0.7)
        legend("topright", legend = names(cols), col = cols, pch = 19, cex = 0.8)
    })

    annot <- data.frame(row.names = samples, condition = groups)
    # Loaded lazily: the placeholder path above never needs pheatmap.
    suppressMessages(library(pheatmap))
    write_svg(out_heatmap, 8, 8, {
        pheatmap(tmat, annotation_col = annot, show_colnames = TRUE,
                 show_rownames = FALSE, main = "Chimera sample-QC: clustering")
    })
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
