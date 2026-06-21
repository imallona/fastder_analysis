#!/usr/bin/env Rscript
# Landscape coverage track for one cryptic-exon locus (STMN2 by default) as a
# tunable SVG and PDF, so Main Figure 2 panel B is vector rather than a raster
# PNG. Track logic mirrors the coverage-plots chunk of recount3.Rmd. Inputs are
# the recount3 manifest (group, sample, bigwig, per-tool GTF columns) and the
# Ensembl reference GTF; the BigWigs the manifest points to must be reachable.
# Usage:
#   Rscript make_stmn2_track.R \
#     --manifest results/config_klim_2019_tdp43_recount3/recount3_manifest.csv \
#     --reference-gtf .../Homo_sapiens.GRCh38.115.chr.gtf \
#     --out fig_tdp43_stmn2.svg --gene STMN2 --width 12 --height 5

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(Gviz)
})

parse_args <- function(args) {
  defaults <- list(
    manifest = "",
    `reference-gtf` = "",
    out = "fig_tdp43_stmn2.svg",
    gene = "STMN2",
    tools = "fastder,derfinder,megadepth_baseline,grohmm",
    width = "4.5",
    height = "6.5",
    pdf = "true")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% names(defaults)) stop("unknown argument: ", args[[i]])
    defaults[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (!nzchar(opt$manifest) || !file.exists(opt$manifest))
  stop("missing or unreadable --manifest: ", opt$manifest)

width <- as.numeric(opt$width)
height <- as.numeric(opt$height)
want_tools <- strsplit(opt$tools, ",")[[1]]

# Cryptic exon panel for the TDP-43 knockdown, hg38. Same table as recount3.Rmd:
# the plot window is the gene body plus a small flank; ce_start/ce_end mark the
# cryptic exon read from the knockdown-specific novel splice junctions.
regions_all <- data.frame(
  gene = c("STMN2", "HDGFL2", "ELAVL3", "CELF5", "KCNQ2"),
  chrom = c("chr8", "chr19", "chr19", "chr19", "chr20"),
  start = c(79606000, 4467000, 11446000, 3219000, 63395000),
  end = c(79671000, 4507000, 11486000, 3302000, 63478000),
  ce_start = c(79616821, 4492015, 11463496, 3278200, 63444558),
  ce_end = c(79617200, 4492152, 11463662, 3278400, 63444758),
  stringsAsFactors = FALSE)
r <- regions_all[regions_all$gene == opt$gene, ]
if (nrow(r) != 1) stop("gene not in panel: ", opt$gene)

manifest <- read.csv(opt$manifest, stringsAsFactors = FALSE)
manifest <- manifest[order(manifest$group, manifest$sample), ]
groups <- unique(manifest$group)
# WT blue, KD red, matching panels A and C-E. Short titles keep the track compact;
# the legend and scheme spell out TDP-43 WT / TDP-43 KD.
group_colors <- c(control = "#4575b4", knockdown = "#d73027")[groups]
group_short <- c(control = "WT", knockdown = "KD")
manifest$idx <- ave(seq_len(nrow(manifest)), manifest$group, FUN = seq_along)

# Tool columns, palette and labels, identical to the report.
tools <- c(fastder = "fastder_gtf",
           derfinder = "derfinder_gtf",
           megadepth_baseline = "megadepth_gtf",
           grohmm = "grohmm_gtf")
tools <- tools[names(tools) %in% want_tools]
tool_colors <- c(derfinder = "#66C2A5", fastder = "#FC8D62",
                 grohmm = "#8DA0CB", megadepth_baseline = "#E78AC3")
tool_labels <- c(fastder = "fastder", derfinder = "derfinder",
                 megadepth_baseline = "megadepth", grohmm = "grohmm")

# Ensembl models: the GTF drops the "chr" prefix the BigWigs and tool GTFs use.
gene_models <- NULL
if (nzchar(opt$`reference-gtf`) && file.exists(opt$`reference-gtf`)) {
  gene_models <- tryCatch({
    ann <- import(opt$`reference-gtf`)
    seqlevels(ann) <- paste0("chr", sub("^chr", "", seqlevels(ann)))
    ann
  }, error = function(e) NULL)
}

# Keep exon features: fastder emits one transcript per stitched ER spanning its
# introns, so plotting transcripts hides whether the cryptic exon was its own ER.
load_regions <- function(gtf) {
  if (length(gtf) == 0 || !nzchar(gtf) || !file.exists(gtf)) return(GRanges())
  gr <- tryCatch(import(gtf), error = function(e) GRanges())
  if (length(gr) > 0 && "type" %in% colnames(mcols(gr))) {
    ex <- gr[!is.na(gr$type) & gr$type == "exon"]
    if (length(ex) > 0) gr <- ex
  }
  gr
}
regions_by <- list()
for (tool in names(tools)) {
  regions_by[[tool]] <- list()
  for (group in groups) {
    gtf <- unique(manifest[[tools[[tool]]]][manifest$group == group])
    regions_by[[tool]][[group]] <- load_regions(if (length(gtf)) gtf[1] else "")
  }
}

# Per-sample CPM scaling: fastder's library size is the area under the coverage
# curve, so 1e6 / sum(width * coverage) puts the Y axis in the same CPM units as
# the --min-coverage threshold.
cpm_scale <- function(bigwig) {
  if (!nzchar(bigwig) || !file.exists(bigwig)) return(1.0)
  gr <- tryCatch(rtracklayer::import(bigwig, as = "GRanges"),
                 error = function(e) NULL)
  if (is.null(gr) || length(gr) == 0) return(1.0)
  lib <- sum(as.numeric(width(gr)) * as.numeric(score(gr)))
  if (!is.finite(lib) || lib <= 0) 1.0 else 1e6 / lib
}
manifest$cpm_factor <- vapply(manifest$bigwig, cpm_scale, numeric(1))

# Coverage threshold fastder used, parsed from the GTF parameter id (mc<value>).
extract_min_cov <- function() {
  for (group in groups) {
    gtf <- unique(manifest$fastder_gtf[manifest$group == group])
    if (length(gtf) && nzchar(gtf[1]) && file.exists(gtf[1])) {
      m <- regmatches(gtf[1], regexpr("mc[0-9.]+", gtf[1]))
      if (length(m) && nchar(m)) return(as.numeric(sub("mc", "", m)))
    }
  }
  0.05
}
min_cov_cpm <- extract_min_cov()

build_tracks <- function() {
  # Labels below the scale line so they are not clipped at the top of the canvas.
  axis_track <- GenomeAxisTrack(labelPos = "below", fontsize = 9)
  body <- list()
  body_sizes <- numeric(0)

  if (!is.null(gene_models)) {
    in_window <- gene_models[seqnames(gene_models) == r$chrom &
                             start(gene_models) <= r$end &
                             end(gene_models) >= r$start &
                             gene_models$type == "exon"]
    if (length(in_window) > 0) {
      # Merge all transcripts into one clear exon model rather than stacking many
      # faint transcript rows: the marker track only needs the gene's exon layout.
      merged <- reduce(granges(in_window), ignore.strand = TRUE)
      body <- c(body, AnnotationTrack(
        merged, name = "Ensembl", stacking = "dense",
        fill = "grey45", col = "grey45"))
      body_sizes <- c(body_sizes, 1.0)
    }
  }

  for (j in seq_len(nrow(manifest))) {
    m <- manifest[j, ]
    f <- m$cpm_factor
    body <- c(body, DataTrack(
      range = m$bigwig, genome = "hg38", chromosome = r$chrom,
      name = paste0(group_short[m$group], " ", m$idx),
      type = c("histogram", "h"),
      transformation = local({ factor <- f; function(x) x * factor }),
      baseline = min_cov_cpm,
      col.baseline = "#d73027", lty.baseline = 2, lwd.baseline = 0.7,
      col = group_colors[m$group],
      col.histogram = group_colors[m$group],
      fill.histogram = group_colors[m$group]))
    body_sizes <- c(body_sizes, 1.1)
  }

  for (tool in names(tools)) {
    for (group in groups) {
      gr <- regions_by[[tool]][[group]]
      if (length(gr) > 0) {
        gr <- gr[seqnames(gr) == r$chrom & start(gr) <= r$end & end(gr) >= r$start]
      }
      body <- c(body, AnnotationTrack(
        gr, name = paste0(tool_labels[tool], " ", group_short[group]),
        chromosome = r$chrom, fill = tool_colors[tool], col = NA,
        stacking = "squish", shape = "box"))
      body_sizes <- c(body_sizes, 0.7)
    }
  }

  ce_w <- r$ce_end - r$ce_start
  highlight <- HighlightTrack(
    trackList = body, chromosome = r$chrom,
    start = r$ce_start - 1.5 * ce_w, end = r$ce_end + 1.5 * ce_w,
    col = "#d73027", fill = "#fde0dd")
  list(axis = axis_track, highlight = highlight, sizes = body_sizes)
}

draw <- function(tr) {
  plotTracks(list(tr$axis, tr$highlight),
             from = r$start, to = r$end, chromosome = r$chrom,
             background.title = "white", col.title = "black",
             col.axis = "black", fontsize = 10, cex.title = 0.85, cex.axis = 0.85,
             rotation.title = 0, title.width = 1.6, sizes = c(1.0, tr$sizes))
}

tr <- build_tracks()

svglite::svglite(opt$out, width = width, height = height)
draw(tr)
invisible(dev.off())
cat("wrote", opt$out, "\n")

if (identical(opt$pdf, "true")) {
  pdf_out <- sub("\\.svg$", ".pdf", opt$out)
  grDevices::cairo_pdf(pdf_out, width = width, height = height)
  draw(tr)
  invisible(dev.off())
  cat("wrote", pdf_out, "\n")
}
