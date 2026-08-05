#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(magrittr)
  library(sangerseqR)
  library(editR)
})

args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) == 0) return(default)
  args[idx + 1]
}

has_flag <- function(name) {
  name %in% args
}

ab1_file <- get_opt("--ab1")
folder <- get_opt("--folder")
guide <- get_opt("--guide")
pval <- as.numeric(get_opt("--pval", 0.01))
trim5 <- as.numeric(get_opt("--trim5", NA))
trim3 <- as.numeric(get_opt("--trim3", NA))
guide_is_rc <- has_flag("--rc")
example <- has_flag("--example")
output_dir <- get_opt("--out", "results")

if (example) {
  ab1_file <- system.file("example.ab1", package = "editR")
  guide <- "CACTGGAATGACACACGCCC"
}

if (is.null(guide)) {
  stop("Please provide a guide RNA sequence with --guide")
}

files <- if (!is.null(folder)) {
  list.files(folder, pattern = "\\.ab1$", full.names = TRUE, ignore.case = TRUE)
} else if (!is.null(ab1_file)) {
  ab1_file
} else {
  stop("Provide --ab1 <file.ab1> or --folder <dir>")
}

if (length(files) == 0) {
  stop("No .ab1 files found", if (!is.null(folder)) paste0(" in ", folder))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

make_quad_plot <- function(editing.df, null.m.params, p.val.cutoff) {
  edit.long <- editing.df %>% gather(key = focal.base, value = value,
                                     A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
    separate(col = focal.base, into = c("focal.base", "measure"))

  p.a <- makeEditingBarPlot(edit.long = edit.long, null.m.params = null.m.params$a,
                            base = "A", pval = p.val.cutoff, editing.df)
  p.c <- makeEditingBarPlot(edit.long, null.m.params$c,
                            base = "C", pval = p.val.cutoff, editing.df)
  p.g <- makeEditingBarPlot(edit.long, null.m.params$g,
                            base = "G", pval = p.val.cutoff, editing.df)
  p.t <- makeEditingBarPlot(edit.long, null.m.params$t,
                            base = "T", pval = p.val.cutoff, editing.df)
  grid.arrange(p.a, p.c, p.g, p.t)
}

make_tile_plot <- function(editing.df, sangs.filt, null.m.params, p.val.cutoff) {
  avg.base <- sangs.filt %>% gather(key = focal.base, value = value,
                                    A.area:T.area, A.perc:T.perc) %>%
    separate(col = focal.base, into = c("focal.base", "measure")) %>%
    spread(key = measure, value = value) %>%
    filter(base.call == focal.base) %>%
    group_by(focal.base) %>%
    summarize(avg.percsignal = mean(perc),
              avg.areasignal = mean(area))

  mul <- lapply(null.m.params, FUN = function(x) {x$mu})
  mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)

  edit.long <- editing.df %>% gather(key = focal.base, value = value,
                                     A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
    separate(col = focal.base, into = c("focal.base", "measure"))

  edit.spread <- edit.long %>% spread(key = measure, value = value)

  color.cutoff <- min(avg.base$avg.percsignal - mulvec)
  edit.color <- edit.spread %>%
    mutate(adj.perc = {ifelse(perc >= color.cutoff, 100, perc)} %>% as.numeric) %>%
    filter(pval < p.val.cutoff)

  tile_theme <- theme(axis.ticks = element_blank(),
                      axis.text = element_text(size = 16),
                      plot.title = element_text(hjust = 0, size = 16),
                      plot.margin = unit(c(0, 0, 0, 2), "cm"),
                      panel.background = element_rect(fill = "transparent", colour = NA),
                      plot.background = element_rect(fill = "transparent", colour = NA))

  if (any(edit.color$adj.perc != 100)) {
    p <- edit.spread %>%
      ggplot(aes(x = as.factor(index), y = focal.base)) +
      geom_tile(data = edit.color, aes(fill = adj.perc)) +
      geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +
      guides(fill = "none") +
      scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") +
      scale_x_discrete(position = "top", labels = editing.df$guide.seq) +
      labs(x = NULL, y = NULL) +
      tile_theme +
      coord_fixed(1)
  } else {
    p <- edit.spread %>%
      ggplot(aes(x = as.factor(index), y = focal.base)) +
      geom_tile(data = edit.color, fill = "#9acdee") +
      geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +
      guides(fill = "none") +
      scale_x_discrete(position = "top", labels = editing.df$guide.seq) +
      labs(x = NULL, y = NULL) +
      tile_theme +
      coord_fixed(1)
  }
  p
}

make_chromatogram <- function(res, output_dir, base_name) {
  png(file.path(output_dir, paste0(base_name, "_chromatogram.png")),
      width = 1500, height = 600, res = 150)
  chromatogram(obj = res$input.basecalls,
               showcalls = "none",
               showhets = FALSE,
               trim5 = res$guide.coord$start - 1,
               trim3 = length(res$input.basecalls@primarySeq) - res$guide.coord$end,
               width = res$guide.coord$end - res$guide.coord$start + 1)
  dev.off()
}

analyze_one <- function(ab1_file, guide, pval, trim5, trim3, guide_is_rc, output_dir) {
  message("Analyzing ", ab1_file)
  res <- editR::run_editr(
    ab1_file = ab1_file,
    guide = guide,
    pval_cutoff = pval,
    trim5 = trim5,
    trim3 = trim3,
    guide_is_rc = guide_is_rc
  )

  base_name <- tools::file_path_sans_ext(basename(ab1_file))

  quad <- make_quad_plot(res$editing.df, res$null.m.params, pval)
  ggsave(file.path(output_dir, paste0(base_name, "_quad.png")), quad,
         width = 10, height = 8, dpi = 150)

  tile <- make_tile_plot(res$editing.df, res$sangs.filt, res$null.m.params, pval)
  ggsave(file.path(output_dir, paste0(base_name, "_tile.png")), tile,
         width = 10, height = 4, dpi = 150)

  make_chromatogram(res, output_dir, base_name)

  message("  saved: ", base_name, "_quad.png, ", base_name, "_tile.png, ",
          base_name, "_chromatogram.png")
}

failures <- character(0)
for (f in files) {
  ok <- tryCatch({
    analyze_one(f, guide, pval, trim5, trim3, guide_is_rc, output_dir)
    TRUE
  }, error = function(e) {
    message("  ERROR on ", f, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) failures <- c(failures, f)
}

cat("\nDone. Plots saved to ", output_dir, "\n", sep = "")
if (length(failures) > 0) {
  cat("Failed files:\n", paste0("  ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
