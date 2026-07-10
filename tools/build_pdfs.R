#!/usr/bin/env Rscript
# Regenerate every rendered document for sliR.
#
#   Rscript tools/build_pdfs.R                 # write PDFs next to the repo, in ./Derived/
#   Rscript tools/build_pdfs.R /some/other/dir # write them somewhere else
#
# Produces:
#   README.md                    (regenerated in place from README.Rmd)
#   sliR_reference_manual.pdf    every help page, alphabetical -- the R convention
#   sliR_overview.pdf            README.Rmd typeset, with executed code
#   sliR_design_notes.pdf        Project_notes.md typeset
#
# There is no single source file for the reference manual. It is generated from
# the roxygen `#'` comments in R/*.R, via man/*.Rd. Edit the roxygen, run
# devtools::document(), then run this script. Never edit man/*.Rd by hand.

# Setup ----

repo <- normalizePath(getwd())
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from the sliR repository root: Rscript tools/build_pdfs.R", call. = FALSE)
}

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else file.path(repo, "Derived")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## MacTeX lives outside the PATH that Rscript inherits from a GUI session.
tex <- "/Library/TeX/texbin"
if (dir.exists(tex) && !grepl(tex, Sys.getenv("PATH"), fixed = TRUE)) {
  Sys.setenv(PATH = paste(tex, Sys.getenv("PATH"), sep = ":"))
}

scratch <- tempfile("sliR-docs-")
dir.create(scratch)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE)


# 1. README.md, regenerated from README.Rmd ----

message("Rendering README.md ...")
rmarkdown::render(file.path(repo, "README.Rmd"), quiet = TRUE)
unlink(file.path(repo, "README.html"))


# 2. Reference manual: every help page, alphabetical ----

message("Building sliR_reference_manual.pdf ...")
manual <- file.path(out_dir, "sliR_reference_manual.pdf")
status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "Rd2pdf", shQuote(repo),
    paste0("--output=", shQuote(manual)), "--force", "--no-preview"),
  stdout = FALSE, stderr = FALSE
)
if (status != 0) stop("Rd2pdf failed. Run devtools::document() first.")


# 3. Overview: README.Rmd typeset as a standalone PDF ----

## Rebuilt from README.Rmd rather than kept as a separate document, so the two cannot drift apart.
message("Building sliR_overview.pdf ...")
src  <- readLines(file.path(repo, "README.Rmd"), warn = FALSE)
ends <- which(trimws(src) == "---")
body <- src[(ends[2] + 1):length(src)]
body <- body[!grepl("^# sliR\\s*$", body)]        # duplicates the title page
body <- body[!grepl('^<div id="refs">', body)]    # citeproc places these itself

## Drop the badge block entirely. The badges are remote SVGs, which LaTeX cannot fetch and would not want on a title page anyway.
b_start <- grep("<!-- badges: start -->", body, fixed = TRUE)
b_end   <- grep("<!-- badges: end -->",   body, fixed = TRUE)
if (length(b_start) && length(b_end)) body <- body[-(b_start[1]:b_end[1])]

overview <- file.path(scratch, "sliR_overview.Rmd")
writeLines(c(
  "---",
  'title: "sliR"',
  'subtitle: "Standardized Length Index and Allometric Data Simulation"',
  'author: "Aaron Skinner"',
  "date: \"`r format(Sys.Date(), '%d %B %Y')`\"",
  paste0('bibliography: "', file.path(repo, "references.bib"), '"'),
  "link-citations: true",
  "output:",
  "  pdf_document:",
  "    toc: true",
  "    number_sections: true",
  "    highlight: tango",
  "    latex_engine: xelatex",
  "geometry: margin=1in",
  "urlcolor: blue",
  "linkcolor: blue",
  "---",
  "",
  "```{r, include = FALSE}",
  'knitr::opts_chunk$set(collapse = TRUE, comment = "#>")',
  "options(width = 76)",
  "```",
  body
), overview)
rmarkdown::render(overview, quiet = TRUE)
invisible(file.copy(file.path(scratch, "sliR_overview.pdf"), file.path(out_dir, "sliR_overview.pdf"), overwrite = TRUE))


# 4. Design notes: Project_notes.md typeset ----

## xelatex, not pdflatex: the notes contain Greek and mathematical symbols that latin1 inputenc cannot encode.
message("Building sliR_design_notes.pdf ...")
notes_rmd <- file.path(scratch, "sliR_design_notes.Rmd")
writeLines(c(
  "---",
  'title: "sliR: design notes and decisions"',
  'author: "Aaron Skinner"',
  "date: \"`r format(Sys.Date(), '%d %B %Y')`\"",
  "output:",
  "  pdf_document:",
  "    toc: true",
  "    number_sections: false",
  "    latex_engine: xelatex",
  "geometry: margin=1in",
  "urlcolor: blue",
  "linkcolor: blue",
  "---",
  "",
  '```{r, echo = FALSE, results = "asis"}',
  paste0('notes <- readLines("', file.path(repo, "Project_notes.md"), '", warn = FALSE)'),
  'cat(notes[!grepl("^# sliR", notes)], sep = "\\n")',
  "```"
), notes_rmd)
rmarkdown::render(notes_rmd, quiet = TRUE)
invisible(file.copy(file.path(scratch, "sliR_design_notes.pdf"), file.path(out_dir, "sliR_design_notes.pdf"), overwrite = TRUE))

message("\nWrote to ", out_dir, ":")
message(paste(" -", basename(list.files(out_dir, pattern = "\\.pdf$", full.names = TRUE)), collapse = "\n"))
