# EditR

Predicts where base editing occurred from a single Sanger sequencing run. You provide a guide RNA protospacer sequence (~20 bp) and a `.ab1` Sanger file of the edited region (~300–700 bp). EditR fits zero-adjusted gamma distributions to the background noise of your sample and assigns a *P*-value to each background peak in the guide region — the probability that it is a product of base editing as opposed to noise.

## Dependencies

### R packages

Required at runtime (`Imports` in `DESCRIPTION`):

- Biostrings, pwalign, sangerseqR
- gamlss, gamlss.dist
- magrittr, dplyr, tidyr
- ggplot2, cowplot, gridExtra
- rmarkdown (HTML reports), plotly, yaml

Used for tests: testthat. Used by the CLI: pandoc (via `rmarkdown`/HTML reports).

### Via Nix (recommended)

A Nix flake (`flake.nix` + `flake.lock`) installs R, all R packages, and `pandoc` in project scope, pinned to a reproducible nixpkgs revision.

```bash
nix develop          # shell with R + editR + all deps
nix build .#default  # build the editR R package artifact
```

### Manual install

```r
install.packages(c("gamlss", "magrittr", "dplyr", "tidyr", "ggplot2",
                   "cowplot", "gridExtra", "rmarkdown", "plotly", "yaml"))
BiocManager::install(c("Biostrings", "pwalign", "sangerseqR"))
install.packages("path/to/editr", repos = NULL, type = "source")
```

## Usage

### Command line

```bash
# Example analysis
nix run .# -- --example --out results

# Analyze your own sample
nix run .# -- --ab1 sample.ab1 --guide CACTGGAATGACACACGCCC --out results

# Batch: analyze every .ab1 in a directory
nix run .# -- --folder data/ --guide CACTGGAATGACACACGCCC --out results

# Optional parameters
nix run .# -- --ab1 sample.ab1 --guide CACTGGAATGACACACGCCC \
  --pval 0.01 --trim5 20 --trim3 318 --rc --out results
```

Flags:

| Flag | Description |
| ---- | ----------- |
| `--ab1 <file>` | Path to the Sanger `.ab1` file |
| `--folder <dir>` | Analyze all `.ab1` files in a directory (batch mode) |
| `--guide <seq>` | gRNA protospacer sequence, 5' to 3' |
| `--pval <num>` | P-value cutoff (default `0.01`) |
| `--trim5 <n>` | 5' trim position (optional) |
| `--trim3 <n>` | 3' trim position (optional) |
| `--rc` | Guide sequence is the reverse complement |
| `--example` | Use the bundled example data |
| `--out <dir>` | Output directory for plots (default `results`) |

For each `.ab1` file, three PNG plots are written to the output directory: `<name>_quad.png`, `<name>_tile.png`, and `<name>_chromatogram.png`. In batch mode the run continues past individual failures and exits non-zero if any file failed.

### R API

```r
library(editR)

result <- run_editr(
  ab1_file = "path/to/sample.ab1",
  guide = "CACTGGAATGACACACGCCC",   # ~20 bp gRNA, 5' to 3'
  pval_cutoff = 0.01,
  trim5 = NA,                        # optional 5' trim position
  trim3 = NA,                        # optional 3' trim position
  guide_is_rc = FALSE,               # TRUE if the guide is antisense
  output_dir = "results"             # optional: also write an HTML report here
)

# Use the bundled example data
res <- run_editr(example = TRUE, output_dir = "results")
```

`run_editr()` returns a list with:

- `editing.df` — guide region with per-base peak areas and *P*-values
- `base.info` — average percent signal, critical percent value, model mu, and Filliben's correlation per base
- `null.m.params` — fitted zero-adjusted gamma parameters per base
- `sangs` / `sangs.filt` — peak area data before and after filtering
- `guide.coord` — start/end indices of the guide match
- `input.seq` / `input.basecalls` — the raw and basecalled Sanger objects

Bases with a *P*-value below the cutoff are called as significantly different from noise — i.e. editing has occurred.

## Example data

Bundled `.ab1` files for testing:

```r
system.file("example.ab1", package = "editR")
system.file("testfiles", package = "editR")
```

## License

GPL-3
