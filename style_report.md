# Plain-language styling pass: R Markdown reports

Scope: markdown prose only. YAML front matter, code chunks, inline `r ...`
expressions, kable captions inside code, and all technical terms were left
untouched. No facts, numbers, or dataset names were changed.

## Files edited

- workflow/reports/benchmarks.Rmd
- workflow/reports/gtex_concordance.Rmd
- workflow/reports/recount3.Rmd
- workflow/reports/summary.Rmd

## Files left unchanged

- workflow/reports/meta.Rmd: prose already short and plain; nothing to fix.
- workflow/reports/summary_custom.Rmd: prose already plain and concise.

## Main kinds of change

- Split long sentences (over about 25 words) into one idea per sentence.
- Replaced semicolon sentence joins with full stops where they marked two
  separate ideas.
- Minor word swaps: "a number of" to "many", "apples-to-apples" to
  "like-for-like".
- Added a missing verb ("max_rss from ..." to "max_rss comes from ...") and
  used active phrasing where the agent was clear.

No em-dashes were present or introduced. No AI-writing tics (empty openers,
vague intensifiers) were found in the prose.

## Before/after examples

benchmarks.Rmd:
- Before: "`max_rss` from snakemake's benchmark TSVs. It is the per-process
  peak resident set size in MiB sampled by snakemake every 10 s;
  multi-process rules report the parent only, so high-fan-out runs ..."
- After: "`max_rss` comes from snakemake's benchmark TSVs. It is the
  per-process peak resident set size in MiB, sampled by snakemake every 10 s.
  Multi-process rules report the parent only, so high-fan-out runs ..."

summary.Rmd:
- Before: "... so an exon that the variant skips still receives reads from
  the template and shows up as a coverage dip rather than a hard zero. In
  *variant only* the FASTQ is post-filtered ... and the truth GFF is
  rewritten to remove templates, so the alternative-splicing signal is
  unambiguous: ..."
- After: "... both the canonical transcript and the alternative form are
  simulated. An exon that the variant skips then still receives reads from
  the template and shows up as a coverage dip rather than a hard zero. In
  *variant only* the FASTQ is post-filtered ... and the truth GFF is
  rewritten to remove templates. The alternative-splicing signal is then
  unambiguous: ..."

recount3.Rmd:
- Before: "Loss of TDP-43 lets cryptic exons appear in a number of genes.
  ... STMN2 is the strongest case; the others are weaker ..."
- After: "Loss of TDP-43 lets cryptic exons appear in many genes. ... STMN2
  is the strongest case. The others are weaker ..."
