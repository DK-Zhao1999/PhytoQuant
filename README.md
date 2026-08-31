# PhytoQuant

PhytoQuant is a cohort-resolved, direction-aware transcriptome-signature discovery engine specifically designed for plant functional genomics. Unlike conventional differential expression workflows that reduce multi-cohort evidence to a flat list of genes, PhytoQuant triangulates three orthogonal statistical frameworks—DESeq2, edgeR, and limma-voom—under distinct distributional assumptions, fusing their evidence into a signed gene-level importance score. This statistical triangulation, rather than a simple consensus vote, attenuates method-specific biases while preserving biological directionality. Through cohort-level Robust Rank Aggregation, PhytoQuant separates reproducible regulatory signals from cohort-specific technical variation, enabling robust cross-study meta-analysis across ecotypes, stress treatments, developmental stages, or genetic backgrounds. The framework further models activation and inhibition as separate evidence streams, avoiding the information loss inherent in unsigned enrichment scores, and combines prior pathway annotations with data-driven de novo candidates through set-theoretic concordance filtering. Finally, redundancy-resolving consolidation via Sorensen-Dice similarity and hierarchical clustering delivers compact, directionally labeled gene modules—designated as ags* (activation) and igs* (inhibition)—that are immediately reusable for enrichment analysis, signature scoring, validation cohort testing, and wet-lab candidate prioritization. The result is a transparent, auditable, and experimentally tractable pathway signature rather than a static list of differentially expressed genes, positioning PhytoQuant as a principled solution for hypothesis-driven discovery in plant multi-cohort studies.

## Highlights

- **Statistical triangulation, not consensus voting.** DESeq2, edgeR, and
  limma-voom operate in parallel under distinct distributional assumptions.
  Their evidence is fused into one signed importance score, attenuating
  method-specific biases while preserving biological direction.
- **Cohort-level rank meta-analysis.** Robust Rank Aggregation treats every
  independent plant experiment as a ranking source, separating reproducible
  regulatory signal from cohort-specific technical variation.
- **Direction-resolved signature algebra.** Activation and inhibition are
  modeled as separate evidence streams, avoiding the information loss that
  occurs when unsigned enrichment scores collapse opposing regulators into one
  value.
- **Knowledge-aware de novo discovery.** Prior pathway annotations and
  data-driven de novo candidates are combined through set-theoretic concordance
  filtering, connecting hypothesis-driven biology with discovery-driven
  statistics.
- **Redundancy-resolving consolidation.** Sorensen-Dice similarity, hierarchical
  clustering, and minimum-cardinality filtering jointly remove overlapping
  modules while retaining the most parsimonious interpretable gene sets.
- **Translation-first output.** Final `ags*` activation and `igs*` inhibition
  panels are compact, named, and directly reusable for enrichment analysis,
  signature scoring, validation cohorts, and candidate prioritization.

## Pipeline overview

```mermaid
flowchart TB
    subgraph Input
        A["Plant count matrices<br/>(<dataset>.Rdata)"]
        B["Group annotations<br/>(<dataset>_G.Rdata)"]
        C["Prior knowledge gene sets"]
    end

    subgraph DE["Differential-expression engines"]
        D1["DESeq2"]
        D2["edgeR"]
        D3["limma-voom"]
    end

    subgraph Evidence["Evidence synthesis"]
        E1["Gene-wise merge"]
        E2["Signed importance score<br/>= logFC x -log10(p)"]
        E3["De novo up/down gene sets"]
    end

    subgraph RRA["Robust Rank Aggregation"]
        R1["Cross-cohort rank integration"]
        R2["Significance threshold"]
    end

    subgraph Refine["Signature refinement"]
        K1["Knowledge + de novo integration"]
        K2["Concordance filtering"]
        K3["Sorensen-Dice clustering"]
        K4["Minimum-size filter"]
    end

    subgraph Output
        O1["Activation gene sets (ags*)"]
        O2["Inhibition gene sets (igs*)"]
        O3["Importance matrices"]
    end

    A --> D1
    A --> D2
    A --> D3
    B --> D1
    B --> D2
    B --> D3
    D1 --> E1
    D2 --> E1
    D3 --> E1
    E1 --> E2 --> E3
    E3 --> R1 --> R2
    C --> K1
    E3 --> K1
    R2 --> K1
    K1 --> K2 --> K3 --> K4
    K4 --> O1
    K4 --> O2
    R2 --> O3
```

## Detailed plant-oriented workflow

```mermaid
flowchart TB

    %% ================= PLANT MULTI-COHORT INPUT =================
    subgraph INP["Plant multi-cohort input"]
        direction TB
        IN_COUNT["Raw RNA-seq count matrices<br/>genes x samples"]
        IN_GROUP["Group annotations<br/>Tag / group (H vs L)"]
        IN_PRIOR["Plant knowledge gene sets<br/>auxin, ethylene, GA, JA, SA, and others"]
        IN_META["Cohort context<br/>ecotype, tissue, stress, development, genotype"]
    end

    %% ================= DATA CONTRACT AND QUALITY CONTROL =================
    subgraph QC["Data contract and quality control"]
        direction TB
        QC_NAME["Resolve object names<br/><dataset>.Rdata and <dataset>_G.Rdata"]
        QC_MATCH["Align samples between counts and annotations"]
        QC_CONTRAST["Fix contrast orientation<br/>levels = c('L', 'H')"]
        QC_MISSING{"Missing cohort files?"}
        QC_SKIP["Warn and skip cohort"]
    end

    %% ================= TRIANGULATED DIFFERENTIAL EXPRESSION =================
    subgraph DE["Triangulated differential expression"]
        direction TB
        DE_DES2["DESeq2<br/>size factors + NB dispersion + Wald test"]
        DE_EDGER["edgeR<br/>norm factors + estimateDisp + QL F-test"]
        DE_VOOM["limma-voom<br/>precision weights + lmFit + eBayes"]
        DE_TRY{"Method error handling"}
        DE_FAIL["Warn and retain successful methods"]
    end

    %% ================= SIGNED EVIDENCE FUSION =================
    subgraph FUSE["Signed evidence fusion"]
        direction TB
        FUSE_MERGE["Gene-wise merge of method tables"]
        FUSE_NA["Repair missing p-values and logFC"]
        FUSE_P["Geometric-mean combined p-value"]
        FUSE_LOGFC["Mean log-fold change"]
        FUSE_IMP["Importance = logFC x -log10(p)"]
        FUSE_ORDER["Rank genes by signed importance"]
        FUSE_UP["De novo activation candidates"]
        FUSE_DN["De novo inhibition candidates"]
        FUSE_NO{"No successful method?"}
        FUSE_NEXT["Skip cohort"]
    end

    %% ================= CROSS-COHORT ROBUST RANK AGGREGATION =================
    subgraph RRA["Cross-cohort robust rank aggregation"]
        direction TB
        RRA_UP_RANK["Rank genes by decreasing importance"]
        RRA_DN_RANK["Rank genes by increasing importance"]
        RRA_AGG_UP["aggregateRanks: activation evidence"]
        RRA_AGG_DN["aggregateRanks: inhibition evidence"]
        RRA_THRESH{"Score < rra_p_threshold?"}
        RRA_UP_SIG["Reproducible up-regulated genes"]
        RRA_DN_SIG["Reproducible down-regulated genes"]
    end

    %% ================= KNOWLEDGE-GUIDED REFINEMENT =================
    subgraph KNOW["Knowledge-guided refinement"]
        direction TB
        KNOW_PRIOR["Phytohormone and pathway gene sets"]
        KNOW_DENOVO["Per-cohort de novo sets"]
        KNOW_UNION["Candidate set pool"]
        KNOW_CONC_UP["Intersect activation candidates with RRA-up"]
        KNOW_CONC_DN["Intersect inhibition candidates with RRA-down"]
        KNOW_FILTER{"Non-empty after filtering?"}
    end

    %% ================= REDUNDANCY-RESOLVING CONSOLIDATION =================
    subgraph DEDUP["Redundancy-resolving consolidation"]
        direction TB
        DEDUP_SIM["Pairwise Sorensen-Dice similarity"]
        DEDUP_DIST["Convert similarity to distance"]
        DEDUP_HC["Hierarchical clustering"]
        DEDUP_CUT["cutree at 1 - sorenson_threshold"]
        DEDUP_MERGE["Merge cluster members by set union"]
        DEDUP_UNIQUE["Remove duplicated modules"]
        DEDUP_SIZE["Enforce min_geneset_size"]
    end

    %% ================= PARAMETER CONTROLS =================
    subgraph PAR["Parameter controls"]
        direction TB
        PAR_VECTOR["vector: cohort names"]
        PAR_DIR["data_dir: input directory"]
        PAR_METHODS["methods: DESeq2 / edgeR / limma"]
        PAR_IMP["importance_threshold"]
        PAR_NA["max_na_ratio"]
        PAR_RRA["rra_p_threshold"]
        PAR_SDS["sorenson_threshold"]
        PAR_LINK["linkage_method"]
        PAR_MIN["min_geneset_size"]
        PAR_SAVE["save_intermediate"]
    end

    %% ================= PLANT TRANSLATION-READY OUTPUT =================
    subgraph OUT["Plant translation-ready output"]
        direction TB
        OUT_AGS["Activation modules<br/>(ags*)"]
        OUT_IGS["Inhibition modules<br/>(igs*)"]
        OUT_IMP["Gene-by-cohort importance matrix"]
        OUT_RRA["RRA evidence tables"]
        OUT_DOWN["Enrichment, scoring, validation, prioritization"]
    end

    %% ================= EDGES =================
    IN_COUNT --> QC_NAME
    IN_GROUP --> QC_NAME
    PAR_VECTOR --> QC_NAME
    PAR_DIR --> QC_NAME

    QC_NAME --> QC_MATCH --> QC_CONTRAST --> QC_MISSING
    QC_MISSING -->|yes| QC_SKIP
    QC_MISSING -->|no| DE_DES2
    QC_MISSING -->|no| DE_EDGER
    QC_MISSING -->|no| DE_VOOM
    PAR_METHODS --> DE_DES2
    PAR_METHODS --> DE_EDGER
    PAR_METHODS --> DE_VOOM

    DE_DES2 --> DE_TRY
    DE_EDGER --> DE_TRY
    DE_VOOM --> DE_TRY
    DE_TRY -->|success| FUSE_MERGE
    DE_TRY -->|error| DE_FAIL
    DE_FAIL --> FUSE_MERGE

    FUSE_MERGE --> FUSE_NA --> FUSE_P
    FUSE_MERGE --> FUSE_NA --> FUSE_LOGFC
    FUSE_P --> FUSE_IMP
    FUSE_LOGFC --> FUSE_IMP
    FUSE_IMP --> FUSE_ORDER
    FUSE_ORDER --> FUSE_UP
    FUSE_ORDER --> FUSE_DN
    PAR_IMP --> FUSE_UP
    PAR_IMP --> FUSE_DN
    PAR_NA --> FUSE_NA

    FUSE_UP --> FUSE_NO
    FUSE_DN --> FUSE_NO
    FUSE_NO -->|no result| FUSE_NEXT
    FUSE_NO -->|continue| RRA_UP_RANK
    FUSE_NO -->|continue| RRA_DN_RANK
    FUSE_UP --> RRA_UP_RANK
    FUSE_DN --> RRA_DN_RANK

    RRA_UP_RANK --> RRA_AGG_UP
    RRA_DN_RANK --> RRA_AGG_DN
    RRA_AGG_UP --> RRA_THRESH
    RRA_AGG_DN --> RRA_THRESH
    PAR_RRA --> RRA_THRESH
    RRA_THRESH -->|up| RRA_UP_SIG
    RRA_THRESH -->|down| RRA_DN_SIG

    IN_PRIOR --> KNOW_PRIOR
    FUSE_UP --> KNOW_DENOVO
    FUSE_DN --> KNOW_DENOVO
    KNOW_PRIOR --> KNOW_UNION
    KNOW_DENOVO --> KNOW_UNION
    RRA_UP_SIG --> KNOW_CONC_UP
    RRA_DN_SIG --> KNOW_CONC_DN
    KNOW_UNION --> KNOW_CONC_UP
    KNOW_UNION --> KNOW_CONC_DN
    KNOW_CONC_UP --> KNOW_FILTER
    KNOW_CONC_DN --> KNOW_FILTER

    KNOW_FILTER -->|activation| DEDUP_SIM
    KNOW_FILTER -->|inhibition| DEDUP_SIM
    DEDUP_SIM --> DEDUP_DIST --> DEDUP_HC --> DEDUP_CUT --> DEDUP_MERGE --> DEDUP_UNIQUE --> DEDUP_SIZE
    PAR_SDS --> DEDUP_CUT
    PAR_LINK --> DEDUP_HC
    PAR_MIN --> DEDUP_SIZE

    DEDUP_SIZE --> OUT_AGS
    DEDUP_SIZE --> OUT_IGS
    FUSE_IMP --> OUT_IMP
    RRA_UP_SIG --> OUT_RRA
    RRA_DN_SIG --> OUT_RRA
    PAR_SAVE --> OUT_IMP

    OUT_AGS --> OUT_DOWN
    OUT_IGS --> OUT_DOWN
    OUT_IMP --> OUT_DOWN
    OUT_RRA --> OUT_DOWN
```

## Methodological novelty and publication value

PhytoQuant addresses a recurring bottleneck in plant functional genomics:
multi-cohort studies are abundant, but conventional differential-expression
workflows rarely propagate the uncertainty and direction of evidence across
cohorts. PhytoQuant reframes the problem as a hierarchical evidence-integration
task, making it especially suitable for comparative studies across ecotypes,
stress treatments, developmental stages, tissues, or genetic backgrounds.

- **A reproducible analytical object.** Every result retains the underlying
  importance matrix and rank-aggregation statistics, so reviewers can trace how
  a final signature was derived instead of receiving only a gene list.
- **Methodological heterogeneity as a feature.** By deliberately combining
  negative-binomial and linear-model frameworks, PhytoQuant rewards genes whose
  signal is robust across statistical paradigms rather than overfitted to one
  package.
- **Directional biological interpretation.** Separate activation and inhibition
  panels make it possible to propose regulatory hypotheses directly from the
  output, including candidate activators, repressors, and pathway-level switches.
- **Immediate translational utility.** The compact `ags*` and `igs*` modules
  are well suited to signature scoring, cross-study validation, co-expression
  network annotation, and prioritization of targets for wet-lab validation.

## Installation

Start from a clean R session and install the development dependencies:

```r
install.packages(c("pak", "remotes"))
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "RobustRankAggreg"))
```

Then install PhytoQuant directly from GitHub:

```r
pak::pkg_install("DK-Zhao1999/PhytoQuant")
```

The `remotes` equivalent is:

```r
remotes::install_github("DK-Zhao1999/PhytoQuant", upgrade = FALSE)
```

Load the package:

```r
library(PhytoQuant)
```

## Quick start with simulated plant data

The complete example below is designed to run in a fresh R session and to mimic
the heterogeneity of a multi-cohort plant study. It first constructs three
independent RNA-seq-like count cohorts with planted directional signal, then
uses `diffexp_integrate()` to recover the planted activation and inhibition
modules through the full evidence-integration workflow.

```r
library(PhytoQuant)

# Create a temporary directory for the example data.
data_dir <- tempfile("phyto_demo_")
dir.create(data_dir)

set.seed(42)
n_genes <- 100
n_h <- 6
n_l <- 6

genes <- paste0("PHYTO", sprintf("%03d", seq_len(n_genes)))
up_genes <- paste0("PHYTO", sprintf("%03d", 1:10))
dn_genes <- paste0("PHYTO", sprintf("%03d", 11:20))

# Simulate a negative-binomial count matrix with planted directional signal.
make_counts <- function(seed, genes, up_genes, dn_genes, n_h, n_l) {
  set.seed(seed)
  n_samples <- n_h + n_l
  mu <- matrix(100, nrow = length(genes), ncol = n_samples)
  mu[match(up_genes, genes), seq_len(n_h)] <- 1000
  mu[match(dn_genes, genes), n_h + seq_len(n_l)] <- 1000
  matrix(
    rnbinom(length(genes) * n_samples, size = 20, mu = mu),
    nrow = length(genes),
    ncol = n_samples,
    dimnames = list(genes, paste0("S", seq_len(n_samples)))
  )
}

# Save each simulated cohort in the layout expected by PhytoQuant.
for (dataset_name in c("PhytoA", "PhytoB", "PhytoC")) {
  seed <- switch(dataset_name, PhytoA = 1, PhytoB = 2, PhytoC = 3)
  expr <- make_counts(seed, genes, up_genes, dn_genes, n_h, n_l)
  group <- data.frame(
    Tag = colnames(expr),
    group = c(rep("H", n_h), rep("L", n_l)),
    stringsAsFactors = FALSE
  )

  assign(dataset_name, expr)
  assign(paste0(dataset_name, "_G"), group)
  save(list = dataset_name, file = file.path(data_dir, paste0(dataset_name, ".Rdata")))
  save(list = paste0(dataset_name, "_G"), file = file.path(data_dir, paste0(dataset_name, "_G.Rdata")))
}

# Prior knowledge gene sets used as plant pathway hypotheses.
knowledge_genesets <- list(
  auxin_signaling = up_genes,
  ethylene_signaling = dn_genes,
  gibberellin_biosynthesis = paste0("PHYTO", sprintf("%03d", 21:30)),
  jasmonate_response = paste0("PHYTO", sprintf("%03d", 31:40)),
  salicylic_acid_pathway = paste0("PHYTO", sprintf("%03d", 41:50))
)

result <- diffexp_integrate(
  vector = c("PhytoA", "PhytoB", "PhytoC"),
  geneSets_list = knowledge_genesets,
  data_dir = data_dir,
  importance_threshold = 1,
  rra_p_threshold = 0.01,
  sorenson_threshold = 0.7,
  min_geneset_size = 3,
  save_intermediate = FALSE,
  methods = c("DESeq2", "edgeR", "limma"),
  linkage_method = "complete",
  verbose = TRUE
)

str(result$activation_gene_sets)
str(result$inhibition_gene_sets)
head(result$importance_all)
```

The planted activation genes are `PHYTO001`-`PHYTO010`, and the planted
inhibition genes are `PHYTO011`-`PHYTO020`. When the three cohorts are analyzed
and integrated, PhytoQuant should return these genes inside the compact final
`ags*` and `igs*` panels, while the returned importance matrices provide the
full quantitative evidence trail for downstream reuse.

## Input data format

PhytoQuant uses a minimal, cohort-oriented data contract. For every name in
`vector`, two RData files must be present in `data_dir`:

```text
<name>.Rdata
<name>_G.Rdata
```

The first file must contain a raw count matrix assigned to an object named
`<name>`, with genes as rows and samples as columns. The second file must
contain a sample-level annotation data frame assigned to `<name>_G`:

```text
Tag    group
S1     H
S2     H
S3     L
S4     L
```

The `group` values are interpreted relative to `levels = c("L", "H")`. This
fixed contrast orientation ensures that every downstream signed statistic,
including the final importance score, has a consistent biological meaning: a
positive log-fold change corresponds to higher expression in `H`.

## Function reference

### `diffexp_integrate()`

```r
diffexp_integrate(
  vector,
  geneSets_list,
  data_dir = getwd(),
  importance_threshold = 10,
  rra_p_threshold = 0.00001,
  sorenson_threshold = 0.6,
  min_geneset_size = 5,
  max_na_ratio = 0.3,
  save_intermediate = FALSE,
  methods = c("DESeq2", "edgeR", "limma"),
  linkage_method = "complete",
  verbose = TRUE
)
```

| Argument | Description | Default |
| --- | --- | --- |
| `vector` | Cohort identifiers whose rank evidence will be integrated. | required |
| `geneSets_list` | Named prior-knowledge gene sets used to guide and refine discovery. | required |
| `data_dir` | Directory containing the count and annotation RData files. | `getwd()` |
| `importance_threshold` | Signed-evidence cutoff used to partition de novo up- and down-regulated genes. | `10` |
| `rra_p_threshold` | Cross-cohort significance cutoff for robust rank aggregation. | `0.00001` |
| `sorenson_threshold` | Similarity cutoff controlling redundancy removal among final modules. | `0.6` |
| `min_geneset_size` | Minimum cardinality required for a reproducible final module. | `5` |
| `max_na_ratio` | Maximum fraction of cohorts in which a gene may lack evidence. | `0.3` |
| `save_intermediate` | Persist the per-cohort merged evidence tables for external audit. | `FALSE` |
| `methods` | Orthogonal differential-expression backends used for triangulation. | `c("DESeq2", "edgeR", "limma")` |
| `linkage_method` | Agglomerative linkage rule used during similarity clustering. | `"complete"` |
| `verbose` | Emit trace-level progress diagnostics during the workflow. | `TRUE` |

The function returns:

```r
list(
  activation_gene_sets,
  inhibition_gene_sets,
  importance_all,
  importance_up_gene,
  importance_dn_gene
)
```

## Notes

1. PhytoQuant is designed around raw count matrices rather than normalized
   expression values. This preserves the distributional assumptions of DESeq2
   and edgeR and lets the workflow maintain a coherent, auditable evidence chain
   from counts to final modules.
2. Each cohort is analyzed independently before rank-level integration. Missing
   genes are therefore handled explicitly by `max_na_ratio`, avoiding the
   complete-case bias introduced by premature sample-level intersection.
3. For exploratory runs on small or deliberately noisy simulations, relaxed
   thresholds such as `importance_threshold = 1` and `rra_p_threshold = 0.01`
   expose the full pipeline output without altering the statistical method.

## Author

Dingkang Zhao

Email: <dingkang.25@intl.zju.edu.cn>
