# ASD Biomarker Discovery: FA_gene & DMN_miRNA

A Python reimplementation of the two-step, network-based biomarker discovery
framework from:

> Rastegari, M., Salehi, N. & Zare-Mirakabad, F. **Biomarker prediction in
> autism spectrum disorder using a network-based approach.** *BMC Medical
> Genomics* **16**, 12 (2023).
> [https://doi.org/10.1186/s12920-023-01439-5](https://doi.org/10.1186/s12920-023-01439-5)
> (Open Access, CC BY 4.0)

The framework identifies candidate ASD gene and miRNA biomarkers from blood
gene-expression data in two stages:

| Step | Algorithm | Method | Output |
|---|---|---|---|
| 1 | **FA_gene** | Control co-expression network (WGCNA-style) → module preservation testing (control vs. autism) → PPI hub selection | 20 candidate genes, $G^A$ |
| 2 | **DMN_miRNA** | mRNA–miRNA bipartite network over $G^A$ → extended greedy set cover | Minimum miRNA regulator set, $\mathfrak{R}$ |

This repo runs on the paper's **own supplementary data** (bundled in
`data/`) and its results are validated directly against the paper's
published Tables 2 and 3 — see [Validation](#validation) below.

---

## Contents

```
.
├── implementation.ipynb        # full pipeline: data loading, FA_gene, DMN_miRNA, validation
├── README.md
└── data/
    ├── 12920_2023_1439_MOESM1_ESM.csv   # autism expression matrix (Additional file 1)
    ├── 12920_2023_1439_MOESM2_ESM.csv   # control expression matrix (Additional file 2)
    ├── 12920_2023_1439_MOESM3_ESM.csv   # PPI network, non-preserved module (Additional file 3)
    └── 12920_2023_1439_MOESM4_ESM.csv   # mRNA–miRNA interactions (Additional file 4)
```

All four files are the paper's own supplementary material (CC BY 4.0),
downloaded from the [article page](https://doi.org/10.1186/s12920-023-01439-5).

## Installation

```bash
pip install pandas numpy scipy networkx scikit-learn matplotlib dynamicTreeCut requests
jupyter notebook implementation.ipynb
```

`dynamicTreeCut` is a Python port of WGCNA's `cutreeDynamic` hybrid
tree-cutting algorithm, used to reproduce the paper's module-detection step
without requiring R. (The notebook applies a one-line `numpy>=2.0`
compatibility shim — `np.in1d = np.isin` — since the package still calls the
removed `np.in1d`.)

## Quick start

With `data/` present as above, the notebook runs the real pipeline
end-to-end in about 60–90 seconds: real GSE18123 expression data (4,701
genes × 21 samples per group) → real co-expression network → real module
preservation testing → the paper's own precomputed PPI network → the
paper's own mRNA–miRNA network → set cover.

If `data/` is empty, the notebook falls back to simulated expression data
with a planted, known co-expression module, so the pipeline still runs
end-to-end and every step remains checkable — just without the real-data
validation.

## Validation

Two cells in the notebook check the reimplementation directly against the
paper's published numbers, using the paper's own supplementary data as
input:

**FA_gene → PPI hub selection**, against Table 2:
```
computed top-20 degree sequence: [333, 253, 248, 213, 158, 155, 148, 148, 137, 136, 134, 129, 127, 127, 124, 123, 122, 118, 117, 117]
paper Table 2 degree sequence  : [333, 253, 248, 213, 158, 155, 148, 148, 137, 136, 134, 129, 127, 127, 124, 123, 122, 118, 117, 117]
EXACT MATCH: True
```
The module-preservation step independently discovers a **1,173-gene**
non-preserved module — exactly the size of the paper's reported
"Turquoise" module — with Zsummary well below the cutoff, consistent with
(if more conservative than) the paper's own reported value of -0.94.

**DMN_miRNA → extended greedy set cover**, against Table 3:
```
recovered R: ['hsa-mir-155-5p', 'hsa-mir-17-5p', 'hsa-mir-181a-5p', 'hsa-mir-18a-5p', 'hsa-mir-92a-1-5p']
EXACT MATCH to paper's published answer: True
```
Run on the real 17-gene / 35-miRNA / 119-edge network from Additional file
4, this recovers the authors' exact published 5-miRNA set, including the
same 3 unreachable genes they report (`CD19`, `LCK`, `ZAP70`).

## Methodology

### FA_gene

1. **Variance filtering** — keep genes with $\mathrm{Var}(\text{control} \cup \text{autism}) \geq 0.75$.
   (Skipped for the bundled data, which is already the paper's post-filtering
   output — see `pre_filtered` in the notebook.)
2. **Co-expression network** — biweight midcorrelation (`bicor`) → signed-hybrid
   adjacency ($a_{ij} = [\mathrm{cor}_{ij}]_+^{power}$, `power=8`) → topological
   overlap matrix (TOM).
3. **Module detection** — average-linkage clustering on `1 - TOM`, cut with
   `dynamicTreeCut` (`minModuleSize=30`, `deepSplit=2`).
4. **Module preservation** — permutation-based $Z_{summary}$
   ($= (Z_{density} + Z_{connectivity})/2$, `nPermutations=200`), following
   Langfelder & Horvath (2011), *PLoS Comput Biol* 7:e1001057. Modules with
   $Z_{summary} < 2$ are flagged non-preserved between control and autism.
5. **PPI + hub selection** — the bundled PPI network (`MOESM3`, STRING
   protein IDs) is the paper's own precomputed network for its non-preserved
   module; the top 20 nodes by degree give $G^A$.

### DMN_miRNA

1. Query mRNA–miRNA interaction databases (TarBase, miRTarBase, miRecords /
   miRNet) for regulators of $G^A$ — the bundled `MOESM4` is the paper's raw
   miRNet export, filtered here to edges targeting $G^A$, forming the
   bipartite network $G_{miR}$.
2. **Extended greedy set cover** — sort miRNAs by degree in $G_{miR}$
   *once*, then make a single pass, adding a miRNA to the result set if it
   covers at least one previously-uncovered gene. This is a deliberate,
   documented deviation from textbook greedy set cover (which re-sorts by
   marginal gain each iteration) and is implemented exactly as specified in
   the paper's Fig. 2 pseudocode.

## Configuration

All parameters mirror the paper and live in a single `Config` class:

```python
class Config:
    VAR_THRESHOLD    = 0.75   # Step 2 variance cutoff
    SOFT_POWER       = 8      # adjacency power
    MIN_MODULE_SIZE  = 30
    DEEP_SPLIT       = 2
    ZSUMMARY_CUTOFF  = 2.0
    N_PERMUTATIONS   = 200
    PPI_SCORE_CUTOFF = 400    # only used for non-bundled, generic STRING exports
    TOP_N_GENES      = 20     # |G^A|
```

## Notes & limitations

- **PPI node labels**: `MOESM3` is keyed by STRING protein ID
  (`9606.ENSP...`), not gene symbol. 14 of the top-20 hub nodes map
  unambiguously to a gene symbol by exact degree-value matching against
  Table 2; the other 6 form 3 tied-degree pairs (148, 127, 117) that degree
  alone can't resolve to an individual gene without a STRING ID mapping
  service (not reachable offline). The notebook leaves these as raw IDs
  rather than guessing — the provable claim is the degree *sequence* match
  shown above, not a claim about every individual tied label.
- **Chained end-to-end run**: the "Full pipeline" cell feeds FA_gene's own
  output straight into DMN_miRNA, so its miRNA coverage is naturally partial
  (3 of 5) for the 6 tied-degree genes above — this is expected, not a bug;
  the two dedicated validation cells are the exact, unambiguous checks.
- **Module preservation**: this implementation averages the two components
  the paper names ($Z_{density}$, $Z_{connectivity}$); WGCNA's full
  `modulePreservation` blends several more sub-statistics. On the real data
  this simplified version is more conservative — it flags 3 modules
  (including the exact 1,173-gene module the paper calls "Turquoise") where
  the paper's fuller statistic flags only 1. For byte-for-byte replication,
  use R's `WGCNA::modulePreservation` (e.g. via `rpy2`).
- This is a **biomarker discovery** pipeline over a sample cohort, not a
  per-sample diagnostic classifier — it outputs candidate genes and miRNAs,
  not case/control predictions for new individuals.

## Citation

If you use this code, please cite the original paper:

```bibtex
@article{rastegari2023biomarker,
  title   = {Biomarker prediction in autism spectrum disorder using a network-based approach},
  author  = {Rastegari, Maryam and Salehi, Najmeh and Zare-Mirakabad, Fatemeh},
  journal = {BMC Medical Genomics},
  volume  = {16},
  number  = {1},
  pages   = {12},
  year    = {2023},
  doi     = {10.1186/s12920-023-01439-5}
}
```

## License

Code: add a license of your choice (e.g. MIT) before publishing — none is
specified here. Data in `data/`: CC BY 4.0, per the original article.
