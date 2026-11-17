# Identification and Prioritization of Plastic Additives in Chinese Agricultural Soils Using High-Resolution Mass Spectrometry and Network Annotation Propagation

<p align="center">
  <img src="https://img.shields.io/badge/R-HRMS%20processing-276DC3?logo=r&logoColor=white" alt="R">
  <img src="https://img.shields.io/badge/Python-Benchmark%20%26%20mNAP-3776AB?logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/Jupyter-Notebook-F37626?logo=jupyter&logoColor=white" alt="Jupyter">
  <img src="https://img.shields.io/badge/HRMS-Orbitrap-5A67D8" alt="HRMS">
  <img src="https://img.shields.io/badge/GNPS-Molecular%20Networking-2F855A" alt="GNPS">
  <img src="https://img.shields.io/badge/MetFrag-In%20silico%20fragmentation-7B61A8" alt="MetFrag">
</p>

---

## 🧭 Overview

This repository contains the data-processing and computational workflows supporting the study:

**Identification and Prioritization of Plastic Additives in Chinese Agricultural Soils Using High-Resolution Mass Spectrometry and Network Annotation Propagation**

The repository is organized into three main modules:

| Module | Description | Language |
|---|---|---|
| 🔬 `1_HRMS_data_processing.R` | HRMS feature processing, suspect screening, transformation-product generation, diagnostic MS/MS filtering, and structural annotation | ![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white) |
| 📊 `2_Benchmark_evaluation/` | Benchmarking of MetFrag, ChemWalker, and modified Network Annotation Propagation (mNAP) | ![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white) |
| 🕸️ `3_mnap_application/` | Application of mNAP to a GNPS molecular network for candidate re-ranking | ![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white) |

---

## 📁 Repository Structure

```text
.
├── 🔬 1_HRMS_data_processing.R
│
├── 📊 2_Benchmark_evaluation/
│   ├── benchmark_nist.py
│   ├── chemwalker_official_worker.py
│   └── run_benchmark.ipynb
│
└── 🕸️ 3_mnap_application/
    └── mNAP_GNPS_application.ipynb
```

---

# 🔬 1. HRMS_data_processing

![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white)
![patRoon](https://img.shields.io/badge/patRoon-HRMS%20workflow-4C78A8)
![ProteoWizard](https://img.shields.io/badge/ProteoWizard-mzML%20conversion-6B7280)

`HRMS_data_processing.R` provides the HRMS data-processing workflow implemented using the **patRoon** R package.

### Main workflow

```text
Raw HRMS data
      ↓
Feature detection and alignment
      ↓
Feature filtering
      ↓
Isotope/adduct componentization
      ↓
Suspect screening
      ↓
Transformation-product generation
      ↓
MS/MS processing
      ↓
Diagnostic MS/MS filtering
      ↓
Formula and structural annotation
```

The workflow includes:

- conversion of vendor raw files to centroided mzML using ProteoWizard;
- feature detection, grouping, alignment, and blank filtering;
- isotope and adduct componentization;
- suspect screening and transformation-product generation;
- MS/MS peak-list processing;
- diagnostic MS/MS filtering based on characteristic fragment ions and neutral losses;
- molecular formula generation using **GenForm**; and
- structural candidate generation using **MetFrag** and PubChemLite.

The diagnostic MS/MS filtering section uses organophosphorus compounds as an example. The diagnostic fragment ions and neutral losses can be replaced with class-specific features for other chemical groups.

**Main script**

```text
1_HRMS_data_processing.R
```

---

# 📊 2. Benchmark_evaluation

![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![Jupyter](https://img.shields.io/badge/Jupyter-F37626?logo=jupyter&logoColor=white)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-00A6D6)
![MetFrag](https://img.shields.io/badge/MetFrag-Candidate%20ranking-7B61A8)

This module provides a reproducible comparison of:

- 🧩 **MetFrag**
- 🚶 **ChemWalker**
- 🧠 **mNAP**

using the NIST tandem mass spectra selected for the published network-annotation benchmark.

The three methods are evaluated using the **same target spectra, seed assignments, and MetFrag candidate pools**, allowing direct comparison of candidate re-ranking performance.

### Main files

- **`benchmark_nist.py`** – main benchmark implementation.
- **`run_benchmark.ipynb`** – stepwise Jupyter interface for benchmark execution.
- **`chemwalker_official_worker.py`** – subprocess wrapper for running the original ChemWalker `cand_pair` and `random_walk` functions safely under Windows/Jupyter.

Installing the ChemWalker package does not replace `chemwalker_official_worker.py`. Keep the worker beside `benchmark_nist.py`. During the benchmark, the main script automatically starts this worker using the Python interpreter of the active environment, passes the network, candidates, seeds, and parameters through an input file, and reads back the results. The worker calls the installed ChemWalker functions in a separate process to provide a safe multiprocessing entry point on Windows. It does not need to be started manually and is not called by the mNAP application notebook.

```text
run_benchmark.ipynb
    -> benchmark_nist.py
    -> chemwalker_official_worker.py
    -> installed ChemWalker cand_pair() and random_walk()
```

### Benchmark workflow

```text
NIST MS/MS spectra
        ↓
Molecular network
        ↓
MetFrag candidate generation
        ↓
Seed allocation
        ↓
┌───────────┬────────────┬──────────┐
│  MetFrag  │ ChemWalker │   mNAP   │
└───────────┴────────────┴──────────┘
        ↓
Candidate ranking evaluation
```

The molecular network is constructed using modified cosine similarity with:

- cosine similarity ≥ **0.60**;
- at least **2 matched product ions**; and
- mutual Top-K filtering.

Approximately **10% of nodes in each connected component** are assigned as seeds, with repeated seed allocations used to evaluate robustness. This controlled three-method benchmark still requires seeds because the original ChemWalker comparator is seed-based. The default seed allocation and benchmark reports are unchanged; setting `SEED_FRACTION = 0` is not a supported shortcut for a seed-free three-method comparison.

For direct seed-free mNAP scoring with already prepared candidates, `compute_mnap_for_component(spectral_edges, full_tlid)` also supports candidate tables with `is_seed = False` for every row. Use the existing candidate schema (`node_id`, `uid`, `Identifier`, `InChI`, `Score`, `score_norm`, and `is_seed`) and spectral-edge schema (`node1`, `node2`, and `cosine`). This returns mNAP scores and ranks without running ChemWalker or changing the controlled benchmark.

For the controlled benchmark, MetFrag candidate lists are standardized to enable consistent comparison among methods.

### ▶️ Recommended execution

Start JupyterLab from `E:\plastic-additives\2_Benchmark_evaluation`, open `run_benchmark.ipynb`, and select **Python (mNAP)**. The notebook uses `PROJECT_DIR = Path.cwd()` to locate `benchmark_nist.py`; if the working directory differs, set `PROJECT_DIR` explicitly to the benchmark folder in the first code cell.

```python
import benchmark_nist as bm

# Set the actual locations of the benchmark inputs and MetFrag JAR first.
# The configuration cell in run_benchmark.ipynb calls bm.configure_paths(...).

# Validate inputs and software
bm.run_step(0)

# Build the molecular network
bm.run_step(1, force=True)

# Validate MetFrag output parsing
bm.run_step(2, force=True, max_nodes=5)

# Generate complete candidate lists
bm.run_step(2, force=False)

# Run the benchmark
bm.run_step(3, force=False)

# Evaluate ranking performance
outputs = bm.run_step(4)
```

The benchmark outputs include candidate coverage, Top-k accuracy, mean reciprocal rank, pairwise comparisons, bootstrap analyses, and publication-quality figures.

The benchmark requires `nist.msp`, `pcbi.1006089.s007.tsv`, a PubChemLite database, and the MetFrag JAR. Set their locations in the notebook's `bm.configure_paths(...)` cell before Step 0. In the supplied local folder, the JAR is `E:\plastic-additives\2_Benchmark_evaluation\MetFrag-CL.jar`; adjust `metfrag_path` if the configuration cell still points to a `software/` subfolder.

---

# 🕸️ 3. mnap_application

![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![GNPS](https://img.shields.io/badge/GNPS-Molecular%20Networking-2F855A)
![RDKit](https://img.shields.io/badge/RDKit-Tanimoto%20similarity-00A6D6)
![MetFrag](https://img.shields.io/badge/MetFrag-Candidates-7B61A8)

`mNAP_GNPS_application.ipynb` demonstrates the application of mNAP to an experimental **GNPS feature-based molecular network**.

### Required inputs

```text
clusterinfo_summary.tsv
networking_pairs.tsv
spectra.mgf
seed_nodes.csv  (optional)
PubChemLite.csv
MetFrag-CL.jar
```

The network files and PubChemLite database are read from `DATA_DIR`. The MetFrag JAR can be stored separately; set `METFRAG_PATH` to its actual location. The notebook reads locally downloaded GNPS files and does not download them automatically.

Verified structures supplied in `seed_nodes.csv` are used as **seed nodes**, whereas the remaining nodes in the selected molecular-network component are treated as unknown features. Seed structures are optional: set `SEED_PATH = None` in the input-path cell to run without any seed file. An empty seed table or a selected component with no seed annotations is also supported. A configured but missing seed file still raises a file error, so that seed annotations are not silently ignored.

When no verified seeds are present in the selected component, every spectrum retains its MetFrag candidate list. Each spectrum receives equal total restart probability, distributed among its candidates according to their normalized MetFrag scores (with a minimum score of `score_eps = 0.05`). A random walk with restart then propagates support over the existing weighted candidate graph. Candidate scores can therefore change when neighboring spectra have structurally compatible candidates supported by their MetFrag scores. Spectra with no usable candidate connections retain their original MetFrag ordering.

When verified seeds are available, each seed retains its confirmed structure and the existing seed-based random-walk implementation is used. This update does not change the original seeded propagation weights, candidate filtering, or ranking normalization. The default restart probability remains 0.5 for the application notebook and 0.1 for the controlled benchmark; it has not been retuned.

After changing `SEED_PATH` or selecting another component, restart the notebook kernel and execute the cells from top to bottom. Set `component_id` to the component to analyze and `target_node` to an unknown node in that component. In a component without seeds, all of its nodes are unknown targets.

### Application workflow

```text
GNPS molecular network
        ↓
Select network component
        ↓
Identify verified seeds, if available
        ↓
MetFrag candidates for unknown nodes
        ↓
Spectral similarity
        +
Structural similarity
        +
MetFrag score
        ↓
mNAP propagation
        ↓
Candidate re-ranking
```

The notebook first summarizes all available GNPS network components before a component is selected for analysis.

For the example dataset:

```python
component_id = 6164
```

The verified seeds and unknown nodes within this component are then identified automatically.

One unknown node can subsequently be selected for detailed inspection:

```python
target_node = 88874
```

The target node is selected **after** the complete component has been processed and does not alter the mNAP propagation network.

### 🧠 mNAP candidate re-ranking

mNAP integrates:

- MS/MS spectral similarity between connected features;
- structural similarity between candidate molecules;
- normalized MetFrag candidate scores; and
- structural information propagated from verified seed nodes, when available.

Structural similarity is calculated using the **Tanimoto coefficient** based on RDKit molecular fingerprints.

Unlike the controlled benchmark, the experimental application retains **all valid MetFrag candidates** without applying a Top-50 cutoff.

### 📤 Outputs

Two result modes are provided:

**Single target node**

```text
mNAP_node_88874.csv
```

**All unknown nodes in the selected component**

```text
mNAP_all_unknown_nodes.csv
```

Candidate structures before and after mNAP re-ranking can additionally be exported as publication-quality **SVG** and **PDF** figures.

---

# ⚙️ Software Requirements

### ![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white) R workflow

- R
- patRoon
- ProteoWizard
- external software required by the selected patRoon annotation workflow

### ![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white) Python workflow

- pandas
- numpy
- networkx
- RDKit
- matchms *(benchmark only; not required by the application notebook)*
- pyteomics
- scipy
- matplotlib
- ChemWalker
- Jupyter

### ☕ Additional software

- Java
- MetFrag command-line JAR
- CairoSVG *(optional; PDF structure export)*

### Local installation with conda on Windows

The following commands create a dedicated environment for the supplied application notebook. Run them in **Anaconda Prompt** or **Miniforge Prompt**. The version constraints are recommended for this code; install the dependencies first, then install ChemWalker without replacing them.

```bat
conda create -n mnap -c conda-forge --override-channels python=3.11 numpy=1.26 pandas=2.2 scipy=1.14 networkx=3.4 rdkit=2024.03 scikit-learn=1.5 matplotlib=3.8 jupyterlab ipykernel requests xmltodict click pillow cairosvg openjdk=19 git pip -y
conda activate mnap
python -m pip install "pyteomics==4.7.2"
python -m pip install --no-deps git+https://github.com/computational-chemical-biology/ChemWalker.git
python -m pip check
```

ChemWalker is installed from its [official repository](https://github.com/computational-chemical-biology/ChemWalker). The local mNAP implementation is supplied by this repository's scripts and notebook; installing ChemWalker alone does not install mNAP.

Check the imports, seeded random walk, and Java executable:

```bat
python -c "from rdkit import Chem; from chemwalker.gnps import Proteosafe; from chemwalker.utils import run_metfrag; from chemwalker.rwalker import random_walk; import networkx as nx; g=nx.Graph(); g.add_edge('a','b',weight=1.0); print(random_walk(g,['a'],restart_prob=0.5)); print('mNAP dependencies OK')"
java -version
```

The first command should finish with `mNAP dependencies OK`, and the second should display the Java version. These checks verify the main imports and random walk; MetFrag execution is checked when the notebook runs candidate generation.

If running the NIST benchmark as well, install its additional dependency in the same environment:

```bat
python -m pip install "matchms==0.29.0"
python -m pip check
```

### Start the application notebook

Register the environment as a Jupyter kernel and launch JupyterLab from the application folder:

```bat
python -m ipykernel install --user --name mnap --display-name "Python (mNAP)"
cd /d E:\plastic-additives\3_mnap_application
python -m jupyterlab
```

Open `mNAP_GNPS_application.ipynb` and select **Python (mNAP)** as the kernel. Update the first code cell to use the local paths, for example:

```python
DATA_DIR = r"E:\plastic-additives\3_mnap_application\data"
PUBCHEMLITE_PATH = os.path.join(DATA_DIR, "PubChemLite.csv")
METFRAG_PATH = r"E:\plastic-additives\2_Benchmark_evaluation\MetFrag-CL.jar"
SEED_PATH = None
```

For a run with verified seeds, replace the last line with:

```python
SEED_PATH = os.path.join(DATA_DIR, "seed_nodes.csv")
```

Use the component summary to choose `component_id`, and select `target_node` from the unknown nodes in that component. The supplied example uses `component_id = 6164` and `target_node = 88874`. Restart the kernel and execute all cells in order after changing the seed settings or component. To confirm that the correct environment is selected, run `import sys; print(sys.executable)` in a notebook cell; the path should point to the `mnap` environment.

On subsequent runs, activate `mnap`, change to the application folder, and run `python -m jupyterlab` again. Results are written to the notebook's working directory.

---

# 🚀 Workflow Summary

```text
🔬 HRMS data processing
        │
        ▼
HRMS_data_processing.R
        │
        ▼
📊 Independent mNAP benchmark
        │
        ▼
run_benchmark.ipynb
        │
        ▼
🕸️ Application to GNPS molecular networks
        │
        ▼
mNAP_GNPS_application.ipynb
```

---

# 📖 Citation

If you use this repository, please cite the associated manuscript once it becomes available:

> **Su, W. et al. Identification and Prioritization of Plastic Additives in Chinese Agricultural Soils Using High-Resolution Mass Spectrometry and Network Annotation Propagation.**

**DOI:** *To be added after publication.*

### Related methods and software

1. **Helmus, R.; ter Laak, T. L.; van Wezel, A. P.; de Voogt, P.; Schymanski, E. L.**  
   patRoon: open source software platform for environmental mass spectrometry based non-target screening. *J. Cheminform.* **2021**, *13*, 1.  
   DOI: `10.1186/s13321-020-00477-w`

2. **Ruttkies, C.; Schymanski, E. L.; Wolf, S.; Hollender, J.; Neumann, S.**  
   MetFrag relaunched: incorporating strategies beyond in silico fragmentation. *J. Cheminform.* **2016**, *8*, 3.  
   DOI: `10.1186/s13321-016-0115-9`

3. **da Silva, R. R. et al.**  
   Propagating annotations of molecular networks using in silico fragmentation. *PLoS Comput. Biol.* **2018**, *14*, e1006089.  
   DOI: `10.1371/journal.pcbi.1006089`

4. **Borelli, T. C. et al.**  
   Improving annotation propagation on molecular networks through random walks: introducing ChemWalker. *Bioinformatics* **2023**, *39*, btad078.  
   DOI: `10.1093/bioinformatics/btad078`

5. **Nothias, L.-F. et al.**  
   Feature-based molecular networking in the GNPS analysis environment. *Nat. Methods* **2020**, *17*, 905–908.  
   DOI: `10.1038/s41592-020-0933-6`

---

## 📌 Notes

- Local file paths should be adjusted according to the user's environment.
- Large raw MS files, PubChemLite databases, NIST spectra, and third-party software may need to be obtained separately.
- The controlled benchmark and experimental application use different candidate-retention strategies and should be treated as separate workflows.
- Parameter changes should be documented when adapting the workflow to other datasets.
