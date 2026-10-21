# mPHGR: microbiome personalized hyper-graph ranking for identifying critical states of complex biological processes at an individual level
## Overview

A sudden and non-linear transition often occurs just before a system state change in complex biological processes. Detecting such key transitions or critical states is essential for understanding intricate biological mechanisms and for timely preventing or delaying destructive effects, particularly in the identification of critical states during disease progression. In this study, from an individual-based perspective, we propose microbiome personalized hyper-graph ranking (mPHGR), a single-sample method based on high-order networks, to detect critical transitions of complex biological systems using microbiome data. Specifically, the sample-specific microbial hyper-graph is constructed based on microbiome data from both a set of reference samples and an individual case sample. Then, a personalized PageRank model is proposed to calculate the local mPHGR score for each node/species based on a state transition matrix derived from the sample-specific microbial hyper-graph. A sudden increase in the mPHGR score indicates the imminent critical stages or tipping points in complex biological processes, as our proposed method quantifies dynamic changes induced by individual case samples in the higher-order interactions among microbial species.
<img width="4500" height="4249" alt="Fig1" src="https://github.com/user-attachments/assets/bb28fb4a-8822-44f5-ac86-157aa413d891" />

## Usage

mPHGR can be used directly in Python or through an R package. Both approaches use the same computational algorithm. 

### 1. Python implementation

Download `mPHGR.py` from this repository and install the required Python packages:

```bash
pip install numpy pandas matplotlib
```

Import the `run_mphgr()` function to calculate mPHGR scores:

```python
from mPHGR import run_mphgr
```

The function requires the input data directory (`data_dir`), analysis mode (`data_type`), parameters (`d` and `Q`), and input file prefix (`prefix_data`). For stage-based analysis, the `stages` argument is additionally required.

### 2. R package

The mPHGR R package (v1.0.0) provides the `mphgr_run()` function, which calls the Python implementation through `reticulate`.

Install the R package directly from GitHub:

```r
install.packages(c("remotes", "reticulate"))

remotes::install_github(
  "WeQiao/mPHGR",
  subdir = "R_package"
)

library(mPHGR)
```

A Python environment with `numpy`, `pandas`, and `matplotlib` is also required.

For detailed installation instructions, Python environment configuration, and package documentation, please refer to the [R package README](R_package/README.md).

## Examples

We use an individual (10013) from the preterm birth (PTB) dataset and a gastric cancer (GC) dataset to demonstrate the application of mPHGR. Users can also apply the method to their own microbiome datasets using the same workflow.

### 1. Example datasets

**Gastric cancer (GC)**

The GC dataset contains gastric mucosal microbiome samples collected from 47 participants, divided into four groups according to disease progression: superficial gastritis (SG, 17 participants), atrophic gastritis with intestinal metaplasia (AG, 10 participants), gastric intraepithelial neoplasia (GIN, 5 participants), and intestinal-type gastric cancer (GC, 15 participants).

The first 10 subjects from the SG group were selected as the control group. The raw sequencing data are available under BioProject accession [PRJNA634837](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA634837).

**Preterm birth (PTB)**

The PTB dataset was obtained from a longitudinal study of 40 pregnant women, 11 of whom experienced preterm delivery. Vaginal, stool, saliva, and tooth/gum samples were self-collected weekly during pregnancy and monthly after delivery for up to 12 months.

Since the original study indicated that vaginal microbiota were more closely associated with preterm birth risk, we analyzed vaginal microbiome data. The first three time points from the remaining 29 women were used as the control group. The raw sequencing data are available in the Sequence Read Archive under accession SRP288562.

### 2. Data preprocessing

Convert the microbial abundance data into tab-separated TXT files, with microbial species in rows and samples in columns. Normalize the data by column before applying mPHGR.

As a single-sample method, mPHGR constructs an individual test dataset by adding one case sample column to the reference samples.

The example input data are provided in the `GC` and `PTB10013` archives in this repository. After extraction, each data directory should contain `control.txt` and the corresponding test files.

### 3. Run mPHGR

#### Example 1: Preterm birth (PTB)

For longitudinal microbiome data, use `data_type="1"`. The mPHGR scores are normalized to the score of the first time point.

**Python:**

```python
from mPHGR import run_mphgr

run_mphgr(
    data_dir="C:/mPHGR/PTB10013",
    data_type="1",
    d=30,
    Q=0.5,
    prefix_data="test2_data"
)
```

**R:**

```r
library(mPHGR)

mphgr_run(
  data_dir = "C:/mPHGR/PTB10013",
  data_type = "1",
  d = 30,
  Q = 0.5,
  prefix_data = "test2_data"
)
```

#### Example 2: Gastric cancer (GC)

For stage-based microbiome data, use `data_type="2"`. 

**Python:**

```python
from mPHGR import run_mphgr

run_mphgr(
    data_dir="C:/mPHGR/GC",
    data_type="2",
    d=145,
    Q=0.7,
    prefix_data="test3_data",
    stages={
        "SG": list(range(1, 8)),
        "AG": list(range(8, 18)),
        "GIN": list(range(18, 23)),
        "GC": list(range(23, 38))
    }
)
```

**R:**

```r
library(mPHGR)

mphgr_run(
  data_dir = "C:/mPHGR/GC",
  data_type = "2",
  d = 145,
  Q = 0.7,
  prefix_data = "test3_data",
  stages = list(
    SG = 1:7,
    AG = 8:17,
    GIN = 18:22,
    GC = 23:37
  )
)
```
For both implementations, `d` and `Q` should be specified according to the dataset. The `prefix_data` argument determines the input file naming pattern. For `data_type="2"`, `stages` specifies the file IDs belonging to each stage, and the first specified stage serves as the normalization reference.

The paths shown above are examples and should be replaced with the actual local data directories.

### 4. Outputs

Both the Python implementation and the R package save two files in the specified input directory:

- `Score.csv`: Normalized mPHGR scores.
- `d.csv`: Node-selection vector.


## Citation
Qiao Wei, Jiayuan Zhong, et al. mPHGR: microbiome personalized hyper-graph ranking for identifying critical states of complex biological processes
## Contact
Qiao Wei: scut_wq@163.com

Jiayuan Zhong: Zjiayuan@foshan.edu.cn
