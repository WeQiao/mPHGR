# mPHGR R package (v1.0.0)

The `mPHGR` R package provides one function, `mphgr_run()`, which calls the **unchanged Python mPHGR algorithm** through `reticulate`.

## Install

After uploading this `R_package/` directory to the existing GitHub repository:

```r
install.packages(c("remotes", "reticulate"))
remotes::install_github("WeQiao/mPHGR", subdir = "R_package", ref = "v1.0.0")
library(mPHGR)
```


Python >= 3.8 with `numpy`, `pandas`, and `matplotlib` is required:

```bash
python -m pip install numpy pandas matplotlib
```

If necessary, select the same Python environment in R **before** running the package:

```r
reticulate::use_python("C:/path/to/python.exe", required = TRUE)
```

## Data format

Put `control.txt` and the individual test files in the same folder. Inputs must already be preprocessed; they are tab-separated matrices with **taxa in rows and samples in columns**. Files are selected using the supplied prefix followed by a numeric ID, e.g. `test2_data1.txt`, `test2_data2.txt`, etc.

## Run

```r
library(mPHGR)

# Mode 1: longitudinal samples, scores normalized to the first sample.
mphgr_run(
  data_dir = "C:/mPHGR/PTB10013", data_type = "1",
  d = 30, Q = 0.5, prefix_data = "test2_data"
)

# Mode 2: disease stages, stage means normalized to the first stage.
mphgr_run(
  data_dir = "C:/mPHGR/GC", data_type = "2",
  d = 145, Q = 0.7, prefix_data = "test3_data",
  stages = list(SG = 1:7, AG = 8:17, GIN = 18:22, GC = 23:37)
)
```

## Outputs

The algorithm saves `Score.csv` and `d.csv` in `data_dir`. Temporary transition matrices are deleted. 


## Python code

The R package includes `inst/python/mPHGR.py`; it is byte-for-byte identical to `mPHGR.py` at the top level of this distribution. No Python computational formulas or execution steps were modified.
