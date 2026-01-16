# mPHGR: microbiome personalized hyper-graph ranking for identifying critical states of complex biological processes at an individual level
## Overview

A sudden and non-linear transition often occurs just before a system state change in complex biological processes. Detecting such key transitions or critical states is essential for understanding intricate biological mechanisms and for timely preventing or delaying destructive effects, particularly in the identification of critical states during disease progression. In this study, from an individual-based perspective, we propose microbiome personalized hyper-graph ranking (mPHGR), a single-sample method based on high-order networks, to detect critical transitions of complex biological systems using microbiome data. Specifically, the sample-specific microbial hyper-graph is constructed based on microbiome data from both a set of reference samples and an individual case sample. Then, a personalized PageRank model is proposed to calculate the local mPHGR score for each node/species based on a state transition matrix derived from the sample-specific microbial hyper-graph. A sudden increase in the mPHGR score indicates the imminent critical stages or tipping points in complex biological processes, as our proposed method quantifies dynamic changes induced by individual case samples in the higher-order interactions among microbial species.
<img width="4500" height="4249" alt="Fig1" src="https://github.com/user-attachments/assets/bb28fb4a-8822-44f5-ac86-157aa413d891" />
## Usage
Download the source codes and upzip the data.zip.
## Examples
This project takes an individual 10013 in the preterm births (PTB) dataset and gastric cancer (GC) as examples to illustrate the specific application of the mPHGR method. The input data can be changed with the other datasets if necessary.
The gastric cancer (GC) dataset is a gastric mucosal microbiome dataset collected from 47 participants. In the original study, all subjects were divided into four groups based on the disease progression: 17 patients of superficial gastritis (SG group), 10 patients of atrophic gastritis with intestinal metaplasia (AG group), five patients of gastric intraepithelial neoplasia (GIN group) and 15 patients of intestinal-type gastric cancer (GC group). In this study, selected the first 10 subjects from the SG group as the control group. The Raw sequence data of GC can be found at the following URL: https://www.ncbi.nlm.nih. gov/bioproject/PRJNA634837.
The preterm dataset is derived from a study including 40 pregnant women throughout gestation and 11 of whom experienced preterm delivery. Where specimens from the vagina, stool, saliva and tooth/gum were self-collected by participants weekly from the time of study enrollment until delivery and monthly from the time of delivery for up to 12 mo. Since original study indicates that the vaginal microbiota is more associated with the risk of preterm birth and has greater predictive value, in this study, we analyzed the vaginal microbiota data and used the first three time points from the remaining 29 women as the control group for analysis. The Raw sequence data of PTB have been deposited at the Sequence Read Archive (SRP no. 288562).
### Step1. Data preprocessing
Convert the original microbial abundance matrix data into a txt file. The rows represent species and the columns represent samples. Normalize the data by column as the input file. As a single-sample method, the experimental group data can be obtained by adding one column of samples to be tested on the basis of the control group data.
In the given example, the input data are contained within the folder "PTB10013" or "GC".
### Step2. Calculate mPHGR score and identify the pre-disease stage
Running pipeline: mPHGR_PTB.py (or mPHGR_GC.py), which has been tested in Python 3.8.
Output data: scores_PTB_10013.txt (or score_GC.txt)
## Citation
Qiao Wei, Jiayuan Zhong, et al. mPHGR: microbiome personalized hyper-graph ranking for identifying critical states of complex biological processes at an individual level
## Contact
Qiao Wei: scut_wq@163.com

Jiayuan Zhong: Zjiayuan@foshan.edu.cn
