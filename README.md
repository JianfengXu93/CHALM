# Cell Heterogeneity Accounted cLonal Methylation (CHALM)

Different from calculating the traditional mean methylation level of a predefined region (e.g., promoter CpG island), CHALM directly uses the aligned sequencing reads as input and quantifies the cell heterogeneity accounted clonal methylation level of this region, which are more powerful in predicting gene expression. CHALM can also calculate cell heterogeneity accounted clonal methylation ratio for single CpG sites which are often required for DMR/UMR detection analysis. Furthermore, in order to illustrate the importance of clonal information, CHALM provides a CNN deep learning framework to predict expression directly by aligned sequencing reads produced by high throughput methylation profiling technologies like whole genome bisulfite sequencing (WGBS). 
## Authors
- Jianfeng Xu (jianfenx@bcm.edu)
- Jianzhong Su (sujz@wibe.ac.cn)
- Wei Li (wl1@bcm.edu)
## Availability
All documents can be downloaded from GitHub link: [https://github.com/JR0202/CHALM][1]
## Dependencies
- samtools/0.1.19
- anaconda/2.5.0 

_note: CHALM depends on above packages, but some tools included in CHALM may have different dependencies._
## Installation
No installations are needed. Simply run the python scripts as:  `python <python scripts>`

CHALM is tested on python 2.7.
## Code Example (download the example data by this link: http://dldcc-web.brc.bcm.edu/lilab/JX/CHALM\_example.tar.gz)
### 1. Calculate the traditional methylation level for promoter CGIs
#### a. CpG meth ratio
```bash
python ../src/CHALM.py trad -d hg19.fa -x CG -i no-action -p -r -m 4 -o output_examples/CD3_primary_chr13_trad_CpG_methratio.txt CD3_primary_CGI_chr13.sam
# time cost: ~15 min
```
#### b. Traditional mean methylation
```bash
python ../src/Trad_meth_mean.py -r CGI_Gene_match_human_sorted_header.txt -m output_examples/CD3_primary_chr13_trad_CpG_methratio.txt -o output_examples/CD3_primary_chr13_trad_meth_mean_promoter_CGI.txt
# time cost: ~3 min
```
### 2. Quantify the heterogeneity-accounted methylation ratio of CpG sites (for downstream DMR/UMR detection)
```bash
python ../src/CHALM.py trad -d hg19.fa -x CG -i no-action -p -r -m 4 -l 1 -o output_examples/CD3_primary_chr13_CHALM_CpG_methratio.txt CD3_primary_CGI_chr13.sam
# time cost: ~15 min
```
### 3. Calculate the CHALM methylation level for predefined region (e.g., promoter CGIs)
```bash
python ../src/CHALM.py CHALM -d hg19.fa -x CG -R Human_CGI_bedfile.txt -L 99 -l 1 -p -r -o output_examples/CD3_primary_chr13_CHALM.txt CD3_primary_CGI_chr13.sam
# time cost: ~12 min
```
### 4. Calculate differential CHALM methylation level
```bash
python ../src/CHALM_dif.py -f1 output_examples/CD3_primary_chr13_CHALM.txt -f2 output_examples/CD14_primary_chr13_CHALM.txt -o output_examples/CD3_CD14_chr13_CHALM_dif.txt
# time cost: ~10 sec
```
_note: for replicates, separate the files by "," (e.g., Condition\_1\_replicate1.txt,Condition\_1\_replicate2.txt,Condition\_1\_replicate3.txt)_
### 5. Imputation to extend the sequencing read length
#### a. Dependencies
- anaconda2/4.3.1
- R/3.3.0 
- samtools/0.1.19
#### b. Command example
```bash
python ../src/CHALM_SVD_imputation.py -d hg19.fa -e 100 -x CG -R Human_CGI_bedfile.txt -L 99 -l 1 -p -r -o output_examples/CD3_primary_chr13_CHALM_extend_100.txt CD3_primary_CGI_chr13.sam
# time cost: ~11 min
```
### 6. Deep learning prediction of expression by CHALM
#### a. Dependencies
- samtools/0.1.19
- anaconda/2.5.0 
- bedtools/2.17.0
- pytorch (install from [http://pytorch.org/][2])
#### b. Command example
##### (1) Process aligned reads for deep learning
```bash
python ../src/Deep_learning_read_process.py -d hg19.fa -x CG -p -r -o output_examples -n CD3_primary --region Gene_CGI_match_TSS_sorted.txt --depth_cut 50 --read_bins 200 CD3_primary_CGI.sam
# time cost: ~15min
```
As control, add '-S' to disrupt the clonal information
```bash
python ../src/Deep_learning_read_process.py -d hg19.fa -x CG -p -r -S -o output_examples -n CD3_primary --region Gene_CGI_match_TSS_sorted.txt --depth_cut 50 --read_bins 200 CD3_primary_CGI.sam
# time cost: ~18min
```
##### (2) Train deep learning model and do expression prediction
```bash
python ../src/Deep_learning_prediction.py -f1 output_examples/CD3_primary_meth_2D_code.txt -f2 output_examples/CD3_primary_distance_2_TSS.txt -m output_examples/CD3_primary_trad_meth_mean_promoter_CGI.txt -e CD3_primary_RSEM.genes.results -s CD3_primary -d -o output_examples/
# time cost: ~5min
```
![][image-1]

Train the control data (with disrupted clonal information)
```bash
python ../src/Deep_learning_prediction.py -f1 output_examples/CD3_primary_meth_2D_code_control.txt -f2 output_examples/CD3_primary_distance_2_TSS_control.txt -m output_examples/CD3_primary_trad_meth_mean_promoter_CGI.txt -e CD3_primary_RSEM.genes.results -s CD3_primary_control -d -o output_examples/
# time cost: ~5min
```
![][image-2]

_note: CD3\_primary\_RSEM.genes.results contains the expression level calculated by RSEM (rsem-calculate-expression)_
### 7. Deep learning prediction of expression by CHALM (pre-trained model)
#### a. Command example
##### (1) Expression prediction of CD3 primary cell by pre-trained model
```bash
python ../src/Deep_learning_prediction_pretrained.py -f1 output_examples/CD3_primary_meth_2D_code.txt -f2 output_examples/CD3_primary_distance_2_TSS.txt -m output_examples/CD3_primary_trad_meth_mean_promoter_CGI.txt -e CD3_primary_RSEM.genes.results -s CD3_primary_pretrained --model pretrained_model.pt -d -o output_examples/ 
# time cost: ~5min
```
![][image-3]

Train the control data (with disrupted clonal information)
```bash
python ../src/Deep_learning_prediction_pretrained.py -f1 output_examples/CD3_primary_meth_2D_code_control.txt -f2 output_examples/CD3_primary_distance_2_TSS_control.txt -m output_examples/CD3_primary_trad_meth_mean_promoter_CGI.txt -e CD3_primary_RSEM.genes.results -s CD3_primary_pretrained_control --model pretrained_model.pt -d -o output_examples/
# time cost: ~5min
```
![][image-4]


[1]:	https://github.com/JR0202/CHALM
[2]:	http://pytorch.org/

[image-1]:	https://github.com/JR0202/CHALM/blob/master/example_figures/CD3_primary_deep_learning_prediction.png
[image-2]:	https://github.com/JR0202/CHALM/blob/master/example_figures/CD3_primary_control_deep_learning_prediction.png
[image-3]:	https://github.com/JR0202/CHALM/blob/master/example_figures/CD3_primary_pretrained_deep_learning_prediction.png
[image-4]:	https://github.com/JR0202/CHALM/blob/master/example_figures/CD3_primary_pretrained_control_deep_learning_prediction.png