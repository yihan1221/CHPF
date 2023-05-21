Description
===========
This package performs cellular hypoxia status predicting by integrating single cell transcriptome data and hypoxic gene sets. It will:
* identify high-confidence hypoxic and normoxic cells in your data.
* infer hypoxia status of other cells.

System requirements and dependency
==================================
This package runs on Python 3.7.3.
It also requires R/4.2.3 to run and has dependency on the R packages:
	corrplot and mclust.

Installation
============
Please download and copy the distribution to your specific location. If you are cloning from github, ensure that you have git-lfs installed.

For example, if the downloaded distribuition is CHPF.tar.gz.
	Type 'tar zxvf CHPF.tar.gz'

Then, run CHPF.py in the resulting folder.

Usage
=====
```
Options:
 
  -h, --help            Show this help message and exit.
  -g GMT, --gmt GMT     The input file of gene sets. 
  -e EXPR,--expr EXPR   The input file of expression profile from scRNA-seq.Columns correspond to cells.Rows correspond to genes.                          
  -o OUTPATH, --outpath OUTPATH
                        The output path.
  -n N_SPLITS, --n_splits N_SPLITS
                        The fold of cross-validation.Default value is 5.
  -t TREE, --tree TREE  The number of decision trees bulit in one cross-validation.Default value is 100.
  -T THEO, --THEO THEO  The threshold value of the Px.Default value is 0.5. If the Px > THEO, the predict outcome is 1, meaning the cell will
                        be predicted as a hypoxic cell; otherwise, the predict outcome is 0, meaning the cell will be predicted as a normoxic cell.

```
Input files
===========

  *Gene sets input file:*

  	geneset1	describe1	gene1	gene2	gene3	......
  	geneset2	describe2	gene4	gene5	gene6	......

  *scRNA-seq expression input file:*

	Two types of input files are allowed in CHPF:
		(1).RData file contains only one object.
		(2).txt or .csv file

          cell1  cell2  cell3  ......
    gene1  0      1.5    3.2   ......
    gene2  1.4    0      0     ......

>For gene sets file input,it contains at least two gene sets and the type of file should be .gmt.
>For scRNA-seq expression file input, the value can be count or tpm.


Run CHPF package
==================
    python CHPF.py [-n <n fold cross-validation>] [-t <decision tree number>] [–T <threshold value>] -g <Gene sets> -e <Expression profile> -o <output path> 
    [...] contains optional parameters.
    The mandatory arguments are -g, -e and -o.

Examples
========
Try CHPF in the package directory on the different example datasets

    python CHPF.py -g ./example/Hypoxia_geneset.gmt -e ./example/expr.RData -o ./example/result

Output files
============

    cell_status.csv which includes the cells'id and the corresponding hypoxia status.(1-hypoxia and 0-normoxia)

    *temp output:*
    (1)exp_highconfi.csv which is the expression profile of high-confidence cells in 500 feature genes.
    (2)exp_others.csv which is the expression profile of non high-confidence cells in 500 feature genes.
    (3)label_highconfi.csv which is the cell status of high-confidence cells predicted by GMM.
    (4)label_others.csv which is the cell status of non high-confidence cells predicted by lightgbm.

    Two figures:
    (1)GSVA_corr.pdf which is a visualization of correlation of hypoxia activity scores between any two gene sets.
    (2)GMM_classification.pdf which is a heatmap of cell classification using GMM based on each gene set.

    *feature output:*
     The number of files depends on  the product of the fold of cross-validation and the number of trees built in each cross-validation.
     In each file, it includes the important of features of tree. 
