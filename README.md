# scRNA-seq Cell-Type Identification via Probabilistic Barcodes

Single-cell RNA-sequencing cell-type identification that leverages reference data to combine information across thousands of genes, learn probabilistic barcode representations of cell-types, and probabilistically assign cell-type identity to unknown cells.

# System Requirements and Installation

Our code requires the ```sads``` package, version ```0.4.2```. We tested our code in ```R``` version ```3.6.2``` and ```4.0.0```. We are currently working on a package. In the meantime, our code can be used by downloading the files ```cell_type_identification.R``` and ```params_EM_81020.rda``` from this repository into the same directory, and calling ```source('cell_type_identification.R')```. This should only take a few seconds. 

# Usage

To learn a barcode representation for cell-types in reference data, given an annotated dataset in numeric, genes by cells format with genes represented by Ensembl IDs:

```
fit <- trainAllReference(ref,labels)
barcodes <- getBarcode(fit)
```

To identify unknown cells in a numeric, genes by cells format from a set of possible reference cell-types with genes represented by Ensembl IDs:

```
target_labels <- classifyTarget(target,fit)
```

We currently support human, UMI count data. 

# Demo

The Rmarkdown file ```PBMCs.Rmd``` provides links to example data from 10X Genomics and demonstrates both training and classification. This file assumes that the filtered genes-by-cell count data have been downloaded from the provided links. The directories in the ```Read10X``` commands in the Rmarkdown file may need to be updated depending on where and with what name the data files were saved. Running this file should take approximately thirty minutes or so. The expected output is shown in ```PBMCs.pdf```, which resulted from running this Rmarkdown file in ```R``` version ```4.0.0```.
