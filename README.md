<p align="center"><img src="misc/cargo.jpg" alt="Exploring gene context" width="450"></p>

EARGC is a tool to discover the gene organization around specific genes in annotation files. For now, solely the the Genbank format is supported.

It was designed to find the genetic context around antibiotic resistance gene which are annotated very specifically in Genbank files. The main drawback of this tool is therefore that you have to follow a specific procedure to annotate your antibiotic resistance genes files. Briefly, the procedure require to use PROKKA with a custom database (namely the reference CARD database).
The produre is well described here (put link here).




## Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Example commands (quick)](#example-commands-quick)
* [Acknowledgements](#acknowledgements)



## Requirements

* Linux or macOS
* Biopython
* RangeFeature



## Installation

EARGC is a stand-alone executable:


```
git clone https://github.com/soda460/exploring_genetic_context.git
```

## The class GeneCluster.py

The idea behind this class was to represent a gene cluster object that can be easily manipulated. Gene cluster objects were designed to be able to growth (when parsing annotation file) , to be reversed, or to be compared with other gene cluster objects for identity (with operator overloading) or inclusion.


### Quick overview of the class

Initialization of a geneCluster object

```
a0 = geneCluster("a0")

c1 = geneCluster("c1")
```


*  #### load method

This method allow to load genes from scratch. Note that in practice genes are more often load one at a time with the add method when parsing an annotation file. 

```
geneList = ['A', 'B', 'C']
strandList = ['+','+','+']
a0.load(geneList, strandList)
```

* #### add method

This method allow to add a new gene cluster in a geneCluster object.  The geneName and the strand orientation are required.

```
c1.add('E', '+')
```


* #### read method

Used to read the content of a gene cluster. This method add the proper 5'- 3'- stuff around the genes.
```
print(c1.name + '\t' + c1.read())
```

* #### Testing the equality/difference between two gene clusters




* #### The isin method - Testing the equality between two gene clusters






## The dock Dock.py 





## Example commands (quick)



## Acknowledgements

I owe many thanks to DPL and GT from AAFC.

