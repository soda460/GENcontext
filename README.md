<p align="center"><img src="misc/logo_github_2.png" alt="Exploring gene context" width="1000"></p>

GENcontext is a tool designed to discover the genetic organization around specific loci (e.g. antibiotic resistance genes, IS elements) in annotation files in the Genbank format. You need to define your loci of interest and the requested number of adjacent features. The tool outputs text files and schematic representations of the genetic organization around the loci of interest.

In the current version of GENcontext,you have to follow [a specific procedure](https://github.com/soda460/prokka_tag_ARGs_with--proteins) to annotate ARGs in your assemblies files. Briefly, the procedure requires to use [PROKKA](https://github.com/tseemann/prokka) with the [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/ "The Comprehensive Antibiotic Resistance Database").


## Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [The class GeneCluster.py](#The-class-GeneClusterpy)
* [The class Dock.py](#The-class-Dockpy)
* [Example commands (quick)](#example-commands-quick)
* [Contributors](#contributors)




## Requirements

* Linux or macOS
* Biopython
* GenomeDiagram (For graphics)



## Installation

GENcontext is a python stand-alone executable:

A good idea to test the program would be to clone this repository.

```shell
git clone https://github.com/soda460/GENcontext.git
```

and install the dependencies in a conda environment

```shell
conda create -n GeneContextEnv -c bioconda biopython reportlab
```



## The class GeneCluster.py

The idea behind this class was to represent a gene cluster object that can be easily manipulated. Gene cluster objects were designed to be able to growth (e.g. when parsing annotation file) , to be reversed, or to be compared with other gene cluster objects for identity or inclusion. Here we present a quick overview of the class

#### Initializing a geneCluster object

```python
a0 = geneCluster("a0")
c1 = geneCluster("c1")
```


#### Reading a gene cluster

Use the read method to read the content of a gene cluster. This method add the proper 5'- 3'- stuff around the genes.

```python
# copy the object a0 and rename the new object a1 accordingly
a1 = a0
a1.name = 'a1_is_a_copy_of_a0'
print(a1.name + '\t' + a1.read())
```

#### Defining the gene content (and strand orientation) of a gene cluster

Gene clusters are basically defined by their gene content and their relative orientation. The  load method take two lists as arguments to load genes from scratch. Alternatively, genes can be load one at a time with the add method.

```python
geneList = ['A', 'B', 'C']
strandList = ['+','+','+']
a0.load(geneList, strandList)
```

#### Reversing a gene cluster
GeneCluster objects can be reversed with the reversed method. 
```python
r1 = a0.reverse()
print('r1: ' + r1.name + '\t' + r1.read())
```






#### Adding one gene to a gene cluster

To add a new gene in a geneCluster object, simply use the add method. This is useful when parsing annotation files. The name and the strand orientation of the new gene are required.

```python
c1 = a0
c1.name='c1'
c1.add('D', '+')
c1.add('E', '+')
print(c1.name + '\t' + c1.read() + ' after the addition')
```



#### Testing the equality between two gene clusters

With gene clusters objects, you can use the overloaded operators == or !=  to check if two gene clusters are identical. For example, 5'-A-B-C-3' will be equivalent to 3'-C-B-A-5'.


```python
print ('--Testing operator overloading--')
	if a1 == b1:
		print('c1 equal a1')
	else:
		print('c1 not equal a1')
	if a0 != a1:
		print('a0 not equal a1')
	else:
		print('a0 equal a1')
```



#### The isin method - Testing the inclusion of a gene cluster into another gene cluster

The isin method is usefull to known whether a gene cluster is entirely included in another gene cluster. Obviously, testing the inclusion of a1 into r1 is different from testing the inclusion of r1 into a1.

```python
print ('--Testing isin method--')
# a1 vs r1
if a1.isin(r1):
	print ('yes a1 is included into r1')
else:
	print ('no a1 is not included into r1')

# d1 vs a1
if d1.isin(a1):
	print ('yes d1 is included into in a1')
else:
	print ('no d1 is not included into a1')
```

#### Ignoring hypothetical genes

Sometimes, hypothetical genes (labelled with ?) can blur the genetic context around a specific genes. The pretty_read method of this class ignore these genes facilitating the creation of altered geneCluster objects devoid of hypothical genes.

```python
u0 = geneCluster("u0")
geneList = ['C', '?', 'A']
strandList = ['-','-','-']
u0.load(geneList, strandList)
print('read method: ' + u0.name + '\t' + u0.read())
print('pretty_read method: ' + u0.name + '\t' + u0.pretty_read())
```

#### Other methods

There a few other methods in this class (remove, get_size and reverse_read). One can look at geneCluster.py for more details.




## The class Dock.py 

This dock class has been defined to regroup similar geneCluster objects.

Each dock object is defined by :

* a numeric id
* a name
* elements, i.e. a list of similar geneCluster objects
* its size
* an HEAD, i.e. an arbitrary geneCluster having the greatest number of genes among all elements of the dock object

 By similar geneCluster objects, understand gene clusters (devoid of hypothetical genes) that can be entirely included in other same-size or larger elements of the same dock object.
 
For example, after a classification, one will want that 5'-A-B-?-?-C-3', 5'-A-B-C-?-D-3' and 5'-A-B-C-?-D-E-3' be part to the same dock object because, when excluding hypothetical genes, the first cluster is included in the two larger elements of the dock and because the second one is included in the third one.

If interested, have a look at geneCluster.py for more details.






#### Setting up your annotation files when input is the MOB-suite output



However, to benefit from good metadata labeling in output files, the best is to have your genbank files organized in subfolders. At first level, create folders with strain_names and at second level, create folders with molecule names in which your genbank files will be located.

```shell
tree ./annotation folder/ -P *.gbk


├── strain1
│   └── a_molecule
│       └── your_genbank_file_here.gbk
├── strain2
│   └── a molecule
│       └── another_genbank_file.gbk
├── PC29
│   └── plasmid_476
│       └── PC-29-I_plasmid_476.gbk
├── Res13-Lact-PER04-34
│   ├── chromosome
│   │   └── Res13-Lact-PER04-34_chromosome.gbk
│   ├── plasmid_1009
│   │   └── Res13-Lact-PER04-34_plasmid_1009.gbk
│   ├── plasmid_1068
│   │   └── Res13-Lact-PER04-34_plasmid_1068.gbk
etc
```


## Example commands (quick)

```
conda activate GeneContextEnv
```

```python
cd scripts/
./expl_gen_context.py -i mob-suite -t 'CTX-M-1' -c 'card' -p '../gbk/mini_prokka' -n 6 -e chromosome plasmid_973

# For IS
./expl_gen_context.py -i mob-suite -t 'IS26' -c 'IS' -p '../gbk/mini_prokka/' -n 10

# To get a infos, notably the distance in bp,  about the other ARG presents
./general_context.py -t 'CTX-M-1' -c 'card' -p '../gbk/mini_prokka/'



```

The optional 'i' argument refer to the kind of input given to the programs. When specifying 'mob-suite' the program assume that indivual gbk files exist for each molecule (for each plasmid and for the chromomosme) and are organized in subfolders. This is important because the program will get the strain and molecule from the PATH of the gbk files. Otherwise, the default value is 'plain' and in that case, the strain and molecule name are obtained directly from specific fields in the genbank file.


The 't' argument refer to the name of the targeted antibiotic resistance gene, for example, CTX-M-1 as labeled in the CARD antibiotic resistance gene database.

Note that you can search for multiple alleles of an ARG. If you type CTX, the program will find the context around all CTX resistance genes.

The 'p' argument is a single path where are located annotation files. The program will retrieved all genbank files nested in this folder.

The 'n' argument is the number of genes (CDS features actually) on both sides of the targeted genes that will be considered by the program.


## Remote repositories

Exact copies of this project exist at two locations :

old-origin	https://github.com/soda460/GENcontext (fetch)  
old-origin	https://github.com/soda460/GENcontext (push)  
origin	https://gccode.ssc-spc.gc.ca/ac_/dpl/genecontext.git (fetch)  
origin	https://gccode.ssc-spc.gc.ca/ac_/dpl/genecontext.git (push)  




## Contributors

  * Jean-Simon Brouard, Ph.D.  
Biologist in bioinformatics, Science and Technology Branch  
Agriculture and Agri-Food Canada / Government of Canada  
jean-simon.brouard@canada.ca


  * Dominic Poulin-Laprade, Ph.D.  
Research Scientist, Science and Technology Branch  
Agriculture and Agri-Food Canada / Government of Canada  
Dominic.Poulin-Laprade@Canada.ca


