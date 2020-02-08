# YATAK
YATAK is a tool for detecting somatic structural variants from a pair of whole genome sequencing datasets of a germline genome and a somatic genome. The main feature of YATAK is that it directly compares germline reads and somatic reads to call structural variants.
In the first step of YATAK, it a creates a succinct index for the assembly graph and store all paths that correspond to important germline reads together with a reference genome after filtering out the reads that are well aligned to the reference genome. In the next step, it conducts string-to-graph search for all the remaining somatic reads to find the reads that potentially support somatic breakpoints. Finally, such reads are assembled into contigs to detect somatic structural variants by aligning the contigs to the reference genome.

<img src="https://github.com/cBioLab/YATAK/blob/master/figure2.png"
width=800/>

# Installation requirements
* [Succinct Data Structure Library 2.0](https://github.com/simongog/sdsl-lite)

# Installation

      clone git@github.com:cBioLab/YATAK.git
      cd YATAK
      mkdir bin
      make

* You need to edit Makefile to set correct paths for "SDSLLIBS".
# Running

      pipeline.sh
      
# Examples
