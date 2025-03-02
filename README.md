# TS-PBWT
Contact Author: Yuan Wei (yuan.wei@ucf.edu)

## Description
TS-PBWT is a program designed to infer identity by descent (IBD) segments using a recombination-based definition from ancestral recombination graphs (ARGs) in a [tree sequence](https://tskit.dev/tutorials/what_is.html) format.

## Run TS-PBWT Program
TS-PBWT program uses [Meson](https://mesonbuild.com/) as its build system. By default, Meson build system uses [Ninja](https://ninja-build.org/) in the backend to build C/C++ applications. Below is an example command of building the TS-PBWT program using Meson:
```
git clone https://github.com/ucfcbb/TS-PBWT.git
cd TS-PBWT
git submodule update --init
meson setup build
meson configure --optimization 3 build
cd build
meson compile
```

TS-PBWT program has the following parameters:
- t, or trees `<INPUT TREE SEQUENCE FILE>`, where `<INPUT TREE SEQUENCE FILE>` is the input tree sequence path and file name (required).
- m, or map `<INPUT GENETIC MAP FILE>`, where `<INPUT GENETIC MAP FILE>` is the input genetic map path and file name (optional; required if unit of cutoff length is set to genetic).
- u, or unit `<UNIT OF CUTOFF LENGTH>`, where `<UNIT OF CUTOFF LENGTH>` is the unit of cutoff length of an IBD (optional; 0 if physical distance is used; 1 if genetic distance is used (input genetic map is required); default is 0).
- l, or length `<MINIMUM CUTOFF LENGTH>`, where `<MINIMUM CUTOFF LENGTH>` is the minimum cutoff length of an IBD (optional; default is 0).
- o, or output `<OUTPUT IBD FILE>`, where `<OUTPUT IBD FILE>` is the output directory path for all files (optional; default is "./output.ibd").

An example command of running the TS-PBWT program using physical distance 1 megabases (Mb) as the minimum cutoff length of an IBD:
```
./build/TSPBWT_v0.1 -t ./test/sample_10.ts -l 1000000
```

Another example command of running the TS-PBWT program using genetic distance 1 centiMorgan (cM) as the minimum cutoff length of an IBD:
```
./build/TSPBWT_v0.1 -t ./test/sample_10.ts -u 1 -l 1.0 -m ./test/sample_10.map -o ./output_genetic_1.0.ibd
```

The command to get the help of the program:
```
./build/TSPBWT_v0.1 -h
```

## Input and Output Files
The required input file for TS-PBWT program is a tree sequence file format (*.ts or *.trees). It is supported by the [tskit](https://github.com/tskit-dev/tskit) library and can be simulated using [msprime](https://github.com/tskit-dev/msprime) or [SLiM](https://github.com/MesserLab/SLiM). The tree sequence file can also be inferred from DNA sequence data using an ARG inference method such as [tsinfer](https://github.com/tskit-dev/tsinfer). The optional genetic map file uses the HapMap format, whose description can be found in the first line of the file. The format of each line starting with the second line contains four tab-delimited fields: *Chromosome*, *Position(bp)*, *Rate(cM/Mb)*, and *Map(cM)*. Note that only *Position(bp)* and *Map(cM)* field values are used for genetic position identification. If the genetic position corresponding to the physical position in the tree sequence file is not found in the genetic map file, interpolation is used to estimate the genetic position of such physical position.

The output file contains the inferred IBD segments of sample individual haplotypes in tree sequence in a tab-delimited text format with four fields: *Sample individual haplotype id 1*, *Sample individual haplotype id 2*, *IBD start position (inclusive)*, *IBD end position (exclusive)*. The *Sample individual haplotype id* is a zero-based indexing of sample individual haplotypes (or leaf nodes) in the order of appearance in tree sequence.
