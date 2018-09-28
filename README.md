## NOTE: The latest version of ANGSD-wrapper is at [ANGSD-wrapper/angsd-wrapper](https://github.com/ANGSD-wrapper/angsd-wrapper)

The version in this repository is out of date, please go to the version linked above.

---

# ANGSD-wrapper

<!--
> **NOTE:** Due to changes to Apple's C/C++ compiler, ANGSD-wrapper is unavailable for Mac OS X Version 10.11.4
> This affects setting up dependencies, creating a system-wide installation of ANGSD-wrapper, downloading example data, and performing any computational analyses with ANGSD-wrapper
> The Shiny graphing interface is still available through `./angsd-wrapper shiny graphing`
> Please use an older version of Mac OS X or another POSIX-compliant system to regain full capabilities of ANGSD-wrapper
-->

ANGSD-wrapper is a utility developed to aid in the analysis of next generation sequencing data. Users can do the following with this suite:
- Calculate a site frequency spectrum
- Calculate a 2D site frequency spectrum with corresponding *F*<sub>ST</sub> estimations
- Perform ABBA/BABA tests
- Extract a FASTA sequence from BAM files
- Calculate genotype likelihoods
- Estimate Thetas and various neutrality statistics
- Calculate per-individual inbreeding coefficient
- Find admixture proportions

Likelihood based approaches are used in ANGSD to calculate summary statistics from next generation sequencing data. The wrapper scripts and documentation are designed to make ANGSD user-friendly.

## Installing ANGSD-wrapper

To install ANGSD-wrapper, download from GitHub

```shell
git clone https://github.com/mojaveazure/angsd-wrapper.git
```

Go into the ANGSD-wrapper directory

```shell
cd angsd-wrapper/
```

Run the setup command

```shell
./angsd-wrapper setup dependencies
```

Download the example dataset (optional)

```shell
./angsd-wrapper setup data
```

Finish the installation

```shell
source ~/.bash_profile
```

## A note about BAM files

ANGSD requires BAM files as input, and ANGSD-wrapper passes a list of BAM files to ANGSD. These BAM files have a few requirements:

- The BAM files must have an '@HD' header line
- The BAM files must be indexed (.bai)

To see whether or not the BAM files have an '@HD' header line, run the following on your list of samples:
```shell
for sample in `cat ~/path/to/sample_list.txt`
do
    echo $sample
    samtools view -H $sample | head -1
done
```

If any samples start with '@SQ' instead of '@HD', ANGSD and ANGSD-wrapper will fail. This [Gist](https://gist.github.com/mojaveazure/d194c4705642eecf8437) will add an `@HD` header lines to your BAM files.

The index files must be generated after the BAM files. To index the BAM files using SAMTools, run the following on your sample list:

```shell
for sample in `cat ~/path/to/sample_list.txt`
do
    samtools index $sample
done
```

If you have GNU Parallel installed on your system, this process can be sped up:

```shell
cat ~/path/to/sample_list.txt | parallel samtools index {}
```

## Basic usage

To run ANGSD-wrapper, run

```shell
angsd-wrapper <wrapper> <config>
```

Where `wrapper` is one of the methods that ANGSD-wrapper can run and `config` is the relative path to the corresponding configuration file.

To see a list of available wrappers, run

```shell
angsd-wrapper
```

## Configuration files

There is a configuration (config) file for each method available through `angsd-wrapper.` The configuration files hold variables used by the wrappers. This is where you need to modify and save the variables (i.e., specify filepaths of indexed BAM files/CRAM files,  FASTA files, sample lists, etc.) to suit your samples before running angsd-wrapper with a specified method.

The default config files can be found in the `Configuration_Files` directory. You will need to modify them to suit your samples. Please refer to the config files or the [wiki](https://github.com/mojaveazure/angsd-wrapper/wiki) to see what each variable is used for and how they should be specified. If you run `angsd-wrapper` without any arguments, it will return a usage message.

Example config files can be found in `Example_Data/Configuration_Files` upon running `angsd-wrapper setup data`.

## Futher Information

For more information about ANGSD-wrapper, the methods availble through ANGSD-wrapper, and a comprehensive tutorial, please see the [wiki](https://github.com/mojaveazure/angsd-wrapper/wiki).

## Dependencies
This package requires the following dependencies:
 - [ANGSD](https://github.com/angsd/angsd)
 - [ngsPopGen](https://github.com/mfumagalli/ngsPopGen)
 - [ngsF](https://github.com/fgvieira/ngsF)
 - [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix)

These are downloaded and installed automatically when angsd-wrapper is [installed](https://github.com/mojaveazure/angsd-wrapper#installing-angsd-wrapper)

There are a few other dependencies that are **not** automatically downloaded during the installation:
 - [SAMTools](http://samtools.github.io/)
 - [GNU Scientific Library](http://www.gnu.org/software/gsl/)
 - [Git](http://www.git-scm.com/)
 - [Wget](http://www.gnu.org/software/wget/)

## Supported methods

 - [Site frequency spectrum (SFS)](https://github.com/mojaveazure/angsd-wrapper/wiki/Site-Frequency-Spectrum)
 - [Thetas estimations](https://github.com/mojaveazure/angsd-wrapper/wiki/Thetas)
 - [2D SFS and *F*<sub>ST</sub>](https://github.com/mojaveazure/angsd-wrapper/wiki/2D-Site-Frequency-Spectrum-and-Fst)
 - [ABBA/BABA](https://github.com/mojaveazure/angsd-wrapper/wiki/Abbababa)
 - [Ancestral sequence extractions](https://github.com/mojaveazure/angsd-wrapper/wiki/Ancestral-Sequence)
 - [Genotype likelihood estimations](https://github.com/mojaveazure/angsd-wrapper/wiki/Genotype-Likelihoods)
 - [Inbreeding coefficients calculations](https://github.com/mojaveazure/angsd-wrapper/wiki/Inbreeding-Coefficients)
 - [Principal component analysis](https://github.com/arundurvasula/angsd-wrapper/wiki/Principle-Components-Analysis)
 - [Admixture analysis](https://github.com/mojaveazure/angsd-wrapper/wiki/Admixture-Analysis)
 
## Citing ANGSD-wrapper
 
ANGSD-wrapper was published in [Molecular Ecology Resources](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12578/abstract); if you use this in your work please cite the paper. For BibTeX users, the citation is as follows:

```latex
@article {MEN:MEN12578,
author = {Durvasula, Arun and Hoffman, Paul J. and Kent, Tyler V. and Liu, Chaochih and Kono, Thomas J. Y. and Morrell, Peter L. and Ross-Ibarra, Jeffrey},
title = {angsd-wrapper: utilities for analysing next-generation sequencing data},
journal = {Molecular Ecology Resources},
issn = {1755-0998},
url = {http://dx.doi.org/10.1111/1755-0998.12578},
doi = {10.1111/1755-0998.12578},
pages = {n/a--n/a},
keywords = {domestication, population genetics, software, visualization, Zea},
year = {2016},
}
```
