# VOE: A python code for identifying variant located on epitope of SARS-CoV-2

# Installing VOE (Linux)
* We recommended running all steps on Linux. However, you can perform some steps on Windows.
## Step 1: Installing Miniconda 3

This installation requires the 64-bit operating system and install the miniconda via the terminal. Simply type

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Please answer "yes" for the question "Do you wish the installer to initialize Miniconda3 by running conda init? [yes|no]".
Then relogin.
For more information of conda installation, please see https://docs.conda.io/en/latest/miniconda.html#linux-installers

## Step 2: Installing Pandas

Install python and pip before installing Pandas via conda for linux version.

```bash
conda install -c anaconda python=3.6
conda install -c anaconda pip
conda activate base
conda install pandas
```
Install python and pip before installing Pandas via pip for window version.
```bash
pip install pandas
```

## Step 3: Installing Blast

Install python and pip before installing Pandas via conda.

```bash
conda install -c bioconda blast
```

## Step 4: Preparing your Variant database

```bash
Example Variant database
BCFtool_variant_database.ann.vcf
Freebayes_variant_database.ann.vcf
```

## Step 5: Downloading workflow, software, and prepare working directory

Download the VOE and an example of input file. 
```bash
git clone https://github.com/lee99dn/SARS-CoV-2-VOE.git
cd SARS-CoV-2-VOE

```

Then the VOE is ready to use!

# For VOE Windows installation

Do all steps in windows version.

# Running VOE

Put an input file in VCF format to your working directory.
The VOE is developed based on the python. Therefore, to run VOE, simply type

Windows version
```bash
path/python path/VOE_window.py [-e, --epitope <epitope sequence> -v, --variantdatabase <variantdatabase> -b, --blastdatabase <blastdatabase> -o, --output <output TSV format>]
```
Linux version
```bash
path/./VOE_linux.py [-e, --epitope <epitope sequence> -v, --variantdatabase <variantdatabase> -b, --blastdatabase <blastdatabase> -o, --output <output TSV format>]
```
For example
```bash
path/./VOE_linux.py -e KLNDLCFTNV -v BCFtool_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_1_output.tsv
```
```bash
path/python VOE_window.py -e KLNDLCFTNV -v BCFtool_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_1_output.tsv
```

Output:

Sensitivity:100 %

https://github.com/lee99dn/SARS-CoV-2-VOE/blob/main/Ex_1_output.tsv

```bash
./VOE_linux.py -e APGQTGK -v Freebayes_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_2_output.tsv
```
```bash
python VOE_window.py -e APGQTGK -v Freebayes_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_2_output.tsv
```

Output:

Sensitivity:6.8249 %

https://github.com/lee99dn/SARS-CoV-2-VOE/blob/main/Ex_2_output.tsv

Then the workflow will process automatically.
The output will be show on the screen and create in TSV format.

# License

VOE was be released under GPL license and be widely available on GitHub: https://github.com/lee99dn/SARS-CoV-2-VOE.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11197886.svg)](https://doi.org/10.5281/zenodo.11197886)
