# VOE: A python code for identifying variant located on epitope of SARS-CoV-2

# Installing VOE
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

# Running VOE

Put an input file in VCF format to your working directory.
The VOE is developed based on the python. Therefore, to run VOE, simply type

Window version
```bash
python VOE_window.py [-e, --epitope <epitope sequence> -v, --variantdatabase <variantdatabase> -b, --blastdatabase <blastdatabase> -o, --output <output TSV format>]
```
Linux version
```bash
./VOE_linux.py [-e, --epitope <epitope sequence> -v, --variantdatabase <variantdatabase> -b, --blastdatabase <blastdatabase> -o, --output <output TSV format>]
```
For example
```bash
./VOE_linux.py -e KLNDLCFTNV -v BCFtool_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_1_output.tsv
```
```bash
python VOE_window.py -e KLNDLCFTNV -v BCFtool_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_1_output.tsv
```

[UploadinEpitope sequence: KLNDLCFTNV Gene: S
Sensitivity = 100.00%
Variant database: BCFtool_variant_database.ann.vcf
g Ex_1_output.tsv…]()

```bash
./VOE_linux.py -e APGQTGK -v Freebayes_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_2_output.tsv
```
```bash
python VOE_window.py -e APGQTGK -v Freebayes_variant_database.ann.vcf -b db_cds_nucl_covid.fasta -o Ex_2_output.tsv
```

[Uploading Ex_2_ouEpitope sequence: APGQTGK Gene: S
POS_Genome	Type	ALT	AA_change	Allele_Count(AC)	Sample_Count(NS)	Allele_frequency(AF)	Chance(%)
22802	missense_variant	c.1240C>A	p.Gln414Lys	1.0	1011	0.001	0.0989
22811	missense_variant	c.1251G>T	p.Lys417Asn	1	1011	0.001	0.0989
22811	missense_variant	c.1249_1251delAAGinsCAT	p.Lys417His	1	1011	0.001	0.0989
22811	missense_variant	c.1251G>T	p.Lys417Asn	1	1011	0.001	0.0989
22812	missense_variant	c.1250_1251delAGinsCT	p.Lys417Thr	1.0	1011	0.001	0.0989
22812	missense_variant	c.1250_1251delAGinsTT	p.Lys417Ile	1.0	1011	0.001	0.0989
22813	missense_variant	c.1251G>T	p.Lys417Asn	938.0	1011	0.9278	92.7794
FN:942
TP:69
Sensitivity:6.8249 %
Variant database: Freebayes_variant_database.ann.vcf
tput.tsv…]()


Then the workflow will process automatically.
The output will be show on the screen and create in TSV format.
