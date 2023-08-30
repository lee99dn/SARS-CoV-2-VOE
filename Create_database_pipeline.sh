#!/bin/bash

#Read argument -l -r -t
function usage()
{
  echo "Automate_Database Usage: $0 Automate_Database.sh [-h] [ -l <Accession_list.tabular> ] [ -r <Reference_Genome.fasta> ] [ -t <Threads number> ]
			-h Show help text
			-l File or path to SRA accession list in tabular format
			-r File or path to Reference_Genome in fasta format
			-t Threads number"
  exit 1
}

while getopts "hl:r:t:" option; do
  case "${option}" in
	h) usage; exit;;
	l) l=$OPTARG;;
	r) r=$OPTARG;;
	t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; usage >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$l" ] || [ ! "$r" ] || [ ! "$t" ]; then
  echo "arguments -l, -r, and -t must be provided"
  usage >&2; exit 1
fi


begin=`date +%s`
# Read working directory as variable "dir"
dir=$(pwd)
mkdir all_Fastq

##1 Download SRA files from the given list
~/sratoolkit.3.0.5-centos_linux64/bin/prefetch --option-file $l

#1.1 Move .sra to working directory and remove directory
echo "##############################################################"
echo "				Coverting sra files to fastq files				"
echo "##############################################################"
while IFS=: read -r List; do
	#Reformat name to without \n
	fix_List=${List::-1}
	mv ${fix_List}/* $dir
	~/sratoolkit.3.0.5-centos_linux64/bin/fasterq-dump ${fix_List}.sra --split-3 --threads $t
	#Rename files or _1 and _2 to _R1 and _R2
	mv ${fix_List}_1.fastq ${fix_List}_R1.fastq
	mv ${fix_List}_2.fastq ${fix_List}_R2.fastq
	#Remove directory
	rm -d ${fix_List}/
	#Move files fastq to dir all_Fastq
	mv *.fastq  ${dir}/all_Fastq/
done < $dir/$l
echo "Download SRA and coverting to fastq..."
echo "Complete"

echo "##############################################################"
echo "				Assembling with HAVoC pipeline"
echo "##############################################################"
#2. Excuting HAVoC pipeline and Merging > files_consensus.fa

#2.1 Creating directory for _consensus.fa 
mkdir $dir/consensus
mkdir $dir/consensus/merged_fa
mv $dir/consensus/*_consensus.fa $dir/consensus/merged_fa

#HAVoC pipeline
cd $dir/havoc/
echo "Processig in" $dir/havoc/
bash HAVoC.sh $dir/all_Fastq/
echo "HAVoC pipeline Done"

cd $dir/
mv $dir/all_Fastq/*/*_consensus.fa $dir/consensus
#2.2 Merging all files_consensus.fa > Merged_fa.fasta
echo "...Merging files consensus.fa > Merged_fa.fasta"
cat $dir/consensus/*_consensus.fa > $dir/consensus/Merged_fa.fasta

#2.3 creating directory for pre-trim fastq
mkdir $dir/before_trimming
mv $dir/all_Fastq/*/*.fastq $dir/before_trimming/
echo "Done"

echo "##############################################################"
echo "				Identifying lineage with Pangolin"
echo "##############################################################"
#3. Pangolin
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pangolin
~/miniconda3/envs/pangolin/bin/pangolin --update
~/miniconda3/envs/pangolin/bin/pangolin $dir/consensus/Merged_*.fasta -o $dir/consensus/ -t $t
mv -f $dir/consensus/Merged_fa.fasta $dir/consensus/merged_fa/
rm -r all_Fastq/
conda deactivate
mkdir $dir/Pre_SRA_Omicron
echo "Identifying lineage done"

#4 Trimming with FastP
#4.1 Moving AutoRead_Omicron.py and lineage_report to working directory
mv $dir/consensus/lineage_report.csv $dir/before_trimming/
mkdir $dir/Trim
cd $dir/before_trimming/
mv $dir/AutoRead_Omicron.py $dir/before_trimming/
#Excuting Filter Omicron Accession_number
echo "Reading omicron accession numbers...."
./AutoRead_Omicron.py
echo "Omicron accession numbers"
cat $dir/before_trimming/Omicron_Lineage.txt
echo "complete"

#4.2 Read Accession_number in Omicron_Lineage.txt
echo "Selecting omicron sra files to dir Pre_SRA_Omicron"
while IFS=: read -r name; do
	mv $dir/${name}.sra $dir/Pre_SRA_Omicron/
	gzip ${name}_R1.fastq ${name}_R2.fastq
	fastp --in1 ${name}_R1.fastq.gz --in2 ${name}_R2.fastq.gz --out1 ${name}_1.fastq.gz.fastp.gz --out2 ${name}_2.fastq.gz.fastp.gz -w $t --compression 4
done < Omicron_Lineage.txt
#Remove non Omicron SRA
echo "removing non-omicron"
rm $dir/*.sra
echo "complete"

mv $dir/before_trimming/AutoRead_Omicron.py $dir/
mv $dir/before_trimming/lineage_report.csv $dir/consensus/

#5 Quality checking with FastQC & MultiQC
echo "Quality checking with FastQC....."
rm $dir/before_trimming/fastp.html $dir/before_trimming/fastp.json
mkdir $dir/fastqc_report
mv $dir/before_trimming/*.fastp.gz $dir/fastqc_report/
cd $dir/fastqc_report/
fastqc *fastp.gz -t $t
echo "Merging fastqc results with MultiQC....."
multiqc . -f -k tsv
echo "Removing unrevelant files"
echo "complete"

#5.1 Keep Trim fastq .fastq.gz.fastp.gz in Trim/ directory
echo "Moving Omicron fastq to dir Trim"
mv $dir/fastqc_report/*.fastp.gz $dir/Trim/
rm -r $dir/before_trimming/
echo "complete"

#6. Filter Samples for SARS-CoV-2-Freebayes and VCFtools pipeline
#6.1 Variant calling with Sam_tools pipeline
mv $dir/fastqc_report/multiqc_data/multiqc_fastqc.txt  $dir/Trim/
mv $dir/Filter_QC.py $dir/Trim/
cd  $dir/Trim/

# Executing Filter_QC.py for choosing Pass QC Sample
echo "Selecting PassQC omicron fastq"
./Filter_QC.py
mv $dir/Trim/Filter_QC.py $dir/
cd $dir/SAM_tool_pipeline/
cat $dir/Trim/QC_Pass.txt
echo "complete"

#6.1 Variant calling with Sam_tools pipeline
echo "################## Calling Sam_tools pipeline ##################"
bowtie2-build -f $r ref
while IFS=: read -r sample; do
	#Move PassQC_Omicron SRA to SARS-CoV-2-Freebayes pipeline
	mv $dir/Pre_SRA_Omicron/${sample}.sra $dir/SARS-CoV-2-Freebayes/
	
	mv $dir/Trim/${sample}_* $dir/SAM_tool_pipeline/
	#6.1.1 Aligning trimmed Reads output .sam
	bowtie2 -p $t -x ref -1 ${sample}_1.fastq.gz.fastp.gz -2 ${sample}_2.fastq.gz.fastp.gz -S ${sample}.sam
	rm ${sample}*.fastq.gz.fastp.gz
	#6.1.2 Converting SAM to BAM
	samtools view -b -S -o ${sample}.bam ${sample}.sam
	#Remove .sam
	rm ${sample}.sam
	#6.1.3 Sorting and Indexing
	#6.1.3.1 Sorting output .sorted.bam
	samtools sort ${sample}.bam -T temp -o ${sample}.sorted.bam
	#Remove .bam
	rm ${sample}.bam
	#6.1.3.2 Indexing output .sorted.bam.bai
	samtools index ${sample}.sorted.bam
	#6.1.4 Identifying Genomic Variants
	echo "variant calling" ${sample}.sorted.bam
	bcftools mpileup -f $r ${sample}.sorted.bam  --threads $t | bcftools call --ploidy 1 -mv -Oz --threads $t | bcftools view -i 'DP>100' > ${sample}.sorted.bam.vcf --threads $t 
	#Remove .bam .sorted.bam .bam.bai
	rm ${sample}.sorted.bam
	rm ${sample}.sorted.bam.bai
done < $dir/Trim/QC_Pass.txt

#Remove non Omicron and passQC directory
rm -r $dir/Pre_SRA_Omicron/
rm -r $dir/Trim/
echo "Variant calling Done."
 
#Activate SARS-CoV-2-Freebayes
source ~/miniconda3/etc/profile.d/conda.sh
conda activate SARS-CoV-2-freebayes

#######################################
### Merging variants using jacquard ###
#######################################
echo "Merging variants using jacquard"
echo ""
echo "for information, please see: https://jacquard.readthedocs.io/en/v0.42/overview.html#why-would-i-use-jacquard"
echo ""
# Removing 0-byte files in folder
find . -size 0 -delete
echo "Merge VCFs using jacquard"
echo ""
jacquard merge --include_all ./ merged.vcf

echo "Left only genotypes in merged VCF"
echo ""
vcfkeepgeno merged.vcf GT > merged.GT.vcf

echo "Splitting variants and header from merged.GT.vcf"
echo ""
grep "#" merged.GT.vcf > header
grep -v "#" merged.GT.vcf > variants.vcf

sed -i 's|1/1|1|'g variants.vcf   # convert diploid to haploid
sed -i 's|0/1|1|'g variants.vcf   # convert diploid to haploid
sed -i 's|1/0|1|'g variants.vcf   # convert diploid to haploid
sed -i 's/[.]/0/'g variants.vcf   # convert points to zeros

# Reconstitute vcf file
echo "Reconstituting vcf file"
echo ""
cat header variants.vcf > merged.fixed.vcf
rm header variants.vcf
sed -i 's/NC_04551202/NC_045512.2/'g merged.fixed.vcf

echo "left-aligning vcf file and fixing names"
echo ""
vcfleftalign -r ${r} merged.fixed.vcf > merged.left.vcf
sed -i 's/|unknown//'g merged.left.vcf

# Calculating Viral Frequencies
echo "Calculating viral frequencies with vcflib"
echo ""
vcffixup merged.left.vcf > merged.AF.raw.vcf
rm merged.fixed.vcf merged.left.vcf merged.GT.vcf merged.vcf

echo "Merging Done."


#snpEff for merged.AF.raw.vcf by SAM_tool_pipeline
conda deactivate
mkdir merged_vcf_Samtool/
~/jdk-17.0.1/bin/java -Xmx4g -jar ~/snpEff/snpEff.jar NC_045512.2 $dir/SAM_tool_pipeline/merged.AF.raw.vcf > $dir/SAM_tool_pipeline/merged_vcf_Samtool/merged_vcf_Samtool.ann.vcf
mv -f $dir/SAM_tool_pipeline/merged.AF.raw.vcf $dir/SAM_tool_pipeline/merged_vcf_Samtool/

#6.2 Variant calling with SARS-CoV-2-Freebayes_pipeline
cd $dir/SARS-CoV-2-Freebayes/
./SARS-CoV-2-NGS-freebayes.sh -l $l -g $r -a 0.1 -t $t

#snpEff for merged.AF.raw.vcf by SARS-CoV-2-Freebayes
conda deactivate
mkdir merged_freebayes/
~/jdk-17.0.1/bin/java -Xmx4g -jar ~/snpEff/snpEff.jar NC_045512.2 $dir/SARS-CoV-2-Freebayes/merged.AF.raw.vcf > $dir/SARS-CoV-2-Freebayes/merged_freebayes/merged_freebayes.ann.vcf
mv -f $dir/SARS-CoV-2-Freebayes/merged.AF.raw.vcf $dir/SARS-CoV-2-Freebayes/merged_freebayes/



rm -r $dir/fastqc_report/

end=`date +%s`
elapsed=`expr $end - $begin`
echo "Time taken: $elapsed seconds"