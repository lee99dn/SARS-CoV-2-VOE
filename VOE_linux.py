#!/usr/bin/env python3
#Developer: Danusorn Lee
import pandas as pd
import subprocess
import os
import getopt
import sys

Epitope = ""
vd = ""
bd = ""
o = ""
# get all arguments except file name
argv = sys.argv[1:]
try:
    args, opts = getopt.getopt(argv, ":he:v:b:o:",
                               ["help", "epitope=", "variantdatabase=", "blastdatabase=", "output="])
    for arg, value in args:
        if arg in ('-e', '--epitope'):
            Epitope = value.upper().strip()
        elif arg in ('-v', '--variantdatabase'):
            vd = value
        elif arg in ('-b', '--blastdatabase'):
            bd = value
        elif arg in ('-o', '--output'):
            o = value
        elif arg == '-h' or arg == '--help':
            print(
                "Usage:\nVOE_linux.py [-e, --epitope <epitope sequence> -v, --variantdatabase <variantdatabase> -b, --blastdatabase <blastdatabase> -o, --output <output> TSV format>]")
except getopt.GetoptError as err:
    print("Error try again")
    print(
        "Usage:\nVOE_linux.py [-h, --help] [-e, --epitope <epitope sequence> -v, --variantdatabase <variantdatabase> -b, --blastdatabase <blastdatabase> -o, --output <output> TSV format>]")
if Epitope != "" and vd != "" and bd != "" and o != "":
    print("\n*******************************************\n"
          "         Variant on epitope analysis"
          "\n*******************************************\n")
    file = open(vd, "r")
    new_database = open("variant_database.tsv", "a+")
    Read_file = file.readlines()
    VCF = ""
    # Clean_Header and changing Merged VCF format
    for line in Read_file:
        if "##" in line:
            pass
        else:
            data = line.strip().split("\t")
            CHROM = data[0].strip()
            POS = data[1].strip()
            REF = data[3].strip()
            ALT = data[4].strip()
            INFO = data[7].strip()
            FORMAT = data[8].strip()
            for i in range(9, len(data)):
                VCF = VCF + str(data[i]).strip() + "\t"
                if i == (len(data) - 1):
                    VCF.strip("\t")
            Read = POS + "\t" + CHROM + "\t" + REF + "\t" + ALT + "\t" + INFO + "\t" + FORMAT + "\t" + VCF + "\n"
            VCF = ""
            new_database.write(Read)
    file.close()
    new_database.close()
    print()
    print("Epitope sequence:", Epitope)
    print("variant database:", vd)
    print("blast database:", bd)
    print("output:", o)
    print("\n*******************************************\n"
          "                Process blast"
          "\n*******************************************\n")
    # Creating blast result
    file_Epitope = open("Epitope.fasta", "a")
    file_Epitope.write(">Epitope1")
    file_Epitope.write("\n" + Epitope)
    file_Epitope.close()
    Format = "7 qacc sacc evalue pident qstart qend sstart send"
    # MakeCDSdatabase for TBLASTN
    Title = bd.strip(".fasta") + " Database nucleotide"
    path_db_in = bd
    path_db_out = bd.strip(".fasta")
    command_CDS = ['makeblastdb', '-in', path_db_in, '-dbtype', 'nucl',
                   '-parse_seqids', '-out', path_db_out, '-title', Title]
    subprocess.call(command_CDS)
    # TBLASTN
    in_blast = "Epitope.fasta"
    out_blast = "BlastResults.txt"
    command_TBLASTN = ['tblastn', '-query', in_blast, '-db', path_db_out, '-seg', 'no', 
                       '-outfmt',
                       Format, '-out', out_blast]
    subprocess.call(command_TBLASTN)
    # Use pandas to read .tsv
    df = pd.read_table("variant_database.tsv", low_memory=False)
    file = open(out_blast, "r+")
    Read_file = file.readlines()
    Start = 0
    Stop = 0
    # For fiding the lowest E_value
    Gene_Name = ""
    Real_Gene_Name = ""
    List_Min = list()
    E_value = 0
    First = 0
    End = 0
    s_start = 0
    s_end = 0
    Sensitivity = 0
    Dict_Gene = {"ORF1ab": ["NC_045512.2_cds_YP_009724389.1_1", 266, 21555],
                 "ORF1ab_prot_a": ["NC_045512.2_cds_YP_009725295.1_2", 266, 13483],
                 "S": ["NC_045512.2_cds_YP_009724390.1_3", 21563, 25384],
                 "ORF3a": ["NC_045512.2_cds_YP_009724391.1_4", 25393, 26220],
                 "E": ["NC_045512.2_cds_YP_009724392.1_5", 26245, 26472],
                 "M": ["NC_045512.2_cds_YP_009724393.1_6", 26523, 27191],
                 "ORF6": ["NC_045512.2_cds_YP_009724394.1_7", 27202, 27387],
                 "ORF7a": ["NC_045512.2_cds_YP_009724395.1_8", 27394, 27759],
                 "ORF7b": ["NC_045512.2_cds_YP_009725318.1_9", 27756, 27887],
                 "ORF8": ["NC_045512.2_cds_YP_009724396.1_10", 27894, 28259],
                 "N": ["NC_045512.2_cds_YP_009724397.2_11", 28274, 29533],
                 "ORF10": ["NC_045512.2_cds_YP_009725255.1_12", 29558, 29674]}
    for line in Read_file:
        if "# " not in line:
            List_line = line.split("\t")
            # Read E-value
            E_value = float(List_line[2])
            List_Min.append(E_value)
    E_value = min(List_Min)
    for line in Read_file:
        if "# " not in line:
            List_line = line.split("\t")
            # Comparing %Identidy with 100, Finding the lowest E_value and query.stop = length of Epitope then Read gene_name
            if float(List_line[3]) == float(100) and float(List_line[2]) == Min and int(List_line[5]) == len(Epitope):
                Gene_Name = List_line[1]
                # Read nucleotide position from start to stop of query amino acid sequences
                s_start = int(List_line[6].strip())
                s_end = int(List_line[7].strip("\n"))
    Blast_file = open(out_blast)
    Blast_result = Blast_file.read()
    Blast_file.close()
    print(Blast_result)
    if all("#" in check_line for check_line in Read_file):
        print("***Blast results are not found***")
    for gene in Dict_Gene:
        # Calculating nucloetide position on SARS-CoV-2 genome from query amino acid sequences
        if Gene_Name == Dict_Gene[gene][0]:
            Real_Gene_Name = gene
            First = s_start + Dict_Gene[gene][1] - 1
            End = s_end + Dict_Gene[gene][1] - 1
    file.close()
    POS = df["POS"]
    INFO = df["INFO"]
    Output_Code = ""
    Dict_Sensitivity = {}
    List_Sense = []
    # List of SRA accession number
    List_Column = list(df)[6:len(list(df)) - 1]
    for name in List_Column:
        Output_Code = ""
        for i in range(First, End + 1):
            # List of INFO at POS in Merged VCF
            Check_Epitope = list(df.loc[df["POS"] == i, "INFO"])
            # List of ALT at POS in Merged VCF
            List_ALT = list(df.loc[df["POS"] == i, "ALT"])
            if i in list(POS):
                for j in range(0, len(List_ALT)):
                    Check_sense = list(df.loc[df["POS"] == i, str(name)])
                    if "," not in List_ALT[j]:
                        List_AF = Check_Epitope[j].split(";")
                        List_INFO = Check_Epitope[j].split("|")
                        AF = float(List_AF[1][3:])
                        AC = float(List_AF[0][3:])
                        NS = int(List_AF[3][3:])
                        Chance = round(float(AF * 100), 4)
                        if List_INFO[2] == "HIGH" and AC != 0:
                            Output_Code = Output_Code + str(i) + "\t" + List_INFO[1] + "\t" + List_INFO[9] + "\t" + \
                                          List_INFO[10] + "\t" + str(AC) + "\t" + str(NS) + "\t" + str(
                                round(AF, 4)) + "\t" + str(
                                Chance) + "\n"
                            List_Sense.append(str(Check_sense[j]))
                        elif List_INFO[2] == "MODERATE" and AC != 0:
                            Output_Code = Output_Code + str(i) + "\t" + List_INFO[1] + "\t" + List_INFO[9] + "\t" + \
                                          List_INFO[10] + "\t" + str(AC) + "\t" + str(NS) + "\t" + str(
                                round(AF, 4)) + "\t" + str(
                                Chance) + "\n"
                            List_Sense.append(str(Check_sense[j]))
                        else:
                            pass
                    else:
                        List_AF = Check_Epitope[j].split(";")
                        New_List_AF = List_AF[1].split(",")
                        List_AC = List_AF[0].split(",")
                        New_List_ALT = List_ALT[j].split(",")
                        NS = int(List_AF[3][3:])
                        AF = 0
                        Chance = 0
                        AC = 0
                        count = 0
                        New_Check_sense = str(Check_sense[j]).split("/")
                        for k in range(0, len(New_List_ALT)):
                            List = [x + "|" for x in Check_Epitope[j].split("|,") if x]
                            List_INFO = List[k].split("|")
                            if k == 0:
                                AC = int(List_AC[k][3:])
                                AF = AC / NS
                            elif k > 0:
                                AC = int(List_AC[k][0:])
                                AF = AC / NS
                            Chance = round((AF * 100), 4)
                            # Read Variant type, Variant frequency and Allele frequency
                            if List_INFO[2] == "HIGH" and AC != 0:
                                Output_Code = Output_Code + str(i) + "\t" + List_INFO[1] + "\t" + List_INFO[9] + "\t" + \
                                              List_INFO[10] + "\t" + str(AC) + "\t" + str(NS) + "\t" + str(
                                    round(AF, 4)) + "\t" + str(Chance) + "\n"
                                if len(New_Check_sense) > 1:
                                    count = 2
                                else:
                                    count = 1
                            elif List_INFO[2] == "MODERATE" and AC != 0:
                                Output_Code = Output_Code + str(i) + "\t" + List_INFO[1] + "\t" + List_INFO[9] + "\t" + \
                                              List_INFO[10] + "\t" + str(AC) + "\t" + str(NS) + "\t" + str(
                                    round(AF, 4)) + "\t" + str(Chance) + "\n"
                                if len(New_Check_sense) > 1:
                                    count = 2
                                else:
                                    count = 1
                        if count == 1:
                            List_Sense.append(str(Check_sense[j]))
                        if count == 2:
                            for x in New_Check_sense:
                                List_Sense.append(str(x))
        # Are variants found on epitope from each SRA accession
        Dict_Sensitivity[name] = List_Sense
        List_Sense = []
    for key in Dict_Sensitivity:
        New_Data = []
        for data in Dict_Sensitivity[key]:
            if "/" not in data:
                New_Data = New_Data + [data]
            else:
                New_Data = New_Data + data.split("/")
        Dict_Sensitivity[key] = New_Data
    FN = 0
    TP = 0
    for i in Dict_Sensitivity:
        # If variant found on epitope from each SRA accession > False negative (Value in SRA accession column != 0)
        if any(int(c) > 0 for c in Dict_Sensitivity[i]):
            FN += 1
        # variant not found > True positive (Value in SRA accession column == 0)
        else:
            TP += 1
    Sensitivity = (1 - (FN / (FN + TP))) * 100
    Output = open(o, "w+")
    Output.write("Epitope sequence: " + Epitope + " Gene: " + Real_Gene_Name + "\n")
    if Sensitivity == 100:
        Output.write("Sensitivity = 100.00%")
        print("\n*******************************************\n"
              "                 VOE result"
              "\n*******************************************\n")
        print("Epitope sequence:", Epitope, " Gene:", Real_Gene_Name)
        print("Sensitivity =", Sensitivity, "%")
        print("Done")
    else:
        print("\n*******************************************\n"
              "                 VOE result"
              "\n*******************************************\n")
        Output.write(
            "POS_Genome" + "\t" + "Type" + "\t" + "ALT" + "\t" + "AA_change" + "\t" + "Allele_Count(AC)" + "\t" + "Sample_Count(NS)" + "\t" + "Allele_frequency(AF)" + "\t" + "Chance(%)\n")
        Output.write(Output_Code)
        Output.write("FN:" + str(FN) + "\n")
        Output.write("TP:" + str(TP) + "\n")
        Output.write("Sensitivity:" + str(round(Sensitivity, 4)) + " %")
        print("Epitope sequence:", Epitope, " Gene:", Real_Gene_Name)
        print(
            "POS_Genome" + "\t" + "Type" + "\t" + "ALT" + "\t" + "AA_change" + "\t" + "Allele_Count(AC)" + "\t" + "Sample_Count(NS)" + "\t" + "Allele_frequency(AF)" + "\t" + "Chance(%)")
        print(Output_Code.strip("\n"))
        print("FN=", FN)
        print("TP=", TP)
        print("Sensitivity=", round(Sensitivity, 4), "%")
        print("Done")
    Output.write("\nVariant database: " + vd) 
    os.remove("BlastResults.txt")
    os.remove("Epitope.fasta")
    os.remove("variant_database.tsv")
    Output.close()
else:
    pass
