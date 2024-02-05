# Discover
import argparse
import subprocess
from pathlib import Path
import os
import glob
import operator
import csv

def getStxSubType():
    subprocess.call("blastn -query stx.fasta -db ../../discover/data/stx -task blastn -evalue "+shigatoxin_v[0]+" -out shigatoxin -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -num_threads 8 -strand both -dust yes -max_target_seqs 1 -perc_identity "+shigatoxin_v[1],shell=True)
    stx = open("shigatoxin_fc", "w")
    with open("shigatoxin", "r") as csvfile:
        read = list(csv.reader(csvfile, delimiter="\t"))
        for line in read:
            coverage=(int(line[2])/int(line[-1])*100)
            if coverage >= float(shigatoxin_v[2]) and float(line[8]) >= float(shigatoxin_v[3]):
                stx.write(line[1] + "\t")
                stx.write(str(line[3]) + "\t")
                stx.write(line[8] + "\t")
                stx.write("\n")
    stx.close()
    
    with open("shigatoxin_fc", "r") as csvfile:
        shigatoxin_csv = list(csv.reader(csvfile, delimiter="\t"))
        shigatoxin_typing=sorted(shigatoxin_csv, key=lambda x: float(x[2]), reverse=True)

    with open("../../discover/data/stx_subtypes.txt", "r") as csvfile:
        shigatoxin_types = list(csv.reader(csvfile, delimiter="\t"))

    if len(shigatoxin_typing) == 0:
        str_shigatoxin_subtype = "No subtype match found"
    else:
        # get corresponding subtypes
        str_shigatoxin_subtype = ""
        shigatoxin_subtypes = []
        shigatoxin_subtypes_raw = []

        for subtype in shigatoxin_typing:
            if float(subtype[2]) == 100.0:
                for item in shigatoxin_types:
                    if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                        shigatoxin_subtypes.append(item[1])
                        shigatoxin_subtypes_raw.append(item[1])

        # partial matches
        for subtype in shigatoxin_typing:
            for item in shigatoxin_types:
                if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                        shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[2])) + ")")
                        shigatoxin_subtypes_raw.append(item[1])

        shigatoxin_subtypes.sort()
        str_shigatoxin_subtype = "; ".join(shigatoxin_subtypes)
    return str_shigatoxin_subtype

def getSeroGroup(sero_file):
    antigens=[]
    antigen=sero_file.split("_")[0]
    with open(sero_file, "r") as csvfile:
        read = list(csv.reader(csvfile, delimiter="\t"))
        
    for line in read:
        coverage=(int(line[2])/int(line[-1])*100)
        if coverage >= float(serotyper_v[2]) and float(line[8]) >= float(serotyper_v[3]):
            results=[line[1], coverage, float(line[8])]
            antigens.append(results)
    list_antigen_sorted=sorted(antigens, key=operator.itemgetter(1, 2), reverse=True)
    sero=open("serotype_results_"+antigen, "w")
    serogroup=''
    if list_antigen_sorted:
        serogroup = list_antigen_sorted[0][0].split("_")[-1]
        for line in list_antigen_sorted:
            sero.write(line[0]+"\t")
            sero.write(str(line[1])+"\t")
            sero.write(str(line[2])+"\n")
    else:
        serogroup = antigen + "?"
    return serogroup

def trimm_contigs(list_ngs_data, directory):
    log = open(directory.replace("/", "") + '_logError.txt', 'a')
    for file in list_ngs_data:
        if  file.find("_R")!=-1:
            sample = file.split('_')[0]
        else:
            sample = file.split('.')[0]
        print(sample)
        if file.find('_R1') != -1 and os.path.isfile(file.replace('_R1', '_R2')) == True:
            p1_file = file
            p2_file = file.replace('_R1', '_R2')
        # TRIMMING
            subprocess.call(
            "trimmomatic PE -phred33 " + p1_file + " " + p2_file + " trimmed_" + sample + "_1.fq " + sample + "_1unpaired trimmed_" + sample + "_2.fq " + sample + "_2unpaired" + " SLIDINGWINDOW:"+trimmomatic_v[0]+":"+trimmomatic_v[1] + " LEADING:"+trimmomatic_v[2] + " TRAILING:"+trimmomatic_v[3] + " MINLEN:"+trimmomatic_v[4],
            shell=True)
        # ASSEMBLY
            subprocess.call(
            "perl " + "../discover/scripts/spades.pl " + sample + "_spades_contigs " + sample + "_spades_contig_stats " + sample + "_spades_scaffolds " + sample + "_spades_scaffold_stats " + sample + "_spades_log NODE spades.py --disable-gzip-output --isolate --pe1-ff --pe1-1" + " trimmed_" + sample + "_1.fq " + "--pe1-2 trimmed_" + sample + "_2.fq",
            shell=True)
        # SCAFFOLD FILTERING
            subprocess.call(
            "perl " + "../discover/scripts/filter_spades_repeats.pl -i " + sample + "_spades_contigs -t " + sample + "_spades_contig_stats -c "+spades_v[0] + " -r "+spades_v[1]+ " -l "+ spades_v[2]+ " -o " + sample + "_output_with_repeats -u " + sample + "_scaffolds -n " + sample + "_repeat_sequences_only -e "+spades_v[3]+" -f " + sample + "_discarded_sequences -s " + sample + "_summary",
            shell=True)
            if len(open(sample + "_scaffolds").readlines()) != 0:
                Path(sample + "_disc").mkdir(exist_ok=True)
                subprocess.call("mv " + "trimmed_" + sample + "_1.fq " + sample + "_disc/", shell=True)
                subprocess.call("mv " + "trimmed_" + sample + "_2.fq " + sample + "_disc/", shell=True)
                subprocess.call("cp " + sample + "_scaffolds Scaffolds/", shell=True)
                subprocess.call("mv " + sample + "_scaffolds " + sample + "_disc/", shell=True)
                subprocess.call("rm -rf output_dir", shell=True)
            else:
                log.write(sample + ': No contigs' + "\n")
                subprocess.call("rm " + sample + "_scaffolds", shell=True)

        elif file.find('_R1') != -1 and os.path.isfile(file.replace('_R1', '_R2')) == False:
            log.write(sample + "_R2: NOT FOUND" + "\n")
        elif file.find('_R2') != -1 and os.path.isfile(file.replace('_R2', '_R1')) == False:
            log.write(sample + "_R1: NOT FOUND" + "\n")

        elif file.find('_R1') == -1 and file.find('_R2') == -1:
        # TRIMMING
            subprocess.call(
            "trimmomatic SE -phred33 " + file + " trimmed_" + sample + ".fq" + " SLIDINGWINDOW:"+trimmomatic_v[0]+":"+trimmomatic_v[1] + " LEADING:"+trimmomatic_v[2] + " TRAILING:"+trimmomatic_v[3] + " MINLEN:"+trimmomatic_v[4],
            shell=True)
        # ASSEMBLY
            subprocess.call(
            "perl " + "../discover/scripts/spades.pl " + sample + "_spades_contigs " + sample + "_spades_contig_stats " + sample + "_spades_scaffolds " + sample + "_spades_scaffold_stats " + sample + "_spades_log NODE spades.py --disable-gzip-output --isolate --iontorrent -s trimmed_" + sample + ".fq",
            shell=True)
        # SCAFFOLD FILTERING
            subprocess.call(
            "perl " + "../discover/scripts/filter_spades_repeats.pl -i " + sample + "_spades_contigs -t " + sample + "_spades_contig_stats -c "+spades_v[0] + " -r "+spades_v[1]+ " -l "+ spades_v[2]+ " -o " + sample + "_output_with_repeats -u " + sample + "_scaffolds -n " + sample + "_repeat_sequences_only -e "+ spades_v[3] +" -f " + sample + "_discarded_sequences -s " + sample + "_summary",
            shell=True)
            if len(open(sample + "_scaffolds").readlines()) != 0:
                Path(sample + "_disc").mkdir(exist_ok=True)
                subprocess.call("mv " + " trimmed_" + sample + ".fq " + sample + "_disc/", shell=True)
                subprocess.call("cp " + sample + "_scaffolds Scaffolds/", shell=True)
                subprocess.call("mv " + sample + "_scaffolds " + sample + "_disc/", shell=True)
                subprocess.call("rm -rf output_dir", shell=True)
            else:
                log.write(sample + ': No contigs' + "\n")
                subprocess.call("rm " + sample + "_scaffolds", shell=True)
    log.close()
    return True

def rm_files():
    subprocess.call("rm *spades*", shell=True)
    subprocess.call("rm *unpaired", shell=True)
    subprocess.call("rm *summary", shell=True)
    subprocess.call("rm *sequences", shell=True)
    subprocess.call("rm *sequences_only", shell=True)
    subprocess.call("rm *repeats", shell=True)
    return True

def discover(list_ngs_directories):
    for file in list_ngs_directories:
        sample = file.split('_')[0]
        os.chdir(file)
        
        # VIRULOTYPER
        subprocess.call("abricate --db virulotyper " + sample + "_scaffolds > virulotyper_results.txt", shell=True)

        # MLST
        subprocess.call("mlst --scheme ecoli_achtman_4 " + sample + "_scaffolds > mlst_results.txt", shell=True)

        # AMRGENES
        subprocess.call("abricate -db resfinder " + sample + "_scaffolds > amr_abricate_results.txt", shell=True)

        # CHEwBBACCA
        chewbbca = os.popen("which chewBBACA.py").read().strip()
        subprocess.call("mkdir chewbbca", shell=True)
        subprocess.call("cp " + sample + "_scaffolds ./chewbbca/" + sample + "_scaffolds.fasta", shell=True)
        subprocess.call(
        "python " + chewbbca + " AlleleCall -i ./chewbbca/ -g  ../../discover/chewBBACA_db/schema_chewBBACA_cgMLST_V4 -o chewbacca_results/ --cpu 6 --bsr 0.6 --ptf ../../discover/chewBBACA_db/trained_eColi.trn --fr",
        shell=True)
        subprocess.call("cp chewbacca_results/results_*/results_statistics.tsv ./chewbbca_results.tsv", shell=True)
        subprocess.call("cp chewbacca_results/results_*/results_alleles.tsv ./results_alleles.tsv", shell=True)
        subprocess.call("cp chewbacca_results/results_*/RepeatedLoci.txt ./RepeatedLoci.txt", shell=True)

        if os.path.isfile("trimmed_" + sample + "_1.fq") == True and os.path.isfile("trimmed_" + sample + "_2.fq") == True:

        # SHIGATOXIN TYPER
            subprocess.call("../../discover/scripts/stx_subtype_pe.sh ../.. " + "trimmed_" + sample + "_1.fq " + "trimmed_" + sample + "_2.fq " + sample + "_scaffolds "+shigatoxin_v[0]+" "+shigatoxin_v[1],shell=True)
            results_stx = getStxSubType()
            stx = open("shigatoxin_results.txt", "w")
            stx.write("SHIGATOXIN TYPER: " + "\n")
            if results_stx == "No subtype match found":
                stx.write("No match found")
            else:
                stx.write(results_stx)
            stx.close()

        # SEROTYPER O&H
            subprocess.call("../../discover/scripts/serotype.sh ../.. y " + "trimmed_" + sample + "_1.fq " + "trimmed_" + sample + "_2.fq " + sample + "_scaffolds "+serotyper_v[0]+" "+serotyper_v[1],shell=True)
            results_o = getSeroGroup("serogroup_O")
            results_h = getSeroGroup("serogroup_H")
            serotyper = open("serotyper_results.txt", "w")
            serotyper.write("SEROTYPER O&H " + "\n")
            serotyper.write("SEROTYPER O: " + results_o + "\n")
            serotyper.write("SEROTYPER H: " + results_h)
            serotyper.close()
        else:

        # SHIGATOXIN TYPER
            subprocess.call("../../discover/scripts/stx_subtype_se.sh ../.. " + "trimmed_" + sample + ".fq " + sample + "_scaffolds "+shigatoxin_v[0]+" "+shigatoxin_v[1], shell=True)
            results_stx = getStxSubType()
            stx = open("shigatoxin_results.txt", "w")
            stx.write("SHIGATOXIN TYPER: " + "\n")
            if results_stx == "No subtype match found":
                stx.write("No match found")
            else:
                stx.write(results_stx)
            stx.close()

        # SEROTYPER O&H
            subprocess.call("../../discover/scripts/serotype.sh ../.. n trimmed_" + sample + ".fq xxx " + sample + "_scaffolds "+serotyper_v[0]+" "+serotyper_v[1],shell=True)
            results_o = getSeroGroup('serogroup_O')
            results_h = getSeroGroup('serogroup_H')
            serotyper = open("serotyper_results.txt", "w")
            serotyper.write("SEROTYPER O&H " + "\n")
            serotyper.write("SEROTYPER O: " + results_o + "\n")
            serotyper.write("SEROTYPER H: " + results_h)
            serotyper.close()

        os.chdir('../')
    return True

def read_conf_file(configuration):
    trimmomatic_conf=[]
    spades_conf=[]
    serotyper_conf=[]
    shigatoxin_conf=[]
    csv_file = open(configuration)
    read_csv = list(csv.reader(csv_file, delimiter="\t"))
    for line in read_csv:
        if line[0]=="Trimmomatic":
            if line[2]==line[3]:
                trimmomatic_conf.append(line[2])
            else:
                trimmomatic_conf.append(line[3])
        elif line[0]=="Spades":
            if line[2]==line[3]:
                spades_conf.append(line[2])
            else:
                spades_conf.append(line[3])
        elif line[0]=="Serotyper":
            if line[2]==line[3]:
                serotyper_conf.append(line[2])
            else:
                serotyper_conf.append(line[3])
        elif line[0]=="Shigatoxin typer":
            if line[2]==line[3]:
                shigatoxin_conf.append(line[2])
            else:
                shigatoxin_conf.append(line[3])
    return trimmomatic_conf,spades_conf,serotyper_conf,shigatoxin_conf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='directory', help='relative path to directory with FASTQ data')
    args = parser.parse_args()
    os.chdir('./' + args.directory)
    os.system("ln -s " + os.popen("which trimmomatic").read().strip() + " trimmomatic.jar")
    subprocess.call("mkdir Scaffolds/", shell=True)
    list_ngs_data = glob.glob("*.fastq") + glob.glob("*.fastqsanger") + glob.glob("*.fastq.gz") + glob.glob("*.fastqsanger.gz") + glob.glob("*.fq")+ glob.glob("*.fq.gz")
    trimm_contigs(list_ngs_data, args.directory)
    rm_files()
    list_ngs_directories=glob.glob("*_disc/")
    discover(list_ngs_directories)
    return True

if __name__ == "__main__":
    conf=read_conf_file("conf.txt")
    trimmomatic_v, spades_v, serotyper_v, shigatoxin_v = conf[0], conf[1], conf[2], conf[3]
    main()
    subprocess.call("python ../discover_results2.0.py", shell=True)
