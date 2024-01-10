#Discover
import os, sys
import glob
import csv
import operator
from numpy import mean
import subprocess
from Bio import SeqIO

#extract gene with coverage >80% from abricate output and return two lists:
#1) geni_ED: list of gene with information of coverage (lenght), % identity and coverage
#2) list_geni: list of found gene
def extract_abricate(csv_abricate):
    geni_ED=[]
    list_geni=[]
    if csv_abricate[0][4]!="STRAND":
        for line in csv_abricate[1:]:
            if float(line[8])>=float(abricate_v[0]) and float(line[9])>=float(abricate_v[1]):
                in_line=[]
                in_line.append(line[4])
                list_geni.append(line[4].split("_")[0])
                if float(line[8])<99.5:
                    in_line.append(float(line[8]))
                else:
                    in_line.append(round(float(line[8])))
                in_line.append(float(line[9]))
                in_line.append(float(line[1].split('_')[-1]))
                geni_ED.append(in_line)
    else:
        for line in csv_abricate[1:]:
            if float(line[9])>=float(abricate_v[0]) and float(line[10])>=float(abricate_v[1]):
                in_line=[]
                in_line.append(line[5])
                list_geni.append(line[5].split("_")[0])
                if float(line[9])<99.5:
                    in_line.append(float(line[9]))
                else:
                    in_line.append(round(float(line[9])))
                in_line.append(float(line[10]))
                in_line.append(float(line[1].split('_')[-1]))
                geni_ED.append(in_line)
    list_geni=list(set(list_geni))
    return geni_ED

#order genes by coverage, identity and collect the mean coverage. Return the list of best gene and the mean K-mer scaffold coverage
def select_best_abricate (gene_to_select, signal):
    ED_best_gene=[]
    tot_mean=[]
    if signal=="virulotyper":
        genes=[line[0].split("_")[0].split("-")[0] for line in gene_to_select]
    else:
        genes=[line[0].split("_")[0] for line in gene_to_select]
    for geni in genes:
        list_gene=[]
        for line in gene_to_select:
            if line[0].find(geni)!=-1:
                list_gene.append(line)
                tot_mean.append(line[3])
        list_gene_sorted=sorted(list_gene, key=operator.itemgetter(1, 2, 3), reverse=True)
        ED_best_gene.append(list_gene_sorted[0])
    #coverage mean
    mean_gene=round(mean(tot_mean))
    return [ED_best_gene, mean_gene]

#list of Virulence gene - OLD
#list_gene_fine=['toxb', 'neuc', 'saa', 'agga', 'espp', 'hlye', 'epea', 'bfpa', 'agg3c', 'aggr', 'f17g', 'capu', 'tir', 'agg3a', 'eata', 'ipah9', 
#'foccsfae', 'papc', 'irea', 'sta1', 'sfad', 'aaic', 'iucc', 'focg', 'espj', 'etpd', 'fyua', 'fim41a', 'vat', 'cia', 'kpsmii', 'orf4', 'mcma', 'chua', 
#'usp', 'afac', 'nlec', 'stb', 'terc', 'katp', 'aafd', 'eila', 'espa', 'ibea', 'espf', 'celb', 'virf', 'afae', 'tsh', 'aggc', 'aar', 'pic', 'cea', 'agg4a', 
#'aggd', 'iuta', 'mchf', 'stx1b', 'agg4c', 'nfae', 'ehxa', 'tccp', 'aata', 'aafc', 'mchc', 'cba', 'aap', 'agg4d', 'iha', 'mchb', 'stx2a', 'afad', 'eae', 
#'pera', 'sfas', 'siga', 'nlea', 'sepa', 'agg3b', 'aafb', 'espb', 'cib', 'tia', 'cif', 'afaa', 'cvac', 'aggb', 'agg5a', 'pet', 'efa1', 'ipad', 'cma', 
#'cdtb', 'gad', 'afab', 'irp2', 'cci', 'afae8', 'rpea', 'yfcv', 'agg4b', 'clbb', 'f17a', 'trat', 'cnf1', 'sita', 'tcpc', 'sth-sta3', 'kpse', 'stx1a', 
#'asta', 'foci', 'sat', 'fedf', 'hlyf', 'espc', 'fasa', 'cofa', 'cfac', 'lpfa', 'orf3', 'papa', 'fana', 'stx2b', 'sfae', 'ompt', 'agg3d', 'lnga', 'hra', 
#'espi', 'aafa', 'suba', 'focc', 'iron', 'nleb', 'senb', 'iss', 'feda', 'etsc', 'kpsmiii', 'kpsm', 'k88ab', 'ltca', 'air', 'mcba']

#list of Virulence gene - NEW
list_gene_fine=[]
fasta_sequences = SeqIO.parse(open("../discover/data/virulence_ecoli.fsa"),'fasta')
for fasta in fasta_sequences:
    gene_vir= str(fasta.id).split(":")[0].split("_")[0]
    if gene_vir.split("-")[0] not in list_gene_fine:
        list_gene_fine.append(gene_vir.split("-")[0])

#Loci from chewBBACA allelecall
loci=list(csv.reader(open("../discover/chewBBACA_db/results_alleles_example.tsv"),delimiter="\t"))[0]

#reading the configuration file to extract minimum coverage and perc identity values to filter abricate results
csv_file = open("../conf.txt")
read_csv = list(csv.reader(csv_file, delimiter="\t"))
abricate_v=[]
for line in read_csv:
    if line[0]=="Abricate":
        if line[2]==line[3]:
            abricate_v.append(line[2])
        else:
            abricate_v.append(line[3])

#Table creation
entries = os.listdir('./')
if 'Tab_results.txt' not in entries:
    tab1=open('Tab_results.txt', 'w')
    prima_riga = 'Sample' + '\tAvg Scaffold coverage' + '\tBurst size' + '\tMLST' + '\tadk' + '\tfumC' + '\tgyrB' + '\ticd' + '\tmdh' + '\tpurA' + '\trecA'+'\tSTX subtype' + '\tSerotype O' + '\tSerotype H\t' 
    tab1.write(prima_riga)
    [tab1.write(gene + '\t') for gene in list_gene_fine[:-1]]
    tab1.write(list_gene_fine[-1]+'\n')
else:
    tab1 = open('Tab_results.txt', 'a')

if 'Tab_AMR.txt' not in entries:
    tab2 = open('Tab_AMR.txt', 'w')
    tab2.write('Sample' + '\t' + 'AMR genes' + "\n")
else:
    tab2 = open('Tab_AMR.txt', 'a')

if 'Tab_cgMLST.txt' not in entries:
    tab3 = open('Tab_cgMLST.txt', 'w')
    # extract info for MLST gene results
    tab3.write("Sample\t")
    [tab3.write(gene + '\t') for gene in loci]
    tab3.write("\n")
else:
    tab3 = open('Tab_cgMLST.txt', 'a')

if 'Contamination_sheet.txt' not in entries:
    tab4 = open('Contamination_sheet.txt', 'w')
    tab4.write('Sample' + '\t' + 'Species'+ '\t' + 'Gene' + '\t' + 'PC' + '\t' + 'NDC' + '\n')
else:
    tab4 = open('Contamination_sheet.txt', 'a')

#gathering results from discover analysis
for file in glob.glob("*_disc/"):
    namesample = file.split('_')[0]
    os.chdir(file)
    row_sample=namesample+'\t'
    print(namesample)

    #VIRULOTYPER RESULTS: read the file virulotyper_results.txt
    with open('virulotyper_results.txt') as file1:
        read_csv = list(csv.reader(file1, delimiter="\t"))

    #empty list from list_gene_fine
    find_gene=['0']*len(list_gene_fine)
    mean_gene = "ND"

    if len(read_csv)>1:
        #exctract best virulence gene
        geni_abricate=extract_abricate(read_csv)
        if geni_abricate:
            extract_geni_abricate=select_best_abricate(geni_abricate, "virulotyper")
            best_geni_abricate=extract_geni_abricate[0]
            mean_gene=str(extract_geni_abricate[1])
    
            for x in best_geni_abricate:
                gen = x[0].split('_')[0]
                for c in list_gene_fine:
                    if gen.find(c)!=-1:
                        find_gene[list_gene_fine.index(c)]=x[0]

    #MLST
    mlst=''
    species=''
    with open('mlst_results.txt') as file2:
        read_csv = list(csv.reader(file2, delimiter="\t"))
    for line in read_csv:
        mlst+='ST'+line[2]+'; '
        st_allele=line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]
        species+=line[1]+'; '

    if os.path.isfile("./chewbbca_results.tsv") == True:
        #extract EXC+INF from chewBBACA
        with open('chewbbca_results.tsv') as file3:
            read_csv = list(csv.reader(file3, delimiter="\t"))
        if len(read_csv)>1:
            exc=str(int(read_csv[1][1])+int(read_csv[1][2]))
    else:
        exc="ND"

    # SHIGATOXIN TYPER
    file4 = open("shigatoxin_results.txt")
    file4_lines=file4.readlines()
    list_of_stx = ''
    if file4_lines[1].find("No match found")==-1:
        for c in file4_lines[1:]:
            list_of_stx += c[:5] + '; '
    else:
        list_of_stx+="ND  "
    file4.close()        

    # SEROTYPER- STX O&H
    file5 = open("serotyper_results.txt")
    file5_lines=file5.readlines()
    sero_o = file5_lines[1].strip('\n').split(': ')[1]
    sero_h = file5_lines[2].strip('\n').split(': ')[1]
    file5.close()

    #add info to sample row and create tab1 Tab_results.txt
    row_sample=namesample+'\t'
    row_sample+=mean_gene+'\t'
    row_sample+=exc+'\t'
    row_sample+=mlst[:-2]+'\t'
    row_sample+=st_allele+'\t'
    row_sample += list_of_stx[:-2] + '\t'
    row_sample += sero_o + '\t'
    row_sample +=sero_h + '\t'
    tab1.write(row_sample)
    [tab1.write(gene + '\t') for gene in find_gene]
    tab1.write('\n')

    #extract info for AMR results
    with open('amr_abricate_results.txt') as file6:
        read_csv = list(csv.reader(file6, delimiter="\t"))
    tab2.write(namesample + ' :' + '\t')
    if len(read_csv)>1:
        geni_amr_abricate=extract_abricate(read_csv)
        if geni_amr_abricate:
            best_geniAmr_abricate=select_best_abricate(geni_amr_abricate, "amr")[0]
            #mean_gene=extract_geniAmr_abricate[1]
            [tab2.write(gene[0] + "; ") for gene in best_geniAmr_abricate]
            tab2.write('\n')
        else:
            tab2.write(' ND' + "\n")
    else:
        tab2.write(' ND' + "\n")

    #extract info for cgMLST gene results
    if os.path.isfile("results_alleles.tsv") == True:
        tab3.write(namesample + '\t')
        with open('results_alleles.tsv') as file7:
            read_csv = list(csv.reader(file7, delimiter="\t"))
        [tab3.write(value + '\t') for value in read_csv[1][1:]]
        tab3.write("\n")
        csv_file.close()
    else:
        [tab3.write('ND' + '\t') for locus in loci]
        tab3.write("\n")

    if os.path.isfile("RepeatedLoci.txt") == True:
        #MLST

        with open('RepeatedLoci.txt') as file7:
            read_csv = list(csv.reader(file7, delimiter="\t"))
        row_sample = namesample + '\t' + species[:-2] + '\t'
        if len(read_csv)>1:
            list_loci = []
            for line in read_csv[1:]:
                tab4.write(row_sample + line[0] + '\t' + line[1] + '\t' + line[2]+ '\n')
        else:
            tab4.write(row_sample + '\t' + '\t' + '\t' + '\n')
    os.chdir('../')

tab1.close()
tab2.close()
tab3.close()
tab4.close()
