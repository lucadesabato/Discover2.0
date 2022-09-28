tooldir="$1";
paired="$2";
fastqfile1="$3";
fastqfile2="$4";
fastafile="$5";
evalue="$6"
perc_id="$7"

ln -s $fastqfile1 fastqfile1;
ln -s $fastqfile2 fastqfile2;
# FILTER + ASSEMBLE + BLAST FASTQ
if [ $paired = "y" ]
then
  $tooldir/discover/scripts/duk/duk -m filteredO1.fq -k 23 $tooldir/discover/data/HO/O_type.fasta $fastqfile1;
  $tooldir/discover/scripts/duk/duk -m filteredH1.fq -k 23 $tooldir/discover/data/HO/H_type.fasta $fastqfile1;
  cat filteredO1.fq > filteredOH1.fq;
  cat filteredH1.fq >> filteredOH1.fq;
  $tooldir/discover/scripts/duk/duk -m filteredO2.fq -k 23 $tooldir/discover/data/HO/O_type.asta $fastqfile2;
  $tooldir/discover/scripts/duk/duk -m filteredH2.fq -k 23 $tooldir/discover/data/HO/H_type.fasta $fastqfile2;
  cat filteredO2.fq > filteredOH2.fq;
  cat filteredH2.fq >> filteredOH2.fq;
  $tooldir/discover/scripts/fastq-pair-master/build/fastq_pair filteredOH1.fq filteredOH2.fq;
  $tooldir/discover/scripts/fastq-pair-master/build/fastq_pair filteredOH1.fq.single.fq fastqfile2;
  $tooldir/discover/scripts/fastq-pair-master/build/fastq_pair filteredOH2.fq.single.fq fastqfile1;
  cat filteredOH1.fq.paired.fq > filteredOH1_paired.fq;
  cat filteredOH1.fq.single.fq.paired.fq >> filteredOH1_paired.fq;
  cat fastqfile1.paired.fq >> filteredOH1_paired.fq;
  cat filteredOH2.fq.paired.fq > filteredOH2_paired.fq;
  cat fastqfile2.paired.fq >> filteredOH2_paired.fq;
  cat filteredOH2.fq.single.fq.paired.fq >> filteredOH2_paired.fq;
  dukst1filesize=$(wc -c "filteredOH1_paired.fq" | awk '{print $1}');
  dukst2filesize=$(wc -c "filteredOH2_paired.fq" | awk '{print $1}');
  if [ $dukst1filesize -gt 0 ] && [ $dukst2filesize -gt 0 ]
  then
    perl $tooldir/discover/scripts/spades.pl duk_spades.fasta duk_spades_contig_stats duk_spades_scaffolds duk_spades_scaffold_stats duk_spades_log NODE spades.py --disable-gzip-output --isolate --pe1-ff --pe1-1 fastq:filteredOH1_paired.fq --pe1-2 fastq:filteredOH2_paired.fq
    rm -r output_dir;    
    blastn -query duk_spades.fasta -subject $tooldir/discover/data/HO/.fasta -task blastn -evalue $evalue -out duk_O_seqs -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
    blastn -query duk_spades.fasta -subject $tooldir/discover/data/HO/H_type.fasta -task blastn -evalue $evalue -out duk_H_seqs -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
  else
    touch duk_O_seqs;
    touch duk_H_seqs;
  fi
else
  $tooldir/discover/scripts/duk/duk -m filteredO1.fq -k 23 $tooldir/discover/data/HO/O_type.fasta $fastqfile1;
  $tooldir/discover/scripts/duk/duk -m filteredH1.fq -k 23 $tooldir/discover/data/HO/H_type.fasta $fastqfile1;
  cat filteredO1.fq > filteredOH1.fq;
  cat filteredH1.fq >> filteredOH1.fq;
  dukstx1filesize=$(wc -c "filteredOH1.fq" | awk '{print $1}');
  if [ $dukstx1filesize -gt 0 ]
  then
    perl $tooldir/discover/scripts/spades.pl duk_spades.fasta duk_spades_contig_stats duk_spades_scaffolds duk_spades_scaffold_stats duk_spades_log NODE spades.py --disable-gzip-output --isolate --iontorrent -s fastq:filteredOH1.fq;
    rm -r output_dir;
    blastn -query duk_spades.fasta -subject $tooldir/discover/data/HO/O_type.fasta -task blastn -evalue $evalue -out duk_O_seqs -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
    blastn -query duk_spades.fasta -subject $tooldir/discover/data/HO/H_type.fasta -task blastn -evalue $evalue -out duk_H_seqs -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
  else
    touch duk_O_seqs;
    touch duk_H_seqs;
  fi	
fi
# BLAST FASTA
blastn -query $fastafile -subject $tooldir/discover/data/HO/O_type.fasta -task blastn -evalue $evalue -out fasta_O_seqs -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
blastn -query $fastafile -subject $tooldir/discover/data/HO/H_type.fasta -task blastn -evalue $evalue -out fasta_H_seqs -outfmt '6 qseqid sseqid length qcovs score nident positive gaps ppos qframe sframe qseq sseq qlen slen' -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;

# COMBINE
cat duk_O_seqs > serogroup_O;
cat fasta_O_seqs >> serogroup_O;
cat duk_H_seqs > serogroup_H;
cat fasta_H_seqs >> serogroup_H;
