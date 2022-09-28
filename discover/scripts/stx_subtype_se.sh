tooldir="$1";
fastqfile="$2";
fastafile="$3";
evalue="$4"
perc_id="$5"

# ASSEMBLY
mkdir stxdir;
skesa --fastq $fastqfile --contigs_out stxdir/skesa.fasta;
cp $fastafile stxdir/spades.fasta;
rm -r output_dir;

# FILTER + ASSEMBLY
$tooldir/discover/scripts/duk/duk -m stxdir/duk.fq -k 23 $tooldir/discover/data/stx.fasta $fastqfile;
dukfilesize=$(wc -c "stxdir/duk.fq" | awk '{print $1}');
if [ $dukfilesize -gt 0 ]
then
  skesa --fastq stxdir/duk.fq --contigs_out stxdir/duk_skesa.fasta;
  perl $tooldir/discover/scripts/spades.pl duk_spades_contigs duk_spades_contig_stats duk_spades_scaffolds duk_spades_scaffold_stats duk_spades_log NODE spades.py --disable-gzip-output --isolate -t 8 --iontorrent -s stxdir/duk.fq;
  mv duk_spades_contigs stxdir/duk_spades.fasta;
  rm -r output_dir;
  blastn -query stxdir/duk_skesa.fasta -db $tooldir/discover/data/stx -task blastn -evalue $evalue -out stxdir/duk_skesa_seqs -outfmt '6 qseqid sseqid sframe qseq' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
  blastn -query stxdir/duk_spades.fasta -db $tooldir/discover/data/stx -task blastn -evalue $evalue -out stxdir/duk_spades_seqs -outfmt '6 qseqid sseqid sframe qseq' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
else
  touch stxdir/duk_skesa_seqs;
  touch stxdir/duk_spades_seqs;
fi

# SEQUENCE SEARCH
blastn -query stxdir/skesa.fasta -db $tooldir/discover/data/stx -task blastn -evalue $evalue -out stxdir/skesa_seqs -outfmt '6 qseqid sseqid sframe qseq' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
blastn -query stxdir/spades.fasta -db $tooldir/discover/data/stx -task blastn -evalue $evalue -out stxdir/spades_seqs -outfmt '6 qseqid sseqid sframe qseq' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity $perc_id;
# DIVIDE STX1 FROM STX2
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx1") print $1,$3,$4}' stxdir/skesa_seqs > stxdir/stx1_skesa_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx2") print $1,$3,$4}' stxdir/skesa_seqs > stxdir/stx2_skesa_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx1") print $1,$3,$4}' stxdir/duk_skesa_seqs > stxdir/dukstx1_skesa_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx2") print $1,$3,$4}' stxdir/duk_skesa_seqs > stxdir/dukstx2_skesa_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx1") print $1,$3,$4}' stxdir/spades_seqs > stxdir/stx1_spades_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx2") print $1,$3,$4}' stxdir/spades_seqs > stxdir/stx2_spades_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx1") print $1,$3,$4}' stxdir/duk_spades_seqs > stxdir/dukstx1_spades_seqs;
awk 'BEGIN { OFS="\t" }!seen[$1]++ {if (substr($2,1,4)=="stx2") print $1,$3,$4}' stxdir/duk_spades_seqs > stxdir/dukstx2_spades_seqs;
# CREATE COMBINED MULTIFASTA FROM SEQUENCES
perl $tooldir/discover/scripts/MultifastaFromBlast.pl "stxdir/stx1_skesa_seqs,stxdir/dukstx1_skesa_seqs,stxdir/stx1_spades_seqs,stxdir/dukstx1_spades_seqs" "stxdir/multiassembly_stx1.fasta";
perl $tooldir/discover/scripts/MultifastaFromBlast.pl "stxdir/stx2_skesa_seqs,stxdir/dukstx2_skesa_seqs,stxdir/stx2_spades_seqs,stxdir/dukstx2_spades_seqs" "stxdir/multiassembly_stx2.fasta";

# ALIGN AND GET CONSENSUS
stx1filesize=$(wc -c "stxdir/multiassembly_stx1.fasta" | awk '{print $1}');
if [ $stx1filesize -eq 0 ]
then
  touch stxdir/multiassembly_stx1_consensus.fasta;
else
  cat $tooldir/discover/data/stx1.fa >> stxdir/multiassembly_stx1.fasta;
  muscle -in stxdir/multiassembly_stx1.fasta -out stxdir/multiassembly_stx1_aligned.fasta;
  awk 'BEGIN {RS=">" ; ORS=""} substr($1,1,4)!="stx1" {print ">"$0}' stxdir/multiassembly_stx1_aligned.fasta > stxdir/multiassembly_stx1_aligned_clean.fasta;
  awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' stxdir/multiassembly_stx1_aligned_clean.fasta > stxdir/multiassembly_stx1_aligned_linear.fasta;
  python $tooldir/discover/scripts/GetConsensus.py -i stxdir/multiassembly_stx1_aligned_linear.fasta -o stxdir/multiassembly_stx1_consensus.fasta;
fi

stx2filesize=$(wc -c "stxdir/multiassembly_stx2.fasta" | awk '{print $1}');
if [ $stx2filesize -eq 0 ]
then
  touch stxdir/multiassembly_stx2_consensus.fasta;
else
  cat $tooldir/discover/data/stx2.fa >> stxdir/multiassembly_stx2.fasta;
  muscle -in stxdir/multiassembly_stx2.fasta -out stxdir/multiassembly_stx2_aligned.fasta;
  awk 'BEGIN {RS=">" ; ORS=""} substr($1,1,4)!="stx2" {print ">"$0}' stxdir/multiassembly_stx2_aligned.fasta > stxdir/multiassembly_stx2_aligned_clean.fasta;
  awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' stxdir/multiassembly_stx2_aligned_clean.fasta > stxdir/multiassembly_stx2_aligned_linear.fasta;
  python $tooldir/discover/scripts/GetConsensus.py -i stxdir/multiassembly_stx2_aligned_linear.fasta -o stxdir/multiassembly_stx2_consensus.fasta;
fi
cat stxdir/multiassembly_stx1_consensus.fasta > stx.fasta;
cat stxdir/multiassembly_stx2_consensus.fasta >> stx.fasta;
