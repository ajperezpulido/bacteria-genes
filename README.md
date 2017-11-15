# bacteria-genes
grep -Po "[^\s\t:;]+cin-?[^a-z;/{:]" /data/databases/uniprot/bacteria_jul2017/uniprot_sprot_bacteria.annot | sort | uniq -c | sort -n > cin.terms
grep -P "[^\s\t:;]+cin-?[^a-z;/{:]" /data/databases/uniprot/bacteria_jul2017/uniprot_sprot_bacteria.annot > cin.annot
cut -f1 cin.annot > cin.id
getfastaf.pl cin.id /data/databases/uniprot/tax/uniprot_sprot_bacteria.fasta > cin.fasta
makeblastdb -in cin.fasta -dbtype prot

for i in {1..7}; do blastp -query cin.fasta -db ../s$i.fasta -out s$i\_cin.blastp -num_threads 8 -outfmt '6 qseqid sseqid qlen slen qcovs pident evalue' -evalue 1e-5; sed 's/\./,/g' s$i\_cin.blastp --in-place; done

cd-hit -i cin.fasta -aS 0.9 -c 0.5 -g 1 -o cin_cdhit.fasta -n 2 -T 8

for i in {1..7}; do cut -f1 s$i\_cin.blastp | sort | uniq | wc -l; done

BACTERICINE
./searchBactericines.pl uniprot_sprot_bacteria.annot i bactericin | sort | uniq > cin.id
getfastaf.pl bactericin.id uniprot_sprot_bacteria.fasta > bactericin.fasta
for i in {1..7}; do blastp -query bactericin.fasta -db s$i.fasta -out s$i\_bactericin.blastp -num_threads 8 -outfmt '6 qseqid sseqid qlen slen qcovs pident evalue' -evalue 1e-5; done
for i in {1..7}; do ./filterBlastp.pl s$i\_bactericin.blastp | sort | uniq > s$i\_bactericin.id; grep -wf s$i\_bactericin.id ../contigs/annots/s$i\_UniBac_go_goslim_source.tsv > s$i\_bactericin.annot; done

RESISTANCE
./searchBactericines.pl uniprot_sprot_bacteria.annot i | sort | uniq > cin.id
getfastaf.pl resistance.id uniprot_sprot_bacteria.fasta > resistance.fasta
for i in {1..7}; do blastp -query resistance.fasta -db s$i.fasta -out s$i\_resistance.blastp -num_threads 8 -outfmt '6 qseqid sseqid qlen slen qcovs pident evalue' -evalue 1e-5; done
for i in {1..7}; do ./filterBlastp.pl s$i\_resistance.blastp | sort | uniq > s$i\_resistance.id; grep -wf s$i\_resistance.id ../contigs/annots/s$i\_UniBac_go_goslim_source.tsv > s$i\_resistance.annot; done

wc -l s?_bactericin.id s?_resistance.id

term="quorum"; ./searchAdhesion.pl uniprot_sprot_bacteria.annot i $term | sort | uniq > $term.id; ./searchAdhesion.pl uniprot_sprot_bacteria.annot t $term | sort | uniq > $term.terms; getfastaf.pl $term.id uniprot_sprot_bacter
ia.fasta > $term.fasta; for i in {1..7}; do blastp -query $term.fasta -db s$i.fasta -out s$i\_$term.blastp -num_threads 8 -outfmt '6 qseqid sseqid qlen slen qcovs pident evalue' -evalue 1e-5; done; for i in {1..7}; do ./filter
Blastp.pl s$i\_$term.blastp | sort | uniq > s$i\_$term.id; grep -wf s$i\_$term.id ../contigs/annots/s$i\_UniBac_go_goslim_source.tsv > s$i\_$term.annot; done; wc -l s?_$term.id
