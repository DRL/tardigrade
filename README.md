# tardigrade
Scripts and relevant processed data files for Boothby et al 2015 and Koutsovoulos et al 2015 tardigrade genome papers

## rDNA screening using SILVA
```
blastn -task megablast -query assembly/tg.genome.fsa -db silva/SILVA_123_SSUParc_tax_silva.rDNA.fasta -outfmt '6 qseqid std sscinames sskingdoms stitle' -num_threads 48 -evalue 1e-65 -out silva/unc.vs.SSUParc.1e65.megablast.out 
blastn -task megablast -query assembly/tg.genome.fsa -db silva/SILVA_123_LSUParc_tax_silva.rDNA.fasta -outfmt '6 qseqid std sscinames sskingdoms stitle' -num_threads 48 -evalue 1e-65 -out silva/unc.vs.LSUPARC.1e65.megablast.out  

blastn -task megablast -query georgios_assembly/nHd.2.3.abv500.fna -db silva/SILVA_123_SSUParc_tax_silva.rDNA.fasta -outfmt '6 qseqid std sscinames sskingdoms stitle' -num_threads 48 -evalue 1e-65 -out silva/nHd.vs.SSUPARC.1e65.megablast.out
blastn -task megablast -query georgios_assembly/nHd.2.3.abv500.fna -db silva/SILVA_123_LSUParc_tax_silva.rDNA.fasta -outfmt '6 qseqid std sscinames sskingdoms stitle' -num_threads 48 -evalue 1e-65 -out silva/nHd.vs.LSUPARC.1e65.megablast.out
````

## calculate TPM for CDS from both assemblies
```
# get CDS for UNC
/exports/software/cufflinks/cufflinks-2.2.1.Linux_x86_64/gffread ../gff/tg.default.final.gff -g tg.genome.fsa -x tg.cds.fna 

# run kallisto index
/exports/software/kallisto/kallisto_linux-v0.42.4/kallisto index -i tg.cds.fna.idx -k 31 tg.cds.fna 
/exports/software/kallisto/kallisto_linux-v0.42.4/kallisto index -i nHd.2.3.1.aug.transcripts.UC.idx -k 31 nHd.2.3.1.aug.transcripts.UC.fna 

# run kallisto quant
/exports/software/kallisto/kallisto_linux-v0.42.4/kallisto quant -i tg.cds.fna.idx -o quant_10BS -b 10 -t 64 Tardi_RNASeq_ATCACG_L007_R1_001.fastq.gz Tardi_RNASeq_ATCACG_L007_R2_001.fastq.gz

/exports/software/kallisto/kallisto_linux-v0.42.4/kallisto quant -i nHd.2.3.1.aug.transcripts.UC.idx -o quant_10BS -b 10 -t 64 Tardi_RNASeq_ATCACG_L007_R1_001.fastq.gz Tardi_RNASeq_ATCACG_L007_R2_001.fastq.gz 
```

## get RNAseq average basecov
```
# build index
/exports/software/gsnap/gmap-2015-11-20/bin/bin/gmap_build --dir . -d tg.genome.fsa -k 16 tg.genome.fsa
/exports/software/gsnap/gmap-2015-11-20/bin/bin/gmap_build --dir . -d -k 16 nHd.2.3.abv500.fna

# map RNAseq
/exports/software/gsnap/gmap-2015-11-20/bin/bin/gsnap -t 64 --nofails --novelsplicing=1 --format=sam -D nHd.1.0.contigs.cov.fna/ -d nHd.1.0.contigs.cov.fna Tardi_RNASeq_ATCACG_L007_R1_001.fastq.gz Tardi_RNASeq_ATCACG_L007_R2_001.fastq.gz --gunzip | samtools view -Sbh - > tardi_RNASeq.vs.nHd.1.0.bam
samtools sort -m 100G -O bam -T temp -@ 8 -o tardi_RNASeq.vs.nHd.1.0.sorted.bam tardi_RNASeq.vs.nHd.1.0.bam 
samtools index tardi_RNASeq.vs.nHd.1.0.bam

/exports/software/gsnap/gmap-2015-11-20/bin/bin/gsnap -t 64 --nofails --novelsplicing=1 --format=sam -D nHd.2.3.abv500.fna/ -d nHd.2.3.abv500.fna Tardi_RNASeq_ATCACG_L007_R1_001.fastq.gz Tardi_RNASeq_ATCACG_L007_R2_001.fastq.gz --gunzip | samtools view -Sbh - > tardi_RNASeq.vs.nHd.2.3.bam
samtools sort -m 100G -O bam -T temp -@ 8 -o tardi_RNASeq.vs.nHd.2.3.sorted.bam tardi_RNASeq.vs.nHd.2.3.bam 
samtools index tardi_RNASeq.vs.nHd.2.3.sorted.bam 

/exports/software/gsnap/gmap-2015-11-20/bin/bin/gsnap -t 32 --nofails --novelsplicing=1 --format=sam -D tg.genome.fsa/ -d tg.genome.fsa Tardi_RNASeq_ATCACG_L007_R1_001.fastq.gz Tardi_RNASeq_ATCACG_L007_R2_001.fastq.gz --gunzip | samtools view -Sbh - > tardi_RNASeq.vs.unc.bam
samtools sort -m 100G -O bam -T temp -@ 8 -o tardi_RNASeq.vs.unc.sorted.bam tardi_RNASeq.vs.unc.bam
samtools index tardi_RNASeq.vs.unc.sorted.bam
 
# compute coverage with samtools depth 
samtools depth -Q 30 tardi_RNASeq.vs.unc.sorted.bam > tardi_RNASeq.vs.unc.bam.sorted.depth.txt
samtools depth -Q 30 tardi_RNASeq.vs.nHd.2.3.sorted.bam > tardi_RNASeq.vs.nHd.2.3.bam.sorted.depth.txt

# calculate base cov
./norm_rnaseq_read_cov.py unc.genome.lengths.txt tardi_RNASeq.vs.unc.bam.sorted.depth.txt > tardi_RNASeq.vs.unc.bam.reads_per_contig.norm_rnaseq_depth.tsv
./norm_rnaseq_read_cov.py nHd.2.3.abv500.lengths.txt tardi_RNASeq.vs.nHd.2.3.bam.sorted.depth.txt > tardi_RNASeq.vs.nHd.2.3.bam.reads_per_contig.norm_rnaseq_depth.tsv

# get basecov
cut -f1,4 tardi_RNASeq.vs.unc.bam.reads_per_contig.norm_rnaseq_depth.tsv > tardi_RNASeq.vs.unc.bam.reads_per_contig.norm_rnaseq_depth.name.basecov.txt
cut -f1,4 tardi_RNASeq.vs.nHd.2.3.bam.reads_per_contig.norm_rnaseq_depth.tsv > tardi_RNASeq.vs.nHd.2.3.bam.reads_per_contig.norm_rnaseq_depth.name.basecov.txt

# write cat colour file
perl -ne 'chomp; @temp = split("\t"); if ($temp[1] > 100){print $temp[0].",>100 RNAseq basecov\n"}elsif($temp[1] >= 10 && $temp[1] < 99){print $temp[0].",10-99 RNAseq basecov\n"}elsif($temp[1] >= 1 && $temp[1] < 9){print $temp[0].",1-9 RNAseq basecov\n"}else{print $temp[0].",0 RNAseq basecov\n"}' tardi_RNASeq.vs.unc.bam.reads_per_contig.norm_rnaseq_depth.name.basecov.txt > tardi_RNASeq.vs.unc.bam.reads_cov.catcolour.txt
~/git/blobtools/blobtools plot -i unc.TG-cov.BlobDB.json --catcolour tardi_RNASeq.vs.unc.bam.reads_cov.catcolour.txt

perl -ne 'chomp; @temp = split("\t"); if ($temp[1] > 100){print $temp[0].",>100 RNAseq basecov\n"}elsif($temp[1] >= 10 && $temp[1] < 99){print $temp[0].",10-99 RNAseq basecov\n"}elsif($temp[1] >= 1 && $temp[1] < 9){print $temp[0].",1-9 RNAseq basecov\n"}else{print $temp[0].",0 RNAseq basecov\n"}' tardi_RNASeq.vs.nHd.2.3.bam.reads_per_contig.norm_rnaseq_depth.name.basecov.txt > tardi_RNASeq.vs.nHd.2.3.bam.reads_cov.catcolour.txt
````

## Make blobplots
```
# create
/exports/software/blobtools/blobtools create -i nHd.2.3.abv500.fna -c nHd_Lib350.sam.cov -t ../blast/nHd.2.3.abv500.vs.uniref90.dmnd.out -t nHd.2.3.abv500.vs.nt.10cul2.1e25.megablast.out -o nHd.2.3.nHd_lib350-cov -x bestsumorder -x bestsum

/exports/software/blobtools/blobtools create -i tg.genome.fsa -c TG-300.vs.tg.genome.bam.cov -c TG-500.vs.tg.genome.bam.cov -c TG-800.vs.tg.genome.bam.cov -t tg.genome.fsa.vs.uniref90.dmnd.out -t tg.genome.fsa.vs.nt.megablast.out -o unc.TG-cov -x bestsumorder -x bestsum

/exports/software/blobtools/blobtools create -i assembly/tg.genome.fsa -c unc.Hd_Lib350.sam.cov -t tg.genome.fsa.vs.uniref90.dmnd.out -t tg.genome.fsa.vs.nt.megablast.out -o unc.nHd-cov -x bestsum -x bestsumorder

/exports/software/blobtools/blobtools create -i nHd.1.0.contigs.cov.fna -y velvet -o nHd.1.0 -t nHd.1.0.abv500.vs.uniref90.dmnd.daa.out -t nHd.1.0.vs.nt.10cul2.1e25.megablast.out -x bestsum -x bestsumorder

# plot 
~/git/blobtools/blobtools plot -i nHd.1.0.BlobDB.json -x bestsumorder --format pdf
~/git/blobtools/blobtools plot -i nHd.2.3.nHd_lib350-cov.BlobDB.json -x bestsumorder --format pdf
~/git/blobtools/blobtools plot -i unc.TG-cov.BlobDB.json -x bestsumorder --format pdf
~/git/blobtools/blobtools plot -i unc.nHd-cov.BlobDB.json -x bestsumorder --format pdf

# plot catcolour
~/git/blobtools/blobtools plot -i unc.TG-cov.BlobDB.json --catcolour tardi_RNASeq.vs.unc.bam.reads_cov.catcolour.txt --format pdf
~/git/blobtools/blobtools plot -i nHd.2.3.nHd_lib350-cov.BlobDB.json --catcolour  tardi_RNASeq.vs.nHd.2.3.bam.reads_cov.catcolour.txt --format pdf

# plot comparecov
~/git/blobtools/blobtools comparecov -i unc.TG-cov.BlobDB.json -c unc.Hd_Lib350.sam.cov -r superkingdom --ylabel "Coverage of nHd reads" --xlabel "Coverage of UNC TG reads" -x bestsumorder --log -s --format pdf
~/git/blobtools/blobtools comparecov -i nHd.2.3.nHd_lib350-cov.BlobDB.json -c TG-SUM.cov -r superkingdom --ylabel "Summed coverage of UNC TG libraries" --xlabel "Coverage of nHd reads" -x bestsumorder --log -s --format pdf
```
