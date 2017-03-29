STRAINS=001-1 001-2 002-1 002-2 003-1 003-2 004-1 004-2 005-1 005-2 006-1 006-2 007-1 007-2 008-1 008-2 009-1 009-2 010-1 010-2 011-1 011-2 012-1 012-2 013-1 013-2 014-1 014-2 015-1 015-2 016-1 016-2 017-1 017-2 018-1 018-2 019-1 019-2 020-1 020-2 021-1 021-2 022-1 022-2 023-1 023-2 024-1 024-2 025-1 025-2 026-1 026-2 027-1 027-2 028-1 028-2 029-1 029-2 030-1 030-2 031-1 031-2 032-1 032-2 033-1 033-2 034-1 034-2 035-1 035-2 036-1 036-2 
PAIRS=001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036

CPUS=1
REF=NC_000962

## file locations
GENOMES    = genomes
ALIGNMENTS = alignments
PHYLOGENY  = phylogeny
DATA       = data
VCF        = vcf
SPOLPRED   = spolpred

## exe locations

## Workflow

indexed = $(addsuffix .bwt, $(GENOMES)/$(REF).fna )
bam     = $(addprefix $(ALIGNMENTS)/$(REF)_, $(addsuffix .bam, $(STRAINS) ) )
vcf     = $(addprefix $(VCF)/$(REF)_, $(addsuffix .vcf, $(STRAINS) ) )
sites   = $(addprefix $(VCF)/$(REF)_, $(addsuffix .all.vcf.gz, $(STRAINS) ) )
pairs   = $(addprefix $(PHYLOGENY)/individuals/$(REF)_, $(addsuffix .infile, $(PAIRS) ) )
mixed   = $(addprefix $(VCF)/$(REF)_, $(addsuffix .all.vcf.gz.mixed.png, $(STRAINS) ) )
fastq   = $(addprefix $(SPOLPRED)/, $(addsuffix .fastq, $(STRAINS) ) )
spolpred = $(addprefix $(SPOLPRED)/, $(addsuffix .txt, $(STRAINS) ) )
indels  = $(addprefix indels/$(REF)_, $(addsuffix .indels.vcf, $(PAIRS) ) )

all: index alignments sites mixed trees pairs figures
index: $(indexed)
alignments: $(bam)
sites: $(sites)
snps: $(vcf)
pairs: $(pairs)
mixed: $(mixed)
fastq: $(fastq)
spolpred: $(spolpred)
indels: $(indels)

versions:
	@samtools 2>&1 | grep Version | perl -p -e 's/Version/samtools/' > versions.txt
	@bwa 2>&1 | grep Version | perl -p -e 's/Version/bwa/' >> versions.txt
	@raxmlHPC-PTHREADS-SSE3 -v | grep "RAxML version" | perl -p -e 's/.+RAxML version (\d+.\d+.\d+).+/RAxML: \1/' >> versions.txt
	@Rscript --version 2>&1 | grep version | perl -pe 's/R scripting front-end version (.+)/Rscript: \1/' >> versions.txt

$(GENOMES)/$(REF).fna.bwt:
	bwa index $(GENOMES)/$(REF).fna

$(VCF)/$(REF)_%.all.vcf.gz: $(ALIGNMENTS)/$(REF)_%.bam
	samtools mpileup -gf $(GENOMES)/$(REF).fna $(ALIGNMENTS)/$(REF)_$*.bam | bcftools view -cg - > $(VCF)/$(REF)_$*.all.vcf
	gzip -S .gz $(VCF)/$(REF)_$*.all.vcf

$(VCF)/$(REF)_%.vcf: $(ALIGNMENTS)/$(REF)_%.bam
	samtools mpileup -L 20000 -d 20000 -F 0.05 -ugf $(GENOMES)/$(REF).fna $(ALIGNMENTS)/$(REF)_$*.bam | bcftools view -bvcg - > $(VCF)/$*.raw.bcf && bcftools view $(VCF)/$*.raw.bcf  | vcfutils.pl varFilter -D 20000 > $@
	rm $(VCF)/$*.raw.bcf

$(ALIGNMENTS)/$(REF)_%.bam: $(GENOMES)/$(REF).fna.bwt
	bwa mem -t $(CPUS) $(GENOMES)/$(REF).fna $(DATA)/$*_1.fastq.gz $(DATA)/$*_2.fastq.gz | samtools view -bS - | samtools sort -o - $(ALIGNMENTS)/temp.$* | samtools rmdup - - > $(ALIGNMENTS)/$(REF)_$*.bam
	samtools index $@

$(VCF)/$(REF)_%.all.vcf.gz.mixed.png: $(VCF)/$(REF)_%.all.vcf.gz
	perl scripts/mixed_genotype.pl --total --qual 30 --dp 4 --dp4 75 --dpmax 5000 --mq 30 --noindels --vcf $(VCF)/$(REF)_$*.all.vcf.gz

$(PHYLOGENY)/NC_000962.b1.infile:
	perl scripts/snp_caller_mixed.pl --vcf vcf/NC_000962_035-1.all.vcf.gz
	perl scripts/snp_caller.pl --chrom NC_000962.3 --qual 30 --dp 4 --dp4 75 --dpmax 5000 --af 0 --mq 30 --noindels --noheader \
                --vcf 001-1,vcf/NC_000962_001-1.all.vcf.gz \
                --vcf 001-2,vcf/NC_000962_001-2.all.vcf.gz \
                --vcf 002-1,vcf/NC_000962_002-1.all.vcf.gz \
                --vcf 002-2,vcf/NC_000962_002-2.all.vcf.gz \
                --vcf 003-1,vcf/NC_000962_003-1.all.vcf.gz \
                --vcf 003-2,vcf/NC_000962_003-2.all.vcf.gz \
                --vcf 004-1,vcf/NC_000962_004-1.all.vcf.gz \
                --vcf 004-2,vcf/NC_000962_004-2.all.vcf.gz \
                --vcf 005-1,vcf/NC_000962_005-1.all.vcf.gz \
                --vcf 005-2,vcf/NC_000962_005-2.all.vcf.gz \
                --vcf 006-1,vcf/NC_000962_006-1.all.vcf.gz \
                --vcf 006-2,vcf/NC_000962_006-2.all.vcf.gz \
                --vcf 007-1,vcf/NC_000962_007-1.all.vcf.gz \
                --vcf 007-2,vcf/NC_000962_007-2.all.vcf.gz \
                --vcf 008-1,vcf/NC_000962_008-1.all.vcf.gz \
                --vcf 008-2,vcf/NC_000962_008-2.all.vcf.gz \
                --vcf 009-1,vcf/NC_000962_009-1.all.vcf.gz \
                --vcf 009-2,vcf/NC_000962_009-2.all.vcf.gz \
                --vcf 010-1,vcf/NC_000962_010-1.all.vcf.gz \
                --vcf 010-2,vcf/NC_000962_010-2.all.vcf.gz \
                --vcf 011-1,vcf/NC_000962_011-1.all.vcf.gz \
                --vcf 011-2,vcf/NC_000962_011-2.all.vcf.gz \
                --vcf 012-1,vcf/NC_000962_012-1.all.vcf.gz \
                --vcf 012-2,vcf/NC_000962_012-2.all.vcf.gz \
                --vcf 013-1,vcf/NC_000962_013-1.all.vcf.gz \
                --vcf 013-2,vcf/NC_000962_013-2.all.vcf.gz \
                --vcf 014-1,vcf/NC_000962_014-1.all.vcf.gz \
                --vcf 014-2,vcf/NC_000962_014-2.all.vcf.gz \
                --vcf 015-1,vcf/NC_000962_015-1.all.vcf.gz \
                --vcf 015-2,vcf/NC_000962_015-2.all.vcf.gz \
                --vcf 016-1,vcf/NC_000962_016-1.all.vcf.gz \
                --vcf 016-2,vcf/NC_000962_016-2.all.vcf.gz \
                --vcf 017-1,vcf/NC_000962_017-1.all.vcf.gz \
                --vcf 017-2,vcf/NC_000962_017-2.all.vcf.gz \
                --vcf 018-1,vcf/NC_000962_018-1.all.vcf.gz \
                --vcf 018-2,vcf/NC_000962_018-2.all.vcf.gz \
                --vcf 019-1,vcf/NC_000962_019-1.all.vcf.gz \
                --vcf 019-2,vcf/NC_000962_019-2.all.vcf.gz \
                --vcf 020-1,vcf/NC_000962_020-1.all.vcf.gz \
                --vcf 020-2,vcf/NC_000962_020-2.all.vcf.gz \
                --vcf 021-1,vcf/NC_000962_021-1.all.vcf.gz \
                --vcf 021-2,vcf/NC_000962_021-2.all.vcf.gz \
                --vcf 022-1,vcf/NC_000962_022-1.all.vcf.gz \
                --vcf 022-2,vcf/NC_000962_022-2.all.vcf.gz \
                --vcf 023-1,vcf/NC_000962_023-1.all.vcf.gz \
                --vcf 023-2,vcf/NC_000962_023-2.all.vcf.gz \
                --vcf 024-1,vcf/NC_000962_024-1.all.vcf.gz \
                --vcf 024-2,vcf/NC_000962_024-2.all.vcf.gz \
                --vcf 025-1,vcf/NC_000962_025-1.all.vcf.gz \
                --vcf 025-2,vcf/NC_000962_025-2.all.vcf.gz \
                --vcf 026-1,vcf/NC_000962_026-1.all.vcf.gz \
                --vcf 026-2,vcf/NC_000962_026-2.all.vcf.gz \
                --vcf 027-1,vcf/NC_000962_027-1.all.vcf.gz \
                --vcf 027-2,vcf/NC_000962_027-2.all.vcf.gz \
                --vcf 028-1,vcf/NC_000962_028-1.all.vcf.gz \
                --vcf 028-2,vcf/NC_000962_028-2.all.vcf.gz \
                --vcf 029-1,vcf/NC_000962_029-1.all.vcf.gz \
                --vcf 029-2,vcf/NC_000962_029-2.all.vcf.gz \
                --vcf 030-1,vcf/NC_000962_030-1.all.vcf.gz \
                --vcf 030-2,vcf/NC_000962_030-2.all.vcf.gz \
                --vcf 031-1,vcf/NC_000962_031-1.all.vcf.gz \
                --vcf 031-2,vcf/NC_000962_031-2.all.vcf.gz \
                --vcf 032-1,vcf/NC_000962_032-1.all.vcf.gz \
                --vcf 032-2,vcf/NC_000962_032-2.all.vcf.gz \
                --vcf 033-1,vcf/NC_000962_033-1.all.vcf.gz \
                --vcf 033-2,vcf/NC_000962_033-2.all.vcf.gz \
                --vcf 034-1,vcf/NC_000962_034-1.all.vcf.gz \
                --vcf 034-2,vcf/NC_000962_034-2.all.vcf.gz \
                --vcf 035-1-min,vcf/NC_000962_035-1.all.vcf.gz-min.vcf \
                --vcf 035-1-maj,vcf/NC_000962_035-1.all.vcf.gz-maj.vcf \
                --vcf 035-2,vcf/NC_000962_035-2.all.vcf.gz \
                --vcf 036-1,vcf/NC_000962_036-1.all.vcf.gz \
                --vcf 036-2,vcf/NC_000962_036-2.all.vcf.gz \
	--dir $(PHYLOGENY)/ --phylip NC_000962.b1.infile --verbose 1 --refilter -b 1 --cpus $(CPUS)

$(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1: $(PHYLOGENY)/NC_000962.b1.infile
	raxmlHPC-PTHREADS-SSE3 -T $(CPUS) -f a -s $(PHYLOGENY)/NC_000962.b1.infile -x 12345 -p 1234 -# 1000 -m GTRGAMMA -w $(shell pwd)/$(PHYLOGENY) -n NC_000962.b1 -o NC_000962.3
	perl -pi -e "s/NC_000962.3/H37Rv/" $(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1

trees: $(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1

$(PHYLOGENY)/individuals/$(REF)_%.infile:
	perl scripts/snp_caller.pl --chrom NC_000962.3 --qual 30 --dp 4 --dp4 75 --dpmax 5000 --af 1 --mq 30 --noindels --noheader \
		--vcf $*-1,$(VCF)/$(REF)_$*-1.all.vcf.gz \
		--vcf $*-2,$(VCF)/$(REF)_$*-2.all.vcf.gz \
	-dir $(PHYLOGENY)/individuals --phylip $(REF)_$*.infile --verbose 0 --refilter -b 1 --cpus $(CPUS)

figures/figure_4a.png: $(pairs)
	head -1 phylogeny/individuals/*.infile | perl -p -e 's/^\n//' | perl -p -e 's/==> phylogeny\/individuals\/NC_000962_(.+?).infile <==\n/\1/' > figures/figure_4a.txt
	Rscript scripts/figure_4a.R

figures/figure_4b.png: $(pairs)
	Rscript scripts/figure_4b.R

figures/figure_4c.png: $(pairs)
	Rscript scripts/figure_4c.R

figures/figure_4d.png: $(pairs)
	Rscript scripts/figure_4d.R

figures/figure_2.png: $(PHYLOGENY)/RAxML_bipartitions.NC_000962.b1
	Rscript scripts/figure_2.R

figures: figures/figure_2.png figures/figure_4a.png figures/figure_4b.png figures/figure_4c.png figures/figure_4d.png


$(SPOLPRED)/%.fastq:
	perl scripts/merge_fastq.pl $(DATA)/$*_1.fastq.gz $(DATA)/$*_2.fastq.gz $(SPOLPRED)/$*.fastq

$(SPOLPRED)/%.txt: $(SPOLPRED)/%.fastq
	spolpred $(SPOLPRED)/$*.fastq -l 100 -o $(SPOLPRED)/$*.txt

lineages:
	for i in vcf/*.gz ; do perl scripts/check_lineages.pl scripts/lineages.txt $$i; done

indels/$(REF)_%.indels.vcf:
	perl scripts/snp_caller-indels.pl --chrom NC_000962.3 --qual 30 --dp 4 --dp4 75 --dpmax 5000 --af 1 --mq 30 \
		--vcf $*-1,$(VCF)/$(REF)_$*-1.all.vcf.gz \
		--vcf $*-2,$(VCF)/$(REF)_$*-2.all.vcf.gz \
		--verbose 1 --gff /homedirs8/share/NGS/Mycobacterium/Genomes/NC_000962.gff

