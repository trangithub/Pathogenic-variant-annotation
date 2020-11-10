Sample_root="/dragennfs/area7/annotate/tran_population/testpopulation/504_samples";

chrom_list="/dragennfs/area7/annotate/tran_population/testpopulation/504_samples/all_chrom_list.csv";

readarray -t myarray < $chrom_list

# vep annotation

Cache_root="...";
Fasta_root=".../Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";

for chrom_name in "${myarray[@]}"
do 
	echo $chrom_name

	# select only snps, assume that multiallelic snps have been decomposed with vt tools
	bcftools view -v snps ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.vcf.gz -Oz > ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.snps.vcf.gz


	#run vep annotation throughout chromosomes
	vep -i ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.snps.vcf.gz --cache --dir "${Cache_root}" \
 	--offline \
 	-a GRCh38 \
 	--everything \
 	--custom /vep/plugins/gnomADc/hg38_gnomad_genome.txt.gz,gnomADg,vcf,exact,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
 	--plugin dbscSNV,/vep/plugins/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
 	--plugin dbNSFP,/vep/plugins/dbNSFP/dbNSFP4.0a.gz,gnomAD_genomes_controls_AF,gnomAD_genomes_controls_AFR_AF,gnomAD_genomes_controls_EAS_AF,gnomAD_genomes_controls_NFE_AF,gnomAD_genomes_controls_POPMAX_AF,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_AFR_AF,gnomAD_exomes_controls_EAS_AF,gnomAD_exomes_controls_NFE_AF  \
 	--fork 1 --buffer_size 50 --vcf --force_overwrite \
 	--stats_text  \
 	-fasta "${Fasta_root}" \
 	-o ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.snps.vep.vcf

 	
	# extract data from VEP annoatted files
	# change name for samples file from "AF" to "AF_1KGP". This allow crrect extracttion. First, identify row number xxxx first. The run this script:
	sed -i 'xxxxs/.*/##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF_1KGP|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|gnomAD_exomes_controls_AF|gnomAD_exomes_controls_AFR_AF|gnomAD_exomes_controls_EAS_AF|gnomAD_exomes_controls_NFE_AF|gnomAD_exomes_controls_POPMAX_AF|gnomAD_genomes_controls_AF|gnomAD_genomes_controls_AFR_AF|gnomAD_genomes_controls_EAS_AF|gnomAD_genomes_controls_NFE_AF|gnomAD_genomes_controls_POPMAX_AF|ada_score|rf_score">/' ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.snps.vep.vcf

	#extract
	bcftools +split-vep -f '%CHROM:%POS %SYMBOL %Allele %Existing_variation %IMPACT %Consequence %CLIN_SIG %INFO/AC %INFO/AN %INFO/AF %AF_1KGP %gnomAD_genomes_controls_AF %gnomAD_genomes_controls_EAS_AF %gnomAD_exomes_controls_AF %gnomAD_exomes_controls_EAS_AF \n' -d  ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.snps.vep.vcf > ${Sample_root}/${chrom_name}/testMASH_504.${chrom_name}.vqsr.pass.decompose.snps.vep.extract.txt

	




