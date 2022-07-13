## VCF-file-handeling-and-processing

## Task 1
  
  Current exome/sequencing variant calls are standardly output in variant call format (VCF) format. We provided you with a small VCF file (Test_annotate.vcf) and ask you to check out the features of this file and process it further.

  #### Task 1a: Use the provided Test_annotate.vcf to select only the SNPs using the Genomic Analysis Toolkit aka GATK (thus removing the INDELS). Provide us with the code you used to do this.
  
    ./gatk SelectVariants -select-type SNP -V Test_annotate.vcf -O snp_only.vcf
    # I confirmed the output with bcftools 
           bcftools view --types snps Test_annotate.vcf > test.vcf

  #### Task 1b: Annotate the SNP VCF file you just created using the Refseq and cytoBand databases with the latest standalone version of ANNOVAR. Use a 12bp splice boundary as an argument instead of the default splice boundary definitions.
  ###### Downloading and unzipping annovar
      wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
      tar xvzf annovar.latest.tar.gz
      
  ###### Downloading RefSeq and cytoBand databases for annotation
      annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
      
  ###### Annotating the VCF file:    
       perl  table_annovar.pl  snp_only.vcf   humandb/ -buildver hg19   -out snp_anno -remove  -protocol refGene  -nastring . -arg '-splicing 12'   -vcfinpu --operation g

      How many unique variants are in this file?
          5323 variants are unique out of 5747 SNPs
     
      How many “splicing” variants do you see in the refseq annotation?
         I found 221 splicing variant which included: 
            194 splicing 
            2 ncRNA_splicing 
            1 ncRNA_exonic;splicing 
            24 exonic;splicing
            
      How many “startloss” variants do you see?
          I found 4 "startloss" variant
