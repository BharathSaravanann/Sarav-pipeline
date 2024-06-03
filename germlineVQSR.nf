/*
    Usage:
    nextflow run <nextflow_script> --ref <ref_file> --fastq_dir <father_fastq_pattern> --output_dir <output_dir> --known_sites_dir <known_sites_file>
    
    
    Example:
    nextflow run germlinecopy.nf --ref /home/bharath/Lifecell/ref/hg38.fa --fastq_dir '/home/bharath/Lifecell/fastq/*_{R1,R2}*' --output_dir /home/bharath/Lifecell/output --known_sites_dir /home/bharath/Lifecell/ref/Homo_sapiens_assembly38.dbsnp138.vcf --data_sources /home/bharath/germline/funcotator_dataSources.v1.8.hg38.20230908g --hapmap /home/bharath/final/vcfdata/hapmap_3.3.hg38.vcf.gz --omni /home/bharath/final/vcfdata/1000G_omni2.5.hg38.vcf.gz --indelsites /home/bharath/final/vcfdata/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --snpsites /home/bharath/final/vcfdata/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    

    Example:
    docker run -it --rm \
    -v /home/bharath/Lifecell/ref:/home/docker/ref \
    -v /home/bharath/Lifecell/fastqnew:/home/docker/fastq \
    -v /home/bharath/Lifecell/output:/home/docker/output \
    -v /home/bharath/somatic:/home/docker/data_sources \
    wxs1


*/

ref_flag = false
fastq_dir_flag = false
output_dir_flag = false
known_sites_dir_flag = false
data_sources=false
hapmap=false
omni=false
indelsites=false
snpsites=false


params.ref = ""
params.fastq_dir = ""
params.output_dir = ""
params.known_sites_dir = ""
params.data_sources=""
params.hapmap=""
params.omni=""
params.indelsites=""
params.snpsites=""


if (params.ref && params.fastq_dir && params.output_dir && params.known_sites_dir && params.data_sources && params.hapmap && params.omni && params.indelsites && params.snpsites ) {
    ref_flag = true
    fastq_dir_flag = true
    output_dir_flag = true
    known_sites_dir_flag = true
    data_sources=true
    hapmap=true
    omni=true
    indelsites=true
    snpsites=true
    
} else {
    println "Missing required argument(s). Please provide all required arguments."
    println "Usage: nextflow run somatic_variant_calling.nf --ref <ref_file> --fastq_dir <father_fastq_pattern> --output_dir <output_dir> --known_sites <known_sites_file> --data_sources <data_source_file> --hapmap <path_hapmap> --omni <path_omni> --indelsites <path_indesites> --snpsites <path_snpsites>"
    exit 1
}

println """\
        =================================================
        |              G E R M L I N E   VC               |
        =================================================

             genome reference    : ${params.ref}
           father FASTQ location : ${params.fastq_dir}
             output directory    : ${params.output_dir}
            known sites for BQSR : ${params.known_sites_dir}
             data for annotation : ${params.data_sources}
             hapmap              : ${params.hapmap}
             omni                : ${params.omni}
             indelknownsites     : ${params.indelsites}
             snpknownsites       : ${params.snpsites}
             
        =================================================
        >               E X E C U T I N G               <
        =================================================
        """
        .stripIndent()
        
        

workflow {                  

    process align {
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        val ref
        tuple val(sample_id), path(fastq)

        output:
        path "${sample_id}.paired.sam"

        script:
        """
        bwa mem -t 4 -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tSM:${sample_id}" ${params.ref} ${fastq} > ${sample_id}.paired.sam
        """
    }

    process markDuplicates {
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        path sam

        output:
        path "${sam.baseName}_sorted_dedup_reads.bam"

        script:
        """
        gatk MarkDuplicatesSpark -I ${sam} -O ${sam.baseName}_sorted_dedup_reads.bam
        """
    }

    process baseRecalibrator {
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        path bam

        output:
        path "${bam.baseName}.recal_data.table"

        script:
        """
        gatk BaseRecalibrator -I ${bam} -R ${params.ref} --known-sites ${params.known_sites_dir} -O ${bam.baseName}.recal_data.table
        """
    }

    process applyBQSR {
        cpus 4 
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        path bam
        path recal_table

        output:
        path "${bam.baseName}.sorted_dedup_bqsr_reads.bam"

        script:
        """
        gatk ApplyBQSR -I ${bam} -R ${params.ref} --bqsr-recal-file ${recal_table} -O ${bam.baseName}.sorted_dedup_bqsr_reads.bam
        """
    }

    process haplotypeCaller {
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')
    
        input:
        path bam

        output:
        path "${bam.baseName}.haplotype_variants.vcf"

        script:
        """
        gatk HaplotypeCaller -R ${params.ref} -I ${bam} -O ${bam.baseName}.haplotype_variants.vcf
        """
    }
   
    process snp {
       cpus 4
       memory '15.5GB'
       publishDir("${params.output_dir}", mode: 'copy')
     
       input:
       val ref
       path vcf
    
       output:
       path "${vcf.baseName}_snp.vcf"
    
       script:
       """
       gatk SelectVariants -R ${params.ref} -V ${vcf} --select-type SNP -O ${vcf.baseName}_snp.vcf
       """
    }

    process indels {
      cpus 4
      memory '15.5GB'
      publishDir("${params.output_dir}", mode: 'copy')
    
      input:
      val ref
      path vcf
    
      output:
      path "${vcf.baseName}_indels.vcf"
    
      script:
      """
      gatk SelectVariants -R ${params.ref} -V ${vcf} --select-type INDEL -O ${vcf.baseName}_indels.vcf
      """
    }

     process VQSRsnp {

        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        val ref
        path vcf

        output:
        path ("${vcf.baseName}_recalibrated_variants1.recal"),emit:recalibratedsnpVcf

        
        

        script:
        """
        gatk VariantRecalibrator \
            -R ${ref} \
            -V ${vcf} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 ${params.omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.snpsites} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.known_sites_dir} \
            -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
            -mode SNP \
            -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
            --output ${vcf.baseName}_recalibrated_variants1.recal \
            --tranches-file ${vcf.baseName}_recalibrated_variants.tranches \

        
        """
    }

     process VQSRsnp1 {

       cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        val ref
        path vcf

        output:
        path("${vcf.baseName}_recalibrated_variants1.tranches"), emit: tranchesFile
        
        

        script:
        """
        gatk VariantRecalibrator \
            -R ${ref} \
            -V ${vcf} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 ${params.omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.snpsites} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.known_sites_dir} \
            -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
            -mode SNP \
            -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
            --output ${vcf.baseName}_recalibrated_variants.recal \
            --tranches-file ${vcf.baseName}_recalibrated_variants1.tranches \

        
        """
    }
     process VQSRindels{

        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        val ref
        path vcf

        output:
        path ("${vcf.baseName}_recalibrated_indels2.recal"),emit:recalibratedindelsVcf
        
           
        

        script:
        """
        gatk VariantRecalibrator \
            -R ${ref} \
            -V ${vcf} \
            -resource:mills,known=true,training=true,truth=true,prior=12.0 ${params.indelsites} \
            -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
            -mode INDEL \
            -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
            --max-gaussians 4 \
            --output ${vcf.baseName}_recalibrated_indels2.recal \
            --tranches-file ${vcf.baseName}_recalibrated_indels.tranches \

        
        """
   
}

    process VQSRindels1 {

    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')

    input:
    val ref
    path vcf

    output:
    path ("${vcf.baseName}_recalibrated_indels.tranches"),emit : tranchesFile

    script:
    """
    gatk VariantRecalibrator \
        -R ${ref} \
        -V ${vcf} \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 ${params.indelsites} \
        -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
        -mode INDEL \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        --max-gaussians 4 \
        --output ${vcf.baseName}_recalibrated_indels.recal \
        --tranches-file ${vcf.baseName}_recalibrated_indels.tranches
    """
}

     process indexsnp{
     cpus 4
     memory '15.5GB'
     publishDir("${params.output_dir}", mode: 'copy')

     input:
     path recal_file


     output:
     path "${recal_file.baseName}.recal.idx"


     script:

     """
      gatk IndexFeatureFile -I ${recal_file} -O ${recal_file.baseName}.recal.idx

     """

   }



     process indexindels{
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        path recal_file


        output:
        path "${recal_file.baseName}.recal.idx"


        script:


        """

        gatk IndexFeatureFile -I ${recal_file} -O ${recal_file.baseName}.recal.idx


        """

      } 





    process applyVQSRsnp{

        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        val ref
        path vcf
        path recalibratedsnpVcf
        path tranchesFile

        output:
        path "${vcf.baseName}_recalibrated_filtered_variants.vcf"

        script:
        """
         gatk ApplyVQSR -R ${ref} -V ${vcf} --truth-sensitivity-filter-level 99.0 \
        --recal-file /home/bharath/final/VQSR/indexfiles/father.paired_sorted_dedup_reads.sorted_dedup_bqsr_reads.haplotype_variants_snp_recalibrated_variants1.recal --tranches-file ${tranchesFile} \
        -mode SNP -O ${vcf.baseName}_recalibrated_filtered_variants.vcf \
        --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
        """

   }


     process applyVQSRindels{


        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')

        input:
        val ref
        path vcf
        path recalibratedindelsVcf
        path tranchesFile

        output:
        path "${vcf.baseName}_recalibrated_filtered_indels.vcf"

        script:
        """
       
            
             gatk ApplyVQSR -R ${ref} -V ${vcf} --truth-sensitivity-filter-level 99.0 \
        --recal-file /home/bharath/final/VQSR/indexfiles/father.paired_sorted_dedup_reads.sorted_dedup_bqsr_reads.haplotype_variants_indels_recalibrated_indels2.recal --tranches-file ${tranchesFile} \
        -mode INDEL -O ${vcf.baseName}_recalibrated_filtered_indels.vcf \
        --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
        """




     }



     process snpannotation{
       cpus 4
       memory '15.5GB'
       publishDir("${params.output_dir}", mode: 'copy')


      input:
      val ref
      path vcf

      output:
      path "${vcf.baseName}analysis-ready-snps-filteredGT-functotated.vcf"


      script:

      """
      gatk Funcotator \
	--variant ${vcf} \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path ${params.data_sources} \
	--output ${vcf.baseName}analysis-ready-snps-filteredGT-functotated.vcf \
	--output-file-format VCF

      """

     }

      process indelsannotation{
         cpus 4
         memory '15.5GB'
         publishDir("${params.output_dir}", mode: 'copy')


         input:
         val ref
         path vcf

         output:
         path "${vcf.baseName}analysis-ready-indels-filteredGT-functotated.vcf"


         script:

         """
         gatk Funcotator \
	--variant ${vcf} \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path ${params.data_sources} \
	--output ${vcf.baseName}analysis-ready-indels-filteredGT-functotated.vcf \
	--output-file-format VCF

         """

      }



    ref_ch = Channel.of(params.ref)
    
    fastq_ch = Channel.fromFilePairs(params.fastq_dir)
    
    mapped_reads = align(ref_ch, fastq_ch)
    
    marked_duplicates = markDuplicates(mapped_reads)
    
    recalibrated = baseRecalibrator(marked_duplicates)
    
    bqsred = applyBQSR(marked_duplicates, recalibrated)
    
    variants = haplotypeCaller(bqsred)
    
    snp_variants = snp(ref_ch,variants)
    
    indels_variants = indels(ref_ch,variants)
    
    snpvariantrecalibration1 = VQSRsnp(ref_ch,snp_variants)
    
    indelsvariantrecalibration1 = VQSRindels(ref_ch,indels_variants)
    
    snpvariantrecalibration2 = VQSRsnp1(ref_ch,snp_variants)
    
    indelsvariantrecalibration2 = VQSRindels1(ref_ch,indels_variants)
    
    indexrecal1 = indexsnp(snpvariantrecalibration1)
    
    indexrecal2 = indexindels(indelsvariantrecalibration1)
    
    VQSRSNPAPPLY = applyVQSRsnp(ref_ch,snp_variants,VQSRsnp.out.recalibratedsnpVcf,VQSRsnp1.out.tranchesFile)
    
    VQSRINDELSAPPLY = applyVQSRindels(ref_ch,indels_variants,VQSRindels.out.recalibratedindelsVcf,VQSRindels1.out.tranchesFile)
    
    snpannotation(ref_ch,VQSRSNPAPPLY)
    
    indelsannotation(ref_ch,VQSRINDELSAPPLY)
    
    
    
    
}



