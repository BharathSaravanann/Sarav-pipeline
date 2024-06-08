/*
    Usage:
    nextflow run <nextflow_script> --ref <ref_file> --fastq_dir <father_fastq_pattern> --output_dir <output_dir> --known_sites_dir <known_sites_file>
    Example:
    nextflow run somatichard.nf --ref /home/bharath/Lifecell/ref/hg38.fa --fastq_dir '/home/bharath/Lifecell/fastqnew/*_{R1,R2}*' --output_dir /home/bharath/Lifecell/output --known_sites_dir /home/bharath/Lifecell/ref/Homo_sapiens_assembly38.dbsnp138.vcf --data_sources /home/bharath/somatic/funcotator_dataSources.v1.8.hg38.20230908g -resume

    Example:
    docker run -it --rm \
    -v /home/bharath/Lifecell/ref:/home/docker/ref \
    -v /home/bharath/Lifecell/fastqnew:/home/docker/fastq \
    -v /home/bharath/Lifecell/output:/home/docker/output \
    -v /home/bharath/somatic:/home/docker/data_sources \
    wxs1
    
    nextflow run somatic.nf --ref /home/docker/ref/hg38.fa --fastq_dir '/home/docker/fastq/*_{R1,R2}*' --output_dir /home/docker/output --known_sites_dir /home/docker/ref/Homo_sapiens_assembly38.dbsnp138.vcf --data_sources /home/docker/data_sources/funcotator_dataSources.v1.8.hg38.20230908g
    
*/

ref_flag = false
fastq_dir_flag = false
output_dir_flag = false
known_sites_dir_flag = false
data_sources=false

params.ref = ""
params.fastq_dir = ""
params.output_dir = ""
params.known_sites_dir = ""
params.data_sources=""

if (params.ref && params.fastq_dir && params.output_dir && params.known_sites_dir && params.data_sources) {
    ref_flag = true
    fastq_dir_flag = true
    output_dir_flag = true
    known_sites_dir_flag = true
    data_sources=true
} else {
    println "Missing required argument(s). Please provide all required arguments."
    println "Usage: nextflow run somatic_variant_calling.nf --ref <ref_file> --fastq_dir <father_fastq_pattern> --output_dir <output_dir> --known_sites <known_sites_file> --data_sources <data_source_file>"
    exit 1
}

println """\
        =================================================
        |              S O M A T I C   VC               |
        =================================================

             genome reference : ${params.ref}
             father FASTQ location : ${params.fastq_dir}
             output directory : ${params.output_dir}
             known sites for BQSR : ${params.known_sites_dir}
             data for annotation : ${params.data_sources}

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
    
    process alignmentmetrics {
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')
          
          
        input:
        val ref
        path bam
        
        
        output:
        
        path "${bam.baseName}.alignment_summary_metrics.txt"
        
        
        script:
        
        """
          gatk CollectAlignmentSummaryMetrics -R ${params.ref} -I ${bam} -O ${bam.baseName}.alignment_summary_metrics.txt
         

        
        """
        
        
        }
    
    process insertmetrics {
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}",mode:'copy')
    
    input:
    path bam
    
    output:
    path "${bam.baseName}insert_size_metrics.txt"
    path "${bam.baseName}insert_size_histogram.pdf"
    
    
    
    script:
    """
    gatk CollectInsertSizeMetrics \
    -I ${bam} \
    -O ${bam.baseName}insert_size_metrics.txt \
    -H ${bam.baseName}insert_size_histogram.pdf \
    --METRIC_ACCUMULATION_LEVEL LIBRARY

    
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
    
    process analyzecovariates{
        cpus 4 
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')
    
    input:
    path recal_table
    
    output:
    
    path "recalibration_plots.pdf"
    
    
    script:
    """
    gatk AnalyzeCovariates -bqsr ${recal_table} -plots recalibration_plots.pdf
  
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

    process mutect2caller {
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')
    
        input:
        path bam

        output:
        path "${bam.baseName}.mutect2_variants.vcf"

        script:
        """
         gatk Mutect2 -R ${params.ref} -I ${bam} -O ${bam.baseName}.mutect2_variants.vcf
         
        """
    }
    
    process variantevaluation{
    
        cpus 4
        memory '15.5GB'
        publishDir("${params.output_dir}", mode: 'copy')   
    
    input:
    val ref
    path vcf
    
    output:
    path "variant_evaluation_metrics.txt"
    
    
    script:
    
    """
    gatk VariantEval -R ${params.ref} -O variant_evaluation_metrics.txt --eval ${vcf}

    
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

process filtersnp {
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')
    input:
    val ref
    path vcf
    
    output:
    path "${vcf.baseName}_filtered_snps.vcf"
    
    script:
    """
    
      gatk VariantFiltration -R ${params.ref} -V ${vcf} -O ${vcf.baseName}_filtered_snps.vcf \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -8.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -genotype-filter-expression "DP < 20" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 20" \
        -genotype-filter-name "GQ_filter"
    """
}

process filterindels {
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')
    input:
    val ref
    path vcf
    
    output:
    path "${vcf.baseName}_filtered_indels.vcf"
    
    script:
    """
     gatk VariantFiltration -R ${params.ref} -V ${vcf} -O ${vcf.baseName}_filtered_indels.vcf \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -8.0" \
        -genotype-filter-expression "DP < 20" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 20" \
        -genotype-filter-name "GQ_filter"
    """
}

process passsnps{
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')

input:
path vcf


output:
path "${vcf.baseName}analysis-ready-snps.vcf"

script:

"""
gatk SelectVariants \
	--exclude-filtered \
	-V ${vcf} \
	-O ${vcf.baseName}analysis-ready-snps.vcf

"""

}

process passindels{
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')

input:
path vcf


output:
path "${vcf.baseName}analysis-ready-indels.vcf"

script:

"""
gatk SelectVariants \
	--exclude-filtered \
	-V ${vcf} \
	-O ${vcf.baseName}analysis-ready-indels.vcf

"""

}

process failedgenotypesnp {
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')

input:
path vcf



output:
path "${vcf.baseName}analysis-ready-snps-filteredGT.vcf"

script:

"""
cat ${vcf}|grep -v -E "DP_filter|GQ_filter" > ${vcf.baseName}analysis-ready-snps-filteredGT.vcf

"""

}

process failedgenotypeindels {
    cpus 4
    memory '15.5GB'
    publishDir("${params.output_dir}", mode: 'copy')

input:
path vcf



output:
path "${vcf.baseName}analysis-ready-indels-filteredGT.vcf"

script:

"""
cat ${vcf}|grep -v -E "DP_filter|GQ_filter" > ${vcf.baseName}analysis-ready-indels-filteredGT.vcf

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
    insertmetrics(marked_duplicates)
    alignmentmetrics(ref_ch,marked_duplicates)
    recalibrated = baseRecalibrator(marked_duplicates)
    bqsred = applyBQSR(marked_duplicates, recalibrated)
    variants=mutect2caller(bqsred)
    variantevaluation(ref_ch,variants)
    snp_variants=snp(ref_ch,variants)
    indels_variants=indels(ref_ch,variants)
    snpsfilter=filtersnp(ref_ch,snp_variants)
    indelsfilter=filterindels(ref_ch,indels_variants)
    snpfail=passsnps(snpsfilter)
    indelsfail=passindels(indelsfilter)
    finalsnp=failedgenotypesnp(snpfail)
    finalindels=failedgenotypeindels(indelsfail)
    snpannotation(ref_ch,finalsnp)
    indelsannotation(ref_ch,finalindels)
    
    
    
    
}


