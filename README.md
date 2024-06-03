# Sarav--Automated and Scalable pipeline for WES data analysis
Sarav is an open-source application for automated sequencing data analysis. The pipeline is written in nextflow, and the running environment is provided in a docker container to make them scalable and reproducible. This README provides information on how to setup and run Sarav pipeline for DNA-seq data analysis.
# Abstract:
The advancement of high-throughput sequencing technologies has propelled genomic research into new frontiers, enabling the rapid generation of vast amounts of sequencing data. However, the analysis of this data presents a significant computational challenge, requiring robust and efficient pipelines. In this study, our team presents an automated pipeline for whole exome sequencing (WES) data analysis, developed using Nextflow and based on the Genome Analysis Toolkit (GATK). The pipeline encompasses two main workflows: germline variant calling using the GATK HaplotypeCaller tool and somatic variant calling employing Mutect2. To enhance reproducibility and scalability, the pipeline was containerized using Docker, ensuring independence and ease of deployment across different computing environments. Furthermore, we have integrated GATK Funcotator into our pipeline for variant annotation, enhancing the interpretability and biological relevance of detected variants. This addition expands the utility of our pipeline by providing detailed functional annotations for genomic variants, facilitating downstream analysis and interpretation. Our pipeline streamlines the analysis process for WES data, from read alignment to variant annotation. The germline workflow identifies genetic variations in the germline genome, while the somatic workflow detects mutations specific to tumor samples, providing valuable insights into cancer genetics and personalized medicine.This thesis presents the design, implementation, and validation of the automated pipeline, highlighting its efficiency, accuracy, and scalability. The pipeline's modular architecture allows for easy customization and adaptation to diverse research needs. Our findings demonstrate the utility of automated pipelines in accelerating genomic data analysis and advancing precision medicine research.
# Workflow used in this pipeline
![image](https://github.com/BharathSaravanann/Sarav-pipeline/assets/167014670/20c44f97-0883-44dc-82e1-5eb650cdc59c)

# Datasource used in this pipeline
1) BWA-MEM (Genome Alignment)
     https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
2) Base Recalibration (Known-sites)
     https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
3) Annotate variants (Germline)
     https://storage.cloud.google.com/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908g.tar.gz
4) Annotate variants (Somatic)
     https://storage.googleapis.com/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz
     
# Setting up pipeline:
# Ensure you have Docker installed on your system. Docker is compatible with various operating systems, including Linux, macOS, and Windows.
# Pull the Docker image containing the pipeline's environment and dependencies
                                            docker pull <image_name>
# Building Docker Images:
# Build the Docker image using the provided Dockerfile:
                                           docker build -t <image_name> .
# Running Docker Container:
# Create and run a Docker container using the built image. Ensure to mount the required directories from your host system to the container for data access and output storage:
                                           docker run -it --rm \
                                           -v /path/to/host/ref:/home/docker/ref \
                                           -v /path/to/host/fastq:/home/docker/fastq \
                                           -v /path/to/host/output:/home/docker/output \
                                           -v /path/to/host/datasource:/home/docker/data_sources \
                                           <image_name>

# Running Nextflow Script:
# Once the Docker container is running, execute the Nextflow script within the container. Use the following command, providing the necessary parameters such as reference file name, input fastq files, output directory, known sites file, and data sources directory:
                                   nextflow run <file.nf> \
                                   --ref /home/docker/ref/<ref_file_name> \
                                   --fastq_dir '/home/docker/fastq/*_{R1,R2}*' \
                                   --output_dir /home/docker/output \
                                   --known_sites_dir /home/docker/ref/<known_sites.vcf> \
                                   --data_sources /home/docker/data_sources/<data_sources_containing_file>
             
Adjust the paths and filenames according to your specific setup.
These instructions outline the process of setting up the Docker environment, building the Docker image, running the Docker container with mounted volumes, and executing the Nextflow script within the container. Ensure that all paths and filenames are accurately specified to avoid any errors during execution.

# Overall architecture of this pipeline:
![image](https://github.com/BharathSaravanann/Sarav-pipeline/assets/167014670/337fff70-0b3a-4521-bb46-6e1048d390b2)

