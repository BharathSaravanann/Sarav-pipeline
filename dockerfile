FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    unzip \
    openjdk-17-jdk-headless \  
    bwa \
    samtools \
    python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip gatk-4.5.0.0.zip && \
    mv gatk-4.5.0.0 /usr/local/bin/

RUN chmod +x /usr/local/bin/gatk-4.5.0.0
ENV PATH="/usr/local/bin/gatk-4.5.0.0:${PATH}"

RUN wget -qO- https://get.nextflow.io | bash -s -- -v 23.10.1 && \
    mv nextflow /usr/local/bin/

COPY germlinehard.nf /
COPY somatichard.nf /

RUN rm -f /usr/bin/python && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    chmod +x /usr/bin/python

RUN mkdir -p /home/docker/ref /home/docker/fastq /home/docker/output /home/docker/datasource

VOLUME /home/docker/ref
VOLUME /home/docker/fastq
VOLUME /home/docker/output
VOLUME /home/docker/datasource

CMD ["/bin/bash"]
