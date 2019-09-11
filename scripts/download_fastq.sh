#!/bin/bash

echo "This script requires SRA tools to work"
prefetch SRX6612940
prefetch SRX6612943
prefetch SRX298000
prefetch SRX2768920
prefetch SRR3048043
prefetch SRR3048044
prefetch SRR3048049
prefetch SRR3048050
prefetch SRR5427862
prefetch SRR5427863
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR301/008/ERR3013408/ERR3013408_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR301/008/ERR3013408/ERR3013408_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR301/009/ERR3013409/ERR3013409_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR301/009/ERR3013409/ERR3013409_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR301/000/ERR3013410/ERR3013410_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR301/000/ERR3013410/ERR3013410_2.fastq.gz