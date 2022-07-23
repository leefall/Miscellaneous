#!/usr/bin/bash


sed -n '1~4s/^@/>/p;2~4p' /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTQ/KMT2B_08-54084.fastq > /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTA/KMT2B_08-54084.fasta
