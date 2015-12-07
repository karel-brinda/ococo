#! /usr/bin/env bash

cmake . && make && samtools view -h BWA-MEM.bam | ./ococo -d

