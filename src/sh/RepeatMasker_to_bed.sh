#!/bin/bash

cd ../../data/annotations/
zcat RepeatMasker_RepeatLibrary20140131_mm10.fa.out.gz | awk -v OFS="\t" 'NR>3{print $5,$6,$7,$10,$1,".",$11}' > RepeatMasker_RepeatLibrary20140131_mm10.bed
cd -
