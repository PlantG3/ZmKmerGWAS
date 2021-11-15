#!/bin/bash
perl homostack \
        --prefix ENF2_3hap \
        --seq Zm00001d002425_homo2Ky21-1.0.4.target.fas \
        --annot Ky21_Zm00031ab074260_T001.adjusted.bed \
        --seq Zm00001d002425_homo2Il14H-1.0.4.target.fas \
        --annot Il14H_Zm00028ab072610_T001.adjusted.bed \
        --seq Zm00001d002425_homo2CML277-1.0.4.target.fas \
        --annot CML277_Zm00024ab072670_T001.adjusted.bed
