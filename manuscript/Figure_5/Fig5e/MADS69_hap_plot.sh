#!/bin/bash
perl homostack \
        --prefix Fig5e.MADS69_4hap \
        --seq MAD69.3.Zm00018ab147610.fasta \
        --annot BZm00018ab147610_T001.adjusted.bed \
        --seq MAD69.3.Zm00041ab149100.fasta \
        --annot Zm00041ab149100_T001.adjusted.bed \
        --seq MAD69.3.Zm00001eb143080.fasta \
        --annot Zm00001eb143080_T001.adjusted.bed \
        --seq MAD69.3.Zm00035ab147480.fasta \
        --annot Zm00035ab147480_T001.adjusted.bed
