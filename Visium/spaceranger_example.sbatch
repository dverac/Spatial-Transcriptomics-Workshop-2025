#!/bin/bash

#SBATCH -J spaceranger_example
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -t 12:00:00
#SBATCH -o /example/file/path/spaceranger_example.log
#SBATCH --partition=tier2q
#SBATCH --mem=48GB

#export MRO_DISK_SPACE_CHECK=disable   #NOTE: some SLURM systems require this option to avoid an error

SPACERANGER=/example/file/path/to/spaceranger

TMPDIR=/file/path/to/tmp
REFDIR=/file/path/to/references
READDIR=/file/path/to/raw/reads
WORKDIR=/file/path/to/workding/directory
IMAGEDIR=$WORKDIR/images

cd $WORKDIR

#Example uses the Human 2020-A reference from 10X
#The "sample_id" in "sample=sample_id" below is the prefix for the sample in the raw reads
#Slide information should be a serial number matching your slide. Possible to run without.
#area will typically be A1, B1, C1, or D1
#jobmode and options are designed for running on the HPC with 8 threads on a single node.

$SPACERANGER count --id=spaceranger_example \
                   --transcriptome=$REFDIR/refdata-gex-GRCh38-2020-A\
                   --probe-set=$REFDIR/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
                   --fastqs=$READDIR \
                   --sample=sample_id \
                   --cytaimage=$IMAGEDIR/imagefile.tif \
                   --slide=slide-number \
                   --area=slide-area \
                   --no-bam \
                   --loupe-alignment=$WORKDIR/alignments/image-alignment.json \
                   --jobmode=local \
                   --localcores=8 \
                   --localmem=48
