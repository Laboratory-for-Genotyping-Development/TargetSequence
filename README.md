# Target Sequence
Pipeline script collection for NGS Target Sequence.

## Brief workflow
fastq -> SampleBamCall -> IndividualBamCall -> ListupVariants -> GroupGenotyping -> DrawVAFPlot -> ReviseGenotyping -> vcf

## Software requirements
- GATK3
- bwa
- samtools
- picard-tools
- vcftools
- cutadapt or fastx_toolkit
- R

This pipeline requires Grid Engine environment.
- Like Sun Grid Engine (Univa Grid Engine, Open Grid Scheduler, etc...)

## Description of SampleTable
- Name: Fastq file name
- Plate_WellNo: Plate well position number (e.g. Plate01_001)
- Plate_Set: Plate name - primer pool name (e.g. Plate01_set1A)
- Class: Case/Control/NTC/NA
- Bed: Path to BED file
- FCno: Path to fastq dir
- Naming: Illumina NGS machine name (e.g. HiSeq, MiSeq, etc...)
- Fasta: Path to primer fasta file

## Description of IndividualTable
- IDY: Individual ID
- Plate_WellNo: Plate well position number (e.g. Plate01_001)
- Plate: Plate name
- Class: Case/Control/NTC/NA
- Bed: Path to BED file

## Requests
When using these scripts, please include the following paper as citations.
- Hum Mol Genet. 2016 Nov 15;25(22):5027-5034