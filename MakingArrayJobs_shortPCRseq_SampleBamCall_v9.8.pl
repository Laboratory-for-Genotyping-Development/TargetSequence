#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;
use TS_Params;

# GATK3.7-0でコール
# fastaファイルの拡張子はfasta。samtools用のは直接書き込まれている。他動物では修正

# SampleCallの実行有無
my $job_groups      = 1;                              # まとめて行うジョブ数
my $sample_table    = $TS_Params::sample_table;
my $ontarget_script = 'OnTargetProportion_v3.4.pl';
my $sample_call     = 0;
my $fastqc          = 0;
my $keep_untrimmed  = 0;
my $help            = 0;
GetOptions(
	'call_sample'    => \$sample_call,
	'n_jobs=i'       => \$job_groups,
	'table=s'        => \$sample_table,
	'fastqc'         => \$fastqc,
	'keep_untrimmed' => \$keep_untrimmed,
	'ontarget=s'     => \$ontarget_script,
	'help'           => \$help
);
if ($help) { &usage(); exit; }

unless ( -f $ontarget_script ) {
	die "File not found: $ontarget_script. Exit.";
}
makeDirs( 0, 2, 3, 4 );

my %count      = (), my $line = 1;
my $job_script = 'SampleBamCall.sh';
open( OUT, ">$job_script" ) or die;
open( IN,  $sample_table )  or die "Can't open file: $TS_Params::sample_table";
while (<IN>) {
	if ( $. == 1 ) { next; }
	s/(\r?\n)\z//;
	my ( $no, $plate_wellno, $plate_set, $class, $bed, $input_fastq, $naming,
		$primer_fasta )
	  = split(/\t/);
	unless ( -f $bed ) { next; }
	$bed =~ s/(\(|\))/\\$1/g;

	# naming systemに応じてファイル名変更
	my $filename1, my $filename2;
	if ( $naming eq 'MiSeq' ) {
		$count{$input_fastq}++;
		$filename1 =
		  "${input_fastq}${no}_S" . $count{$input_fastq} . '_L001_R1_001';
	}
	elsif ( $naming eq 'bbj' ) {
		$filename1 = "${input_fastq}bbj_amplicon_${no}_1";
	}
	elsif ( $naming eq 'HiSeq1' ) {
		$filename1 = "${input_fastq}${no}_R1";
	}
	elsif ( $naming =~ /^(Hi|Nova)Seq$/ ) {
		$count{$input_fastq}++;
		( my $path = $input_fastq ) =~ s/(\(|\)|\{|\})/\\$1/g;
		my @f1 = glob("${path}${no}_S*_R1_001.fastq.gz");
		if ( scalar(@f1) > 0 ) {
			$filename1 = $input_fastq . basename( shift(@f1), '.fastq.gz' );
		}
		else {
			$filename1 =
			  "${input_fastq}${no}_S" . $count{$input_fastq} . '_R1_001';
		}
	}
	else {
		print STDERR "Unknown naming system $naming\n";
		exit;
	}

	# fastq.gzが存在しているか確認。
	unless ( -f "${filename1}.fastq.gz" ) {
		print STDERR "No ${filename1}.fastq.gz\n";
		next;
	}
	( $filename2 = $filename1 ) =~ s/_R1/_R2/g;

	# FastQCの結果を出力する場合
	if ($fastqc) {
		makeDirs(1);
		print OUT "$TS_Params::fastqc ${filename1}.fastq.gz"
		  . " --outdir=$TS_Params::output[1]/;";
		print OUT "$TS_Params::fastqc ${filename2}.fastq.gz"
		  . " --outdir=$TS_Params::output[1]/;";
	}

	# adapter trimming
	if ( -f $primer_fasta ) {

		# primer setのfastaファイルがある場合は、cutadapt
		print OUT $TS_Params::cutadapt
		  . " --untrimmed-output=$TS_Params::output[2]/${no}_F.Untrimmed.fastq.gz"
		  . " --untrimmed-paired-output=$TS_Params::output[2]/${no}_R.Untrimmed.fastq.gz"
		  . ' -u -1 -U -1'
		  . " -g file:$primer_fasta"
		  . " -G file:$primer_fasta"
		  . " -o $TS_Params::output[2]/${no}_F.fastq.gz"
		  . " -p $TS_Params::output[2]/${no}_R.fastq.gz"
		  . " ${filename1}.fastq.gz ${filename2}.fastq.gz;";
	}
	else {

		# primer setのfastaファイルがない場合は、fastx_trimmerで20bp cut
		print OUT "gunzip -c ${filename1}.fastq.gz "
		  . "| $TS_Params::fastx_trimmer -Q33 -f 21 | $TS_Params::fastx_trimmer -Q33 -t 1 -o $TS_Params::output[2]/${no}_F.fastq;";
		print OUT "gzip $TS_Params::output[2]/${no}_F.fastq;";
		print OUT "gunzip -c ${filename2}.fastq.gz "
		  . "| $TS_Params::fastx_trimmer -Q33 -f 21 | $TS_Params::fastx_trimmer -Q33 -t 1 -o $TS_Params::output[2]/${no}_R.fastq;";
		print OUT "gzip $TS_Params::output[2]/${no}_R.fastq;";
	}

	# mapping
	print OUT "$TS_Params::bwa mem -M"
	  . " -R \"\@RG\\tID:$no\\tSM:$no\\tPL:Illumina\""
	  . " $TS_Params::huref $TS_Params::output[2]/${no}_F.fastq.gz"
	  . " $TS_Params::output[2]/${no}_R.fastq.gz > $TS_Params::output[3]/${no}.sam;";

	# サンプル単位で変異のコールまでする場合の処理
	if ($sample_call) {
		makeDirs(5);
		print OUT "$TS_Params::samtools view"
		  . " -bS $TS_Params::output[3]/${no}.sam "
		  . "| $TS_Params::samtools sort -o $TS_Params::output[3]/${no}_sort.bam;"
		  ;    # ここを修正
		print OUT
		  "$TS_Params::samtools index $TS_Params::output[3]/${no}_sort.bam;";

		# Indel Realignment
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -I $TS_Params::output[3]/${no}_sort.bam"
		  . " -R ${TS_Params::huref}.fasta -T RealignerTargetCreator -L $bed"
		  . " -o $TS_Params::output[3]/${no}.intervals -rf BadCigar;";
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -I $TS_Params::output[3]/${no}_sort.bam"
		  . " -R ${TS_Params::huref}.fasta -T IndelRealigner"
		  . " --targetIntervals $TS_Params::output[3]/${no}.intervals -o $TS_Params::output[4]/${no}.bam;";
		print OUT "$TS_Params::samtools index $TS_Params::output[4]/${no}.bam;";

		# Base Recalibration
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -I $TS_Params::output[4]/${no}.bam"
		  . " -R ${TS_Params::huref}.fasta -T BaseRecalibrator -L $bed"
		  . " -knownSites $TS_Params::indel1000GVcf"
		  . " -o $TS_Params::output[3]/${no}.grp -rf BadCigar;";
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -I $TS_Params::output[4]/${no}.bam"
		  . " -R ${TS_Params::huref}.fasta -T PrintReads -BQSR $TS_Params::output[3]/${no}.grp"
		  . " -o $TS_Params::output[3]/${no}_recal.bam -rf BadCigar;";

		# HaplotypeCaller
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T HaplotypeCaller -dfrac 1"
		  . " -R ${TS_Params::huref}.fasta -I $TS_Params::output[3]/${no}_recal.bam"
		  . " -L $bed -o $TS_Params::output[5]/${no}.HC.vcf -rf BadCigar;";

		# UnifiedGenotyper
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T UnifiedGenotyper -dfrac 1 -glm BOTH"
		  . " -R ${TS_Params::huref}.fasta -I $TS_Params::output[3]/${no}_recal.bam"
		  . " -L $bed -o $TS_Params::output[5]/${no}.UG.vcf -rf BadCigar;";

		print OUT "rm $TS_Params::output[3]/${no}.intervals;";
		print OUT "rm $TS_Params::output[3]/${no}_*.bam*;";
		print OUT "rm $TS_Params::output[3]/${no}_*.bai;";
		print OUT "rm $TS_Params::output[3]/${no}.grp;";
	}
	else {
		print OUT "$TS_Params::samtools view"
		  . " -bS $TS_Params::output[3]/${no}.sam "
		  . "| $TS_Params::samtools sort -o $TS_Params::output[4]/${no}.bam;"
		  ;    # ここを修正
		print OUT "$TS_Params::samtools index $TS_Params::output[4]/${no}.bam;";
	}
	print OUT
	  "perl $ontarget_script -n $no -r ${TS_Params::huref}.fasta -b $bed;";

	print OUT "rm $TS_Params::output[3]/${no}.sam;";
	if ($keep_untrimmed) {
		print OUT "rm $TS_Params::output[2]/${no}_*F.fastq.gz;";
		print OUT "rm $TS_Params::output[2]/${no}_*R.fastq.gz;";
	}
	else {
		print OUT "rm $TS_Params::output[2]/${no}_*.fastq.gz;";
	}

	print OUT "echo $no >> Success_SampleBamCall.list;";
	if ( ( ( $line++ % $job_groups ) == 0 ) or eof ) {
		print OUT "\n";
	}
}
close(IN);
close(OUT);

open( GRID, ">ArrayJobs_$job_script" ) or die;
print GRID while (<DATA>);
close(GRID);

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -f|--fastqc: Execute FastQC\n";
	print STDOUT " -c|--call_sample: Execute GATK variant call\n";
	print STDOUT " -n|--n_jobs: N of samples per 1 jobs (default: 1)\n";
	print STDOUT " -t|--table: Path to SampleTable"
	  . " (default: $TS_Params::sample_table)\n";
	print STDOUT " -o|--ontarget: Path to OnTargetProportion.pl"
	  . " (default: OnTargetProportion_v3.4.pl)\n";
	print STDOUT " -k|--keep_untrimmed: Do not delete untrimmed fastq files\n";
	print STDOUT " -h|--help: Show this message\n";
}

__DATA__
#!/bin/bash
#$ -N SampleBam
#$ -l h_vmem=16G
#$ -l mem_free=16G
#$ -e ./Logfiles
#$ -o ./Logfiles
#$ -cwd
#$ -l h_rt=24:00:00
#$ -q all.q

#$-S /bin/bash
command=`head -$SGE_TASK_ID SampleBamCall.sh | tail -1`
eval "$command"
