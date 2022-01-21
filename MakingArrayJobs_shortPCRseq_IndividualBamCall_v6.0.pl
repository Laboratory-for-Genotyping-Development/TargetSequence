#!/usr/bin/perl

use strict;
use Getopt::Long;
use TS_Params;

#
my $sample_table     = $TS_Params::sample_table;
my $individual_table = $TS_Params::individual_table;
my $sample_bam_dir   = $TS_Params::output[4];
my $prefix           = '';
my $help             = 0;
GetOptions(
	'sample_table=s'     => \$sample_table,
	'individual_table=s' => \$individual_table,
	'bam_dir=s'          => \$sample_bam_dir,
	'prefix=s'           => \$prefix,
	'help'               => \$help
);
if ($help) { &usage(); exit; }

makeDirs( 0, 3, 7, 8, 9 );
$sample_bam_dir =~ s/\/+\Z//;

my %sample_name;
open( IN, $sample_table ) or die "Can't open $sample_table";
while (<IN>) {
	s/(\r?\n)\z//;
	my ( $no, $name ) = split(/\t/);
	push( @{ $sample_name{$name} }, $no );
}
close(IN);

if ( $prefix ne '' ) { $prefix .= '_'; }
my $job_script = 'IndividualBamCall.sh';

open( IN,    $individual_table ) or die "Can't open file: $individual_table";
open( OUT,   ">$job_script" )    or die;
open( COUNT, ">${prefix}SampleCountForIndividualBam.list" ) or die;
print COUNT
  join( "\t", ( 'No', 'IDY', 'Plate_Well', 'N_in_SampleTable', 'N_bamfiles' ) )
  . "\n";
my $count = 1;
while (<IN>) {
	if ( $. == 1 ) { next; }
	s/(\r?\n)\z//;
	my ( $idy, $plate_no, $plateName, $class, $bed ) = split(/\t/);
	$bed =~ s/(\(|\))/\\$1/g;
	my $option_bed = ( -f $bed ) ? "-L $bed" : '';

	my $sample_count = my $bam_count =
	  exists( $sample_name{$plate_no} )
	  ? scalar( @{ $sample_name{$plate_no} } )
	  : 0;
	my $merge_cmd = '';
	if ( $sample_count > 0 ) {
		$merge_cmd = "java -Xmx4g -jar $TS_Params::picard MergeSamFiles"
		  . " OUTPUT=$TS_Params::output[3]/${idy}_raw.bam";
		foreach my $sample ( @{ $sample_name{$plate_no} } ) {
			if ( -f "$sample_bam_dir/${sample}.bam" ) {
				$merge_cmd .= " INPUT=${sample_bam_dir}/${sample}.bam";
			}
			else {
				print STDERR
				  "No expected bam file: ${sample}.bam of ${plate_no}\n";
				$bam_count--;
			}
		}
	}
	print COUNT
	  join( "\t", ( $. - 1, $idy, $plate_no, $sample_count, $bam_count ) )
	  . "\n";

	if ( $bam_count > 0 ) {
		print OUT "$merge_cmd;";
		print OUT "$TS_Params::samtools view"
		  . " -h $TS_Params::output[3]/${idy}_raw.bam"
		  . " | perl -pse 'BEGIN{\$IDY=~s/.bam//g};if(\$_ =~ /^@/){s/(?<=SM:)[^(\\t|\\\\t)]+/\$IDY/g}' -- -IDY='$idy' | "
		  . "$TS_Params::samtools view - -hb > $TS_Params::output[3]/${idy}_temp.bam;";
		print OUT
		  "$TS_Params::samtools index $TS_Params::output[3]/${idy}_temp.bam;";

		# Indel Realignment
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T RealignerTargetCreator -R ${TS_Params::huref}.fasta"
		  . " $option_bed -I $TS_Params::output[3]/${idy}_temp.bam"
		  . " -o $TS_Params::output[3]/${idy}.intervals -rf BadCigar;";
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T IndelRealigner -R ${TS_Params::huref}.fasta"
		  . " --targetIntervals $TS_Params::output[3]/${idy}.intervals"
		  . " -I $TS_Params::output[3]/${idy}_temp.bam"
		  . " -o $TS_Params::output[7]/${idy}.bam -rf BadCigar;";
		print OUT
		  "$TS_Params::samtools index $TS_Params::output[7]/${idy}.bam;";

		# Base Recalibraion
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T BaseRecalibrator -R ${TS_Params::huref}.fasta"
		  . " $option_bed -knownSites $TS_Params::indel1000GVcf"
		  . " -I $TS_Params::output[7]/${idy}.bam"
		  . " -o $TS_Params::output[7]/${idy}.grp -rf BadCigar;";

		# HaplotypeCaller
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T HaplotypeCaller -R ${TS_Params::huref}.fasta $option_bed -dfrac 1"
		  . " -BQSR $TS_Params::output[7]/${idy}.grp -I $TS_Params::output[7]/${idy}.bam"
		  . " -o $TS_Params::output[8]/${idy}.HC.vcf.gz -rf BadCigar;";

		# Unified Genotyper
		print OUT "java -Xmx4g -jar $TS_Params::GATK"
		  . " -T UnifiedGenotyper -R ${TS_Params::huref}.fasta $option_bed -dfrac 1 -glm BOTH"
		  . " -BQSR $TS_Params::output[7]/${idy}.grp -I $TS_Params::output[7]/${idy}.bam"
		  . " -o $TS_Params::output[8]/${idy}.UG.vcf.gz -rf BadCigar;";

		# Coverageの計算
		print OUT "java -Xmx4g -jar $TS_Params::GATK -T DepthOfCoverage"
		  . " -R ${TS_Params::huref}.fasta $option_bed"
		  . " -I $TS_Params::output[7]/${idy}.bam -o $TS_Params::output[9]/$idy"
		  . ' -rf BadCigar --omitIntervalStatistics --omitLocusTable --omitPerSampleStats;';
		print OUT "gzip -f $TS_Params::output[7]/${idy}.grp;";
		print OUT "gzip -f $TS_Params::output[9]/$idy;";
		print OUT "echo $count >> Success_IndividualBamCall.list;";

		print OUT "rm $TS_Params::output[3]/${idy}.intervals;";
		print OUT "rm $TS_Params::output[3]/${idy}_*.bam*;";
		print OUT "rm $TS_Params::output[7]/${idy}.bai;";
		print OUT "\n";
	}
	$count++;
}
close(COUNT);
close(OUT);
close(IN);

open( GRID, ">ArrayJobs_$job_script" ) or die;
print GRID while (<DATA>);
close(GRID);

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -s|--sample_table: Path to SampleTable "
	  . "(default: $TS_Params::sample_table)\n";
	print STDOUT " -i|--individual_table: Path to IndividualTable "
	  . "(default: $TS_Params::individual_table)\n";
	print STDOUT " -b|--bam_dir: Path to SampleBam directory "
	  . "(default: $TS_Params::output[4])\n";
	print STDOUT " -p|--prefix: Bam count file name prefix\n";
	print STDOUT " -h|--help: Show this message\n";
}

__DATA__
#!/bin/bash
#$ -N IndividualCall
#$ -l h_vmem=12G
#$ -l mem_free=12G
#$ -e ./Logfiles
#$ -o ./Logfiles
#$ -cwd
#$ -l h_rt=12:00:00
#$ -q all.q

#$-S /bin/bash
command=`head -$SGE_TASK_ID IndividualBamCall.sh | tail -1`
eval "$command"
