#!/usr/bin/perl

use strict;
use Getopt::Long;
use TS_Params;

# v2.0から、MakingArrayJobs_shortPCRseq_SampleBamCall_v5.1.plに取り込んで受けることに。ID、fasta、BED
# v2.1: GATK3対応
my $name, my $ref_fa, my $bed;
my $help = 0;
GetOptions(
	'name=s' => \$name,
	'ref=s'  => \$ref_fa,
	'bed=s'  => \$bed,
	'help'   => \$help
);
if     ($help)                    { &usage(); exit; }
unless ( -f $ref_fa and -f $bed ) { die; }

open( ONTARGET, ">$TS_Params::output[3]/OnTarget_$name" ) or die;
my $trimmed_fastq   = "$TS_Params::output[2]/${name}_F.fastq";
my $untrimmed_fastq = "$TS_Params::output[2]/${name}_F.Untrimmed.fastq.gz";
my $trimmed_reads   = 0, my $untrimmed_reads = 0;
if ( -f $trimmed_fastq ) {
	($trimmed_reads) = split( /\s+/, `wc -l $trimmed_fastq` );
}
elsif ( -f "${trimmed_fastq}.gz" ) {
	($trimmed_reads) = split( /\s+/, `zcat $trimmed_fastq | wc -l` );
}
if ( -f $untrimmed_fastq ) {
	($untrimmed_reads) = split( /\s+/, `zcat $untrimmed_fastq | wc -l` );
}
$trimmed_reads   /= 2;
$untrimmed_reads /= 2;
my $total_reads = $trimmed_reads + $untrimmed_reads;    # トータルリード

if ( -f "$TS_Params::output[4]/${name}.bam" ) {
	my $command =
		"java -Xmx8g -jar $TS_Params::GATK -T CountReads -R $ref_fa -L $bed"
	  . " -I $TS_Params::output[4]/${name}.bam 2>&1 | grep CountReads | cut -d "
	  . q(' ')
	  . ' -f 8 | tail -n 1';
	my $OnTarget = `$command`;
	chomp $OnTarget;

	# -Tとか出力される時はおそらくリード0なので、OnTargetReadも0にする。
	if ( $OnTarget eq '-T' ) { $OnTarget = 0; }

	my $PropTrimmed =
	  ( $total_reads == 0 )
	  ? 0
	  : sprintf( "%.2f", $trimmed_reads / $total_reads * 100 );
	my $PropOnTarget =
	  ( $total_reads == 0 )
	  ? 0
	  : sprintf( "%.2f", $OnTarget / $total_reads * 100 );

	print ONTARGET join(
		"\t",
		(
			$name,        $total_reads, $trimmed_reads,
			$PropTrimmed, $OnTarget,    $PropOnTarget
		)
	) . "\n";
}
close(ONTARGET);

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -n|--name: Sample name (required)\n";
	print STDOUT " -r|--ref: Path to reference fasta file (required)\n";
	print STDOUT " -b|--bed: Path to BED file (required)\n";
	print STDOUT " -h|--help: Show this message\n";
}
