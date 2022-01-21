#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;
use TS_Params;

# GroupGenotypingKnownVariants.plと共に使用。
my $genotype_script = 'GroupGenotypingKnownVariants_v4.2.pl';

# この個体数ずつグループでタイピング。500個は動いたが、2496個はsamtoolsでエラー。
my $group_unit       = 100;
my $bam_dir          = $TS_Params::output[7];
my $column           = 3;
my $variants_list    = $TS_Params::variant_list;
my $individual_table = $TS_Params::individual_table;
my $help             = 0;
GetOptions(
	'num_samples=i' => \$group_unit,
	'column=i'      => \$column,
	'bam_dir=s'     => \$bam_dir,
	'table=s'       => \$individual_table,
	'variants=s'    => \$variants_list,
	'script=s'      => \$genotype_script,
	'help'          => \$help
);

my $n_unit_max = 500;
if     ($help) { &usage(); exit; }
unless ( -f $genotype_script ) {
	die "Script file not found: $genotype_script";
}
unless ( -f $variants_list ) {
	die "Variants list file not found: $variants_list";
}
if ( $group_unit > $n_unit_max ) {
	die "Samples per 1 group is up to $n_unit_max "
	  . "(Your selection: $group_unit)";
}

makeDirs( 0, 3, 10 );
$bam_dir =~ s/\/+\Z//;

my $job_script = 'GroupGenotyping.sh';
open( IND, $individual_table ) or die "Can't open file: $individual_table";
open( OUT, ">$job_script" )    or die;
my $count = 1;
while (<IND>) {
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	if ( $. == 1 ) {
		print STDERR "Selected column is $data[$column]\n";
		next;
	}
	unless ( $data[$column] =~ /Case|Control/ ) { next; }

	if ( -f "${bam_dir}/$data[0].bam" ) {
		if ( $count % $group_unit == 1 ) {
			open( LIST, ">$TS_Params::output[3]/${count}-.bamlist" )
			  or die "Can't open ${count}-.bamlist";
			print OUT "perl $genotype_script -v $variants_list"
			  . " -b $TS_Params::output[3]/${count}-.bamlist"
			  . " -o $TS_Params::output[10]/ || echo ${count}- >> Failure_GroupGenotyping.list;";
			print OUT "gzip -f $TS_Params::output[10]/Group_${count}-;\n";
		}
		print LIST "${bam_dir}/$data[0].bam\n";
		if ( $count % $group_unit == 0 ) {
			close(LIST);
		}
		$count++;
	}
	else {
		print STDERR "No expected bam: $data[0]\n";
	}
}
close(LIST);
close(IND);

open( GRID, ">ArrayJobs_$job_script" ) or die;
print GRID while (<DATA>);
close(GRID);

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -n|--num_samples: N of samples per 1 job"
	  . " (default: 100, max: $n_unit_max)\n";
	print STDOUT " -b|--bam_dir: Path to individual bam data dir"
	  . " (default: $TS_Params::output[7])\n";
	print STDOUT " -c|--column: Case/Control column in IndividualTable"
	  . " (default: 3)\n";
	print STDOUT " -t|--table: Path to IndividualTable"
	  . " (default: $TS_Params::individual_table)\n";
	print STDOUT " -s|--script: Path to GroupGenotypingKnownVariants.pl"
	  . " (default: GroupGenotypingKnownVariants_v4.2.pl)\n";
	print STDOUT " -v|--variants: Path to Variants list"
	  . " (default: $TS_Params::variant_list)\n";
	print STDOUT " -h|--help: Show this message\n";
}

__DATA__
#!/bin/bash
#$ -N Genotyping
#$ -l h_vmem=12G
#$ -l mem_free=12G
#$ -e ./Logfiles
#$ -o ./Logfiles
#$ -cwd
#$ -l h_rt=24:00:00
#$ -q all.q

#$-S /bin/bash
command=`head -$SGE_TASK_ID GroupGenotyping.sh | tail -1`
eval "$command"
