#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;
use Symbol;
use TS_Params;

# VAFの縦軸をDepthにして散布図で描画
# ReviseVariants.listを読み込んで、境界線の位置を調節できるバージョン
my $script = 'DrawVAFDepthPlot_v1.4.R';
my $variants;
my $max_depth        = 10;                             # y軸の最大値
my $column_number    = 3;
my $individual_table = $TS_Params::individual_table;
my $genotype_list    = $TS_Params::genotype_list;
my $revised          = 0;
my $enable_maf0      = 0;
my $help             = 0;
GetOptions(
	'script=s'   => \$script,
	'variants=s' => \$variants,
	'depth=i'    => \$max_depth,
	'column=i'   => \$column_number,
	'table=s'    => \$individual_table,
	'genotype=s' => \$genotype_list,
	'revised'    => \$revised,
	'maf0'       => \$enable_maf0,
	'help'       => \$help
);
if     ($help)        { &usage(); exit; }
unless ( -f $script ) { die "$script not found. exit."; }

makeDirs( 0, 3 );
my $figure_dir =
  ( $revised == 0 ) ? $TS_Params::output[12] : $TS_Params::output[14];
if ( $column_number != 3 ) {
	$figure_dir .= "_$column_number";
}
unless ( -d $figure_dir ) { mkdir($figure_dir); }

my %types;
my ( $basename, $dirname, $ext ) =
  fileparse( $genotype_list, ( '.list', '.list.gz' ) );
if ( -f $variants ) {
	open( LIST, $variants ) or die "Can't open file: $variants";
	while (<LIST>) {
		if ( ( $. == 1 ) or /^#/ ) { next; }
		s/(\r?\n)\z//;
		my @data = split(/\t/);
		my $key  = join( '_', @data[ 0 .. 4 ] );
		if ( scalar(@data) >= 6 ) {
			@{ $types{$key} } = @data[ 5, 6 ];
		}
		else {
			@{ $types{$key} } = ( '', '' );
		}
	}
	close(LIST);
}
else {
	my $genotype_summary = "${dirname}${basename}.summary";
	open( IN,   $genotype_summary ) or die "Can't open file: $genotype_summary";
	open( MAF0, '>MAF0.list' )      or die;
	while (<IN>) {
		if ( ( $. == 1 ) or /^#/ ) { next; }
		s/(\r?\n)\z//;
		my (
			$chr,    $pos,     $type,            $ref,
			$alt,    $chr_pos, $chr_pos_ref_alt, $homo0,
			$hetero, $homo1
		) = split(/\t/);

		my $name = join( '_', ( $chr, $pos, $type, $ref, $alt ) );
		if (   ( $homo0 + $hetero ) == 0
			or ( $hetero + $homo1 ) == 0 )
		{
			print MAF0 "${name}\n";
			unless ($enable_maf0) { next; }
		}    # 頻度が0のものは除外
		@{ $types{$name} } = ( '', '' );
	}
	close(MAF0);
	close(IN);
}
my $VarintNumber = keys(%types);
print "Variant number is ${VarintNumber}\n";

my %case_ctrl;
open( IN, $individual_table ) or die "Can't open file: $individual_table";
while (<IN>) {
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	if ( $. == 1 ) {
		print "Selected column is $data[$column_number]\n";
		next;
	}
	$case_ctrl{ $data[0] } = $data[$column_number];
}
close(IN);

# ArrayJobsで使うファイルを作成
my %file_handles;
open( IN, ( $ext =~ /\.gz\Z/ ) ? "zcat $genotype_list |" : $genotype_list )
  or die "Can't open file $genotype_list";
while (<IN>) {
	if ( ( $. == 1 ) or /^#/ ) { next; }
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	my $name = join( '_', @data[ 1 .. 5 ] );
	unless ( $case_ctrl{ $data[0] } =~ /Case|Control/
		and exists( $types{$name} ) )
	{
		##print "$data[0]\n";
		next;
	}

	# IN/DELが長すぎるとエラーになるので、その場合はファイル名を短縮
	if ( length($name) > 246 ) {
		$name = substr( $name, 0, 246 );
	}

	unless ( exists( $file_handles{$name} ) ) {
		$file_handles{$name} = Symbol::gensym();
		open( $file_handles{$name}, ">$TS_Params::output[3]/Genotype_$name" )
		  or die;
	}

	my $fh = $file_handles{$name};
	if ( $data[$#data] =~ /[01]\/[01]/ ) { print $fh "$_\n"; }
}
foreach my $fh ( values(%file_handles) ) { close($fh); }

my $job_script = 'DrawVAFDepthPlot.sh';
open( OUT, ">$job_script" ) or die;
my $count = 0;
foreach my $type ( sort( keys(%types) ) ) {
	$count++;
	my $number = sprintf( "%05d", $count );
	my @type   = split( '_', $type );
	if ( ( $count % 100 == 1 ) and ( $count > 1 ) ) {
		print OUT "\n";
	}
	elsif ( $count > 1 ) {
		print OUT '; ';
	}

	my ( $low, $high ) = @{ $types{$type} };
	if ( $low eq '' and $high eq '' ) {
		print OUT "$TS_Params::R --vanilla --slave <"
		  . " $script --args $type $number $figure_dir $max_depth";
	}
	else {
		if ( $low eq '' )  { $low  = 0.2; }
		if ( $high eq '' ) { $high = 0.8; }
		print OUT "$TS_Params::R --vanilla --slave <"
		  . " $script --args $type $number $figure_dir $max_depth $low $high";
	}
}
print OUT "\n";
close(OUT);

open( GRID, ">ArrayJobs_$job_script" ) or die;
print GRID while (<DATA>);
close(GRID);

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -s|--script: Path to DrawVAFDepthPlot.R"
	  . " (default: DrawVAFDepthPlot_v1.4.R)\n";
	print STDOUT " -v|--variants: Path to Variants list\n";
	print STDOUT " -d|--depth: Max of Y axis (default: 10)\n";
	print STDOUT " -c|--column: Case/Control column in IndividualTable"
	  . " (default: 3)\n";
	print STDOUT " -t|--table: Path to IndividualTable"
	  . " (default: $TS_Params::individual_table)\n";
	print STDOUT " -g|--genotype: Path to IndividualGenotype.list"
	  . " (default: $TS_Params::genotype_list)\n";
	print STDOUT " -r|--revised: If you want to draw revised VAF figures,"
	  . " spcify this option\n";
	print STDOUT " -m|--maf0: If you want to draw MAF0 VAF figures,"
	  . " spcify this option\n";
	print STDOUT " -h|--help: Show this message\n";
}
__DATA__
#!/bin/bash
#$ -N DrawVAFDepthPlot
#$ -l h_vmem=4G
#$ -l mem_free=4G
#$ -e ./Logfiles
#$ -o ./Logfiles
#$ -cwd
#$ -l h_rt=10:00:00
#$ -q all.q

#$-S /bin/bash
command=`head -$SGE_TASK_ID DrawVAFDepthPlot.sh | tail -1`
eval "$command"
