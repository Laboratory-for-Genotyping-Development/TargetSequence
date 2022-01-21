#!/usr/bin/perl

use strict;
use Getopt::Long;
use List::Util qw/sum/;
use TS_Params;

# 各個体のサマリーにはNTC/NAも出すが、各塩基のサマリーからはNTC/NAを外す。
my $individual_table = $TS_Params::individual_table;
my $column_number    = 3;
my $cov_dir          = $TS_Params::output[9];
my $bed_file;
my $orf       = 0;
my $threshold = 20;    # これ未満のcovはノーデータとする
my $amel      = 0;
my $help      = 0;
my $plate;
GetOptions(
	'individual=s' => \$individual_table,
	'orf=i'        => \$orf,
	'dir=s'        => \$cov_dir,
	'column=i'     => \$column_number,
	'bed=s'        => \$bed_file,
	'th_depth=i'   => \$threshold,
	'plate=s'      => \$plate,
	'amel'         => \$amel,
	'help'         => \$help
);
if ($help) { &usage(); exit; }

unless ( -d $cov_dir ) { die "Coverage output dir not found: $cov_dir"; }

my %bed_target;
if ( defined($bed_file) ) {
	if ( -f $bed_file ) {
		&set_target_region( $bed_file, $orf, $amel, \%bed_target );
	}
	else { die "BED file not found: $bed_file"; }
}

$cov_dir =~ s/\/+\Z//;
my $basename =
  defined($plate)
  ? "IndividualCoverage_${plate}_${column_number}"
  : "IndividualCoverage_${column_number}";

open( OUT, ">${basename}.all" ) or die;
print OUT join( "\t", ( 'IDY', 'chr', 'pos', 'reads' ) ) . "\n";
open( SUM, ">${basename}.summary" ) or die;
print SUM join(
	"\t",
	(
		'IDY',        'Plate_WellNo', 'Plate',    'Class',
		'N_Samples',  'Target_bp',    'Total_bp', 'Ave_depth',
		'Covered_bp', 'Covered_%'
	)
) . "\n";

my $cnt = 0;
my @chrpos_list, my %chk_chrpos;
my %avg, my %enough, my %count;
my $case = 0, my $ctrl = 0;
open( IN, $individual_table ) or die "Can't open file: $individual_table";
while (<IN>) {
	s/(\r?\n)\z//;
	my @data   = split(/\t/);
	my $status = $data[$column_number];
	if ( $. == 1 ) {
		print STDERR "Selected column is ${status}\n";
		next;
	}

	if ( defined($plate) and ( $plate ne $data[2] ) ) { next; }

	if    ( $status =~ /Case/i )    { $case++; }
	elsif ( $status =~ /Control/i ) { $ctrl++; }

	if ( keys(%bed_target) == 0 ) {
		&set_target_region( $data[4], $orf, $amel, \%bed_target );
	}

	my $filename = "${cov_dir}/$data[0].gz";
	unless ( -f $filename ) {
		print STDERR "No expected file: $data[0]\n";
		next;
	}

	my $target = 0, my $total = 0, my $covered = 0;
	my @column, my $N_samples;
	open( COV, "zcat $filename |" ) or die "Can't open file $filename";
	while (<COV>) {
		s/(\r?\n)\z//;
		if ( $. == 1 ) {
			@column    = split(/\t/);
			$N_samples = @column - 3;
		}
		my ( $chr, $position, $reads ) =
		  split(/\t|:/);    # それぞれchr, pos, リード数
		my $chrpos = "${chr}:${position}";
		unless ( exists( $bed_target{$chrpos} ) ) { next; }
		$target++;
		$total += $reads;

		unless ( exists( $chk_chrpos{$chrpos} ) ) {
			push( @chrpos_list, $chrpos );
			$chk_chrpos{$chrpos} = 1;
		}
		if ( $status =~ /Case|Control/i ) {
			print OUT join( "\t", ( $data[0], $chr, $position, $reads ) )
			  . "\n";
			$avg{'all'}{$chrpos} =
			  ( $avg{'all'}{$chrpos} * $count{'all'}{$chrpos} + $reads ) /
			  ( $count{'all'}{$chrpos} + 1 );
			$count{'all'}{$chrpos}++;
			if ( $reads >= $threshold ) { $enough{'all'}{$chrpos}++; }
			if ( $status =~ /Case/i ) {
				$avg{'case'}{$chrpos} =
					( $avg{'case'}{$chrpos} * $count{'case'}{$chrpos} + $reads )
				  / ( $count{'case'}{$chrpos} + 1 );
				$count{'case'}{$chrpos}++;
				if ( $reads >= $threshold ) { $enough{'case'}{$chrpos}++; }
			}
			if ( $status =~ /Control/i ) {
				$avg{'ctrl'}{$chrpos} =
					( $avg{'ctrl'}{$chrpos} * $count{'ctrl'}{$chrpos} + $reads )
				  / ( $count{'ctrl'}{$chrpos} + 1 );
				$count{'ctrl'}{$chrpos}++;
				if ( $reads >= $threshold ) { $enough{'ctrl'}{$chrpos}++; }
			}
		}
		if ( $reads >= $threshold ) { $covered++; }
	}
	close(COV);

	if ( $target <= 0 ) { print STDERR "No coverage info: $data[0]\n"; next; }
	my $average       = sprintf( "%.1f", $total / $target );
	my $covered_ratio = sprintf( "%.2f", $covered / $target * 100 );
	print SUM join(
		"\t",
		(
			$data[0],   $data[1], $data[2], $status,
			$N_samples, $target,  $total,   $average,
			$covered,   $covered_ratio
		)
	) . "\n";
}
close(IN);
close(SUM);
close(OUT);

my @target_header = ( 'Chr', 'Pos', 'Ave_reads', 'N_ind_enough', 'Cover_Rate' );
open( TAR, ">${basename}.target" ) or die;
print TAR join( "\t", @target_header ) . "\n";

if ( $case != 0 and $ctrl != 0 ) {
	open( CASE, ">${basename}_Case.target" )
	  or die;
	open( CTRL, ">${basename}_Control.target" )
	  or die;
	print CASE join( "\t", @target_header ) . "\n";
	print CTRL join( "\t", @target_header ) . "\n";
}

foreach my $chrpos (@chrpos_list) {
	my ( $chr, $pos ) = split( /:/, $chrpos );
	unless ( exists( $enough{'all'}{$chrpos} ) ) {
		$enough{'all'}{$chrpos} = 0;
	}
	my $cover_rate = sprintf( "%.5f",
		( $case + $ctrl == 0 )
		? 0
		: $enough{'all'}{$chrpos} / ( $case + $ctrl ) );
	print TAR join(
		"\t",
		(
			$chr,                 $pos,
			$avg{'all'}{$chrpos}, $enough{'all'}{$chrpos},
			$cover_rate
		)
	) . "\n";

	if ( $case != 0 and $ctrl != 0 ) {
		unless ( exists( $enough{'case'}{$chrpos} ) ) {
			$enough{'case'}{$chrpos} = 0;
		}
		unless ( exists( $enough{'ctrl'}{$chrpos} ) ) {
			$enough{'ctrl'}{$chrpos} = 0;
		}
		$cover_rate = sprintf( "%.5f",
			( $case == 0 ) ? 0 : $enough{'case'}{$chrpos} / $case );
		print CASE join(
			"\t",
			(
				$chr,                  $pos,
				$avg{'case'}{$chrpos}, $enough{'case'}{$chrpos},
				$cover_rate
			)
		) . "\n";
		$cover_rate = sprintf( "%.5f",
			( $ctrl == 0 ) ? 0 : $enough{'ctrl'}{$chrpos} / $ctrl );
		print CTRL join(
			"\t",
			(
				$chr,                  $pos,
				$avg{'ctrl'}{$chrpos}, $enough{'ctrl'}{$chrpos},
				$cover_rate
			)
		) . "\n";
	}
}
if ( $case != 0 and $ctrl != 0 ) {
	close(CTRL);
	close(CASE);
}
close(TAR);
system("gzip -f ${basename}.all");

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -i|--individual: Path to IndividualTable"
	  . " (default: $TS_Params::individual_table)\n";
	print STDOUT " -o|--orf: N of cutoff bp (default: 0)\n";
	print STDOUT " -d|--dir: Path to individual coverage data dir"
	  . " (default: $TS_Params::output[9])\n";
	print STDOUT " -c|--column: Case/Control column number in IndividualTable"
	  . " (default: 3)\n";
	print STDOUT " -b|--bed: Path to another BED file (optional)\n";
	print STDOUT " -t|--threshold: Depth threshold (default: 20)\n";
	print STDOUT " -p|--plate: Plate name (if specified)\n";
	print STDOUT " -a|--amel: Include AMEL gene\n";
	print STDOUT " -h|--help: Show this message\n";
}

sub set_target_region() {
	my ( $file, $buf, $flg, $region ) = @_;
	open( BED, $file ) or die "Can't open file $file";
	while (<BED>) {
		if (/^#+/) { next; }
		s/(\r?\n)\z//;
		my ( $chr, $start, $end, $name ) = split(/\t/);

		# AMELは原則計算から除外
		if ( $flg == 0 and $name =~ /AMEL/ ) { next; }

		# ORFがある場合は指定されたbp分を削除する。
		for ( my $i = $start + $buf ; $i < $end - $buf ; $i++ ) {
			my $pos = $i + 1;
			$$region{"${chr}:${pos}"} = $name;
		}
	}
	close(BED);
	print STDERR "BED file: $file (" . keys(%$region) . "bp).\n";
}
