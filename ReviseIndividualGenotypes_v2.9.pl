#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;
use List::Util qw/sum/;
use TS_Params;

# VAFを目で確認後、ReviseVariants.list (ヘッダーあり)に基づき、遺伝子型を改訂。
my $script = 'HardyWeinberg_v1.0.R';

my $freq_low_set  = 0.2;
my $freq_high_set = 0.8;
my $buffer        = 0.05;    # 遺伝子型決定の際、NAにするbuffer

my $variants      = 'ReviseVariants.list';
my $table         = $TS_Params::individual_table;
my $genotype_list = $TS_Params::genotype_list;
my $column_number = 3;
my $gzip          = 0;
my $help          = 0;
GetOptions(
	'genotype_list=s' => \$genotype_list,
	'variants=s'      => \$variants,
	'column=i'        => \$column_number,
	'table=s'         => \$table,
	'script=s'        => \$script,
	'buffer=f'        => \$buffer,
	'zip'             => \$gzip,
	'help'            => \$help
);
if     ($help)        { &usage(); exit; }
unless ( -f $script ) { die "R script not found: $script. exit."; }

my $head, my %ind, my %ctrl;
open( SAMPLE, $table ) or die "Can't open file: $table";
while (<SAMPLE>) {
	s/(\r?\n)\z//;
	my @line = split(/\t/);
	if ( $. == 1 ) {
		$head = $line[$column_number];
	}
	elsif ( $line[$column_number] =~ /Case|Control/ ) {
		if ( $line[$column_number] eq 'Control' ) {
			$ctrl{ $line[0] }++;
		}
		$ind{ $line[0] }++;
	}
}
close(SAMPLE);

my $totalnumber = keys(%ind);
my $ctrlnumber  = keys(%ctrl);
print STDERR "Selected column is ${head}\n";
print STDERR "Total number is ${totalnumber}\n";
print STDERR "Control number is ${ctrlnumber}\n";

my %low, my %high;
open( VARIANT, $variants ) or die "Can't open file: $variants";
while (<VARIANT>) {
	if ( $. == 1 ) { next; }
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	my $name = join( '_', @data[ 0 .. 4 ] );
	$low{$name}  = ( $data[5] > 0 ) ? $data[5] : $freq_low_set;
	$high{$name} = ( $data[6] > 0 ) ? $data[6] : $freq_high_set;
}
close(VARIANT);

my %genotype, my %geno_ctrl, my %vaf_het;
my $basename = basename( $variants, '.list' );
my $revised_list =
  "RevisedIndividualGenotype_${basename}_${column_number}_$head";
my $input =
  $genotype_list =~ /\.gz\Z/ ? "zcat $genotype_list |" : $genotype_list;
open( IN,  $input )                  or die "Can't open file: $genotype_list";
open( OUT, ">${revised_list}.list" ) or die;

while (<IN>) {
	if ( $. == 1 ) { print OUT; }
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	unless ( exists( $ind{ $data[0] } ) ) { next; }    # 個体レベルの削除
	my $name = join( '_', ( @data[ 1 .. 5 ] ) );
	unless ( exists( $low{$name} ) ) { next; }         # 変異レベルの削除

	if ( $data[9] <= $low{$name} - $buffer ) {
		if ( $data[10] =~ /0\/[01]/ ) {
			$data[10] = '0/0';
			$genotype{$name}[0]++;
			if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[0]++; }
		}
		else {
			$genotype{$name}[4]++;
			if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[4]++; }
		}
	}
	elsif ( ( $low{$name} + $buffer <= $data[9] )
		and ( $data[9] <= $high{$name} - $buffer ) )
	{
		if ( $data[10] =~ /[01]\/[01]/ ) {
			$data[10] = '0/1';
			push( @{ $vaf_het{$name} }, $data[9] );

			#push( @{ $vaf_het{"${name}_st"} }, $data[11] );
			$genotype{$name}[1]++;
			if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[1]++; }
		}
		else {
			$genotype{$name}[4]++;
			if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[4]++; }
		}
	}
	elsif ( $high{$name} + $buffer <= $data[9] ) {
		if ( $data[10] =~ /[01]\/1/ ) {
			$data[10] = '1/1';
			$genotype{$name}[2]++;
			if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[2]++; }
		}
		else {
			$genotype{$name}[4]++;
			if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[4]++; }
		}
	}
	else {
		$data[10] = './.';
		$genotype{$name}[3]++;
		if ( exists( $ctrl{ $data[0] } ) ) { $geno_ctrl{$name}[3]++; }
	}

	print OUT join( "\t", @data ) . "\n";
}
close(OUT);
close(IN);

&output_summary( "${revised_list}.summary", $totalnumber, %genotype );
if ( $ctrlnumber > 0 ) {
	&output_summary( "${revised_list}_ctrl.summary", $ctrlnumber, %geno_ctrl );
}
if ($gzip) { system("gzip ${revised_list}.list;") }

# summaryの出力
sub output_summary {
	my ( $output_filename, $N_individuals, %GT ) = @_;
	open( SUM, ">$output_filename" ) or die;
	print SUM join(
		"\t",
		(
			'chr',             'pos',
			'type',            'ref',
			'alt',             'chr_pos',
			'chr_pos_ref_alt', '0/0',
			'0/1',             '1/1',
			'./.',             'other',
			'call_rate(%)',    'cover_rate(%)',
			'VAF_hetero',      'freq0',
			'HWpvalue'
		)
	) . "\n";
	foreach my $name ( sort( keys(%low) ) ) {
		my @geno = exists( $GT{$name} ) ? @{ $GT{$name} } : ();
		for ( my $i = 0 ; $i < 5 ; $i++ ) {
			unless ( defined( $geno[$i] ) ) { $geno[$i] = 0; }
		}

		#my $name_st = "${name}_st";
		#my $strandbias =
		#  exists( $vaf_het{$name_st} )
		#  ? sum( @{ $vaf_het{$name_st} } ) / @{ $vaf_het{$name_st} }
		#  : 'NA';

		# ヘテロの平均頻度
		my $vaf_hetero =
		  exists( $vaf_het{$name} )
		  ? sprintf( "%.5f",
			sum( @{ $vaf_het{$name} } ) / @{ $vaf_het{$name} } )
		  : 'NA';
		( my $chrpos = $name ) =~ s/_/\t/g;
		my @chrpos          = split( /\t/, $chrpos );
		my $chr_pos         = join( '_', @chrpos[ 0, 1 ] );
		my $chr_pos_ref_alt = join( '_', @chrpos[ 0, 1, 3, 4 ] );
		my $call_rate =
		  sprintf( "%.3f", sum( @geno[ 0 .. 2 ] ) / $N_individuals * 100 );
		my $cover_rate = sprintf( "%.3f", sum(@geno) / $N_individuals * 100 );
		print SUM join( "\t",
			$chrpos,     $chr_pos,   $chr_pos_ref_alt,
			@geno,       $call_rate, $cover_rate,
			$vaf_hetero, '',         '' )
		  . "\n";
	}
	close(SUM);
	system( "$TS_Params::R --vanilla --slave"
		  . " < $script --args $output_filename" );
}

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -g|--genotype_list: Path to IndividualGenotype.list"
	  . " (default: $TS_Params::genotype_list)\n";
	print STDOUT " -v|--variants: Path to call variants list"
	  . " (default: ReviseVariants.list)\n";
	print STDOUT " -c|--column: Case/Control column number in IndividualTable"
	  . " (default: 3)\n";
	print STDOUT " -t|--table: Path to IndividualTable"
	  . " (default: $TS_Params::individual_table)\n";
	print STDOUT " -s|--script: Path to HardyWeinberg.R"
	  . " (default: HardyWeinberg_v1.0.R)\n";
	print STDOUT " -b|--buffer: No call frequency buffer (default: 0.05)\n";
	print STDOUT " -z|--zip: Compress RevisedGenotype.list or not\n";
	print STDOUT " -h|--help: Show this message\n";
}
