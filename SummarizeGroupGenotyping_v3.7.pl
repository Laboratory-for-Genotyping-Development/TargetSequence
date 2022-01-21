#!/usr/bin/perl

use strict;
use Getopt::Long;
use List::Util qw/sum/;
use TS_Params;

# Groupで出力したファイルの行数などを確認？
# IndividualTableの9列目でstatusがNTCは除く
# $th_cov以下はGenotypeに入れない
my $script       = 'HardyWeinberg_v1.0.R';
my $th_cov       = 20;
my $genotype_dir = $TS_Params::output[10];
my $variants_list;
my $sample_list;
my $gzip = 0;
my $help = 0;
GetOptions(
	'threshold=i' => \$th_cov,
	'dir=s'       => \$genotype_dir,
	'script=s'    => \$script,
	'variants=s'  => \$variants_list,
	'id_list=s'   => \$sample_list,
	'zip'         => \$gzip,
	'help'        => \$help
);
if     ($help)              { &usage(); exit; }
unless ( -f $script )       { die "$script not found. exit."; }
unless ( -d $genotype_dir ) { die "Directory not found: $genotype_dir"; }
$genotype_dir =~ s/\/+\Z//;

my %selected_variants;
if ( -f $variants_list ) {
	open( VAR, $variants_list ) or die;
	while (<VAR>) {
		if (/^#/) { next; }
		s/(\r?\n)\z//;
		my @variant = split(/\t/);
		$selected_variants{ join( '_', @variant[ 0 .. 4 ] ) }++;
	}
	close(VAR);
}

my %selected_samples;
if ( -f $sample_list ) {
	open( IDY, $sample_list ) or die;
	while (<IDY>) {
		if (/^#/) { next; }
		s/(\r?\n)\z//;
		my ($sample_id) = split(/\t/);
		$selected_samples{$sample_id}++;
	}
	close(IDY);
}

open( OUT, ">$TS_Params::genotype_list" ) or die;
print OUT join(
	"\t",
	(
		'ind',       'chr',  'pos',         'type',
		'ref',       'alt',  'total_depth', 'ref_depth',
		'alt_depth', 'freq', 'genotype'
	)
) . "\n";

my %type, my %genotype, my %vaf_het;
foreach my $file ( glob("${genotype_dir}/Group_*") ) {
	open( IN, ( $file =~ /\.gz\Z/ ) ? "zcat $file |" : $file )
	  or die "Can't open file: $file";
	while (<IN>) {
		if ( $. == 1 ) { next; }
		s/(\r?\n)\z//;
		my @data = split(/\t/);
		my $name = join( '_', ( @data[ 1 .. 5 ] ) );

		# サンプルのスキップ
		if ( ( %selected_samples > 0 )
			&& !exists( $selected_samples{ $data[0] } ) )
		{
			next;
		}

		# バリアントのスキップ
		if ( ( %selected_variants > 0 )
			&& !exists( $selected_variants{$name} ) )
		{
			next;
		}

		if ( $data[6] < $th_cov ) {
			print STDERR "${name}: $data[6]\n";
			next;
		}
		print OUT "$_\n";

		#$name_st = $name . '_st';
		$type{$name}++;
		if    ( $data[10] eq '0/0' ) { $genotype{$name}[0]++; }
		elsif ( $data[10] eq '0/1' ) {
			push( @{ $vaf_het{$name} }, $data[9] );

			#push( @{ $vaf_het{"${name}_st"} }, $data[11] );
			$genotype{$name}[1]++;
		}
		elsif ( $data[10] eq '1/1' ) { $genotype{$name}[2]++; }
	}
	close(IN);
}
close(OUT);
print STDERR "Finished Reading Group Genotyping List.\n";

open( SUM, ">$TS_Params::genotype_summary" ) or die;
print SUM join(
	"\t",
	(
		'chr', 'pos',     'type',            'ref',
		'alt', 'chr_pos', 'chr_pos_ref_alt', '0/0',
		'0/1', '1/1',     'VAF_hetero',      'freq0',
		'HWpvalue'
	)
) . "\n";
foreach my $name ( sort( keys(%type) ) ) {
	my @geno = exists( $genotype{$name} ) ? @{ $genotype{$name} } : ();
	for ( my $i = 0 ; $i < 3 ; $i++ ) {
		if ( $geno[$i] eq '' ) { $geno[$i] = 0; }
	}

	#my $name_st = "${name}_st";
	#my $strandbias =
	#  exists( $vaf_het{$name_st} )
	#  ? sum( @{ $vaf_het{$name_st} } ) / @{ $vaf_het{$name_st} }
	#  : 'NA';

	# ヘテロの平均頻度
	my $vaf_hetero =
	  exists( $vaf_het{$name} )
	  ? sprintf( "%.5f", sum( @{ $vaf_het{$name} } ) / @{ $vaf_het{$name} } )
	  : 'NA';
	( my $chrpos = $name ) =~ s/_/\t/g;
	my @chrpos          = split( /\t/, $chrpos );
	my $chr_pos         = join( '_', @chrpos[ 0, 1 ] );
	my $chr_pos_ref_alt = join( '_', @chrpos[ 0, 1, 3, 4 ] );
	print SUM join( "\t",
		$chrpos,     $chr_pos, $chr_pos_ref_alt, @geno[ 0 .. 2 ],
		$vaf_hetero, '',       '' )
	  . "\n";
}
close(SUM);
system( "$TS_Params::R --vanilla --slave"
	  . " < $script --args $TS_Params::genotype_summary" );
if ($gzip) { system("gzip $TS_Params::genotype_list;") }

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -t|--threshold: Depth threshold (default: 20)\n";
	print STDOUT " -d|--dir: Path to group genotyping data dir"
	  . " (default: $TS_Params::output[10])\n";
	print STDOUT " -s|--script: Path to HardyWeinberg.R"
	  . " (default: HardyWeinberg_v1.0.R)\n";
	print STDOUT " -v|--variant: Path to selected variant list (optional)\n";
	print STDOUT " -i|--id_list: Path to selected sample list (optional)\n";
	print STDOUT " -z|--zip: Compress $TS_Params::genotype_list or not\n";
	print STDOUT " -h|--help: Show this message\n";
}
