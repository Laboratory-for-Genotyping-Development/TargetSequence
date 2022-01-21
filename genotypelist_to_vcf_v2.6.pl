#!/usr/bin/perl

use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw/sum/;

#
my $genotype_list;
my $filtering_variants;
my $th_call_rate = 0.98;
my $th_depth     = 20;
my $out_target_list;
my @chr     = ();
my $no_call = 1;
my $help    = 0;
GetOptions(
	'genotype_list=s'   => \$genotype_list,
	'filter_variants=s' => \$filtering_variants,
	'th_call_rate=f'    => \$th_call_rate,
	'depth=i'           => \$th_depth,
	'out_of_targets=s'  => \$out_target_list,
	'chr=s{,}'          => \@chr,
	'multigeno!'        => \$no_call,
	'help'              => \$help
);
if ($help) { &usage(); exit; }

unless ( -f $genotype_list ) {
	die "IndividualGenotype list file not found: $genotype_list";
}

# Multipleでフィルターする方のバリアントを読み込む
my %FILTER_VARIANTS;
if ( -f $filtering_variants ) {
	open( FILTER, $filtering_variants )
	  or die "Can't open file: $filtering_variants";
	while (<FILTER>) {
		s/(\r?\n)\z//;
		my @data = split(/\t/);
		my $key  = join( '_', @data[ 0 .. 3 ] );
		$FILTER_VARIANTS{$key} = ();
	}
	close(FILTER);
}

# VCFに出力しないバリアントのリスト(e.g. ORFの外側)を読み込む
my %OUT_TARGET_VARIANTS;
if ( -f $out_target_list ) {
	open( ORF, $out_target_list ) or die "Can't open file $out_target_list";
	while (<ORF>) {
		if (/^#/) { next; }
		s/(\r?\n)\z//;
		my ( $key, @data ) = split(/\t/);
		$OUT_TARGET_VARIANTS{$key} = [@data];
	}
	close(ORF);
}

# Genotype.listの読み込み
my %SAMPLES, my %GENOTYPES;
open( IN,
	$genotype_list =~ /\.gz\Z/ ? "zcat $genotype_list |" : $genotype_list )
  or die "Can't open file: $genotype_list";
while (<IN>) {
	if ( $. == 1 ) { next; }
	s/(\r?\n)\z//;
	my (
		$idy,       $chr,  $pos,         $type,
		$ref,       $alt,  $total_depth, $ref_depth,
		$alt_depth, $freq, $genotype
	) = split(/\t/);

	# 指定した染色体番号以外は出力しない
	unless ( ( @chr == 0 ) or grep { $_ eq $chr } @chr ) {
		next;
	}

	( $idy, my @data ) = split(/\t/);
	my $chr_pos_ref_alt = join( "_", @data[ 0, 1, 3, 4 ] );
	my $chr_pos         = join( "_", @data[ 0, 1 ] );
	$GENOTYPES{$chr_pos_ref_alt}{$idy} = [ @data[ 5 .. $#data ] ];
	$SAMPLES{$idy}++;

	if ( exists( $FILTER_VARIANTS{$chr_pos_ref_alt} )
		and $genotype =~ /0[\/\|]1|1[\/\|]1/ )
	{
		$FILTER_VARIANTS{$chr_pos_ref_alt}{$idy} = $genotype;
	}
}
close(IN);

# サンプルとバリアントの一覧を配列に出力
my @SAMPLES  = sort { $a cmp $b } keys(%SAMPLES);
my @VARIANTS = keys(%GENOTYPES);

# Multiple Variantの確認
my %MultiLoci;
foreach my $chr_pos_ref_alt (@VARIANTS) {
	my @data    = split( /_/, $chr_pos_ref_alt );
	my $chr_pos = join( "_", @data[ 0, 1 ] );
	unless ( exists( $OUT_TARGET_VARIANTS{$chr_pos_ref_alt} ) ) {
		push( @{ $MultiLoci{$chr_pos} }, $chr_pos_ref_alt );
	}

	#$IDY_TABLE->{'CHR_POS'}{'VARIANT'}{$chr_pos}={$chr_pos_ref_alt};
	# カウント1 -> PASS
	# カウント2 -> Multiple Genotype Filter
	# カウント3 -> Multiple Variant Filter
}

# VCFヘッダーの作成
print STDOUT while (<DATA>);
print STDOUT join(
	"\t",
	(
		'#CHROM', 'POS',  'ID',     'REF', 'ALT', 'QUAL',
		'FILTER', 'INFO', 'FORMAT', @SAMPLES
	)
) . "\n";

foreach my $chr_pos_ref_alt (@VARIANTS) {
	my ( $CHR, $POS,  $REF,    $ALT )    = split( /_/, $chr_pos_ref_alt );
	my ( $ID,  $QUAL, $FILTER, $FORMAT ) = ( '.', '.', 'PASS', 'GT:AD:DP' );

	my @loci =
	  exists( $MultiLoci{"${CHR}_${POS}"} )
	  ? @{ $MultiLoci{"${CHR}_${POS}"} }
	  : ($chr_pos_ref_alt);
	my %MissingGenotype;
	my %INFO;

	## Multipleの場合は、MultiAllelicSiteのINFO情報を付ける
	if ( @loci > 1 ) {
		if ( exists( $FILTER_VARIANTS{$chr_pos_ref_alt} ) ) {
			$INFO{'MultiAllelicSite'} = 'Filter';
		}
		else {
			foreach my $variant (@loci) {
				if ( $variant ne $chr_pos_ref_alt ) {
					if ( exists( $FILTER_VARIANTS{$variant} ) ) {
						if ( scalar( $FILTER_VARIANTS{$variant} ) > 0 ) {
							%MissingGenotype = %{ $FILTER_VARIANTS{$variant} };
						}
						$INFO{'MultiAllelicSite'} = 'PASS';
					}
					else {
						## FilteringListに含まれていないMultiAllelicSiteは保留(KEEP)
						$INFO{'MultiAllelicSite'} = 'KEEP';
					}
				}
			}
		}
	}

	# INFO EDIT
	my @INFO;
	while ( my ( $key, $value ) = each(%INFO) ) {
		push( @INFO, "${key}=${value}" );
	}
	my $INFO = ( @INFO > 0 ) ? join( ';', @INFO ) : '.';

	my @FORMAT_PER_SAMPLES, my $call_count = 0;
	foreach my $sample (@SAMPLES) {
		my ( $DP, $ref_depth, $alt_depth, $freq, $genotype ) =
		  exists( $GENOTYPES{$chr_pos_ref_alt}{$sample} )
		  ? @{ $GENOTYPES{$chr_pos_ref_alt}{$sample} }
		  : ( 0, 0, 0, 0, './.' );
		my $GT =
		  ( exists( $MissingGenotype{$sample} ) and $no_call )
		  ? './.'
		  : $genotype;
		if ( $GT !~ /[01][\|\/][01]/ ) {
			$GT = './.';
		}    # No callの条件を厳しく
		my $AD =
		  ( $ref_depth or $alt_depth )
		  ? "$ref_depth,$alt_depth"
		  : '.';
		if ( $DP >= $th_depth ) { $call_count++; }

		## print STDERR $format_per_sample,"\n";
		## print "ここにSAMPLEのVCFをFORMATの仕様にしたがって記入する";
		push( @FORMAT_PER_SAMPLES, join( ':', ( $GT, $AD, $DP ) ) );
	}

	if ( ( $call_count / @SAMPLES ) < $th_call_rate ) {
		$FILTER = 'lowCallRate';
	}
	elsif ( exists( $OUT_TARGET_VARIANTS{$chr_pos_ref_alt} ) ) {
		$FILTER = 'ORF';
	}

	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
	## VCFの一行を出力する。
	print STDOUT join(
		"\t",
		(
			$CHR, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT,
			@FORMAT_PER_SAMPLES
		)
	) . "\n";
}

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -g|--genotype_list: "
	  . "Path to IndividualGenotype.list (required)\n";
	print STDOUT " -f|--filter_variants: "
	  . "Filtering multiple variants list file\n";
	print STDOUT " -t|--th_call_rate: "
	  . "Variant call rate threshold (default: 0.98)\n";
	print STDOUT " -d|--depth: Depth threshold (default: 20)\n";
	print STDOUT " -o|--out_of_targets: Out of target region variants list\n";
	print STDOUT " -nom|--multigeno: Call another mutiple genotypes\n";
	print STDOUT " -h|--help: Show this message\n";
}

__DATA__
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=lowCallRate,Description="Low call rate">
##FILTER=<ID=ORF,Description="Out of target region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
##contig=<ID=GL000207.1,length=4262,assembly=b37>
##contig=<ID=GL000226.1,length=15008,assembly=b37>
##contig=<ID=GL000229.1,length=19913,assembly=b37>
##contig=<ID=GL000231.1,length=27386,assembly=b37>
##contig=<ID=GL000210.1,length=27682,assembly=b37>
##contig=<ID=GL000239.1,length=33824,assembly=b37>
##contig=<ID=GL000235.1,length=34474,assembly=b37>
##contig=<ID=GL000201.1,length=36148,assembly=b37>
##contig=<ID=GL000247.1,length=36422,assembly=b37>
##contig=<ID=GL000245.1,length=36651,assembly=b37>
##contig=<ID=GL000197.1,length=37175,assembly=b37>
##contig=<ID=GL000203.1,length=37498,assembly=b37>
##contig=<ID=GL000246.1,length=38154,assembly=b37>
##contig=<ID=GL000249.1,length=38502,assembly=b37>
##contig=<ID=GL000196.1,length=38914,assembly=b37>
##contig=<ID=GL000248.1,length=39786,assembly=b37>
##contig=<ID=GL000244.1,length=39929,assembly=b37>
##contig=<ID=GL000238.1,length=39939,assembly=b37>
##contig=<ID=GL000202.1,length=40103,assembly=b37>
##contig=<ID=GL000234.1,length=40531,assembly=b37>
##contig=<ID=GL000232.1,length=40652,assembly=b37>
##contig=<ID=GL000206.1,length=41001,assembly=b37>
##contig=<ID=GL000240.1,length=41933,assembly=b37>
##contig=<ID=GL000236.1,length=41934,assembly=b37>
##contig=<ID=GL000241.1,length=42152,assembly=b37>
##contig=<ID=GL000243.1,length=43341,assembly=b37>
##contig=<ID=GL000242.1,length=43523,assembly=b37>
##contig=<ID=GL000230.1,length=43691,assembly=b37>
##contig=<ID=GL000237.1,length=45867,assembly=b37>
##contig=<ID=GL000233.1,length=45941,assembly=b37>
##contig=<ID=GL000204.1,length=81310,assembly=b37>
##contig=<ID=GL000198.1,length=90085,assembly=b37>
##contig=<ID=GL000208.1,length=92689,assembly=b37>
##contig=<ID=GL000191.1,length=106433,assembly=b37>
##contig=<ID=GL000227.1,length=128374,assembly=b37>
##contig=<ID=GL000228.1,length=129120,assembly=b37>
##contig=<ID=GL000214.1,length=137718,assembly=b37>
##contig=<ID=GL000221.1,length=155397,assembly=b37>
##contig=<ID=GL000209.1,length=159169,assembly=b37>
##contig=<ID=GL000218.1,length=161147,assembly=b37>
##contig=<ID=GL000220.1,length=161802,assembly=b37>
##contig=<ID=GL000213.1,length=164239,assembly=b37>
##contig=<ID=GL000211.1,length=166566,assembly=b37>
##contig=<ID=GL000199.1,length=169874,assembly=b37>
##contig=<ID=GL000217.1,length=172149,assembly=b37>
##contig=<ID=GL000216.1,length=172294,assembly=b37>
##contig=<ID=GL000215.1,length=172545,assembly=b37>
##contig=<ID=GL000205.1,length=174588,assembly=b37>
##contig=<ID=GL000219.1,length=179198,assembly=b37>
##contig=<ID=GL000224.1,length=179693,assembly=b37>
##contig=<ID=GL000223.1,length=180455,assembly=b37>
##contig=<ID=GL000195.1,length=182896,assembly=b37>
##contig=<ID=GL000212.1,length=186858,assembly=b37>
##contig=<ID=GL000222.1,length=186861,assembly=b37>
##contig=<ID=GL000200.1,length=187035,assembly=b37>
##contig=<ID=GL000193.1,length=189789,assembly=b37>
##contig=<ID=GL000194.1,length=191469,assembly=b37>
##contig=<ID=GL000225.1,length=211173,assembly=b37>
##contig=<ID=GL000192.1,length=547496,assembly=b37>
##INFO=<ID=MultiAllelicSite,Number=A,Type=String,Description="Multi Allelic Site information.">
