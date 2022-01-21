#!/usr/bin/perl

use strict;
use Getopt::Long;
use TS_Params;

# perl ListUpVariants1.7.pl SampleTable 4.Call # 見に行くテーブル 見に行くフォルダ
my $table         = $TS_Params::individual_table;
my $call_dir      = $TS_Params::output[8];
my $column_number = 3;
my $bed_file;
my $th_frq = 0.3;    # Singleton variantのフィルタリング関連パラメーター
my $th_fs  = 100;
my $help   = 0;
GetOptions(
	'individual=s'  => \$table,
	'dir=s'         => \$call_dir,
	'column=i'      => \$column_number,
	'bed=s'         => \$bed_file,
	'allele_freq=f' => \$th_frq,
	'fs=i'          => \$th_fs,
	'help'          => \$help
);

if ($help) { &usage(); exit; }
$call_dir =~ s/\/+\Z//;

unless ( -f $table )    { die "File not found: Individual(or Sample)Table"; }
unless ( -d $call_dir ) { die "Call dir not found: $call_dir"; }
if     ( defined($bed_file) and !-f $bed_file ) {
	die "BED file not found: $bed_file";
}

# IDYを別変数に入れるので、カラム番号を1引く
$column_number--;

# UGのフィルタリング関連パラメーター
my $col_fs = 'FS';

my $hc_ext = '.HC.vcf';
my $ug_ext = '.UG.vcf';

my %multiple, my %bi_allele, my %filtered;
my $count = 1;
open( TABLE, $table ) or die "Can't open $table";
while (<TABLE>) {
	s/(\r?\n)\z//;
	my ( $idy, @data ) = split(/\t/);
	if ( $. == 1 ) {
		print STDERR "Selected column is $data[$column_number]\n";
		next;
	}

	unless ( defined($bed_file) ) { $bed_file = $data[3]; }
	my $status = $data[$column_number];
	unless ( $status =~ /Case|Control/ ) { next; }

	# HaplotypeCaller
	unless ( -f "${call_dir}/${idy}$hc_ext" ) {
		print "${call_dir}/$idy is not called.\n";
		next;
	}

	open( HC, "${call_dir}/${idy}$hc_ext" ) or die "Can't open ${idy}$hc_ext";
	while (<HC>) {
		if ( $_ =~ /^#/ ) { next; }
		s/(\r?\n)\z//;
		my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $fmt, @id ) =
		  split(/\t/);
		my $type = &get_variant_type( $ref, $alt );
		if ( $type eq '' ) {
			print STDERR "Unexpected variant type at ${chr}:$pos (${ref}, $alt)\n";
			next;
		}
		my $name = join( '_', ( $chr, $pos, $type, $ref, $alt ) );
		if ( $alt =~ /\,/ ) {
			$name = join(
				'_',
				(
					$chr, $pos, 'GATK HaplotypeCaller multipleVariants',
					$ref, $alt
				)
			);
			$multiple{$name}++;
			next;
		}
		else {
			$bi_allele{$name}{$idy} = &calc_allele_freq( $fmt, @id );
		}
	}
	close(HC);

	# UnifiedGenotyper
	open( UG, "${call_dir}/${idy}$ug_ext" ) or die "Can't open ${idy}$ug_ext";
	while (<UG>) {
		if ( $_ =~ /^#/ ) { next; }
		s/(\r?\n)\z//;
		my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $fmt, @id ) =
		  split(/\t/);

		# FSの値をチェック
		my $filter_val = 0;
		my $type       = &get_variant_type( $ref, $alt );
		if ( $type eq '' ) {
			print "Unexpected variant type at ${chr}:$pos (${ref}, $alt)\n";
			next;
		}
		my $name = join( '_', ( $chr, $pos, $type, $ref, $alt ) );
		foreach my $data ( split( /;/, $info ) ) {
			my ( $key, $val ) = split( /=/, $data );
			if ( $key eq $col_fs ) {
				if ( $val > $th_fs ) {
					$filtered{$name} = "FS=$val";
					$filter_val = 1;
					last;
				}
			}
		}
		if    ( $filter_val == 1 ) { next; }
		elsif ( $alt =~ /\,/ ) {
			$name = join(
				'_',
				(
					$chr, $pos, 'GATK UnifiedGenotyper multipleVariants',
					$ref, $alt
				)
			);
			$multiple{$name}++;
			next;
		}
		else {
			my $frq = &calc_allele_freq( $fmt, @id );

			# 頻度が高い方を採用する
			if ( !exists( $bi_allele{$name}{$idy} )
				or $bi_allele{$name}{$idy} < $frq )
			{
				$bi_allele{$name}{$idy} = $frq;
			}
		}
	}
	close(UG);

	# VCMM (if exists)
	if ( -f "${call_dir}/$idy" ) {
		open( VCMM, "${call_dir}/$idy" ) or die "Can't open $idy";
		while (<VCMM>) {
			if ( $. == 1 ) { next; }
			s/(\r?\n)\z//;
			my @data = split(/\t/);
			my $key  = join( '_', @data[ 0 .. 4 ] );

			# VCMMは、InDelのみをチェックする。
			unless ( $data[2] =~ /SNP/ ) {
				if ( $data[2] =~ /multiple/ ) {
					$multiple{$key}++;
				}
				elsif ( !exists( $bi_allele{$key}{$idy} ) ) {
					$bi_allele{$key}{$idy} = $data[$#data];
				}
			}
		}
		close(VCMM);
	}
	if ( ++$count % 100 == 0 ) { print "$count\n"; }
}
close(TABLE);

print STDERR "Finished reading calls\n";

open( MULTI, '>MultiplexVariants.list' ) or die;
print MULTI join( "\t", ( '#chr', 'pos', 'type', 'ref', 'alt' ) ) . "\n";
foreach my $name ( sort( keys(%multiple) ) ) {
	$name =~ s/_/\t/g;
	print MULTI "$name\n";
}
close(MULTI);

# VCF形式にして、bedでフィルター
my $vcf = 'Variants.temp';
open( TEMP, ">${vcf}.vcf" ) or die;
print TEMP while (<DATA>);
foreach my $name ( sort( keys(%bi_allele) ) ) {
	my @idy = keys( %{ $bi_allele{$name} } );

	# Singletonで頻度が閾値未満の変異は除外へ
	if ( scalar(@idy) == 1 and $bi_allele{$name}{ $idy[0] } < $th_frq ) {
		$filtered{$name} = $bi_allele{$name}{ $idy[0] };
	}
	else {
		$name =~ s/_/\t/g;
		print TEMP "$name\t.\t.\t.\t.\n";
	}
}
close(TEMP);

open( FILTER, '>FilteredVariants.list' ) or die;
print FILTER join( "\t", ( '#chr', 'pos', 'type', 'ref', 'alt', 'value' ) )
  . "\n";
foreach my $name ( sort( keys(%filtered) ) ) {
	my $frq = $filtered{$name};
	if ( $frq =~ /^FS/ ) {
		if ( exists( $bi_allele{$name} ) ) { next; }
	}
	else { $frq = 'Freq=' . sprintf( '%.5f', $frq ); }
	$name =~ s/_/\t/g;
	print FILTER "$name\t$frq\n";
}
close(FILTER);

open( VCF,
	"$TS_Params::vcftools --recode --vcf ${vcf}.vcf --bed $bed_file --stdout |"
) or die;
open( OUT, ">$TS_Params::variant_list" ) or die;
print OUT join( "\t", ( '#chr', 'pos', 'type', 'ref', 'alt' ) ) . "\n";
while (<VCF>) {
	if ( $_ =~ /^#/ ) { next; }
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	print OUT join( "\t", @data[ 0 .. 4 ] ) . "\n";
}
close(OUT);
close(VCF);
unlink( glob("${vcf}.*") );

# Alleleの長さからSNP/In/Delを識別する
sub get_variant_type {
	my ( $ref, $alt ) = @_;
	my $type = '';
	if ( length($ref) == 1 and length($alt) == 1 ) {
		$type = 'SNP';
	}
	elsif ( length($ref) == 1 and length($alt) > 1 ) {
		$type = 'IN';
	}
	elsif ( length($ref) > 1 and length($alt) == 1 ) {
		$type = 'DEL';
	}
	elsif ( $ref =~ /\,/ or $alt =~ /\,/ ) {
		$type = 'MULTIPLE';
	}
	return ($type);
}

# AD, DPの値から頻度を計算
sub calc_allele_freq {
	my ( $fmt, @dataset ) = @_;
	my @format = split( /:/, $fmt );
	my $ad_pos = -1, my $dp_pos = -1;
	for ( my $i = 0 ; $i < scalar(@format) ; $i++ ) {
		if ( $format[$i] eq 'AD' ) { $ad_pos = $i; }
		if ( $format[$i] eq 'DP' ) { $dp_pos = $i; }
	}
	if ( $ad_pos == -1 ) { return 0; }

	my $ref = 0, my $alt = 0, my $depth = 0;
	foreach my $data (@dataset) {
		my @val = split( /:/, $data );
		if ( defined( $val[$ad_pos] ) ) {
			my ( $ad1, $ad2 ) = split( /,/, $val[$ad_pos] );
			$ref += $ad1;
			$alt += $ad2;
		}
		if ( defined( $val[$dp_pos] ) ) {
			$depth += $val[$dp_pos];
		}
	}
	if ( $depth == 0 and $ref + $alt > 0 ) { $depth = $ref + $alt; }
	return ( $depth == 0 ? 0 : $alt / $depth );
}

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -i|--individual: Path to IndividualTable"
	  . " (default: $TS_Params::individual_table)\n";
	print STDOUT " -c|--column: Case/Control column in IndividualTable"
	  . " (default: 3)\n";
	print STDOUT " -d|--dir: Path to variant call data dir"
	  . " (default: $TS_Params::output[8])\n";
	print STDOUT " -b|--bed: Path to another BED file (optional)\n";
	print STDOUT " -a|--allele_freq: AF threshold (default: 0.3)\n";
	print STDOUT " -f|--fs: FS threshold (default: 100)\n";
	print STDOUT " -h|--help: Show this message\n";
}
__DATA__
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
