#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;
use TS_Params;

# GroupGenotypingKnownVariants_v3.pl VariantList BamList OutputDir
my $variants_list, my $bam_list, my $output_dir;
GetOptions(
	'variants=s'   => \$variants_list,
	'bamlist=s'    => \$bam_list,
	'output_dir=s' => \$output_dir
);

unless ( -f $bam_list ) { die "File not found: $bam_list"; }
$output_dir =~ s/\/+\Z//;

my $buffer    = 0;
my $freq_low  = 0.2;
my $freq_high = 0.8;

my %double;
open( IN, $variants_list ) or die "Can't open file: $variants_list";
while (<IN>) {
	if ( $. == 1 ) { next; }
	s/(\r?\n)\z//;
	my @data = split(/\t/);
	for ( my $i = 0 ; $i <= $#data ; $i++ ) { $data[$i] = trim( $data[$i] ); }
	my $name = join( '_', @data[ 0, 1 ] );
	push( @{ $double{$name} }, [ @data[ 2 .. $#data ] ] );
}
close(IN);

# bamの個体名を展開
my @inds;
open( BAM, $bam_list ) or die "Can't open file: $bam_list";
while (<BAM>) {
	s/(\r?\n)\z//;
	my $name = basename( $_, '.bam' );
	push( @inds, $name );
}
close(BAM);

my $filename = $bam_list;
if ( $filename =~ /\w+\-/ ) {
	$filename = $&;
}

open( OUT, ">${output_dir}/Group_$filename" )
  or die "Can't create file ${output_dir}/Group_$filename";
print OUT join(
	"\t",
	(
		'ind',       'chr',  'pos',         'type',
		'ref',       'alt',  'total_depth', 'ref_depth',
		'alt_depth', 'freq', 'genotype'
	)
) . "\n";

open( MPILE,
		"$TS_Params::samtools mpileup"
	  . " -l $variants_list -f ${TS_Params::huref}.fasta"
	  . " -A -Q 0 -d 1000000 -L 1000000 -b $bam_list |" )
  or die;
while (<MPILE>) {
	s/(\r?\n)\z//;
	my ( $chr, $pos, $contig, @info ) = split(/\t/);
	my $key = join( '_', ( $chr, $pos ) );
	foreach my $name ( @{ $double{$key} } ) {
		my ( $type, $ref, $alt ) = @{$name};
		my $allele;
		if ( $type eq 'SNP' ) {
			$allele = $alt;
		}
		elsif ( $type eq 'IN' ) {
			my $length = length($alt) - 1;
			$allele = '[\.\,]\+' . $length . substr( $alt, 1 );
		}
		elsif ( $type eq 'DEL' ) {
			my $length = length($ref) - 1;
			$allele = '[\.\,]\-' . $length . substr( $ref, 1 );
		}
		else { next; }

		my $indcount = 0;
		foreach my $ind (@inds) {
			my $ref_depth = 0;
			my $alt_depth = 0;
			my $sequence  = $info[ ( $indcount * 3 + 1 ) ];

			# カウントに不要な文字列を消す
			$sequence =~ s/\^.|\$//g;
			if ( $type eq 'SNP' ) {

				# SNPの場合は、Depthの計算前に除去
				$sequence  = &remove_indel_seq($sequence);
				$alt_depth = ( $sequence =~ s/$allele//gi );
			}
			else {
				$alt_depth = ( $sequence =~ s/$allele//gi );
				$sequence =
				  &remove_indel_seq($sequence);    # In/Delの場合は、Depthの計算後に除去
			}
			if ( $alt_depth eq '' ) { $alt_depth = 0; }

			$ref_depth = ( $sequence =~ s/\.|\,//g );
			if ( $ref_depth eq '' ) { $ref_depth = 0; }

			my $ref_freq =
			  ( $info[ $indcount * 3 ] == 0 )
			  ? 'NA'
			  : sprintf( "%.4f", $ref_depth / $info[ $indcount * 3 ] );
			my $alt_freq =
			  ( $info[ $indcount * 3 ] == 0 )
			  ? 'NA'
			  : sprintf( "%.4f", $alt_depth / $info[ $indcount * 3 ] );

			my $genotype;
			if    ( $info[ $indcount * 3 ] == 0 )       { $genotype = './.'; }
			elsif ( $ref_freq >= $freq_high + $buffer ) { $genotype = '0/0'; }
			elsif ( $alt_freq >= $freq_high + $buffer ) { $genotype = '1/1'; }
			elsif ( $ref_freq + $alt_freq >= $freq_high + $buffer ) {
				$genotype = '0/1';
			}
			elsif ( $ref_freq + $alt_freq <= $freq_low - $buffer ) {
				$genotype = '2/2';
			}
			elsif ( $freq_low + $buffer <= $ref_freq
				and $ref_freq <= $freq_high - $buffer )
			{
				$genotype = '0/2';
			}
			elsif ( $freq_low + $buffer <= $alt_freq
				and $alt_freq <= $freq_high - $buffer )
			{
				$genotype = '1/2';
			}
			else { $genotype = './.'; }
			print OUT join(
				"\t",
				(
					$ind, $chr, $pos, $type, $ref, $alt, $info[ $indcount * 3 ],
					$ref_depth, $alt_depth, $alt_freq, $genotype
				)

			) . "\n";
			$indcount++;
		}
	}
	print "$key finished: " . printTimestamp() . "\n";
}
close(MPILE);
close(OUT);

sub remove_indel_seq() {
	my ($seq) = @_;
	while ( $seq =~ /([ACGTN\.\,])([\+\-]\d+)/i ) {
		my $indel =
		  "\\" . substr( $seq, index( $seq, $2 ), abs($2) + length($2) );
		if ( $1 =~ /([\.\,])/ ) { $indel = "$1$indel"; }
		$seq =~ s/$indel//g;
	}
	return ($seq);
}
