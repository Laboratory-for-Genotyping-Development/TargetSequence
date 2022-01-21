#!/usr/bin/perl

use strict;
use File::Temp qw/tempfile/;
use Getopt::Long;
use TS_Params;

#
my $temp_dir     = $TS_Params::output[3];
my $sample_table = $TS_Params::sample_table;
my $plot         = 1;
my $help         = 0;
GetOptions(
	'dir=s'   => \$temp_dir,
	'table=s' => \$sample_table,
	'plot!'   => \$plot,
	'help'    => \$help
);
if     ($help)              { &usage(); exit; }
unless ( -d $temp_dir )     { die "Directory not found: $temp_dir"; }
unless ( -f $sample_table ) { die "SampleTable file not found: $sample_table"; }
$temp_dir =~ s/\/+\Z//;

my $output_file = 'OnTargetProportion.list';
open( IN,  $sample_table )   or die "Can't open file: $sample_table";
open( OUT, ">$output_file" ) or die;
print OUT join(
	"\t",
	(
		'No',            'TotalReads',
		'TrimmedReads',  'TrimmedReads(%)',
		'OnTargetReads', 'OnTargetReads(%)',
		'Plate_Set',     'Class'
	)
) . "\n";

while (<IN>) {
	if ( $. == 1 ) { next; }
	s/(\r?\n)\z//;
	my ( $no, $plate_wellno, $plate_set, $class, $bed, $input_fastq, $naming )
	  = split(/\t/);
	my $ontarget = "${temp_dir}/OnTarget_$no";
	if ( -z $ontarget ) {
		print STDERR "$no is filesize 0\n";
		print OUT join( "\t", ( $no, 0, 0, 0, 0, 0, $plate_set, $class ) )
		  . "\n";
	}
	elsif ( -f $ontarget ) {
		open( FILE, $ontarget ) or die "Can't open file $ontarget";
		while (<FILE>) {
			s/(\r?\n)\z//;
			print OUT $_ . "\t${plate_set}\t${class}\n";
		}
		close(FILE);
	}
	else {
		print OUT join( "\t", ( $no, 0, 0, 0, 0, 0, $plate_set, $class ) )
		  . "\n";
	}
}
close(OUT);
close(IN);

if ($plot) {
	my ( $fh, $r_script ) = tempfile( SUFFIX => '.R', UNLINK => 1 );
	while (<DATA>) { print $fh $_; }
	system("$TS_Params::R --vanilla --slave < $r_script --args $output_file");
	close($fh);
}

sub usage() {
	print STDOUT "Usage: $0 [OPTION] ...\n";
	print STDOUT "Options:\n";
	print STDOUT " -t|--table: Path to SampleTable "
	  . "(default: $TS_Params::sample_table)\n";
	print STDOUT " -d|--dir: Path to OnTarget file directory "
	  . "(default: $TS_Params::output[3])\n";
	print STDOUT " -nop|--noplot: Do not plot TotalReads vs OnTarget(%)\n";
	print STDOUT " -h|--help: Show this message\n";
}

__DATA__
library(stringr)

title<-"TotalReads vs OnTarget(%)"
args<-commandArgs(trailingOnly=T)
DATA<-read.table(args[1], header=T, sep="\t")

xlabel<-"TotalReads"
ylabel<-"OnTarget(%)"

output_dir<-"OnTargetProportion"
if (!file.exists(output_dir)) {
	dir.create(output_dir)
}

if (nrow(DATA[DATA$Class=="Case",]) > 0) {
	bitmap(paste(output_dir, "OnTargetProportion_Case.png", sep="/"), width=12, height=12)
	par(cex.lab=1.5, oma=c(0, 1, 0, 0), cex.main=1.5)
	plot(DATA[,2], DATA[,6], col=ifelse(DATA$Class=="Case", "red", "transparent"), ylim=c(0, 100), xlab=xlabel, ylab=ylabel, main=title)
}
if (nrow(DATA[DATA$Class=="Control",]) > 0) {
	bitmap(paste(output_dir, "OnTargetProportion_Control.png", sep="/"), width=12, height=12)
	par(cex.lab=1.5, oma=c(0, 1, 0, 0), cex.main=1.5)
	plot(DATA[,2], DATA[,6], col=ifelse(DATA$Class=="Control", "blue", "transparent"), ylim=c(0, 100), xlab=xlabel, ylab=ylabel, main=title)
}
if (nrow(DATA[DATA$Class=="NTC",]) + nrow(DATA[DATA$Class=="GB",]) > 0) {
	bitmap(paste(output_dir, "OnTargetProportion_NTC.png", sep="/"), width=12, height=12)
	par(cex.lab=1.5, oma=c(0, 1, 0, 0), cex.main=1.5)
	plot(DATA[,2], DATA[,6], col=ifelse(DATA$Class=="NTC", "gray", ifelse(DATA$Class=="GB", "green", "transparent")), ylim=c(0, 100), xlab=xlabel, ylab=ylabel, main=title)
}

i<-30
for (set in unique(str_sub(DATA$Plate_Set, regexpr("_?([sS]et*)", DATA$Plate_Set)))) {
	sep <- ifelse(str_sub(set, 1, 1) == "_", "", "_")
	bitmap(paste(output_dir, "/OnTargetProportion", sep, set, ".png", sep=""), width=12, height=12)
	par(cex.lab=1.5, oma=c(0, 1, 0, 0), cex.main=1.5)
	plot(DATA[,2], DATA[,6], col=ifelse((str_sub(DATA$Plate_Set, regexpr("_?([sS]et*)", DATA$Plate_Set)) == set) & (str_count(DATA$Class, "Case|Control") > 0), colors()[i], "transparent"), ylim=c(0, 100), xlab=xlabel, ylab=ylabel, main=title)
	i<-i + 1
}
