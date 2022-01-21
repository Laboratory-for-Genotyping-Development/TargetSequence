#
# 共通関数、定数定義ファイル
# 実行環境によって変更がありそうなもの、各スクリプト共通で使えそうなコードを記述していく。
#
package TS_Params;

use Exporter;

BEGIN {
	@ISA    = (Exporter);
	@EXPORT = qw(printTimestamp makeDirs trim);
}

# Dirs
my $BIN_DIR = '/usr/local/bin/';

our @output = (
	'Logfiles',                  '1.FastQC',
	'2.TrimmedFastq',            '3.TemporaryFiles',
	'4.SampleBam',               '5.SampleCall',
	'6.SampleCoverage',          '7.IndividualBam',
	'8.IndividualCall',          '9.IndividualCoverage',
	'10.IndividualGenotype',     '11.VAFHistogram',
	'12.VAFDepthPlot',           '13.RevisedVAFHistogram',
	'14.RevisedVAFDepthPlot',    '15.IndividualVAFHistogram',
	'16.IndividualVAFDepthPlot', '17.InsertSizeBam',
	'18.InsertSizeMetrics'
);

# Commands
our $fastqc        = $BIN_DIR . 'FastQC/fastqc';
our $fastx_trimmer = $BIN_DIR . 'fastx_toolkit/bin/fastx_trimmer';
our $cutadapt      = '/usr/bin/cutadapt';
our $bwa           = $BIN_DIR . 'bwa-0.7.17/bwa';
our $samtools      = $BIN_DIR . 'samtools-1.6/bin/samtools';
our $samtoolsOld   = $BIN_DIR . 'samtools-0.1.16/samtools';
our $GATK          = $BIN_DIR . 'GATK3.7-0/GenomeAnalysisTK.jar';
our $picard        = $BIN_DIR . 'picard-tools-2.5.0/picard.jar';
our $R             = 'R';
our $VCMM          = $BIN_DIR . 'VCMM-1.0.2/VCMM';
our $MergeSamFiles = $BIN_DIR . 'picard-tools-1.123/MergeSamFiles.jar';
our $vcftools      = 'vcftools';

# Scripts
our $scriptModVCMM    = 'modify_VCMMout_v1.3.pl';
our $scriptUniquePair = 'ExtractUniquePair_v1.2.pl';

# Reference
our $huref         = 'human_g1k_v37';
our $indel1000GVcf = '1000G_phase1.indels.b37.vcf';

# File names
our $sample_table     = 'SampleTable';
our $individual_table = 'IndividualTable';
our $variant_list     = 'Variants.list';
our $genotype_list    = 'IndividualGenotype.list';
our $genotype_summary = 'IndividualGenotype.summary';

# 現在時刻を取得
sub printTimestamp {
	( $sec, $min, $hour, $mday, $mon, $year ) = localtime(time);
	return sprintf(
		"%04d/%02d/%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $mday, $hour, $min, $sec
	);
}

# 解析ファイル出力用ディレクトリ作成
sub makeDirs {
	my $count = 0;
	foreach my $dir (@output) {
		if ( grep { $_ == $count } @_ ) {
			unless ( -d $dir ) { mkdir($dir) }
		}
		$count++;
	}
}

# 前後の空白を削除
sub trim {
	my ($val) = @_;
	$val =~ s/^ *(.+) *$/$1/g;
	return ($val);
}
1;
