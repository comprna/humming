#!/usr/bin/perl -w

use strict;

use libraries::Cluster;
use libraries::FileReader;
use libraries::CompareJunctions;
use libraries::CompareGenes;
use Getopt::Long;

my $verbose  = 0;
my $verbose2 = 0;
my $verbose3 = 0;

my ($file_pred, $file_ann, $format_pred, $format_ann, $help);

my ($compare_junctions, $compare_genes, $compare_transcripts);

# check options
my $result = GetOptions(
    'p|prediction=s' => \$file_pred,
    'a|annotation=s' => \$file_ann,
    'f|format_pred=s'=> \$format_pred,
    'u|format_ann=s' => \$format_ann,
    'j|junctions'    => \$compare_junctions,
    't|transcripts'  => \$compare_transcripts,
    'g|genes'        => \$compare_genes,
    'h|help'         => \$help
    );

print_usage() if( $help);

unless ( $file_pred && $file_ann && $format_pred && $format_ann ){
    print STDERR "ERROR: needs to specify some input information\n";
    print_usage();
    exit(0);
}

unless( ($format_ann eq "GFF" || $format_ann eq "GTF") 
	&& 
	($format_pred eq "GFF" || $format_pred eq "GTF")){
    print STDERR "ERROR: formats supported so far: GFF, GTF\n";
    print_usage();
}

unless ( $compare_junctions || $compare_transcripts || $compare_genes ){
    print STDERR "ERROR: please use at least one of these options: -j -t -g\n";
    print_usage();
}

sub print_usage{
    print STDERR "\n";
    print STDERR "Usage: perl compare_mappings_and_annotations.pl -a <annotations> -p <predictions> -f <format-annotations> -u <format-predictions> < -j | -t | -g >\n";
    print STDERR "\n";
    print STDERR "-a | --annotation:\tFile with the transcript annotations\n";
    print STDERR "-p | --prediction:\tFile with the transcript predictions\n";
    print STDERR "-f | --format_ann:\tFormat of the annotation file: GTF, GFF, PAF\n";
    print STDERR "-u | --format_pred:\tFormat of the prediction file: GTF, GFF, PAF (it can be different from the annotation file)\n";
    print STDERR "-j | --junctions:\tCompare splice sites and exon-exon junctions\n";
    print STDERR "-t | --transcripts:\tCompare transcripts. It matches transcripts within the same loci\n";
    print STDERR "-g | --genes:\t\tCompare gene loci (genome extension in the same strand that produce one or more transcripts)\n";
    print STDERR "-h | --help:\t\tPrint this help\n";
    print STDERR "\n";
    exit(0);
}


my %strand;
my %chr;

# read GFF lines in annotation file 1
# store all transcripts and all exons in each transcript
my $trans_ann = FileReader::read_file($file_ann, $format_ann, $verbose);

# read GFF lines in annotation file 2
# store all transcripts and all exons in each transcript
my $trans_pred = FileReader::read_file($file_pred, $format_pred, $verbose);

# $trans_ann and $trans_pred are two dictionaries (hashes)
# that hold for each transcript_id (from the 9th column of the GTF/GFF file)
# the pointer to a list (arrayref) of exons, where each exon is a pointer to a list (arrayref)
# of exon properties: $exon : [$chr, $start, $end, $strand, $t_id, $g_id];

if ($compare_junctions){
    CompareJunctions::compare_junctions($trans_ann,$trans_pred);
}

if ($compare_genes){
     CompareGenes::compare_genes($trans_ann,$trans_pred);
}

if ($compare_transcripts){
    print STDERR "Operation compare genes not yet in place\n";
    exit(0);
    CompareTranscripts::compare_transcripts($trans_ann,$trans_pred);
}

