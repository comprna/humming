#!/usr/bin/perl -w

use strict;

use libraries::Cluster;
use libraries::FileReader;
use libraries::CompareJunctions;
use libraries::CompareGenes;
use libraries::CompareTranscripts;
use Getopt::Long;

my ($file_pred, $file_ann, $format_pred, $format_ann, $help, $distance, $verbose, $del_size);

# default distance between matches in PAF file to join them into a transcript                                                          
$distance = 100000;

# default size of a deletion in the cigar string to be considered a potential intron                                                       
$del_size = 25;

my ($compare_junctions, $compare_genes, $compare_transcripts);

# check options
my $result = GetOptions(
    'p|prediction=s' => \$file_pred,
    'a|annotation=s' => \$file_ann,
    'f|format_pred=s'=> \$format_pred,
    'u|format_ann=s' => \$format_ann,
    'j|junctions'    => \$compare_junctions,
    't|transcripts'  => \$compare_transcripts,
    'd|distance=i'   => \$distance,
    's|del_size=i'   => \$del_size,
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

# read GFF lines in annotation file 1
# store all transcripts and all exons in each transcript
my $trans_ann = FileReader::read_file($file_ann, $format_ann, $verbose, $distance, $del_size);

# read GFF lines in annotation file 2
# store all transcripts and all exons in each transcript
my $trans_pred = FileReader::read_file($file_pred, $format_pred, $verbose, $distance, $del_size);

# $trans_ann and $trans_pred are two pointers two lists (array-refs)
# that hold the transcripts. Each transcript is a pointer to a list (arrayref) of exons, 
# where each exon is a pointer to a list (arrayref)
# of exon properties read from the GTF/GFF: $exon : [$chr, $start, $end, $strand, $t_id, $g_id];

#print "found ".scalar(@$trans_ann)." transcripts in the annotation\n";
#print "found ".scalar(@$trans_pred)." transcripts in the prediction\n";


if ($compare_junctions){
    CompareJunctions::compare_junctions($trans_ann,$trans_pred);
}

if ($compare_genes){
     CompareGenes::compare_genes($trans_ann,$trans_pred);
}

if ($compare_transcripts){
    CompareTranscripts::compare_transcripts_directional($trans_ann,$trans_pred);
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
    print STDERR "-d | --distance:   Distance between matches in cigar string file to join them into a transcript (default: 100000)\n";
    print STDERR "-s | --del_size:   Smallest size of a deletion in the cigar string to be considered a possible intron (detault: 25)\n";
    print STDERR "-h | --help:\t\tPrint this help\n";
    print STDERR "\n";
    exit(0);
}

