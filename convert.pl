#!/usr/bin/perl -w

use strict;

use libraries::FileReader;
use libraries::FileWriter;
use Getopt::Long;

my ($file, $format_in, $output, $format_out, $help, $distance, $verbose, $del_size);

# default distance between matches in PAF file to join them into a transcript
$distance = 100000;

# default size of a deletion in the cigar string to be considered a potential intron
$del_size = 25;

# check options
my $result = GetOptions(
    'i|input=s'      => \$file,
    'f|format_in=s'  => \$format_in,
    'o|output=s'     => \$output,
    'u|format_out=s' => \$format_out,
    'd|distance=i'   => \$distance,
    's|del_size=i'   => \$del_size,
    'v|verbose'      => \$verbose,
    'h|help'         => \$help
    );

print_usage() if( $help);

unless ( $file && $format_in && $format_out && $output){
    print STDERR "ERROR: need to specify some input information\n";
    print_usage();
    exit(0);
}

unless( ($format_in eq "GFF" || $format_in eq "GTF" || $format_in eq "BED12" || $format_in eq "PAF")
	&&
	($format_in eq "GFF" || $format_in eq "GTF" || $format_in eq "BED12" || $format_in eq "PAF") ){
    print STDERR "ERROR: formats supported so far: GFF, GTF, BED12, PAF\n";
    print_usage();
}


# read input
my $trans = FileReader::read_file($file, $format_in, $verbose, $distance, $del_size);

# write output
FileWriter::write_file($trans, $output, $format_out);


sub print_usage{
    print STDERR "\n";
    print STDERR "Usage: perl convert.pl -i <input> -f <format_in> -o <output> -u <format_out>\n";
    print STDERR "\n";
    print STDERR "-i | --input:      File name with transcripts\n";
    print STDERR "-f | --format_in:  Format of the transcript file: GTF, GFF, BED12, PAF\n";
    print STDERR "-o | --output:     Output file name\n";
    print STDERR "-u | --format_out: Format of the transcript file: GTF, GFF, BED12, PAF\n";
    print STDERR "-d | --distance:   Distance between matches in cigar string file to join them into a transcript (default: 100000)\n";
    print STDERR "-s | --del_size:   Smallest size of a deletion in the cigar string to be considered a possible intron (detault: 25)\n"; 
    print STDERR "-v | --verbose:    Switch on useful prints for debugging\n";
    print STDERR "-h | --help:       Print this help\n";
    print STDERR "\n";
    exit(0);
}

