#!/usr/bin/perl -w

use strict;

use libraries::FileReader;
use libraries::FileWriter;
use libraries::Gene;
use Getopt::Long;

my ($file, $format, $output, $help);
my $constraint = "exon";
my $reuse      = "0";

# check options
my $result = GetOptions(
    'i|input=s'      => \$file,
    'f|format=s'     => \$format,
    'o|output=s'     => \$output,
    'c|constraint=s' => \$constraint,
    'r|reuse'        => \$reuse,
    'h|help'         => \$help
    );

print_usage() if( $help);

unless ( $file && $format && $constraint && $output){
    print STDERR "ERROR: need to specify some input information\n";
    print_usage();
    exit(0);
}

unless( $format eq "GFF" || $format eq "GTF"){
    print STDERR "ERROR: formats supported so far: GFF, GTF\n";
    print_usage();
}

unless ( $constraint eq "exon" || $constraint eq "splice-site" ){
    print STDERR "ERROR: please use at least one of these options: \"-c exon\" or \"-c splice-site\"\n";
    print_usage();
}

sub print_usage{
    print STDERR "\n";
    print STDERR "Usage: perl build_genes.pl -i <input> -f <format> -c <constraint>\n";
    print STDERR "\n";
    print STDERR "-i | --input:      File name with transcripts\n";
    print STDERR "-f | --format:     Format of the transcript file: GTF, GFF\n";
    print STDERR "-o | --output:     Output file name (GTF format)\n";
    print STDERR "-c | --constraint: Label:exon or splice-site. Genes are defined as transcripts overlapping\n";
    print STDERR "                   in the same strand and sharing an exon or a splice-site. Default: exon\n";
    print STDERR "-r | --reuse:      Optional. It will reuse gene IDs if available. Default: not set\n";
    print STDERR "-h | --help:       Print this help\n";
    print STDERR "\n";
    exit(0);
}


# read input
my $trans = FileReader::read_file($file, $format);

# make genes
my ($genes, $label) = Gene::make_genes($trans, $constraint, $reuse);

# output genes
FileWriter::write_genes_in_GTF($genes, $label, $output);


1;
