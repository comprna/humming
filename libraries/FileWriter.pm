package FileWriter;

use strict;


sub write_file{
    my ($trans, $output, $format) = @_;

    if ($format eq "GTF"){
        return write_GTF($trans, $output);
    }
    elsif( $format eq "GFF"){
        return write_GFF($trans, $output);
    }
    elsif( $format eq "BED12"){
        return write_BD12($trans, $output);
    }
    elsif( $format eq "PAF"){
        return write_PAF($trans, $output);
    }
}


sub write_GTF{
    my ($trans, $output) = @_;
    
    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	foreach my $e (get_sorted_exons($t)){
	    my ($chr, $start, $end, $strand, $t_id, $g_id) = @$e;
	    my $source  = "humming";
	    my $feature = "exon";
	    my $frame   = ".";
	    my $score   = 100;
	    my $col9    = "gene_id \"".$g_id."\"; transcript_id \"".$t_id."\";";
	    my $s = join "\t", ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $col9);
	    print OUT $s."\n";
	}
    }
}


sub write_GFF{
    my ($trans, $output) = @_;

    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	foreach my $e (get_sorted_exons($t)){
	    my ($chr, $start, $end, $strand, $t_id, $g_id) = @$e;
	    my $source  = "humming";
	    my $feature = "exon";
	    my $frame   = ".";
	    my $score   = 100;
	    my $col9    = $t_id;
	    my $s = join "\t", ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $col9);
	    print OUT $s."\n";
	}
    }
}

sub write_genes_in_GTF{
    my ($genes, $label, $output) = @_;

    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $g (@$genes){
	foreach my $t (@$g){
	    foreach my $e (get_sorted_exons($t)){
		my ($chr, $start, $end, $strand, $t_id, $g_id) = @$e;
		$g_id = $label->{$g};
		my $source  = "humming";
		my $feature = "exon";
		my $frame   = ".";
		my $score   = 100;
		my $col9    = "gene_id \"".$g_id."\"; transcript_id \"".$t_id."\";";
		my $s = join "\t", ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $col9);
		print OUT $s."\n";
	    }
	}
    }
}



# 1  Query sequence name
# 2  Query sequence length
# 3  Query start (0-based)
# 4  Query end (0-based)
# 5  Relative strand: "+" or "-"
# 6  Target sequence name
# 7  Target sequence length
# 8  Target start on original strand (0-based)
# 9  Target end on original strand (0-based)
# 10 Number of sequence matches
# 11 Alignment block length (total number of sequence matches, mismatches, insertions and deletions in the alignment)
# 12 Mapping quality (0-255; 255 for missing)    
# 13 SAM-like output
sub write_PAF{
    my ($trans, $output) = @_;

    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	my $query_length = 0;
	my @exons = get_sorted_exons($t);
	#exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	# build cigar:
	my @cigar;
	for(my $i=0; $i<scalar(@exons); $i++){
	    my $match = $exons[$i]->[4] - $exons[$i]->[3]+1;
	    $query_length += $match;
	    push( @cigar, $match."M" );
	    if ($i<scalar(@exons)-1){
		my $intron = $exons[$i+1]->[3] - $exons[$i]->[4]-1;
		push( @cigar, $intron."N" );
	    }
	}
	my $cigar = join "", @cigar;
	my $target_start    = $exons[0]->[3];
	my $target_end      = $exons[-1]->[4];
	my $query_name      = $exons[0]->[8];
	my $relative_strand = $exons[0]->[6];
	my $align_length    = $query_length;
	my $query_start     = 1 - 1;
	my $query_end       = $query_length - 1;
	my $target_length   = "NA";
	my $map_quality     = "NA";
	my $target_name     = $exons[0]->[0];

	my $s = "\t", ($query_name,  $query_length,  $query_start,  $query_end, $relative_strand,
		       $target_name, $target_length, $target_start, $target_end,
		       $query_length, $align_length, $map_quality, $cigar);
	print OUT $s."\n";
    }
}

sub get_sorted_exons{
    my ($t) = @_;
    # chromosome, source, feature, start, end, score, strand, frame, transcript ID, gene ID
    my @e = sort {$a->[3] <=> $b->[3]} @{$t};
    return @e;
}



1;
