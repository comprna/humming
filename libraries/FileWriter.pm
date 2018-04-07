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
        return write_BED12($trans, $output);
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
	    my ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $t_id, $g_id) = @$e;
	    $db      = "humming" unless $db;;
	    $feature = "exon" unless $feature;
	    $frame   = "." unless $frame;
	    $score   = 100 unless $score;
	    $score   = 100 if $score eq ".";
	    my $col9  = "gene_id \"".$g_id."\"; transcript_id \"".$t_id."\";";
	    my $s = join "\t", ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $col9);
	    print OUT $s."\n";
	}
    }
}


sub write_GFF{
    my ($trans, $output) = @_;
    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	foreach my $e (get_sorted_exons($t)){
	    my ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $t_id, $g_id) = @$e;
	    $db      = "humming" unless $db;;
	    $feature = "exon" unless $feature;
	    $frame   = "." unless $frame;
	    $score   = 100 unless $score;
	    $score   = 100 if $score eq ".";
	    my $col9    = $t_id;
	    my $s = join "\t", ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $col9);
	    print OUT $s."\n";
	}
    }
}



# BED12                                                                                                                                                      
#
# 1  chr                                                                                                                                        
# 2  chr_start (0 based)                                                                                                                           
# 3  chr_end                                                                                                                                        
# 4  name                                                                                                                                             
# 5  score                                                                                                                                          
# 6  strand                                                                                                                                           
# 7  thick_start                                                                                                                                    
# 8  thick_end                                                                                                                                      
# 9  ItemGBgroup ID / Parent=ID.                                                                                                                      
# 10 block_count                                                                                                                                        
# 11 block_sizes                                                                                                                                        
# 12 block_starts (relative to chr_start)  
sub write_BED12{
    my ($trans, $output) = @_;

    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	my @exons    = get_sorted_exons($t);
	my $t_chr    = $exons[0]->[0];
	my $t_start  = $exons[0]->[3] - 1;
	my $t_end    = $exons[-1]->[4];
	my $t_strand = $exons[0]->[6];
	my $t_score  = 0;
	my $t_id     = $exons[0]->[8];
	my $g_id     = $exons[0]->[9];
	my $name     = $t_id;
	$name = $g_id.":".$t_id unless ($t_id eq $g_id);
	my $block_count = scalar(@exons);
	my @block_starts;
	my @block_sizes;
	my $block_start = 0;
        foreach my $e (get_sorted_exons($t)){
	    my ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $t_id, $g_id) = @$e;
	    $score   = 100 unless $score;
	    $score   = 100 if $score eq ".";
	    my $t_score = $t_score + $score;
            my $size = $end - $start + 1;
	    $block_start = $start - $t_start - 1;
	    push(@block_sizes, $size);
	    push(@block_starts, $block_start);
	}
	my $block_sizes  = join ",",@block_sizes;
	my $block_starts = join ",",@block_starts;
	$t_score = $t_score / $block_count;
	my $s = join "\t", ($t_chr, $t_start, $t_end, $name, $t_score, $t_strand, $t_start, $t_end, "0", $block_count, $block_sizes, $block_starts);
	print OUT $s."\n";
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

	my $s = join "\t", ($query_name,  $query_length,  $query_start,  $query_end, $relative_strand,
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
