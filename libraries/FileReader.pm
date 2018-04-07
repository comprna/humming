package FileReader;

use strict;
use libraries::Cluster;

sub read_file{
    my ($file, $type, $verbose, $distance, $deletion_size) = @_;
    if ($type eq "GTF"){
	return read_GTF($file, $verbose);
    }
    elsif( $type eq "GFF"){
	return read_GFF($file, $verbose);
    }
    elsif( $type eq "BED12"){
	return read_BED12($file, $verbose);
    }
    elsif( $type eq "PAF"){
	return read_PAF($file, $verbose, $distance, $deletion_size);
    }
}


# GTF
#
# 1 chr
# 2 DB
# 3 feature
# 4 start (1-based)
# 5 end
# 6 score
# 7 strand
# 8 frame
# 9 gene_id "<id>"; transcript_id "<id>";

sub read_GTF{
    my ($file, $verbose) = @_;

    my $trans;
    open(IN,"<$file") or die("cannot open $file");
    while(<IN>){
	chomp;
	#Chr1   TAIR10  exon    3631    3913    .       +       .       ...
	my ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $col9) = split /\t/, $_;
	next unless $feature eq "exon";
	my ($t_id, $g_id);
	if ( $col9=~/gene_id\s+\"(\S+)\";\s+transcript_id\s+\"(\S+)\";/){
	    $t_id = $2;
	    $g_id = $1;
	}
	my $exon = [$chr, $db, $feature, $start, $end, $score, $strand, $frame, $t_id, $g_id];
	print "storing exon: $chr, $start, $end, $strand, $t_id, $g_id\n" if $verbose;
	push( @{$trans->{$t_id}}, $exon );
    }
    #print "transcripts: ".scalar(keys %$trans)."\n";

    my @transcripts;
    foreach my $id (keys %$trans ){
	my @e = @{$trans->{$id}};
	my $t = \@e;
	push( @transcripts, $t);
    }
    return \@transcripts;
}

# GFF
#
# 1 chr
# 2 DB
# 3 feature
# 4 start (1-based)
# 5 end
# 6 score
# 7 strand
# 8 frame
# 9 group ID / Parent=ID.
sub read_GFF{
    my ($file, $verbose) = @_;

    my $trans;
    open(IN,"<$file") or die("cannot open $file");
    while(<IN>){
	chomp;
	#Chr1   TAIR10  exon    3631    3913    .       +       .       Parent=AT1G01010.1
	my ($chr, $db, $feature, $start, $end, $score, $strand, $frame, $col9) = split /\t/, $_;
	next unless $feature eq "exon";
	my ($t_id, $g_id);
	if ($col9 =~/Parent=(\S+)/){
	    $t_id = $1;
	    $g_id = $t_id;
	}    
	else{
	    $t_id = $col9;
	    $g_id = $t_id;
	}
	my $exon = [$chr, $db, $feature, $start, $end, $score, $strand, $frame, $t_id, $g_id];
	print "storing exon: $chr, $start, $end, $strand, $t_id, $g_id\n" if $verbose;
	push( @{$trans->{$t_id}}, $exon );
    }
    my @transcripts;
    foreach my $id (keys %$trans ){
	my @e = @{$trans->{$id}};
	my $t = \@e;
	push( @transcripts, $t);
    }
    return \@transcripts;
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
sub read_BED12{
    my ($file, $verbose) = @_;

    my $trans;
    my %count;
    open(IN,"<$file") or die("cannot open $file");
    while(<IN>){
        chomp;
	my ($chr, $chr_start, $chr_end, $name, $score, $strand, $thick_start, $thick_end, $item, 
	    $block_count, $block_sizes, $block_starts) = split /\s+/, $_;
	
        # In case names are duplicated with add a counter "_2, _3..."
	$count{$name}++;
	my $t_id = $name;
	if ( $count{$name} > 1 ){
	    $t_id .= "_".$count{$name};
	}
	my $g_id = $t_id;
	
	my @start = split "\,", $block_starts;
	my @size = split "\,", $block_sizes;
	
	$chr_start++;
	for(my $i=0; $i<$block_count; $i++){
	    my ($exon_start, $exon_end);
	    $exon_start = $chr_start + $start[$i];
	    $exon_end   = $exon_start + $size[$i] - 1;
	    #my $exon = [$chr, $exon_start, $exon_end, $strand, $t_id, $g_id];
	    my $exon = [$chr, "bed12", "exon", $exon_start, $exon_end, $score, $strand, ".", $t_id, $g_id];
	    print "storing exon: $chr, $exon_start, $exon_end, $strand, $t_id, $g_id\n" if $verbose;
	    push( @{$trans->{$t_id}}, $exon );
	}
    }
    my @transcripts;
    foreach my $id (keys %$trans ){
	my @e = @{$trans->{$id}};
	my $t = \@e;
	push( @transcripts, $t);
    }
    return \@transcripts;
}

# PAF
# PAF format:
#5f0ea31a-9cbc-4a00-b113-a76d6e72a742    703     2       433     -       NC_000009.12    138394717       19376508        19378797 2249 2358 60 tp:A:P cm:i:18 s1:i:101 s2:i:78 NM:i:109 ms:i:129 AS:i:65 nn:i:0 ts:A:+ cg:Z:28M2I9M2I5M2I6M1D4M2D7M1D11M2D4M1I3M5I2M5I2M4I5M3D9M1I11M1D25M3I2M1716N2M1I4M1I5M2I2M1D16M2I4M4I5M1I13M3D6M5I9M3I8M3I22M2I2M1I3M1D1M2D8M2I1M1I29M193N13M1I16M1I6M1I11M1I9M7I14M3I6M2I4M1D10M

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
sub read_PAF{
    my ($file, $verbose, $distance, $deletion_size) = @_;
    #$verbose = 1;
    my $trans;
    my %count;    
    open(IN,"<$file") or die("cannot open $file");
    my @ranges;
    while(<IN>){
        chomp;  
	my ($query_name,     $query_length,  $query_start,  $query_end, $relative_strand,
	    $target_name,    $target_length, $target_start, $target_end,
	    $num_of_matches, $align_length,  $map_quality, @sam_fields ) = split /\t/, $_;
	
	$count{$query_name}++;
	print "reading $query_name $target_name\n" if $verbose;
	my $cigar;
	foreach my $f (@sam_fields){
	    #print STDERR $f."\n";
	    # get cigar string
	    if ($f =~/cg:Z:(.*)/){
		$cigar = $1;
		print "getting cigar $cigar\n" if $verbose;
	    }
	}
      	# transform to 1-based
	$query_start++;
	$query_end++;
	$target_start++;
	$target_end++;
	
	# parse cigar:
	my $blocks = parse_cigar($cigar, $verbose);

	# create blocks
	my $new_blocks;
	my $new_block = [];
	push( @$new_blocks, $new_block);
	my $count = 0;	
	foreach my $block ( @$blocks ){
            my ($size, $id) = @$block;
	    if ($count == 0){
		# Create the first block
		push( @$new_block, [$size, "exon"] );
		$count++;
	    }
	    # any cigar block that is a Match, Insertion, or a Deletion smaller than $deletion_size
	    # is considered an exon block
	    elsif ( $id eq "M" || $id eq "I" || ($id eq "D" && $size < $deletion_size) ){
		push( @$new_block, [$size, "exon"] );
		$count++;
	    }
	    # N blocks and Deletions larger (or equal) than $deltion_size are considered as introns
	    elsif( $id eq "N" || ($id eq "D" && $size >= $deletion_size) ){
		$new_block = [];
		push( @$new_blocks, $new_block );
		push( @$new_block, [$size, "intron"] );
		$new_block = [];
		push( @$new_blocks, $new_block );
		$count++;
	    }
	}
	
	# we create the new blocks by putting into a single
	# one the grouped blocks
	my $really_new_blocks;
	foreach my $new_block (@$new_blocks){
	    my $size = 0;
	    my $id;
	    foreach my $block (@$new_block){
		my ($this_size, $this_id) = @$block;
		$size = $size + $this_size;
		$id = $this_id;
	    }
	    push( @$really_new_blocks, [$size, $id] );
	}
	
	# we create exons if there are multiple blocks
	if (scalar(@$really_new_blocks) == 1){
	    # only one exon
	    push( @ranges, [$query_name,  $query_length,  $query_start, $query_end, $relative_strand,
			    $target_name, $target_length, $target_start, $target_end] ); 
	}
	#   ############## 
	#   eeeeeIIIIIeeee
	elsif( scalar(@$really_new_blocks) > 1){
	    my $q_start = $query_start;
	    my $t_start = $target_start;
	    my $q_end;
	    my $t_end;
	    foreach my $really_new_block (@$really_new_blocks){
		my ($block_size, $block_id) = @$really_new_block;
		if ($block_id eq "exon"){
		    $q_end = $q_start + $block_size - 1;
		    $t_end = $t_start + $block_size - 1;
		    my $range = [$query_name,  $query_length,  $q_start, $q_end, $relative_strand,
				 $target_name, $target_length, $t_start, $t_end];
		    push( @ranges, $range );
		    $q_start = $q_end + 1;
		    $t_start = $t_end + 1;
		}
		elsif( $block_id eq "intron"){
		    $t_start = $t_start + $block_size;
		}
	    }
	}
    }
    
    ######################################################
    # now we need to cluster nearby ranges (exons) into "transcripts":
    # range: [$query_name,  $query_length,  $q_start, $q_end, $relative_strand, $target_name, $target_length, $t_start, $t_end];
    # First separate ranges by chromosome and strand;
    my %ranges;
    my %chrs;
    my %strands;
    my %qnames;
    #print "obtained ".scalar(@ranges)." ranges\n";
    foreach my $range (@ranges){
	my $range_chr    = $range->[5];
	my $range_strand = $range->[4]; 
	my $range_qname  = $range->[0];
	#print "storing range in $range_chr $range_strand\n";
	push( @{$ranges{$range_chr}{$range_strand}{$range_qname}}, $range );
	$chrs{$range_chr}++;
	$strands{$range_strand}++;
	$qnames{$range_qname}++;
    }
    # Now cluster by strand and chromosome
    my %seen;
    foreach my $range_chr ( keys %chrs ){
	foreach my $range_strand ( keys %strands ){
	    foreach my $range_qname  ( keys %qnames ){
		next unless $ranges{$range_chr}{$range_strand}{$range_qname};
		my @ranges_list = @{$ranges{$range_chr}{$range_strand}{$range_qname}};
		my ($clusters, $clusters_start, $clusters_end ) = Cluster::cluster_by_proximity( \@ranges_list, 7, 8, $distance, $verbose);
		
		# each cluster is a candidate transcript
		foreach my $cluster (@$clusters){
		    my $append = "";
		    my $this_trans;
		    my @these_ranges = @$cluster;
		    my $this_query_name = $these_ranges[0]->[0];
		    #print "seen $this_query_name\n";
		    $seen{$this_query_name}++;
		    my $q_name = $this_query_name."_".$seen{$this_query_name};
		    
		    #if ( $seen{$this_query_name} ){
		    #    $append = "_".$seen{$this_query_name};
		    #    #print "seen ".$this_query_name."\n";
		    #}
		    foreach my $range (@these_ranges){
			my ($query_name,  $query_length,  $q_start, $q_end, $relative_strand, $target_name, $target_length, $t_start, $t_end) = @$range;
			# build an exon
			# [$chr, db, feature, $exon_start, $exon_end, $score, $strand, ".", $t_id, $g_id];
			my $e = [$target_name, "hum", "exon",$t_start, $t_end, "0", $relative_strand, ".", $q_name, $q_name]; 
			push( @$this_trans, $e );
		    }
		    push( @$trans, $this_trans );
		}
	    }
	}
    }
    return $trans;
}    

#CIGAR  CIGAR strings ("Compact Idiosyncratic Gapped Alignment Report"),
#M   alignment match (can be a sequence match or mismatch)
#I   insertion to the reference
#D   deletion from the reference
#N   skipped region from the reference
#S   soft clipping (clipped sequences present in SEQ)
#H   hard clipping (clipped sequences NOT present in SEQ)
#P   padding (silent deletion from padded reference)
#=   sequence match
#X   sequence mismatch
#
# H can only be present as the first and/or last operation.
# S may only have H operations between them and the ends of the string.
# For mRNA-to-genome alignment, an N operation represents an intron.
# For other types of alignments, the interpretation of N is not defined.
# Sum of the lengths of the M/I/S/=/X operations shall equal the length of SEQ.
#
sub parse_cigar {
    my ($c, $verbose) = @_;
    print "parsing $c\n" if $verbose;
    my @items = ($c =~ /(\d*[MIDNSHP=X])/g);
    my @blocks;
    foreach my $item (@items){
	print $item."\n" if $verbose;
	$item =~/(\d+)([\w=])/;
	my $num = $1;
	my $id  = $2;
	push(@blocks, [$num, $id]);
    }
    return \@blocks;
}


1;
