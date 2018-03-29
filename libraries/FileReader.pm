package FileReader;

use strict;

sub read_file{
    my ($file, $type, $v) = @_;
    if ($type eq "GTF"){
	return read_GTF($file, $v);
    }
    elsif( $type eq "GFF"){
	return read_GFF($file, $v);
    }
    elsif( $type eq "BED12"){
	return read_BED12($file, $v);
    }
    elsif( $type eq "PAF"){
	return read_PAF($file, $v);
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
	my $exon = [$chr, $start, $end, $strand, $t_id, $g_id];
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
	my $exon = [$chr, $start, $end, $strand, $t_id, $g_id];
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
	    
	    my $exon = [$chr, $exon_start, $exon_end, $strand, $t_id, $g_id];
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
sub read_PAF{
    my ($file, $verbose) = @_;

    my $trans;
    my %count;
    
    print STDERR "PAF reader not implemented yet. Apologies.\n";
    exit(0)
}

1;
