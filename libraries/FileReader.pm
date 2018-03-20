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
}


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


1;
