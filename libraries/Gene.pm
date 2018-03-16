package Gene;

use strict;
use libraries::Cluster;

my $verbose = 0;
my $verbose2 = 0;

###########################################
# function to build genes
# We define a gene as a set of transcripts transcribed from the same locus and strand
# and sharing at least one splice-site
#
# input:
#     a pointer to a list of transcripts that overlap in extent: $loci
#     the index position of the start and end coordinate of each transcript: $s, $e

sub build_genes{
    my ($transcripts) = @_;
    
    # use flag to run Depth First afterwards
    my $color; 
    
    # we generate cluster of loci                                                                                                                 
    my $locus_clusters;

    # sort loci by the number of exons in desceding order
    my @sorted_transcripts = sort { transcript_start($a) <=> transcript_end($b) } @$transcripts;
    
    if ($verbose2){
	foreach my $t (@sorted_transcripts){
	    print_Transcript($t);
	}
    }
    
    # clusters transcripts by whether or not any exon overlaps with an exon in                                                     
    # another transcript        
    my %id_to_object;
    foreach my $tran (@sorted_transcripts) {
        $color->{$tran} = 1;
	$id_to_object{$tran} = $tran;
    }
    # build adjacency list:
    my $adj_list;
  TRAN1:
    foreach my $tran1 (@sorted_transcripts) {
      TRAN2:
        foreach my $tran2 (@sorted_transcripts) {
            next if ($tran1 == $tran2);
	    if (test_exon_overlap($tran1, $tran2) ){
		if ($verbose2){
		    print "linking\n";
		    print_Transcript($tran1);
		    print_Transcript($tran2);
		}
		push( @{$adj_list->{$tran1}}, $tran2 );
	    }
	}
    }
    # Now do depth first to reconstruct the connected components                                                                             
    # using the adjacency lists:                                                                                                                 
    my @tran_clusters;
    foreach my $tran (@sorted_transcripts){
	if ( $color->{$tran} == 1 ){
	    my $nodes = depth_first_nonrecursive($adj_list, $tran, $color);
	    my @ids = keys %{$nodes};
	    print "linked nodes @ids\n" if $verbose2;
	    my @this_cluster;
	    foreach my $id (@ids){
		push( @this_cluster, $id_to_object{$id} );
	    }
	    push( @tran_clusters, \@this_cluster) if ( scalar( @this_cluster) );
	}
    }
    return \@tran_clusters;
}

############################
# algorithm to recover the groups of linked transcripts
# input: 
#        adjacency list for a list $adj_list (pointer to a table of transcripts)
#        a transcript $tran (pointer to a list of exons)
#        the color label of the transcripts $color (pointer to a table)
#
sub depth_first_nonrecursive{
    my ($adj_list, $tran, $color) = @_;
    my @s = ();
    push( @s, $tran);
    my %nodes;
    while(scalar(@s)){
	my $this_tran = pop @s;
	if ($color->{$this_tran} == 1 ){
	    $color->{$this_tran} = 2;
	    $nodes{$this_tran}++;
	    foreach my $other_tran ( @{$adj_list->{$this_tran}} ){
		push( @s, $other_tran);
	    }
	}
    }
    return \%nodes;
}


##########################
# returns 1 if the transcripts
# have exons with at least
# a shared splice-site

sub test_exon_overlap{
    my ($tran1, $tran2) = @_;
    my @exons1 = get_exons($tran1);
    my @exons2 = get_exons($tran2);
    
    # each transcript is a pointer to a list (arrayref) of exons
    # each exon is a pointer to a list (arrayref): [$chr, $start, $end, $strand, $frame, $exon_label, $trans_id, $gene_id ]                        
    foreach my $exon1 ( @exons1 ){
	foreach my $exon2 ( @exons2 ){
	    if ( !( $exon1->[2]  < $exon2->[1] || $exon1->[1]  > $exon2->[2] )  # overlap                                   
		 &&
		 ( $exon1->[1] == $exon2->[1] || $exon1->[2] == $exon2->[2] ) ){ # coincide in at least one splice-site
		if($verbose2){
		    print "LINKED\n";
		    #print "link $tran1 $tran2\n";                                                                                              
		    #print "color[$tran1] = ".$color->{$tran1}."\n";                                                                            
		    #print "color[$tran2] = ".$color->{$tran2}."\n";
		}
		return 1;
	    }
	}
    }
}


sub get_exons{
    my ($t) = @_;
    return sort {$a->[1] <=> $b->[1]} @{$t};
}

sub transcript_start{
    my ($t) = @_;
    my @e = sort {$a->[1] <=> $b->[1]} @{$t};
    return $e[0]->[1];
}

sub transcript_end{
    my ($t) = @_;
    my @e = sort {$a->[1] <=> $b->[1]} @{$t};
    return $e[-1]->[2];
}



sub print_Transcript{
    my ($trans, $out) = @_;
    my $s = join "\t", ("Transcript:", @$trans);
    if ($out){
        print STDERR "Transcript:".$s."\n";
    }
    else{
        print "Transcript:".$s."\n";
    }
}

sub print_Exons{
    my ($t) = @_;
    my @exons = get_exons($t);
    foreach my $e (@exons){
        my $p = join "\t", ("INFO:", "Exon:", @$e);
        print $p."\n";
    }
}


1;
