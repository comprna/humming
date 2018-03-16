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
    my ($loci, $s, $e) = @_;
    
    # use flag to run Depth First afterwards
    my $color; 
    
    # we generate cluster of loci                                                                                                                 
    my $locus_clusters;

    # sort loci by the number of exons in desceding order
    #my @sorted_loci = sort { scalar(get_all_Exons($b)) <=> scalar(get_all_Exons($a)) } @$loci;
    my @sorted_loci = sort { $a->[$s] <=> $b->[$s] } @$loci;
    
    if ($verbose2){
	foreach my $t (@sorted_loci){
	    print_Transcript($t);
	}
    }
    
    # clusters transcripts by whether or not any exon overlaps with an exon in                                                     
    # another transcript        
    my %id_to_object;
    foreach my $tran (@sorted_loci) {
        $color->{$tran} = 1;
	$id_to_object{$tran->[4]} = $tran;
    }
    # build adjacency list:
    my $adj_list;
  TRAN1:
    foreach my $tran1 (@sorted_loci) {
      TRAN2:
        foreach my $tran2 (@sorted_loci) {
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
    foreach my $tran (@sorted_loci){
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
	    $nodes{$this_tran->[5]}++;
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


###########################################
# function to get all exons from a transcript
sub get_all_Exons{
    my ($trans) = @_;
    #my $trans = [$chr, $t_start, $t_end, $strand, $coords, $feature_id, $gene_symbol];
    my ($chr, $t_start, $t_end, $strand, $coords, $trans_id, $gene_symbol) = @$trans;
    my @exon_ranges = split ",", $coords;
    my @exons;
    foreach my $range (@exon_ranges){
        my ($start, $end) = split "-", $range;
	push( @exons, [$chr, $start, $end, $strand, $trans_id, $gene_symbol] );
    }
    return @exons;
}


#########################

sub build_genes_sequential{
    my ($loci, $s, $e) = @_;
    
    # we generate cluster of loci                                                                                                                            
    my $locus_clusters;

    # sort loci by the number of exons in desceding order
    #my @sorted_loci = sort { scalar(get_all_Exons($b)) <=> scalar(get_all_Exons($a)) } @$loci;
    my @sorted_loci = sort { $a->[6] cmp $b->[6] } @$loci;

    if ($verbose){
	foreach my $t (@sorted_loci){
	    print_Transcript($t);
	}
    }
    # start and end coordinates of the clusters                                                                                                            
    my @cluster_exons;

    # Create the first locus_cluster with the first locus                                                                                                    
    my @locus_clusters;
    #push( @$locus_clusters, $locus_cluster);

    # we go over all of them in sorted order (from left to right)                                                                                            
    my $cluster_count = 0;
    my $count = 0;

  LOCUS:
    foreach my $locus ( @sorted_loci ){
        if ($count == 0){
            # Create the first locus_cluster with the first locus                                                                                            
            push( @{$locus_clusters[0]}, $locus);
            my @these_exons = get_all_Exons($locus);
            push( @{$cluster_exons[0]}, @these_exons );
            $cluster_count++;
	    $count++;
            next LOCUS;
        }
	my @other_exons = get_all_Exons($locus);
	my $overlap = 0;
	# go through all available clusters:
	for(my $i=0; $i<scalar(@locus_clusters); $i++){
	    foreach my $exon1 (@other_exons){
		foreach my $exon2 (@{$cluster_exons[$i]} ){
		    print "Comparing @$exon1 with @$exon2\n" if $verbose;
		    if ( !( $exon1->[2]  < $exon2->[1] || $exon1->[1]  > $exon2->[2] )  # overlap                                   
			 &&
			 ( $exon1->[1] == $exon2->[1] || $exon1->[2] == $exon2->[2] ) ){ # coincide in at least one splice-site
			$overlap = 1;
		    }
		}
	    }
	    if ($overlap){         	    
		# add locus to cluster
		print "added\n" if $verbose;
		push( @{$locus_clusters[$i]}, $locus );
		
		# update the exons:
		push( @{$cluster_exons[$i]}, @other_exons);
	    }
	}
	my $cluster_index = scalar(@locus_clusters) - 1;
	unless ($overlap){
	    # then we proceed with the next one: 
	    # create new cluster (same variable name, new memory address!!!!)                                     
	    push( @{$locus_clusters[$cluster_index+1]}, $locus );
	    push( @{$cluster_exons[$cluster_index+1]}, @other_exons);
	}
    }
    return  \@locus_clusters;
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
