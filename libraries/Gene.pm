package Gene;

use strict;
use libraries::Cluster;

my $verbose = 0;
my $verbose2 = 0;

###########################################
# function to make genes
# We first cluster transcripts by genomic extent
# then each cluster is analyzed to build genes, where a 
# is defined as a set of transcripts transcribed from the same locus and strand
# and sharing at least one splice-site or exon,
# then genes are labeled. 
#
# input:
#     a pointer to a list of transcripts, 
#     a constrained = "exon" or "splice-site", which defines how to link the transcripts
#     a flag = reuse. It means that the input already had gene Ids, which we reuse for labelling

sub make_genes{
    my ($transcripts, $constraint, $reuse) = @_;

    # cluster transcripts according genomic extension overlap
    my ($t_clusters, $t_clusters_start, $t_clusters_end) = Cluster::cluster_transcripts($transcripts);

    # build genes in each cluster
    my @genes;
    foreach my $cluster (@$t_clusters){
	my $genes = build_genes($cluster, $constraint);
	push(@genes, @$genes);
    }

    # label genes
    my $label = Gene::label_genes(\@genes, $reuse);

    return (\@genes, $label);
}


###########################################
# function to build genes
# We define a gene as a set of transcripts transcribed from the same locus and strand
# and sharing at least one splice-site
#
# input:
#     a pointer to a list of transcripts that overlap in extent: $loci
#     the index position of the start and end coordinate of each transcript: $s, $e

sub build_genes{
    my ($transcripts, $constraint) = @_;
    
    # use flag to run Depth First afterwards
    my $color; 
    
    # we generate cluster of loci                                                                                                                 
    my $locus_clusters;

    # sort loci by the number of exons in desceding order
    my @sorted_transcripts = sort { transcript_start($a) <=> transcript_end($b) } @$transcripts;
    
    if ($verbose2){
	foreach my $t (@sorted_transcripts){
	    print_transcript($t);
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
	    my $test;
	    if ($constraint && $constraint eq "exon"){
		$test = test_exon_overlap($tran1, $tran2);
	    }
	    elsif($constraint && $constraint eq "splice-site"){
		$test = test_splice_site_overlap($tran1, $tran2);
	    }
	    else{
		$test = test_exon_overlap($tran1, $tran2);
	    }
	    if ($test){
		if ($verbose2){
		    print "linking\n";
		    print_transcript($tran1);
		    print_transcript($tran2);
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
    my @exons1 = get_sorted_exons($tran1);
    my @exons2 = get_sorted_exons($tran2);
    
    # each transcript is a pointer to a list (arrayref) of exons
    # each exon is a pointer to a list (arrayref): [$chr, $start, $end, $strand, $frame, $exon_label, $trans_id, $gene_id ]                        
    foreach my $exon1 ( @exons1 ){
	foreach my $exon2 ( @exons2 ){
	    # exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                      
	    if ( $exon1->[3] == $exon2->[3] && $exon1->[4] == $exon2->[4] ){ # coincide in an exon
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
    return 0;
}


sub test_splice_site_overlap{
    my ($tran1, $tran2) = @_;
    my @exons1 = get_sorted_exons($tran1);
    my @exons2 = get_sorted_exons($tran2);
    
    # each transcript is a pointer to a list (arrayref) of exons
    # each exon is a pointer to a list (arrayref): [$chr, $start, $end, $strand, $frame, $exon_label, $trans_id, $gene_id ]                        
    foreach my $exon1 ( @exons1 ){
	foreach my $exon2 ( @exons2 ){
	    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                      
	    if ( !( $exon1->[4]  < $exon2->[3] || $exon1->[3]  > $exon2->[4] )  # overlap                                   
		 &&
		 ( $exon1->[3] == $exon2->[3] || $exon1->[4] == $exon2->[4] ) ){ # coincide in at least one splice-site
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
    return 0;
}

sub label_genes{
    my ($genes, $reuse) = @_;
    my $count = 1;
    my $label1;
    my $label2;
    foreach my $g (@$genes){
	my $labels;
	
	# build simple label;
	my $gene_id = "gene_".$count;
	$label1->{$g} = $gene_id;
	$count++;
	
	# reuse label
	foreach my $t (@$g){
	    my @e = @$t;
	    # exon: [chr, start, end, strand, t_id, g_id]
	    my $t_id = $e[0]->[4]; 
	    my $g_id = $e[0]->[5];
	    $labels->{$g_id}++;
	}
	my @ids = sort {$a cmp $b} keys %$labels;
	my $new_id = join ":", @ids;
	$label2->{$g} = $new_id;
    }
    if ($reuse){
	return $label2;
    }
    else{
	return $label1;
    }
}


sub get_sorted_exons{
    my ($t) = @_;
    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                               
    return sort {$a->[3] <=> $b->[3]} @{$t};
}

sub transcript_start{
    my ($t) = @_;
    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                               
    my @e = sort {$a->[3] <=> $b->[3]} @{$t};
    return $e[0]->[3];
}

sub transcript_end{
    my ($t) = @_;
    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                               
    my @e = sort {$a->[3] <=> $b->[3]} @{$t};
    return $e[-1]->[4];
}

sub print_transcript{
    my ($t) = @_;
    my @e = @$t;
    my @coords;
    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                               
    my $strand = $e[0]->[6];
    my $chr    = $e[0]->[0];
    my $t_id   = $e[0]->[8];
    foreach my $ee (@e){
        push( @coords, $ee->[3]."-".$ee->[4] );
    }
    my $s = join "\t", ($t_id, $chr, $strand, @coords);
    print $s."\n";
}


1;
