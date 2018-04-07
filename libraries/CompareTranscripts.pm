package CompareTranscripts;

use strict;
use libraries::Cluster;

my $verbose  = 1;

####################
# function to compare two lists of transcripts
# It will first create transcript clusters, and
# then will compare transcripts among the two sets
# using an algorithm for the best available reciprocal match
# It will then given the total and proportion of 
# transcripts that match the other set, 
# where "match" means all internal coordinates coincide
# (other definitions to be added)

sub compare_transcripts_bidirectional{
    my ($trans_ann, $trans_pred) = @_;
    
    # count total number of transcripts and overlaps between them:
    my $over_count = 0;
    my $pred_count = 0;
    my $ann_count  = 0;

    # cluster transcripts by genomic extent and strand overlap                                                                                    
    my ($clusters, $label) = Cluster::cluster_transcripts_2sets($trans_ann, $trans_pred);
    
    # from each transcript cluster, we compare transcripts from both labels
    foreach my $cluster (@$clusters){

	# collect all transcripts from each type
	my @ann_trans;
	my @pred_trans;
	
	#print "cluster with ".scalar(@$cluster)." transcripts\n";
	foreach my $t (@$cluster){
	    if ($label->{$t} eq "pred"){
		$pred_count++;
		push(@pred_trans, $t);
	    }
	    if ($label->{$t} eq "ann"){
                $ann_count++;
                push(@ann_trans, $t);
            }
	}
	
	# if we have transcripts from both types, we compare them
	my $some_overlap = 0;
        print "cluster with ".scalar(@ann_trans)." annotations and ".scalar(@pred_trans)." predictions\n";
	if (@ann_trans && @pred_trans){
	    my %potential_partners;
	    my %score; 
	    foreach my $t_ann (@ann_trans){
		foreach my $t_pred (@pred_trans){
		    my ($t_ann_j, $t_pred_j, $over_j ) = compare_transcript_junctions($t_ann, $t_pred);
		    
		    # simple overlap score = 2 overlap/ (sum of junctions with repetitions), in [0,1]
		    if ($t_ann_j + $t_pred_j> 0){
			my $score = 2* $over_j / ( $t_ann_j + $t_pred_j);
			$score{$t_ann}{$t_pred} = $score;
			print "storing score $score\n";
			if ($score > 0 ){
			    push( @{$potential_partners{$t_ann}}, $t_pred );
			    $some_overlap++;
			}
		    }
		}
	    }
	    # Once we have compared them all, we select the best available pairs if there are any matches:	
	    if ($some_overlap){
		my @couples = stable_marriage( \@ann_trans, \@pred_trans, \%potential_partners, \%score);
	    }
	}
    }
    my $s = join "\t", ("TRANS", 
			 "ann", $ann_count, 
			 "pred", $pred_count, 
			 "over",   $over_count,
			 "prop_ann",  $over_count/$ann_count,
			 "prop_pred", $over_count/$pred_count);
    print $s."\n";
}


sub compare_transcripts_directional{
    my ($trans_ann, $trans_pred) = @_;
    
    # count total number of transcripts and overlaps between them:
    my $over_pred_count = 0;
    my $over_ann_count = 0;
    my $pred_count = 0;
    my $ann_count  = 0;

    # cluster transcripts by genomic extent and strand overlap                                                                                    
    my ($clusters, $label) = Cluster::cluster_transcripts_2sets($trans_ann, $trans_pred);
    
    # from each transcript cluster, we compare transcripts from both labels
    foreach my $cluster (@$clusters){

	# collect all transcripts from each type
	my @ann_trans;
	my @pred_trans;
	
	#print "cluster with ".scalar(@$cluster)." transcripts\n";
	foreach my $t (@$cluster){
	    if ($label->{$t} eq "pred"){
		$pred_count++;
		push(@pred_trans, $t);
	    }
	    if ($label->{$t} eq "ann"){
                $ann_count++;
                push(@ann_trans, $t);
            }
	}
	
	# if we have transcripts from both types, we compare them in both directions
	my %found_ann;
	my %found_pred;
	if (@ann_trans && @pred_trans){
	    my %potential_partners;
	    my %score; 
	    foreach my $t_ann (@ann_trans){
		foreach my $t_pred (@pred_trans){
		    my ($t_ann_j, $t_pred_j, $over_j ) = compare_transcript_junctions($t_ann, $t_pred);		    
		    
		    # we keep track of which transcripts have been found entirely (junctions) by the other set
		    $found_ann{$t_ann}++   if ($t_ann_j == $over_j); 
		    $found_pred{$t_pred}++ if ($t_pred_j == $over_j);
		}
	    }
	    $over_pred_count += scalar(keys %found_pred);
	    $over_ann_count  += scalar(keys %found_ann);
	}
    }
    my $s = join "\t", ("TRANS", 
			"ann", $ann_count, 
			"pred", $pred_count, 
			"over_ann",  $over_ann_count,
			"over_pred", $over_pred_count,
			"prop_ann",  $over_ann_count/$ann_count,
			"prop_pred", $over_pred_count/$pred_count);
    print $s."\n";
}



##############################################################
sub compare_transcript_junctions{
    my ($t1, $t2) = @_;

    # count the matching junctions
    my %junctions;
    my $overlap = 0;
    my $junc1 = get_junctions($t1);
    my $junc2 = get_junctions($t2);
    foreach my $junc (@$junc1){
	$junctions{$junc}++;
    }
    foreach my $junc (@$junc2){
	$overlap++ if ($junctions{$junc});
    }
    return (scalar(@$junc1), scalar(@$junc2), $overlap);
}

##############################################################  
                                                                                            
sub match_transcript_introns{
    my ($t1, $t2) = @_;
    
    # count the matching junctions                                                                                                            
    my %junctions;
    my $overlap = 0;
    my $introns1 = Transcript::get_sorted_introns($t1);
    my $introns2 = Transcript::get_sorted_introns($t2);
    
    # which intron starts first
  INT1:
    foreach my $int1 (@$introns1){
	my ($chr1, $source1, $feature1, $start1, $end1, $score1, $strand1, $frame1, $t_id1, $g_id1) = @$int1;   
       INT2:
	foreach my $int2 (@$introns2){
	    my ($chr2, $source2, $feature2, $start2, $end2, $score2, $strand2, $frame2, $t_id2, $g_id2) = @$int2;   
	    
	}
    }
}


###########################################################
sub compare_transcript_exons{
    my ($t1, $t2) = @_;

    # count the matching junctions
    my %exons;
    my $overlap = 0;
    my $e1 = get_exons($t1);
    my $e2 = get_exons($t2);
    foreach my $e (@$e1){
	$exons{exon_id($e)}++;
    }
    foreach my $e (@$e2){
	$overlap++ if ($exons{exon_id($e)});
    }
    return (scalar(@$e1), scalar(@$e2), $overlap);
}


################################################
sub exon_id{
    my ($e) = @_;
    # exon: [chr, start, end, strand, t_id, g_id]
    return $e->[0].":".$e->[1]."-".$e->[2].":".$e->[3];
}

############################################################
# this method takes all the pairs with the corresponding scores
# and finds the best pairs using the 'stable-marriage' algorithm.

# The Algorithm has the following condition for solution:
# there are no two elements that they're not paired-up to each other but
# they have better score with each other ( $score_matrix is higher ) 
# than with their current corresponding partners

# the main different between this optimization algorithm
# and a 'best-reciprocal-pairs' approach is that
# 'stable-marriage' produces rather 'best-available-pairs', so it keeps matches which
# not being maximal are still optimal. It warranties only on pair
# per element, and this is the best available one, i.e. ' you like
# C. Schiffer but she likes D. Copperfield more than she likes you so
# you have to stick to the next one in your priority list if available'.


sub stable_marriage{
    my ($ann_trans, $pred_trans, $potential_partners, $score) = @_;

    my %potential_partners_ann;   # potential partners for proteins in annotation list
    my %potential_partners_pred;  # potential partners for proteins in predicion list
    my %married_ann;              # keeps track of the married ones in annotation list
    my %married_pred;             # keeps track of the married ones in prediction list
    my %partner_in_ann;
    my %partner_in_pred;

    # we iterate over all transcripts in group ann
    my @unmarried_ones = @$ann_trans;
    
  MARRIAGE:
    while ( @unmarried_ones ){
      
	# pick one of them
	my $t_ann = shift @unmarried_ones;
    
	# sort potential partners in descending order according to score:
	my @sorted_partners_in_pred = sort_scores( $t_ann, $potential_partners->{$t_ann}, $score);
  
      # go over the partners until you get married or run out of partners  
    PARTNER:
      while( @sorted_partners_in_pred && !defined($married_ann{$t_ann}) ){
	  
	  print "checking partner list for ". print_transcript($t_ann)."\n" if $verbose;
	  my $potential_partner_in_pred = shift( @sorted_partners_in_pred );
	  
	  print "looking at ". print_transcript($potential_partner_in_pred)."\n" if $verbose;
	  # check whether it is already married
	  if ( $married_pred{ $potential_partner_in_pred } 
	       && $married_pred{ $potential_partner_in_pred } == 1 ){
	      
	      # is it married to another target?
	      if ( $partner_in_ann{$potential_partner_in_pred} 
		   &&  !( $partner_in_ann{ $potential_partner_in_pred } eq $t_ann ) ){
		  
		  # is it a 'worse' marriage?
		  if ( $score->{ $partner_in_ann{$potential_partner_in_pred} }{ $potential_partner_in_pred}
		       < $score->{$t_ann}{$potential_partner_in_pred} ){
		      
		      # put the divorced one back into the pool only if it has more potential partners
		      if ( @{$potential_partners->{ $partner_in_ann{$potential_partner_in_pred}}} ){
			  push ( @unmarried_ones, $partner_in_ann{$potential_partner_in_pred} );
		      }
		      
		      # divorce the 'worse partner'
		      print "divorcing ".print_transcript($partner_in_ann{$potential_partner_in_pred})."\n" if $verbose;
		      delete $married_ann{ $partner_in_ann{ $potential_partner_in_pred } };
		      delete $partner_in_pred{ $partner_in_ann{ $potential_partner_in_pred } };
		      delete $partner_in_ann{ $potential_partner_in_pred };
		      
		      # let be happier marriage
		      $married_ann{ $t_ann } = 1;
		      $married_pred{ $potential_partner_in_pred } = 1;
		      $partner_in_ann{ $potential_partner_in_pred } = $t_ann;
		      $partner_in_pred{ $t_ann } = $potential_partner_in_pred;
		      print "new marriage: ".print_transcript($t_ann)." - ".print_transcript($potential_partner_in_pred)."\n" if $verbose;
		      next MARRIAGE;
		      
		  }
		  else{
		      # look at the next potential partner in list2
		      next PARTNER;
		  }
	      }
	      # hmm, this t_pred is married, to whom?
	      elsif ( $partner_in_ann{ $potential_partner_in_pred } eq $t_ann ) {
		  # hey, we have already a happy couple
		  $partner_in_pred{ $t_ann } = $potential_partner_in_pred;
		  next MARRIAGE;
	      }
	      elsif ( !defined( $partner_in_ann{ $potential_partner_in_pred } ) ){
		  # we have a cheater!
		  $married_pred{ $potential_partner_in_pred } = 0;
		  next PARTNER;
	      }
	  }
	  else{
	      
	      # this prot2 ( in list 2 ) is still single, let be marriage:
	      $married_ann{ $t_ann } = 1;
	      $married_pred{ $potential_partner_in_pred } = 1;
	      $partner_in_ann{ $potential_partner_in_pred } = $t_ann;
	      $partner_in_pred{ $t_ann } = $potential_partner_in_pred;
	      print "setting partner for ".print_transcript($t_ann)." to be ".print_transcript($potential_partner_in_pred)."\n" if $verbose;
	      next MARRIAGE;
	  }
      } # end of PARTNER
    }   # end of MARRIAGE
    
    # return the final list of couples
    my @couples;
    foreach my $t_ann (@$ann_trans){
	if ( $partner_in_pred{$t_ann} ){
	    my $couple = [$t_ann, $partner_in_pred{$t_ann}, $score->{$t_ann}{$partner_in_pred{$t_ann}} ];
	    push( @couples, $couple);
	}
    }
    return @couples;
}


############################################################
# sort the potential partners by score in descending order
sub sort_scores{
    my ($t1, $partners, $score) = @_;
    my @sorted_partners =  map {$_->[1]} sort {$b->[0] <=> $a->[0]}  map { [$score->{$t1}{$_}, $_] } @$partners;
    return @sorted_partners;
}

############################################################

sub print_transcript{
    my ($t) = @_;
    my @e = @$t;
    my @coords;
    my $strand = $e[0]->[3];
    my $chr    = $e[0]->[0];
    my $t_id   = $e[0]->[4];
    foreach my $ee (@e){
	push( @coords, $ee->[1]."-".$ee->[2] );
    } 
    my $s = join "\t", ($t_id, $chr, $strand, @coords);
    print $s."\n";
}

##################################################

sub print_Exons{
    my ($t) = @_;
    my @exons = get_exons($t);
    foreach my $e (@exons){
	my $p = join "\t", ("INFO:", "Exon:", @$e);
	print $p."\n";
    }
}

# function test_overlap
# takes two lists and checks the coincidences 
sub test_overlap{
    my ($ref_list, $pred_list) = @_;
    my %ref_list;
    my %pred_list;
    foreach my $e (@$ref_list){
	
	$ref_list{$e}++;
    }
    foreach my $e (@$pred_list){
	$pred_list{$e}++;
    }
    my $overlap = 0;
    foreach my $e (keys %ref_list){
	$overlap++ if $pred_list{$e};
	print "counting $overlap\n" if $verbose;
    }
    return (scalar(keys %ref_list), scalar(keys %pred_list), $overlap);
}

##################################

sub get_sorted_exons{
    my ($t) = @_;
    my @e = sort {$a->[1] <=> $b->[1]} @{$t};
    return @e;
}




# function: get_junctions
# this function reads a transcript
# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
# and outputs all junctions: chr:start-end:strand
# where start of the junction is the end of the left exon
# and end of the junction is the start of the right exons
sub get_junctions{
    my ($t) = @_;
    
    # we get the exons
    my @exons = get_sorted_exons($t);
    
    # extract junctions
    my @junctions;
    for (my $i=0; $i<scalar(@exons) - 1; $i++){
	# exon : [$chr, $start, $end, $strand, $trans_id, $gene_id]
	my $junc_start  = $exons[$i]->[2];
	my $junc_end    = $exons[$i+1]->[1];
	my $junc_chr    = $exons[$i]->[0];
	my $junc_strand = $exons[$i]->[3];
	my $junc        = $junc_chr.":".$junc_start."-".$junc_end.":".$junc_strand;
	push( @junctions, $junc);
    }
    return \@junctions;
}
    
# function: get_5ss_splice_sites
# this function reads a transcript
# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
# and outputs the 5ss_splice_sites: chr:start-end:strand
# for forward: start of the 5ss is the 3' end position of the exon
# and end is the 5' end position of the intron (the G of the GT)
# for reverse: start is the 5' end of the intron (G in reverse) (right-most position of the intron)
# the end is the left most base of the exon
sub get_5ss_splice_sites{
    my ($t) = @_;
    
    # we get the exons
    my @exons = get_sorted_exons($t);

    # strand
    my $strand = $exons[0]->[3];

    # extract junctions
    my @sites;
    if ($strand eq "+"){
	# #####----#####------#####-----#####
	#     --       --         --
	for (my $i=0; $i<scalar(@exons)-1; $i++){
	    # exon : [$chr, $start, $end, $strand, $trans_id, $gene_id]
	    my $site_start  = $exons[$i]->[2];
	    my $site_end    = $site_start + 1;
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[3];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    if ($strand eq "-"){
	# #####----#####------#####-----#####
	#         --         --        --
	for (my $i=1; $i<scalar(@exons); $i++){
	    # exon : [$chr, $start, $end, $strand, $trans_id, $gene_id]
	    my $site_end    = $exons[$i]->[1];
	    my $site_start  = $site_end - 1;
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[3];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    return \@sites;
}


# function: get_3ss_splice_sites
# this function reads a transcript
# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
# and outputs the 3ss_splice_sites: chr:start-end:strand
# start end are the left/right positions at either side of the splice-site
sub get_3ss_splice_sites{
    my ($t) = @_;
    
    # we get the exons
    my @exons = get_sorted_exons($t);

    # strand
    my $strand = $exons[0]->[3];

    # extract junctions
    my @sites;
    if ($strand eq "+"){
	# #####----#####------#####-----#####
	#         --         --        --
	for (my $i=1; $i<scalar(@exons); $i++){
	    # exon : [$chr, $start, $end, $strand, $trans_id, $gene_id]
	    my $site_start  = $exons[$i]->[1] - 1;
	    my $site_end    = $exons[$i]->[1];
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[3];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    if ($strand eq "-"){
	# #####----#####------#####-----#####
	#     --       --         --
	for (my $i=0; $i<scalar(@exons)-1; $i++){
	    # exon : [$chr, $start, $end, $strand, $trans_id, $gene_id]
	    my $site_start  = $exons[$i]->[2];
	    my $site_end    = $site_start + 1;
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[3];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    return \@sites;
}

###########################################

1;
