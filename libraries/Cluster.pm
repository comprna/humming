package Cluster;

use strict;

my $verbose = 0;
my $verbose2 = 0;
my $verbose3 = 0;

# The funcion cluster() reads a pointer to a list (arrayref) 
# of pointers to lists (arrayrfs). Each of those lists will hold the information 
# about the ranges (e.g. exons, gene loci, transcript genomic extensions). 
# After sorting the ranges, the clustering is performed linearly in the number of elements O(n).
# Additionally it reads the positions in the element lists of the start and end
# to be used for the clustering. This function is agnostic of reference (chromosome)
# and strand (genome). These issues must be sorted out outside the function.

sub cluster{
    my ($loci, $s, $e, $verbose) = @_;
    
    # $loci is an arrayref to arrayrefs
    
    # we generate cluster of loci
    my $locus_clusters;
        
    # sort loci by start coordinate in ascending order
    my @sorted_loci = sort { $a->[$s] <=> $b->[$s] } @$loci;

    # start and end coordinates of the clusters
    my @cluster_start;
    my @cluster_end;
    
    # Create the first locus_cluster with the first locus
    my $locus_cluster = [];
    push( @$locus_clusters, $locus_cluster);
    
    # we go over all of them in sorted order (from left to right)
    my $cluster_count = 0;
    my $count = 0;    

    LOCUS:
    foreach my $locus ( @sorted_loci ){
        if ($count == 0){
            # Create the first locus_cluster with the first locus
            push( @$locus_cluster, $locus);
            $cluster_start[0]  =  $locus->[$s];
            $cluster_end[0]    =  $locus->[$e];
            $count++;
            if ($verbose){
                print  "cluster_count: $cluster_count\n";
                print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
            }
            next LOCUS;
        }
        if ($verbose){
            print  "cluster_count: $cluster_count\n";
            print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
        }
	# test overlap
        if ( !( $locus->[$e] < $cluster_start[$cluster_count] || $locus->[$s] > $cluster_end[$cluster_count]) ){

            # add locus to current cluster
            print "added\n" if $verbose;
            push( @$locus_cluster, $locus );
            
            # update start and end of cluster if necessary
            if ($locus->[$s] < $cluster_start[$cluster_count]) {
                $cluster_start[$cluster_count] = $locus->[$s];
            }
            if ($locus->[$e] > $cluster_end[$cluster_count]) {
                $cluster_end[$cluster_count]   = $locus->[$e];
            }
        }
        else{
            # if there is no overlap, proceed with the next one:
            # create new cluster (same variable name, new memory address!!!!)
            $locus_cluster = [];
            push( @$locus_clusters, $locus_cluster);  # store ref first, fill later

            # add locus in new cluster
            push( @$locus_cluster, $locus );
            $cluster_count++;
            $cluster_start[$cluster_count] = $locus->[$s];
            $cluster_end[$cluster_count]   = $locus->[$e];
        }
    }
    # returns an arrayref of clusters (each cluster is an arrayref of the same elements as in @loci)
    # and the arrayrefs to the start and end of each cluster
    return( $locus_clusters, \@cluster_start, \@cluster_end);
}


sub cluster_transcripts{
    my ($loci) = @_;
    
    # $loci is an arrayref to arrayrefs
    
    # we generate cluster of loci
    my $locus_clusters;
        
    # sort loci by start coordinate in ascending order
    my @sorted_loci = sort { transcript_start($a) <=> transcript_start($b) } @$loci;

    # start and end coordinates of the clusters
    my @cluster_start;
    my @cluster_end;
    
    # Create the first locus_cluster with the first locus
    my $locus_cluster = [];
    push( @$locus_clusters, $locus_cluster);
    
    # we go over all of them in sorted order (from left to right)
    my $cluster_count = 0;
    my $count = 0;    

    LOCUS:
    foreach my $locus ( @sorted_loci ){
        if ($count == 0){
            # Create the first locus_cluster with the first locus
            push( @$locus_cluster, $locus);
            $cluster_start[0]  =  transcript_start($locus);
            $cluster_end[0]    =  transcript_end($locus);
            $count++;
            if ($verbose){
                print  "cluster_count: $cluster_count\n";
                print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
            }
            next LOCUS;
        }
        if ($verbose){
            print  "cluster_count: $cluster_count\n";
            print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
        }
	# test overlap
        if ( !( transcript_end($locus) < $cluster_start[$cluster_count] || transcript_start($locus) > $cluster_end[$cluster_count]) ){

            # add locus to current cluster
            print "added\n" if $verbose;
            push( @$locus_cluster, $locus );
            
            # update start and end of cluster if necessary
            if (transcript_start($locus) < $cluster_start[$cluster_count]) {
                $cluster_start[$cluster_count] = transcript_start($locus);
            }
            if (transcript_end($locus) > $cluster_end[$cluster_count]) {
                $cluster_end[$cluster_count]   = transcript_end($locus);
            }
        }
        else{
            # if there is no overlap, proceed with the next one:
            # create new cluster (same variable name, new memory address!!!!)
            $locus_cluster = [];
            push( @$locus_clusters, $locus_cluster);  # store ref first, fill later

            # add locus in new cluster
            push( @$locus_cluster, $locus );
            $cluster_count++;
            $cluster_start[$cluster_count] = transcript_start($locus);
            $cluster_end[$cluster_count]   = transcript_end($locus);
        }
    }
    # returns an arrayref of clusters (each cluster is an arrayref of the same elements as in @loci)
    # and the arrayrefs to the start and end of each cluster
    return( $locus_clusters, \@cluster_start, \@cluster_end);
}


sub cluster_transcripts_2sets{
    my ($t_ann, $t_pred) = @_;
    my $label;
    my $trans;
    foreach my $t (@$t_ann){
	$label->{$t} = "ann";
        push(@$trans,$t);
    }
    foreach my $t (@$t_pred){
	$label->{$t} = "pred";
        push(@$trans, $t);
    }
    my ( $clusters, $cluster_start, $cluster_end) = cluster_transcripts($trans);
    return ($clusters, $label);
}
###########################################
#
# This method cluster ranges by coordinate proximity
#
sub cluster_by_proximity{
    my ($loci, $s, $e, $distance, $verbose) = @_;
    
    # $loci is an arrayref to arrayrefs
    
    # we generate cluster of loci
    my $locus_clusters;
        
    # sort loci by start coordinate in ascending order
    my @sorted_loci = sort { $a->[$s] <=> $b->[$s] } @$loci;
    
    if ($verbose){
	foreach my $l (@sorted_loci){
	    print "reading: @$l\n";
	}
    }

    # start and end coordinates of the clusters
    my @cluster_start;
    my @cluster_end;
    
    # Create the first locus_cluster with the first locus
    my $locus_cluster = [];
    push( @$locus_clusters, $locus_cluster);
    
    # we go over all of them in sorted order (from left to right)
    my $cluster_count = 0;
    my $count = 0;    

    LOCUS:
    foreach my $locus ( @sorted_loci ){
        if ($count == 0){
            # Create the first locus_cluster with the first locus
            push( @$locus_cluster, $locus);
            $cluster_start[0]  =  $locus->[$s];
            $cluster_end[0]    =  $locus->[$e];
            $count++;
            if ($verbose){
                print  "cluster_count: $cluster_count\n";
                print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
            }
            next LOCUS;
        }
        if ($verbose){
	    print "testing range: ". $locus->[$s]. " - ". $locus->[$e]."\n";
	    print "cluster_count: $cluster_count\n";
            print "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
	    print "distance: ".($locus->[$s] - $cluster_end[$cluster_count])."\n";
        }
	# test proximity
	#     s  e                                         s  e
	#     ####       cs      ce         cs     ce      ####
	#         <----->##########         #########<---->
        #if ( $cluster_start[$cluster_count] - $locus->[$e] <= $distance ||  $locus->[$s] - $cluster_end[$cluster_count] < $distance ){
	if ( abs($locus->[$s] - $cluster_end[$cluster_count]) < $distance ){

            # add locus to current cluster
            print "added\n" if $verbose;
            push( @$locus_cluster, $locus );
            
            # update start and end of cluster if necessary
            if ($locus->[$s] < $cluster_start[$cluster_count]) {
                $cluster_start[$cluster_count] = $locus->[$s];
            }
            if ($locus->[$e] > $cluster_end[$cluster_count]) {
                $cluster_end[$cluster_count]   = $locus->[$e];
            }
        }
        else{
            # if there is no overlap, proceed with the next one:
            # create new cluster (same variable name, new memory address!!!!)
            $locus_cluster = [];
            push( @$locus_clusters, $locus_cluster);  # store ref first, fill later
	    
            # add locus in new cluster
            push( @$locus_cluster, $locus );
            $cluster_count++;
            $cluster_start[$cluster_count] = $locus->[$s];
            $cluster_end[$cluster_count]   = $locus->[$e];
        }
    }
    # returns an arrayref of clusters (each cluster is an arrayref of the same elements as in @loci)
    # and the arrayrefs to the start and end of each cluster
    return( $locus_clusters, \@cluster_start, \@cluster_end);
}

###########################################

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


sub transcript_start{
    my ($t) = @_;
    # exon: chromosome, source, feature, start, end, score, strand, frame, transcript ID, gene ID                                                           
    my @e = sort {$a->[3] <=> $b->[3]} @{$t};
    return $e[0]->[3];
}

sub transcript_end{
    my ($t) = @_;
    # exon: chromosome, source, feature, start, end, score, strand, frame, transcript ID, gene ID                                                          
    my @e = sort {$a->[3] <=> $b->[3]} @{$t};
    return $e[-1]->[4];
}

###########################################


1;
