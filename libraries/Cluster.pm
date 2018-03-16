#!/usr/bin/perl w

use strict;

sub cluster{
    my ($loci, $s, $e, $verbose) = @_;
    
    # $loci is an arrayref to arrayrefs
    my @loci = @$loci;
    
    # we generate cluster of loci
    my $locus_clusters;
        
    # sort loci by start coordinate in ascending order
    # each loci is of the form  [$locus_chr,$locus_start,$locus_end,$locus_id]
    my @sorted_loci = sort { $a->[$s] <=> $b->[$s] } @loci;

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

        
