#!/usr/bin/perl -w

use strict;

use libraries::Cluster;
use libraries::FileReader;
use Getopt::Long;

my $verbose  = 0;
my $verbose2 = 0;
my $verbose3 = 0;

my ($file_pred, $file_ann, $format_pred, $format_ann, $help);

# check options
my $result = GetOptions(
        'p|prediction=s' => \$file_pred,
        'a|annotation=s' => \$file_ann,
        'f|format_pred=s'=> \$format_pred,
        'u|format_ann=s' => \$format_ann,
        'h|help'         => \$help
    );

print_usage() if( $help);

unless ( $file_pred && $file_ann && $format_pred && $format_ann ){
    print STDERR "Error: needs to specify some input information\n";
    print_usage();
    exit(0);
}

unless( ($format_ann eq "GFF" || $format_ann eq "GTF") 
	&& 
	($format_pred eq "GFF" || $format_pred eq "GTF")){
    print STDERR "ERROR: formats supported so far: GFF, GTF\n";
    print_usage();
}

sub print_usage{
    print STDERR "usage: perl compare_mappings_and_annotations.pl -a <annotations> -p <predictions> -f <format-annotations> -u <format-predictions>\n";
    print STDERR "\n";
    print STDERR "-a | --annotation:\tFile with the transcript annotations\n";
    print STDERR "-a | --annotation:\tFile with the transcript annotations\n";
    print STDERR "-f | --format_ann:\tFormat of the annotation file: GTF, GFF, PAF\n";
    print STDERR "-u | --format_pred:\tFormat of the prediction file: GTF, GFF, PAF (it can be different from the annotation file)\n";
    print STDERR "-h | --help:\tPrint this help\n";
    print STDERR "\n";
    exit(0);
}


my %strand;
my %chr;

# read GFF lines in annotation file 1
# store all transcripts and all exons in each transcript
my $trans_ann = FileReader::read_file($file_ann, $format_ann, $verbose);
my %trans_ann = %$trans_ann;

# read GFF lines in annotation file 2
# store all transcripts and all exons in each transcript
my $trans_pred = FileReader::read_file($file_pred, $format_pred, $verbose);
my %trans_pred = %$trans_pred;

#################
# make transcripts out of exons:
my $chr;
my $strand;

my %chr1;
foreach my $t_id (keys %trans_ann){
    # sort exons by the start coordinate in ascending order
    my @exons = sort {$a->[1] <=> $b->[1]} @{$trans_ann{$t_id}};

    # define transcript attributes
    my $t_chr    = $exons[0]->[0];
    my $t_start  = $exons[0]->[1];
    my $t_end    = $exons[-1]->[2];
    my $t_strand = $exons[0]->[3];
    my $g_id     = $exons[0]->[5];

    # these are the "ref"
    my $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "ref"];
    print "storing trans: $t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id\n" if $verbose;
    push( @{$chr1{$t_chr}{$t_strand}}, $tran);
    $chr{$t_chr}++;
    $strand{$t_strand}++;    
}

my %chr2;
foreach my $t_id (keys %trans_pred){
    # sort exons by the start coordinate in ascending order
    my @exons = sort {$a->[1] <=> $b->[1]} @{$trans_pred{$t_id}};

    # define transcript attributes
    my $t_chr    = $exons[0]->[0];
    my $t_start  = $exons[0]->[1];
    my $t_end    = $exons[-1]->[2];
    my $t_strand = $exons[0]->[3];
    my $g_id     = $exons[0]->[5];

    # these are the "pred"
    my $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
    print "storing trans: $t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id\n" if $verbose;
    push( @{$chr2{$t_chr}{$t_strand}}, $tran);
    $chr{$t_chr}++;
    $strand{$t_strand}++;
}
    
# Now we perform the clustering:
my %seen;
my $color;
# for each chr ID (chromosome and haplotypes), and each strand, we cluster the transcripts
# according to their genomic extent
my $overlap_junc_count = 0;
my $pred_junc_count    = 0;
my $ref_junc_count     = 0;

foreach my $chr ( keys %chr ){
    foreach my $strand ( keys %strand ){

	# all my transcripts
	my @all_trans;
	push(@all_trans, @{$chr1{$chr}{$strand}}) if ($chr1{$chr}{$strand});
	push(@all_trans, @{$chr2{$chr}{$strand}}) if ($chr2{$chr}{$strand});

        # cluster transcripts according to their genomic extent
	print "clustering ".scalar(@all_trans)." transcripts\n" if $verbose3;
        my ( $clusters, $cluster_start, $cluster_end) = Cluster::cluster( \@all_trans, 1, 2 );

	# for each cluster check the number of matching splicing sites and junctions
        foreach my $cluster (@$clusters){
	    my @ref_5ss_sites;
	    my @ref_3ss_sites;
	    my @ref_junctions;
	    my @pred_5ss_sites;
	    my @pred_3ss_sites;
	    my @pred_junctions;
	    
	    foreach my $t (@$cluster){
		if ($verbose){
		    # $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
		    my $s = join "\t", @$t;
		    print $s."\n";
		}
		if ( $t->[6] eq "pred"){
		    my $junctions = get_junctions($t);
		    push (@pred_junctions, @$junctions);
		    
		    if ($verbose){
			foreach my $j (@$junctions){
			    print $j."\n";
			}
		    }
		}
		if ( $t->[6] eq "ref"){
		    my $junctions = get_junctions($t);
		    push (@ref_junctions, @$junctions);
		    if ($verbose){
			foreach my $j (@$junctions){
			    print $j."\n";
			}
		    }
		}
	    }
	    my ($ref_j, $pred_j, $overlap) = test_overlap( \@ref_junctions, \@pred_junctions );
	    if ($verbose){
		my $s = join "\t",("counting junctions",$ref_j, $pred_j, $overlap);
		print $s."\n";
	    }

	    $overlap_junc_count += $overlap;
	    $pred_junc_count    += $pred_j;
	    $ref_junc_count     += $ref_j;
	}
    }
}

my $s1 = join "\t", ("Junctions", 
		     "reference", $ref_junc_count, 
		     "predicted", $pred_junc_count, 
		     "overlap",   $overlap_junc_count,
		     "proportion_of_ref",  $overlap_junc_count/$ref_junc_count,
		     "proportion_of_pred", $overlap_junc_count/$pred_junc_count);
print $s1."\n";




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


sub transcript_length{
    my ($t) = @_;
    my @exons = get_exons($t);
    my $l=0;
    foreach my $e (@exons){
	my ($chr, $start, $end, $strand, $trans_id, $gene_id) = @$e;
	$l += ($end - $start + 1);
    }
    return $l;
}

sub get_exons{
    my ($t) = @_;
    my @exons;
    @exons = sort {$a->[1] <=> $b->[1]} @{$trans_ann{$t->[4]}} if ($trans_ann{$t->[4]});
    @exons = sort {$a->[1] <=> $b->[1]} @{$trans_pred{$t->[4]}} if ($trans_pred{$t->[4]});
    return @exons;
}

# function: get_junctions
# this function reads a transcript
# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
# and outputs the junctions: chr:start-end:strand
# where start of the junction is the end of the left exon
# and end of the junction is the start of the right exons
sub get_junctions{
    my ($t) = @_;
    
    # we get the exons
    my @exons = get_exons($t);

    # we extract sorted exons
    my @sorted_exons = sort {$a->[1] <=> $b->[1]} @exons;

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
    my @exons = get_exons($t);

    # strand
    my $strand = $exons[0]->[3];

    # we extract sorted exons
    my @sorted_exons = sort {$a->[1] <=> $b->[1]} @exons;

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
    my @exons = get_exons($t);

    # strand
    my $strand = $exons[0]->[3];

    # we extract sorted exons
    my @sorted_exons = sort {$a->[1] <=> $b->[1]} @exons;

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

###########################################

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

###########################################

sub build_genes{
    my ($loci, $s, $e) = @_;
    
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

##########################
# returns 1 if the transcripts
# have exons with at least
# a shared splice-site

sub test_exon_overlap{
    my ($tran1, $tran2) = @_;
    my @exons1 = get_exons($tran1);
    my @exons2 = get_exons($tran2);
    
    # each exon is a ref to: [$chr, $start, $end, $strand, $frame, $exon_label, $trans_id, $gene_id ]                                           
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


############################
# algorithm to recover the groups of
# linked transcripts

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
