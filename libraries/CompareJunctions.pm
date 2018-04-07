package CompareJunctions;

use strict;
use libraries::Cluster;

my $verbose  = 0;
my $verbose2 = 0;
my $verbose3 = 0;

####################
# function to compare junctions and splice sites
# from two sets of transcripts
# It will give the proprotion of 
# junctions in each set that overlap with the other
#
sub compare_junctions{
    my ($trans_ann, $trans_pred) = @_;

    # We count the total number of unique
    # junctions and splice-sites, and the overlaps between them
    my $over_junc_count = 0;
    my $pred_junc_count = 0;
    my $ann_junc_count  = 0;

    my $over_5ss_count  = 0;
    my $pred_5ss_count  = 0;
    my $ann_5ss_count   = 0;

    my $over_3ss_count  = 0;
    my $pred_3ss_count  = 0;
    my $ann_3ss_count   = 0;
    
    # cluster transcripts by genomic extent and strand overlap
    my ($clusters, $label) = cluster_transcripts_2sets($trans_ann, $trans_pred);
  
    # for each cluster check the number of matching splicing sites and junctions
    foreach my $cluster (@$clusters){
	my @ann_5ss;
	my @ann_3ss;
	my @ann_junc;
	my @pred_5ss;
	my @pred_3ss;
	my @pred_junc;
	
	foreach my $t (@$cluster){
	    if ( $label->{$t} eq "pred"){
		my $junctions = get_junctions($t);
		push (@pred_junc, @$junctions);
			
		my $list_of_5ss = get_5ss_splice_sites($t);
		push(@pred_5ss, @$list_of_5ss);
		
		my $list_of_3ss = get_3ss_splice_sites($t);
		push(@pred_3ss, @$list_of_3ss);
		
		if ($verbose){
		    foreach my $j (@$junctions){
			print $j."\n";
		    }
		}
	    }
	    elsif ( $label->{$t} eq "ann" ){
		my $junctions = get_junctions($t);
		push (@ann_junc, @$junctions);
		
		my $list_of_5ss = get_5ss_splice_sites($t);
		push(@ann_5ss, @$list_of_5ss);
		
		my $list_of_3ss = get_3ss_splice_sites($t);
		push(@ann_3ss, @$list_of_3ss);
		
		if ($verbose){
		    foreach my $j (@$junctions){
			print $j."\n";
		    }
		}
	    }
	}
	# count the number of unique junctions in each set
	# and the number of common junctions
	my ($ann_j, $pred_j, $over_j) = test_overlap( \@ann_junc, \@pred_junc );
	if ($verbose){
	    my $s = join "\t", ("counting junctions",$ann_j, $pred_j, $over_j);
	    print $s."\n";
	}		
	$over_junc_count += $over_j;
	$pred_junc_count += $pred_j;
	$ann_junc_count  += $ann_j;
	my ($ann_5ss, $pred_5ss, $over_5ss) = test_overlap( \@ann_5ss, \@pred_5ss );
	if ($verbose){
	    my $s = join "\t",("counting 5ss",$ann_5ss, $pred_5ss, $over_5ss);
	    print $s."\n";
	}		
	$over_5ss_count += $over_5ss;
	$pred_5ss_count += $pred_5ss;
	$ann_5ss_count  += $ann_5ss;
	
	my ($ann_3ss, $pred_3ss, $over_3ss) = test_overlap( \@ann_3ss, \@pred_3ss );
	if ($verbose){
	    my $s = join "\t",("counting 3ss",$ann_3ss, $pred_3ss, $over_3ss);
	    print $s."\n";
	}		
	$over_3ss_count += $over_3ss;
	$pred_3ss_count += $pred_3ss;
	$ann_3ss_count  += $ann_3ss;
    }
    my $s1 = join "\t", ("JUNC", 
			 "ann", $ann_junc_count, 
			 "pred", $pred_junc_count, 
			 "over",   $over_junc_count,
			 "prop_ann",  $over_junc_count/$ann_junc_count,
			 "prop_pred", $over_junc_count/$pred_junc_count);
    print $s1."\n";
    
    my $s2 = join "\t", ("5SS", 
			 "ann", $ann_5ss_count, 
			 "pred", $pred_5ss_count, 
			 "over",   $over_5ss_count,
			 "prop_ann",  $over_5ss_count/$ann_5ss_count,
			 "prop_pred", $over_5ss_count/$pred_5ss_count);
    print $s2."\n";
    
    my $s3 = join "\t", ("3SS", 
			 "ann", $ann_3ss_count, 
			 "pred", $pred_3ss_count, 
			 "over",   $over_3ss_count,
			 "prop_ann",  $over_3ss_count/$ann_3ss_count,
			 "prop_pred", $over_3ss_count/$pred_3ss_count);
    print $s3."\n";
    
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
    my ( $clusters, $cluster_start, $cluster_end) = Cluster::cluster_transcripts($trans);
    return ($clusters, $label);
}


# function test_overlap
# takes two lists and checks the coincidences 
sub test_overlap{
    my ($ann_list, $pred_list) = @_;
    my %ann_list;
    my %pred_list;
    foreach my $e (@$ann_list){
	
	$ann_list{$e}++;
    }
    foreach my $e (@$pred_list){
	$pred_list{$e}++;
    }
    my $overlap = 0;
    foreach my $e (keys %ann_list){
	$overlap++ if $pred_list{$e};
	print "counting $overlap\n" if $verbose;
    }
    return (scalar(keys %ann_list), scalar(keys %pred_list), $overlap);
}


sub transcript_length{
    my ($t) = @_;
    my @e = @{$t};
    my $l=0;
    foreach my $e (@e){
	# chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $trans_id, $gene_id) = @$e;
	$l += ($end - $start + 1);
    }
    return $l;
}

sub get_sorted_exons{
    my ($t) = @_;
    # chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
    my @e = sort {$a->[3] <=> $b->[3]} @{$t};
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
	# exon : chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	my $junc_start  = $exons[$i]->[4];
	my $junc_end    = $exons[$i+1]->[3];
	my $junc_chr    = $exons[$i]->[0];
	my $junc_strand = $exons[$i]->[6];
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
    my $strand = $exons[0]->[6];

    # extract junctions
    my @sites;
    if ($strand eq "+"){
	# #####----#####------#####-----#####
	#     --       --         --
	for (my $i=0; $i<scalar(@exons)-1; $i++){
	    # exon : chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	    my $site_start  = $exons[$i]->[4];
	    my $site_end    = $site_start + 1;
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[6];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    if ($strand eq "-"){
	# #####----#####------#####-----#####
	#         --         --        --
	for (my $i=1; $i<scalar(@exons); $i++){
	    # exon : chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	    my $site_end    = $exons[$i]->[3];
	    my $site_start  = $site_end - 1;
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[6];
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
    my $strand = $exons[0]->[6];

    # extract junctions
    my @sites;
    if ($strand eq "+"){
	# #####----#####------#####-----#####
	#         --         --        --
	for (my $i=1; $i<scalar(@exons); $i++){
	    # exon : chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	    my $site_start  = $exons[$i]->[3] - 1;
	    my $site_end    = $exons[$i]->[3];
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[6];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    if ($strand eq "-"){
	# #####----#####------#####-----#####
	#     --       --         --
	for (my $i=0; $i<scalar(@exons)-1; $i++){
	    # exon : chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
	    my $site_start  = $exons[$i]->[4];
	    my $site_end    = $site_start + 1;
	    my $site_chr    = $exons[$i]->[0];
	    my $site_strand = $exons[$i]->[6];
	    my $site        = $site_chr.":".$site_start."-".$site_end.":".$site_strand;
	    push( @sites, $site);
	}
    }
    return \@sites;
}

###########################################


1;
