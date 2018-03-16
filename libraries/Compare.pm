package Compare;

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
    my ($trans_ref, $trans_pred) = @_;

    # label transcripts according to their type
    my $label;
    
    # make transcripts out of exons:
    my $chr;
    my $strand;

    my %chr_ref;
    foreach my $t_id (keys %$trans_ref){
    # sort exons by the start coordinate in ascending order
	my @exons = sort {$a->[1] <=> $b->[1]} @{$trans_ref->{$t_id}};
	
	# define transcript attributes
	my $t_chr    = $exons[0]->[0];
	my $t_start  = $exons[0]->[1];
	my $t_end    = $exons[-1]->[2];
	my $t_strand = $exons[0]->[3];
	my $g_id     = $exons[0]->[5];
	
	# each transcript is a pointer to a list:
	my $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id];
	print "storing trans: $t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id\n" if $verbose;
	push( @{$chr_ref{$t_chr}{$t_strand}}, $tran);
	$label->{$tran} = "ref";
	$chr->{$t_chr}++;
	$strand->{$t_strand}++;    
    }
    
    my %chr_pred;
    foreach my $t_id (keys %$trans_pred){
	# sort exons by the start coordinate in ascending order
	my @exons = sort {$a->[1] <=> $b->[1]} @{$trans_pred->{$t_id}};
	
	# define transcript attributes
	my $t_chr    = $exons[0]->[0];
	my $t_start  = $exons[0]->[1];
	my $t_end    = $exons[-1]->[2];
	my $t_strand = $exons[0]->[3];
	my $g_id     = $exons[0]->[5];
	
	# each transcript is a pointer to a list:
	my $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id];
	print "storing trans: $t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id\n" if $verbose;
	push( @{$chr_pred{$t_chr}{$t_strand}}, $tran);
	$label->{$tran} = "pred";
	$chr->{$t_chr}++;
	$strand->{$t_strand}++;
    }

    # We count the total number of unique
    # junctions and splice-sites, and the overlaps between them
    my $over_junc_count = 0;
    my $pred_junc_count    = 0;
    my $ref_junc_count     = 0;

    my $over_5ss_count  = 0;
    my $pred_5ss_count     = 0;
    my $ref_5ss_count      = 0;

    my $over_3ss_count  = 0;
    my $pred_3ss_count     = 0;
    my $ref_3ss_count      = 0;
    
    # For each chromosome and each strand, we cluster the transcripts
    # according to their genomic extent
    foreach my $chr ( keys %$chr ){
	foreach my $strand ( keys %$strand ){
	    
	    # all my transcripts
	    my @all_trans;
	    push(@all_trans, @{$chr_ref{$chr}{$strand}}) if ($chr_ref{$chr}{$strand});
	    push(@all_trans, @{$chr_pred{$chr}{$strand}}) if ($chr_pred{$chr}{$strand});
	    
	    # cluster transcripts according to their genomic extent
	    print "clustering ".scalar(@all_trans)." transcripts\n" if $verbose3;
	    my ( $clusters, $cluster_start, $cluster_end) = Cluster::cluster( \@all_trans, 1, 2 );
	    
	    # for each cluster check the number of matching splicing sites and junctions
	    foreach my $cluster (@$clusters){
		my @ref_5ss;
		my @ref_3ss;
		my @ref_junc;
		my @pred_5ss;
		my @pred_3ss;
		my @pred_junc;
		
		foreach my $t (@$cluster){
		    if ($verbose){
			# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id];
			my $s = join "\t", @$t;
			print $s."\n";
		    }
		    if ( $label->{$t} eq "pred"){
			my $junctions = get_junctions($t, $trans_pred);
			push (@pred_junc, @$junctions);
			
			my $list_of_5ss = get_5ss_splice_sites($t, $trans_pred);
			push(@pred_5ss, @$list_of_5ss);
			
			my $list_of_3ss = get_3ss_splice_sites($t, $trans_pred);
			push(@pred_3ss, @$list_of_3ss);
			
			if ($verbose){
			    foreach my $j (@$junctions){
				print $j."\n";
			    }
			}
		    }
		    elsif ( $label->{$t} eq "ref" ){
			my $junctions = get_junctions($t, $trans_ref);
			push (@ref_junc, @$junctions);
			
			my $list_of_5ss = get_5ss_splice_sites($t, $trans_ref);
			push(@ref_5ss, @$list_of_5ss);
			
			my $list_of_3ss = get_3ss_splice_sites($t, $trans_ref);
			push(@ref_3ss, @$list_of_3ss);
			
			if ($verbose){
			    foreach my $j (@$junctions){
				print $j."\n";
			    }
			}
		    }
		}
		# count the number of unique junctions in each set
		# and the number of common junctions
		my ($ref_j, $pred_j, $over_j) = test_overlap( \@ref_junc, \@pred_junc );
		if ($verbose){
		    my $s = join "\t", ("counting junctions",$ref_j, $pred_j, $over_j);
		    print $s."\n";
		}		
		$over_junc_count += $over_j;
		$pred_junc_count    += $pred_j;
		$ref_junc_count     += $ref_j;
		my ($ref_5ss, $pred_5ss, $over_5ss) = test_overlap( \@ref_5ss, \@pred_5ss );
		if ($verbose){
		    my $s = join "\t",("counting 5ss",$ref_5ss, $pred_5ss, $over_5ss);
		    print $s."\n";
		}		
		$over_5ss_count += $over_5ss;
		$pred_5ss_count    += $pred_5ss;
		$ref_5ss_count     += $ref_5ss;
		
		my ($ref_3ss, $pred_3ss, $over_3ss) = test_overlap( \@ref_3ss, \@pred_3ss );
		if ($verbose){
		    my $s = join "\t",("counting 3ss",$ref_3ss, $pred_3ss, $over_3ss);
		    print $s."\n";
		}		
		$over_3ss_count += $over_3ss;
		$pred_3ss_count    += $pred_3ss;
		$ref_3ss_count     += $ref_3ss;
	    }
	}
    }
    my $s1 = join "\t", ("JUNC", 
			 "ref", $ref_junc_count, 
			 "pred", $pred_junc_count, 
			 "over",   $over_junc_count,
			 "prop_ref",  $over_junc_count/$ref_junc_count,
			 "prop_pred", $over_junc_count/$pred_junc_count);
    print $s1."\n";

    my $s2 = join "\t", ("5SS", 
			 "ref", $ref_5ss_count, 
			 "pred", $pred_5ss_count, 
			 "over",   $over_5ss_count,
			 "prop_ref",  $over_5ss_count/$ref_5ss_count,
			 "prop_pred", $over_5ss_count/$pred_5ss_count);
    print $s2."\n";

    my $s3 = join "\t", ("3SS", 
			 "ref", $ref_3ss_count, 
			 "pred", $pred_3ss_count, 
			 "over",   $over_3ss_count,
			 "prop_ref",  $over_3ss_count/$ref_3ss_count,
			 "prop_pred", $over_3ss_count/$pred_3ss_count);
    print $s3."\n";

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
    my ($t, $trans_ref) = @_;
    my @exons;
    @exons = sort {$a->[1] <=> $b->[1]} @{$trans_ref->{$t->[4]}} if ($trans_ref->{$t->[4]});
    return @exons;
}


# function: get_junctions
# this function reads a transcript
# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id, "pred"];
# and outputs all junctions: chr:start-end:strand
# where start of the junction is the end of the left exon
# and end of the junction is the start of the right exons
sub get_junctions{
    my ($t, $trans) = @_;
    
    # we get the exons
    my @exons = get_exons($t, $trans);

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
    my ($t, $trans) = @_;
    
    # we get the exons
    my @exons = get_exons($t,$trans);

    # strand
    my $strand = $exons[0]->[3];

    # we extract exons sorted from left to right
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
    my ($t, $trans) = @_;
    
    # we get the exons
    my @exons = get_exons($t, $trans);

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

1;
