package CompareGenes;

use strict;
use libraries::Cluster;
use libraries::Gene;

my $verbose  = 0;
my $verbose2 = 0;
my $verbose3 = 0;

####################
# function to compare genes
# from two sets of transcripts
# One gene is considered "found" if prediction and annotation
# overlap in strand and genomic range by any amount
#
# We first build clusters with all transcripts by genomic extension
# then build genes from the transcript clusters 
# (this still allows genes in the introns of other genes and nested genes)
# and then compute how many genes have create gene loci (extensions where a gene lives)
# and compare gene loci directly

sub compare_genes{
    my ($trans_ann, $trans_pred) = @_;

    # We count the total number of genes and the overlaps between them
    my $over_count = 0;
    my $pred_count    = 0;
    my $ann_count     = 0;

    # cluster transcripts by genomic extent and strand overlap                                                                                    
    my ($clusters, $label) = cluster_transcripts_2sets($trans_ann, $trans_pred);

    # from each transcript cluster, we build genes
    foreach my $cluster (@$clusters){
	print "Here: obtained cluster with ".scalar(@$cluster)." transcripts\n" if $verbose;
	
	# genes are defined as sets of transcripts that share at least one splice-site
	# this will split clusters into multiple genes if transcripts do not overlap
	my $genes = Gene::build_genes($cluster);
	print "obtained ".scalar(@$genes)." genes\n" if $verbose;
	
	# count in each gene, whether there are transcripts from different types
	foreach my $gene (@$genes){
	    
	    # each gene is a cluster of transcripts
	    my %these_labels;
	    foreach my $t (@$gene){
		$these_labels{$label->{$t}}++;
	    }
	    if ($these_labels{"ann"}){
		$ann_count++;
	    }
	    if($these_labels{"pred"}){
		$pred_count++;
	    }
	    if ($these_labels{"ann"} && $these_labels{"pred"}){
		$over_count++;
	    }
	}
    }
    my $s1 = join "\t", ("GENE", 
			 "ann", $ann_count, 
			 "pred", $pred_count, 
			 "over",   $over_count,
			 "prop_ann",  $over_count/$ann_count,
			 "prop_pred", $over_count/$pred_count);
    print $s1."\n";
}



sub transcript_length{
    my ($t) = @_;
    my @exons = get_exons($t);
    my $l=0;
    foreach my $e (@exons){
	#exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                                 
	my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $trans_id, $gene_id) = @$e;
	$l += ($end - $start + 1);
    }
    return $l;
}

sub get_exons{
    my ($t) = @_;
    #chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                                         
    return sort {$a->[3] <=> $b->[3]} @{$t};
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
