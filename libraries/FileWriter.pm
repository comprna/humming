package FileWriter;

use strict;

sub write_GTF{
    my ($genes, $label, $output) = @_;

    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $g (@$genes){
	foreach my $t (@$g){
	    foreach my $e (get_sorted_exons($t)){
		my ($chr, $start, $end, $strand, $t_id, $g_id) = @$e;
		$g_id = $label->{$g};
		my $source  = "humming";
		my $feature = "exon";
		my $frame   = ".";
		my $score   = 100;
		my $col9    = "gene_id \"".$g_id."\"; transcript_id \"".$t_id."\";";
		my $s = join "\t", ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $col9);
		print OUT $s."\n";
	    }
	}
    }
}

sub get_sorted_exons{
    my ($t) = @_;
    my @e = sort {$a->[1] <=> $b->[1]} @{$t};
    return @e;
}



1;
