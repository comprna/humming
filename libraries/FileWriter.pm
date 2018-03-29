package FileWriter;

use strict;


sub write_file{
    my ($trans, $format, $output) = @_;

    if ($format eq "GTF"){
        return write_GTF($trans, $output);
    }
    elsif( $format eq "GFF"){
        return write_GFF($trans, $output);
    }
    elsif( $format eq "BED12"){
        return write_BD12($trans, $output);
    }
    elsif( $format eq "PAF"){
        return write_PAF($trans, $output);
    }
}


sub write_GTF{
    my ($trans, $output) = @_;
    
    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	foreach my $e (get_sorted_exons($t)){
	    my ($chr, $start, $end, $strand, $t_id, $g_id) = @$e;
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


sub write_GFF{
    my ($trans, $output) = @_;

    open(OUT,">$output") or die("cannot open file $output for writing");
    foreach my $t (@$trans){
	foreach my $e (get_sorted_exons($t)){
	    my ($chr, $start, $end, $strand, $t_id, $g_id) = @$e;
	    my $source  = "humming";
	    my $feature = "exon";
	    my $frame   = ".";
	    my $score   = 100;
	    my $col9    = $t_id;
	    my $s = join "\t", ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $col9);
	    print OUT $s."\n";
	}
    }
}

sub write_genes_in_GTF{
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
