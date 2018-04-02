package Transcript;

use strict;
use libraries::Cluster;

my $verbose = 0;
my $verbose2 = 0;


# Function: merge_transcripts                                                                                                                          
# It reads a reference to a list of transcripts                                                                                                          
#                                                                                                                                                       
# Each trascript is a reference to a list with the following data:                                                                                      
# $tran = [$t_chr, $t_start, $t_end, $t_strand, $t_id, $g_id];                                                                                          
#    

sub merge_transcripts{
    my ($cluster, $merge_info) = @_;
    my %t_id;
    my %redundant;
    foreach my $t_a ( @$cluster ){
        foreach my $t_b (@$cluster ){
            next if ($t_a == $t_b);
            next if ($redundant{$t_a} || $redundant{$t_b});
            my $merge = test_merge($t_a, $t_b);

            if($merge){
                
		if ($verbose){
                    print "INFO:\twe will merge:\n";
                    print_Exons($t_a);
                    print "INFO:\tand\n";
                    print_Exons($t_b);
                }
                $redundant{$t_b}++ if (transcript_length($t_a) > transcript_length($t_b) );
                $redundant{$t_a}++ if (transcript_length($t_a) <= transcript_length($t_b) );
		
		# we link all the IDs to the transcript that can be merged to each other
		my @e_a = @$t_a;
		my @e_b = @$t_b;
		# chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
		my $t_id_a = $e_a[8];
		my $t_id_b = $e_b[8];
		push( @{$t_id{$t_a}}, $t_id_a);
		push( @{$t_id{$t_a}}, $t_id_b);
		push( @{$t_id{$t_b}}, $t_id_a);
		push( @{$t_id{$t_b}}, $t_id_b);
            }
        }
    }
    my $new_cluster;
    foreach my $t ( @$cluster ){
        unless ($redundant{$t}){
            push( @$new_cluster, $t);
            if ($merge_info){
                print "INFO:\tkeeping\n";
                print_Exons($t);
            }
        }
    }
    
    my %new_t_id;
    foreach my $t (keys %t_id){
	my @ids =  @{$t_id{$t}};
	my %tmp;
	foreach my $id (@ids){
	    $tmp{$id}++;
	}
	my @uniq_ids = sort {$a cmp $b} keys %tmp;
	my $uniq_id  = join ":",@uniq_ids;
	$new_t_id{$t} = $uniq_id;
    }
    return ($new_cluster, \%new_t_id);
}


sub transcript_length{
    my ($t) = @_;
    my $l=0;
    foreach my $e (@$t){
	#exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
        my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $trans_id, $gene_id) = @$e;
        $l += ($end - $start + 1);
    }
    return $l;
}

sub get_sorted_exons{
    my ($t) = @_;
    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id
    return sort {$a->[3] <=> $b->[3]} @{$t};
}


sub get_sorted_introns{
    my ($t) = @_;
    #exon: chromosome, source, feature, start, end, score, strand, frame, t_id, g_id                                                                       
    my @exons = sort {$a->[3] <=> $b->[3]} @{$t};
    my @introns;
    for(my $i=0; $i<scalar(@exons)-1; $i++){
	my $intron_start = $exons[$i]->[4] + 1;
	my $intron_end   = $exons[$i+1]->[3] - 1;
	my $intron       = [ $exons[$i]->[0], $exons[$i]->[1], "intron", $intron_start, $intron_end, 
			     $exons[$i]->[5], $exons[$i]->[6], $exons[$i]->[7], $exons[$i]->[8], $exons[$i]->[9] ];
	push( @introns, $intron);
    }
    return \@introns;
}

sub test_merge{
    my ($t1, $t2) = @_;
    my @exons1 = get_sorted_exons($t1);
    my @exons2 = get_sorted_exons($t2);

    # let's test first the trivial case:                                                                                                                
    if (scalar(@exons1)==1 && scalar(@exons2)==1){
        my ($chr1, $start1, $end1, $strand1, $trans_id1, $gene_id1) = @{$exons1[0]};
        my ($chr2, $start2, $end2, $strand2, $trans_id2, $gene_id2) = @{$exons2[0]};
        if ( !($start2>$end1) && !($start1>$end2) ){
            return 1;
        }
    }
    elsif (scalar(@exons1)==1 || scalar(@exons2)==1){
        my @starts1;
        my @starts2;
        my @ends1;
        my @ends2;
        for(my $k=0; $k<scalar(@exons1); $k++){
            my ($chr, $start, $end, $strand, $trans_id, $gene_id) = @{$exons1[$k]};
            push(@starts1, $start);
            push(@ends1,   $end);
        }
        for(my $k=0; $k<scalar(@exons2); $k++){
            my ($chr, $start, $end, $strand, $trans_id, $gene_id) = @{$exons2[$k]};
            push(@starts2, $start);
            push(@ends2,   $end);
        }
        my $merge = 1;
        my $init  = 0;
        for(my $i=0; $i<scalar(@starts1); $i++){
            for(my $j = $init; $j<scalar(@starts2); $j++){
                # only need to worry about the single exon                                                                                                
                if ( !($starts2[$j]>$ends1[$i]) && !($starts1[$i]>$ends2[$j]) ){
                    if ( scalar(@exons2) == 1 ){
                        $merge = 0 if ( $starts2[$j] < $starts1[$i] || $ends1[$i] < $ends2[$j] );
                    }
                    if ( scalar(@exons1) == 1 ){
                        $merge = 0 if ( $starts1[$i] < $starts2[$j] || $ends2[$j] < $ends1[$i] );
                    }
                }
            }
        }
        if ($merge){
            return 1;
        }
        else{
            return 0;
        }
    }
    else{
        # store all the ends except the last one                                                                                                          
        # store all the starts except the first one                                                                                                       
        my @starts1;
        my @starts2;
        my @ends1;
        my @ends2;
        my @pstarts1;
        my @pstarts2;
        my @pends1;
        my @pends2;
        for(my $k=0; $k<scalar(@exons1); $k++){
            my ($chr, $start, $end, $strand, $trans_id, $gene_id) = @{$exons1[$k]};
            push(@pstarts1, $start) unless $k == 0;
            push(@pends1,   $end) unless $k == scalar(@exons1)-1;
            push(@starts1, $start);
            push(@ends1,   $end);
        }
        for(my $k=0; $k<scalar(@exons2); $k++){
            my ($chr, $start, $end, $strand, $trans_id, $gene_id) = @{$exons2[$k]};
            push(@pstarts2, $start) unless $k == 0;
            push(@pends2,   $end) unless $k == scalar(@exons2)-1;
            push(@starts2, $start);
            push(@ends2,   $end);
        }
	my $pstarts1 = join "_", @pstarts1;
        my $pends1   = join "_", @pends1;
        my $pstarts2 = join "_", @pstarts2;
        my $pends2   = join "_", @pends2;

        # if one is inside the othe, at least all starts and all ends that we stored are completely contained:                                           
        if ( ($pstarts1=~/$pstarts2/ && $pends1 =~/$pends2/)
             ||
             ($pstarts2=~/$pstarts1/ && $pends2 =~/$pends1/)
            ){
            # now check the overlap:                                                                                                                      
            my $merge = 1;
            my $init  = 0;
            for(my $i=0; $i<scalar(@starts1); $i++){
                for(my $j = $init; $j<scalar(@starts2); $j++){
                    # only need to worry about the first and last exon                                                                                    
                    if ( !($starts2[$j]>$ends1[$i]) && !($starts1[$i]>$ends2[$j]) ){
                        if ( $i==0 && $j>0 ){
                            $merge = 0 if ( $starts1[$i] < $starts2[$j] );
                        }
                        if ( $i>0 && $j==0 ){
                            $merge = 0 if ( $starts2[$j] < $starts1[$i] );
                        }
                        if ($i == scalar(@starts1)-1 && $j<scalar(@starts2)-1 ){
                            $merge = 0 if ( $ends1[$i] > $ends2[$j]);
                        }
			if ($i < scalar(@starts1)-1 && $j==scalar(@starts2)-1 ){
			    $merge = 0 if ( $ends2[$j] > $ends1[$i]);
			}
		    }
		}
	    }
	    if ($merge){
		return 1;
	    }
	    else{
		return 0;
	    }
	}
	else{
	    return 0;
	}
    }
}


1;
