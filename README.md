----------------------------
# Humming
----------------------------

Humming is a a set of tools for the analysis of long-read mappings to a genome reference. 

----------------------------
## Set up
----------------------------

Humming is developed in Perl.

Copy the code in a directory /pathto/humming/

and then do:
```
export PERL5LIB=/pathto/humming/
```

Once downloaded, the scripts can be used directly from the command line by specifying the absolute path to the programs. These use functions
defined in the libraries/ directory. 

----------------------------------------
## Comparison of mappings and annotations
----------------------------------------

This program compares a set of mapped long-reads or predicted transcripts (in GFF/GTF/PAF format) with a set of
annotations (GFF/GTF/PAF format). It then calculates the number and proportions of various elements in each set
that appear in the other set: splice-sites, exon-exon junctions, transcripts (it matches all the exon-exon junctions), 
and genes. For genes the comparison is performed by first making genes (see below) and then testing the genomic overlap. 

Usage: 

```
perl compare_mappings_and_annotations.pl -a <annotations> -p <predictions> -f <format-annotations> -u <format-predictions> < -j | -t | -g >
```

- **-a** | **--annotation**:	File with the transcript annotations
- **-p** | **--prediction**:	File with the transcript predictions
- **-f** | **--format_ann**: 	Format of the annotation file: GTF, GFF
- **-u** | **--format_pred**:  	Format of the prediction file: GTF, GFF (it can be different from the annotation file)
- **-j** | **--junctions**:	Compare splice sites and exon-exon junctions
- **-t** | **--transcripts**:	Compare transcripts. It matches transcripts within the same loci
- **-g** | **--genes**:		Compare gene loci (genome extension in the same strand that produce one or more transcripts)
- **-h** | **--help**:		Print this help

Example:

```
perl compare_mappings_and_annotations.pl -a mapped_reads.gtf -p gencode.gtf -f GTF -u GTF -j -t -g
```

----------------------------
## Make genes
----------------------------

From a set of transcripts in GFF format (or GTF format) you can build genes. A gene is defined
as a set of transcripts that overlap in genomic extention in the same strand and share
at least one exon or one splice-site (option to be selected). This definition allows to separate 
genes within genes, i.e. a gene in the intron of another gene in the same strand, or nested genes.
On the other hand, some genes that are labeled as different by Gencode are classified as the same gene with this definition, as
it only uses the exon-intron structure of transcripts without considering any other information (e.g. the encoded protein). 

Usage: 

```
perl build_genes.pl -i <input> -f <format> -c <constraint>
```

- **-i** | **--input**:      File name with transcripts
- **-f** | **--format**:     Format of the transcript file: GTF, GFF
- **-o** | **--output**:     Output file name (GTF format)
- **-c** | **--constraint**: Label:exon or splice-site. Genes are defined as transcripts overlapping
                             in the same strand and sharing an exon or a splice-site. Default: exon
- **-r** | **--reuse**:      Optional. It will reuse gene IDs if available. Default: not set
- **-h** | **--help**:       Print this help

Example: 

```
perl build_genes.pl -i transcripts.gff -f GTF -c exon -o output_test -r
```

The method to build genes performs a fast clustering of transcripts (linear in the number of transcripts)
by genomic extent (and strand) first. In a second step, it analyzes each transcript cluster, links 
transcripts according to the defined constraint, and recovers the sets of transcripts that form a gene 
(It recovers the connected components of the graph built using depth first search). Gene IDs are either numbered clusters,
or made up from pre-existing gene IDs in the input file (if the option `reuse` is set).


----------------------------
## Format conversion
----------------------------

Script to perform conversion between formats generally used to describe transcript annotations and/or mappings to a genome reference.
Format supported so far are: GTF, GFF, BED12, PAF.

Usage:

```
perl convert.pl -i <input> -f <format_in> -o <output> -u <format_out>
```

- **-i** | **--input**:      File name with transcripts
- **-f** | **--format_in**:  Format of the transcript file: GTF, GFF, BED12, PAF
- **-o** | **--output**:     Output file name
- **-u** | **--format_out**: Format of the transcript file: GTF, GFF, BED12, PAF
- **-h** | **--help**:       Print this help

Example:

```
perl convert.pl -i RefSeq.bed -f BED12 -o RefSeq.gff -u GFF
```

Important: RefSeq annotations in BED12 format usually contains multiple lines for the same transcript ID, indicating
paralogs. In the consersion to GFF or GTF, these will be separated by adding to the transcript id a subindex "_i", so that
they are not confused in downstream analyses. 

----------------------------
## Formats
----------------------------

#### GTF

The GTF format required only needs to contain the exon lines, e.g.:

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";

```
The columns correspond to: `chromosome, source, feature, start, end, score, strand, frame`, and the composite `column 9`; where
column 9 requires the format as shown in the example.


#### GFF

The GFF format, as with GTF, only requires exon lines, and the 9th column defines the transcript (group ID):

```
chr14 Ensembl exon  73741918  73744001  0.0 - . ENST00000000001.1
chr14 Ensembl exon  73749067  73749213  0.0 - . ENST00000000001.1
chr14 Ensembl exon  73750789  73751082  0.0 - . ENST00000000001.1
chr14 Ensembl exon  73753818  73754022  0.0 - . ENST00000000001.1

```


#### BED12


The BED12 is an extension of the BED format to describe a transcript (one or more exons) per line:
```
chr21   37085436        37105240        ENSG00000230794.1:ENST00000412240.1     0       -       37085436        37105240        0       2       309,40  0,19764
chr21   35445848        35515334        ENSG00000243927.1:ENST00000477091.1     0       +       35445848        35515334        0       5       197,74,23,140,627       0,29285,29947,51792,68859
chr21   38309513        38362524        ENSG00000159267.10:ENST00000427746.1    0       -       38309513        38362524        0       5       179,163,122,135,67      0,1618,10952,24822,52944
```

The columns correspond to
`chromosome, chr_start (0 based), chr_end, name, score, strand, thick_start, thick_end, ItemGBGroup ID /Parent ID, block_count, block_sizes, block_starts (relative to chr_start)`


#### PAF

The Pairwise Alignment Format (PAF) describes the mapping of a query sequence to a target sequence (https://github.com/lh3/miniasm/blob/master/PAF.md).

This is an example of PAF output from Minimap2: 

```
d63d83e8-630a-4ebb-bc98-b9a058b81a77    1190    47      1151    -       NC_000017.11    83257441        43525694        43526760        939     1163    60      tp:A:P  cm:i:48 s1:i:334        s2:i:0  NM:i:224        ms:i:477        AS:i:477        nn:i:0  cg:Z:18M1I14M5I36M2I11M1D14M2I6M1I4M2I2M4I4M2I21M2I7M2I4M3D11M1I8M3I7M1I9M1D4M2D2M1D10M3D1M2D10M1I12M7I34M2I16M1I5M1D16M1D5M1D5M1I5M1I2M1I14M1I7M2I11M1D6M1I24M2D17M1I4M1I32M1D3M2I10M2I20M1D3M4I8M1D27M2I3M1D13M1D11M1D31M2D53M4I4M2I4M1I13M1D25M8D4M1D19M6I8M1D4M1D18M2I4M1I18M1D6M1D21M1I1M1I12M1I15M2I15M4I12M2D10M1D20M1I14M1D25M1I7M1D15M2D12M5D6M3D8M2I8M2I12M2I11M1D7M1D4M2I3M1I4M1I8M1D7M2I23M
```

The columns correspond to

1  Query sequence name
2  Query sequence length                                                                                                                                
3  Query start (0-based)                                                                                                                                
4  Query end (0-based)                                                                   
5  Relative strand: "+" or "-"                                                                                                               
6  Target sequence name                                                                                                                                 
7  Target sequence length                                                                                                                               
8  Target start on original strand (0-based)                                                                                                            
9  Target end on original strand (0-based)                                                                                                              
10 Number of sequence matches                                                                                                                           
11 Alignment block length (total number of sequence matches, mismatches, insertions and deletions in the alignment)                                     
12 Mapping quality (0-255; 255 for missing)                                                                                                             
13 SAM-like output (key-value pairs) cg:Z: corresponds to the cigar string (see https://samtools.github.io/hts-specs/SAMv1.pdf)                 


----------------------------
## References
----------------------------

Humming code is based on the algorithms used in this publication

Eyras E, Caccamo M, Curwen V, Clamp M. ESTGenes: alternative splicing from
ESTs in Ensembl. Genome Res. 2004 May;14(5):976-87. PMID: 15123595;


----------------------------
## License
----------------------------

Humming code is released under the MIT license. 





