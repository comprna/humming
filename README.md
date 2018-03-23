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

The columns correspond to: `chromosome, source, feature, start, end, score, strand, frame, and transcript ID`.

#### PAF

The PAF format is not yet supported - `under development`


----------------------------
## Citation
----------------------------

Humming code is based on the code from 

Eyras E, Caccamo M, Curwen V, Clamp M. ESTGenes: alternative splicing from
ESTs in Ensembl. Genome Res. 2004 May;14(5):976-87. PMID: 15123595;


----------------------------
## License
----------------------------

Humming code is released under the MIT license. 





