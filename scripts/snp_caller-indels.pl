#! /usr/local/bin/perl -w

use strict;
use warnings;

use lib 'scripts';

use POSIX;
use Getopt::Long;
use Pod::Usage;

use Bio::DB::SeqFeature::Store;
use Bio::Tools::CodonTable;

###############################################################################
# Process command line options
#

my ( $qual_cutoff, $dp4_cutoff, $dp_cutoff, $dpmax_cutoff, $af_cutoff, $mq_cutoff, @vcffiles, $chrom, 
		$gfffile, $cdsonly, $site, $gt, $verbose, $help );

GetOptions(
           'qual|q=i'     => \$qual_cutoff,
           'dp4=i'        => \$dp4_cutoff,
           'dp=i'         => \$dp_cutoff,
           'dpmax=i'      => \$dpmax_cutoff,
           'af=f'         => \$af_cutoff,
           'mq=f'         => \$mq_cutoff,
           'site=f'       => \$site,
           'gt=f'         => \$gt,
           'vcf=s'        => \@vcffiles,
           'gff=s'        => \$gfffile,
		   'cdsonly'      => \$cdsonly,
           'chrom=s'      => \$chrom,
           'verbose=i'    => \$verbose,
           'help|?'       => \$help,
          ) or pod2usage(2);

# Defaults
 $qual_cutoff  ||= 30;
 $dp_cutoff    ||=  0;
 $dpmax_cutoff ||=  100000;
 $dp4_cutoff   ||=  0;
 $af_cutoff    ||=  0;
 $mq_cutoff    ||=  0;
 $cdsonly      ||=  0;
 
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

if ( $help || !@vcffiles || !$chrom  ) { pod2usage(1) }

my $gffstore = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -dsn => $gfffile );

###############################################################################
## Process VCF files
###############################################################################

my $indels;
my $indel_positions;

my $filtered_sites;

foreach my $vcffile ( @vcffiles ) {
	print  "Reading: ".$vcffile."\n" if $verbose;
	
	if ( $vcffile =~ m/\.gz$/ ) {
		open(VCF, "gunzip -c $vcffile | ") || die "\nCannot open gzipped $vcffile: $!\n";
	} else {
		open(VCF, $vcffile) || die "\nCannot open $vcffile: $!\n";
	}

	my $flag = 0;
	my $indel_count = 0;
	
	while (<VCF>) {
		my $line = $_;
		chomp $line;
	
		next if $line eq '';
	
		if ( $flag == 0 ) {
			if ( $line =~ m/CHROM\tPOS\tID/ ) {
#				print $line."\n";
				$flag++;
				next;
			} else {
#				print $line."\n";
				next;
			}
		}
	
		my @filtered;
		my ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info, $format1, $format2) = split("\t", $line);
		
		# Only look at sites from this chromosome
		next if $seqname ne $chrom;

		my $type = $info =~ m/INDEL/ ? 'INDEL' : 'SNP';

		if ( $type eq 'INDEL' ) {
			if ( $info =~ m/IS=(\d+),([\.\d]+);/ ) {
				my $isdp = $1;
				my $is = $2;
			
				push(@filtered, 'IS') if $is < 0.5;
			}

		}

		push(@filtered, 'HET') if $alt =~ m/,/;
		
		# QUAL filtering
		push(@filtered, 'QUAL') if $qual < $qual_cutoff;

		# MQ filtering
		if ( $info =~ m/MQ=([\.\d]+);/ ) {
			my $mq = $1;
			
			push(@filtered, 'MQ') if $mq < $mq_cutoff; 
		} else {
			push(@filtered, 'MQ')
		}

		# Site depth of coverage filtering
		if ( $info =~ m/DP=(\d+);/ ) {
			my $dp = $1;
			
			push(@filtered, 'DP') if $dp < $dp_cutoff; 
			push(@filtered, 'DPMAX') if $dp > $dpmax_cutoff; 
		} else {
			push(@filtered, 'DP')
		}
		
		if ( $af_cutoff ) {
			if ( $info =~ m/AF1=([\.\d]+);/ ) {
				my $af = $1;

				push(@filtered, 'AF') if ( $alt eq '.' && $af > 0 ); 
				push(@filtered, 'AF') if ( $alt ne '.' && $af < $af_cutoff ); 
			}
		}

		# Strand bias filtering
		if ( $info =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/ ) {
			my $reff = $1;
			my $refr = $2;
			my $altf = $3;
			my $altr = $4;
		
			my $total = $reff + $refr + $altf + $altr;

			my $strand_min_depth = $dp_cutoff / 2;
		
			if ( $alt eq '.' ) {
			
				# strand bias for all bases
				push(@filtered, 'SBREF') if ( $reff < $strand_min_depth || $refr < $strand_min_depth );

				push(@filtered, 'DP4') if ( ( $reff + $refr ) / $total < ( $dp4_cutoff / 100 ) );

			} else {
			
				# strand bias for variants
				push(@filtered, 'SBALT') if ( $altf < $strand_min_depth || $altr < $strand_min_depth );				

				push(@filtered, 'DP4') if ( ( $altf + $altr ) / $total < ( $dp4_cutoff / 100 ) ); 
			}
			
		} else {
			push(@filtered, 'SB');
		}
		
		if ( $site ) {
			if ( $site == $pos ) {
				print $vcffile."\t".$line."\t[".join(", ", @filtered)."]\n";
			}
		}
		
		push(@{$filtered_sites->{$pos}} , { vcf => $vcffile, line => $line , filtered => \@filtered }) if @filtered;

		if ( $type eq 'INDEL' && $alt ne '.' ) {
			$indel_positions->{$vcffile}->[$pos] = 1;
			push(@{$indels->{$pos}} , { ref => $ref, alt => $alt, vcf => $vcffile, line => $line , filtered => \@filtered });
		}
	}
}

my $count = 0;

foreach my $pos ( keys %$indels ) {	
	my $flag = 0;
	my $filtered = 0;

	$filtered++ if $filtered_sites->{$pos};

	for (my $i = $pos - 5; $i < $pos + 5; $i++) {
		foreach my $pos2 ( @{$indels->{$i}} ) {
			$filtered++ if @{$pos2->{filtered}};
			$flag++;
		}
	}

	if ( $flag < 2 && $filtered == 0 ) {
		my $info = geneinfo($indels->{$pos}->[0]->{vcf}, $pos, $indels->{$pos}->[0]->{ref}, $indels->{$pos}->[0]->{alt});
		$count++;
	}
}

print join(", ", @vcffiles)."\t".$count." indel differences found\n" if $count == 0;

print "\n";

##############################################################################################################################
#
sub geneinfo {
	my ($vcf, $pos, $ref, $alt) = @_;

	# Look for genes affected by SNPs
	my $start = $pos - 2;
	my $end   = $pos + 2;

#next unless ( $pos == 1673425 || $pos == 2155168 || $pos == 7582 || $pos == 7585 );
#next unless ( $pos == 3884906 );

 	my @features = $gffstore->features( -seq_id => $chrom, -start => $pos, -end => $pos, -type => 'gene');
	
	foreach my $feature ( @features ) {
		
		my $type = $feature->type;
		$type =~ s/:GenBank//g;
		
		my ($name, $snp, $seqobj);
		
		$name = $feature->display_name;
		$name =~ s/CDS:GenBank//g;
		$name =~ s/(\(|\))//g;

		$seqobj = $feature->seq;

		if ( $feature->strand == 1 ) {
			$snp = $pos - $feature->start + 1;    #  +1 as we need to keep first base of start
		} else {
			$snp = $feature->end - ( $pos - 1 );
			
			# complement single base
			$ref = $ref eq 'G' ? 'C' : $ref eq 'C' ? 'G' : $ref eq 'A' ? 'T' : $ref eq 'T' ? 'A' : '';
			$alt = $alt eq 'G' ? 'C' : $alt eq 'C' ? 'G' : $alt eq 'A' ? 'T' : $alt eq 'T' ? 'A' : '';
		}

		my $aaloc = ceil($snp / 3);
		my $codonloc = ($aaloc * 3) - 2;  # -2 as we want first base of a 3 base codon		

		my $protobj = $seqobj->translate;
		my $aaref = $protobj->subseq($aaloc, $aaloc);
			
		my $codon1 = $seqobj->subseq($codonloc, $codonloc+2);
		my $codon2 = $codon1;
		substr($codon2, $snp - $codonloc, 1, $alt );

		my $myCodonTable = Bio::Tools::CodonTable->new();
		my $aaalt = $myCodonTable->translate($codon2);

		my $snptype = $aaref eq $aaalt ? 'S' : 'NS';

		my $resist = 'DELETE';
#		print $name."\t".$snp."\t".$pos."\t".$type."\t".$feature->strand."\t".$ref."-".$alt."\t".$aaref.'/'.$aaalt."\t".$aaloc."\t".$codon1.' / '.$codon2."\t".$snptype."\t".$resist."\n";
		print $vcf."\t".$name."\t".$snp."\t".$pos."\t".$feature->strand."\t".$ref."-".$alt."\t".$aaref.'/'.$aaalt."\t".$aaloc."\t".$codon1.' / '.$codon2."\t".$snptype."\n";

	}
	
	
	# Check if the site was intergenic
	if ( $cdsonly && @features == 0 ) {

 		@features = $gffstore->features( -seq_id => $chrom, -start => ($pos-1000), -end => ($pos+1000), -type => 'gene');

#		print OUT "\nFeatures found: ".join(", ", @features)."\n";
		
		my $coords;

		foreach my $feature ( @features ) {
			my $name = $feature->display_name;
#				print OUT $name."\n";
			
			my $start = $feature->strand > 0 ? $feature->start : $feature->end; 
			my $end   = $feature->strand > 0 ? $feature->end : $feature->start;
				
			$coords->{$chrom}->{$start} = $feature;
		}

		my ($downstream, $upstream);
		foreach my $start ( sort { $a <=> $b } keys %{$coords->{$chrom}} ) {		
			if ( $pos > $start ) {
				$downstream = $coords->{$chrom}->{$start};
			} elsif ( $pos < $start ) {
				$upstream = $coords->{$chrom}->{$start};
			}	
			
			last if $upstream;
		}

		my ( $down_name, $down_dist, $up_name, $up_dist );
		if ( $downstream ) {
			$down_name = $downstream->display_name;
			$down_dist = $pos - $downstream->end;
			
#			print STDERR "SS ".$downstream->display_name."\t".$downstream->strand."\t".$downstream->start."\t".$downstream->end."\n";
			
		} else {
			$down_name = 'NONE';
			$down_dist = $pos;
		}

		if ( $upstream ) {
			$up_name = $upstream->display_name;
			$up_dist = $upstream->start - $pos;
		} else {
			$up_name = 'NONE';
			$up_dist = $pos;
		}
		
#		print OUT "Down: ".$down_name." [".$down_dist."]\tUp: ".$up_name." [-".$up_dist."]\n";
		
		print "INTERG:\t\t".$pos."\tDown: ".$down_name." [".$down_dist."]\tUp: ".$up_name." [-".$up_dist."]\n";
		
# 		my ( $up_match, $down_match ) = ('', '');
# 		$up_match = ($upstream ? ' ('.find_gene_match($upstream->seq).')' : '') if $genomefile;
# 		$down_match = ($downstream ? ' ('.find_gene_match($downstream->seq).')' : '') if $genomefile;
# 		
# 		if ( $filtered_flag ) {
# 			$filtered .= $chrom."\t".$pos."\t".$type."\t".$ref." -> ".$alt."\tDown: ".$down_dist." : ".$down_name.$down_match."; UP: ".$up_dist." : ".$up_name.$up_match."\n\n";		
# 		} else {
# 			$unfiltered .= $chrom."\t".$pos."\t".$type."\t".$ref." -> ".$alt."\tDown: ".$down_dist." : ".$down_name.$down_match."; UP: ".$up_dist." : ".$up_name.$up_match."\n\n";		
# 		}
		
	}
}

###############################################################################

__END__

=head1 NAME

filter_vcf.pl - Filters VCF file based on various criteria

=head1 SYNOPSIS

filter_vcf.pl --vcf vcffile --chrom NC_009777

filter_vcf.pl --vcf vcffile1 --vcf vcffile2 --chrom NC_009777 --qual 30 --dp4 75 --dp 10 --dpmax 100 -af 0.75 --showfiltered --noindels --noheader --phylip 
	              --ignorefile Mobiles.txt --verbose 1

=head1 OPTIONS

  --vcf             Input VCF file and label. Format --vcf label,filename (no spaces between , and text)
  --chrom           Reference sequence name to analyse (In case VCF is mapped against multiple sequences)
  --qual            Minimum QUAL score [30]
  --dp4             % reads supoprting SNP [0]
  --dp              Minimum read depth at site [0]
  --dpmax           Maximum read depth at site [100000]
  --af              Allele frequency cutoff e.g. 0.75 [0]
  --showfiltered    Also show all the sites removed (in a Filtered subsection at end),
                       useful to see how well filtering is performing
  --noindels        Filter out INDEL variants
  --phylip          Export data in phylip alignement format (may have to adjust sequence names)
                    (forces --noindels)
  --ignorefile      Specify a text file containing coordinates of regions to ignore e.g. mobile elements
  --noheader        Do not print VCF header
  
  --verbose         print more information (not much at the moment!) 
                    1 - logging
                    2 - debugging
  
=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=head1 AUTHOR

Adam Witney <awitney@sgul.ac.uk>
BuG@S group, Deptartment of Cellular and Molecular Medicine,
St George's, University of London,
London, UK

=head1 COPYRIGHT

bugasbase_pars.pl is Copyright (c) 2009 Adam Witney. UK. All rights reserved.

You may distribute under the terms of either the GNU General Public License or the Artistic License, as specified in the Perl README file.

=cut

=head1 DISCLAIMER

This software comes with no warranty and you use it at your own risk. There may be bugs!

=cut
