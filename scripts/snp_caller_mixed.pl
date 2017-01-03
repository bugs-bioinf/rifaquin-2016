#! /usr/local/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

###############################################################################
# Process command line options
#

my ( $vcffile, $verbose, $help );

GetOptions(
           'vcf=s'        => \$vcffile,
           'verbose'    => \$verbose,
           'help|?'       => \$help,
          ) or pod2usage(2);

# Defaults
 
pod2usage(-exitstatus => 0 ) if $help;

if ( $help || !$vcffile ) { pod2usage(1) }


###############################################################################

if ( $vcffile =~ m/\.gz$/ ) {
	open(VCF, "gunzip -c $vcffile | ") || die "\nCannot open gzipped $vcffile: $!\n";
} else {
	open(VCF, $vcffile) || die "\nCannot open $vcffile: $!\n";
}

# define major and minor output VCF files
my $vcf_maj = $vcffile.'-maj.vcf';
my $vcf_min = $vcffile.'-min.vcf';

open (VCFMAJ, ">".$vcf_maj) || die "\nCannot open $vcf_maj: $!\n";
open (VCFMIN, ">".$vcf_min) || die "\nCannot open $vcf_min: $!\n";

my $flag = 0;
my $counters;

while (<VCF>) {
	my $line = $_;
	chomp $line;

	# skip the VCF header, only printing to maj/min VCF's
	if ( $flag == 0 ) {
		print VCFMAJ $line."\n";
		print VCFMIN $line."\n";

		if ( $line =~ m/CHROM\tPOS\tID/ ) {
			$flag++;
			next;
		} else {
			next;
		}
	}
	
	my ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info, $format1, $format2) = split("\t", $line);
	
	if ( $info =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/ ) {
		my $reff = $1;
		my $refr = $2;
		my $altf = $3;
		my $altr = $4;
			
		my ($line_maj, $line_min);

		# create variables to manipulate DP4 field in output VCF's
		my $info_ref = $info;
		my $info_alt = $info;

		my $dp4_ref = 'DP4='.$reff.','.$refr.',0,0';
		my $dp4_alt = 'DP4=0,0,'.$altf.','.$altr;
		$info_ref =~ s/DP4=(\d+),(\d+),(\d+),(\d+)/$dp4_ref/;
		$info_alt =~ s/DP4=(\d+),(\d+),(\d+),(\d+)/$dp4_alt/;

		# original VCF calls a REF base, but there are significant reads suggesting ALT
		if ( $alt eq '.' && ($altf + $altr) > 4 && ($altf + $altr) / ($reff + $refr + $altf + $altr ) > 0.1 ) {

#			$line_maj = join("\t", ($seqname, $pos, $id, $ref, '.', $qual, $filter, $info_ref, $format1, $format2, 'HET-ALT'));
#			$line_min = join("\t", ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info_alt, $format1, $format2, 'HET-ALT'));

			# we can't change anything here, as we don't have access to the ALT base (could check pileup?)
			$line_maj = $line."\tHET-ALT";
			$line_min = $line."\tHET-ALT";
	
			$counters->{'ALT'}++;

		# original VCF calls an ALT base, but there are significant reads suggesting REF
		} elsif ( $alt ne '.' && ($reff + $refr) > 4 && ($reff + $refr) / ($reff + $refr + $altf + $altr ) > 0.1 ) {

			if ( ($reff + $refr) > ($altf + $altr) ) {
				$line_maj = join("\t", ($seqname, $pos, $id, $ref, '.', $qual, $filter, $info_ref, $format1, $format2, 'HET-REF-1'));
				$line_min = join("\t", ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info_alt, $format1, $format2, 'HET-REF-1'));
			} else {
				$line_maj = join("\t", ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info_alt, $format1, $format2, 'HET-REF-2'));
				$line_min = join("\t", ($seqname, $pos, $id, $ref, '.', $qual, $filter, $info_ref, $format1, $format2, 'HET-REF-2'));
			}

			$counters->{'REF'}++;

		# oterwise just print original VCF record
		} else {
			$line_maj = $line;
			$line_min = $line;

			$counters->{'NONE'}++;			
		}
		
		# output for logging
		print $line."\n" if $verbose;
		print $line_maj."\n" if $verbose;
		print $line_min."\n\n" if $verbose;
		
		# print the major/minor output VCF files	
		print VCFMAJ $line_maj."\n";
		print VCFMIN $line_min."\n";
	}
}

if ( $verbose ) {
	print "\n\n";
	print "Counts ALT: ".$counters->{'ALT'}."\n";
	print "Counts REF: ".$counters->{'REF'}."\n";
	print "Counts NONE: ".$counters->{'NONE'}."\n";
}

###############################################################################

__END__

=head1 NAME

snp_caller_mixed.pl - Reads mixed infection VCF file and splits into major and minor isolate vcf files

=head1 SYNOPSIS

snp_caller_mixed.pl --vcf vcffile

snp_caller_mixed.pl --vcf vcffile --verbose 1

=head1 OPTIONS

  --vcf             Input VCF file
  
  --verbose         print more information (not much at the moment!) 
  
=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given VCF file and split the sequence into
major and minor isolate VCF's.

=head1 AUTHOR

Adam Witney <awitney@sgul.ac.uk>
BuG@S group, Deptartment of Cellular and Molecular Medicine,
St George's, University of London,
London, UK

=head1 COPYRIGHT

snp_caller_mixed.pl is Copyright (c) 2016 Adam Witney. UK. All rights reserved.

You may distribute under the terms of either the GNU General Public License or the Artistic License, as specified in the Perl README file.

=cut

=head1 DISCLAIMER

This software comes with no warranty and you use it at your own risk. There may be bugs!

=cut
