#! /usr/local/bin/perl -w

use strict;
use warnings;

my $lineage_data = $ARGV[0];
my $vcffile = $ARGV[1];

open(LINEAGES, $lineage_data) || die "\nCannot open file: $!\n";

my $lineages;

while ( <LINEAGES> ) {
	my $line = $_;
	chomp $line;

	next if $. == 1;
	
	my ($lineage, $pos, $gene, $mutation, $therest) = split("\t", $line);
	my ($from, $to) = split("/", $mutation);
	

	my $match = $to;
	
	if ( $lineage =~ m/\*\*/ ) {
		$match = $from;
	}	

#	print $lineage."\t".$pos."\t".$mutation."\t".$from."\t".$to."\t".$match."\n";
	$lineage =~ s/\*\*//;
	
	$lineages->{$pos} = {
		lineage => $lineage,
		from    => $from,
		to      => $to,
		match   => $match
	}
}

if ( $vcffile =~ m/\.gz$/ ) {
	open(VCF, "gunzip -c $vcffile | ") || die "\nCannot open gzipped $vcffile: $!\n";
	} else {
	open(VCF, $vcffile) || die "\nCannot open $vcffile: $!\n";
}

my $flag = 0;
my @matches;

while (<VCF>) {
	my $line = $_;
	chomp $line;

	next if $line eq '';

	if ( $flag == 0 ) {
		if ( $line =~ m/CHROM\tPOS\tID/ ) {
			$flag++;
			next;
		} else {
			next;
		}
	}

	my ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info, $format1, $format2) = split("\t", $line);

	next unless $lineages->{$pos};
	
	if ( ($alt eq '.' && $ref eq $lineages->{$pos}->{match}) || $alt eq $lineages->{$pos}->{match} ) {
		push(@matches, $lineages->{$pos}->{lineage});
	}
}

print $vcffile."\t".join("\t", sort { length($b) <=> length($a) } @matches)."\n";

