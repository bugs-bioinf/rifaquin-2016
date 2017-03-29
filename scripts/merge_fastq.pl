#!/usr/bin/perl

use strict;
use warnings;

my $fastq_1 = $ARGV[0];
my $fastq_2 = $ARGV[1];
my $fastq_out = $ARGV[2];

if ( $fastq_1 =~ m/\.gz$/ ) {
	open(FASTQ1, "gunzip -c $fastq_1 | ") || die "\nCannot open gzipped $fastq_1: $!\n";
	open(FASTQ2, "gunzip -c $fastq_2 | ") || die "\nCannot open gzipped $fastq_2: $!\n";
} else {
	open(FASTQ1, $fastq_1) || die "\nCannot open $fastq_1: $!\n";
	open(FASTQ2, $fastq_2) || die "\nCannot open $fastq_2: $!\n";
}

open(FASTQ, "> $fastq_out") || die "\nCannot open gzipped $fastq_out: $!\n";

while(<FASTQ1>) {
	print FASTQ $_;
	$_ = <FASTQ1>;
	print FASTQ $_;
	$_ = <FASTQ1>;
	print FASTQ $_;
	$_ = <FASTQ1>;
	print FASTQ $_;

	$_ = <FASTQ2>;
	print FASTQ $_;
	$_ = <FASTQ2>;
	print FASTQ $_;
	$_ = <FASTQ2>;
	print FASTQ $_;
	$_ = <FASTQ2>;
	print FASTQ $_;
}
