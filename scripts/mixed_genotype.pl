#! /usr/local/bin/perl -w

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;

use Statistics::R;

###############################################################################
# Process command line options
#

my ( $vcffile, $qual_cutoff, $dp4_cutoff, $dp_cutoff, $dpmax_cutoff, $af_cutoff, $mq_cutoff, $mixed_cutoff, $noindels, $only_total, $verbose, $help );

GetOptions(
	'qual|q=i'     => \$qual_cutoff,
	'dp4=i'        => \$dp4_cutoff,
	'dp=i'         => \$dp_cutoff,
	'dpmax=i'      => \$dpmax_cutoff,
	'af=f'         => \$af_cutoff,
	'mq=f'         => \$mq_cutoff,
	'mc=f'         => \$mixed_cutoff,
	'vcf=s'        => \$vcffile,
	'noindels'     => \$noindels,
	'total'        => \$only_total,
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
 $mixed_cutoff ||=  0.05;
 $verbose      ||=  0;
 $only_total   ||=  0;
  
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

if ( $help || !$vcffile ) { pod2usage(1) }

###############################################################################

if ( $vcffile =~ m/\.gz$/ ) {
	open(VCF, "gunzip -c $vcffile | ") || die "\nCannot open gzipped $vcffile: $!\n";
} else {
	open(VCF, $vcffile) || die "\nCannot open $vcffile: $!\n";
}
print STDERR "Reading ".$vcffile."\n";

my $flag = 0;
my $data;
my $fails = 0;

my $total_percent = 0;
my $het_counter = 0;
my $af1_counter = 0;
my $fq_counter = 0;
my $uncertain_total = 0;
my $het_percents;

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
	
	my @filtered;
	my $percent = 0;
	my $het_percent1 = 0;
	my $het_percent2 = 0;
	my ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info, $format1, $format2) = split("\t", $line);
	
	my $type = $info =~ m/INDEL/ ? 'INDEL' : 'SNP';
	
	push(@filtered, 'INDEL') if $noindels && $type eq 'INDEL'; # have to skip these for behaviour 2
	next if $noindels && $type eq 'INDEL';
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
			
		push(@filtered, 'DP') if $dp <= $dp_cutoff; 
		push(@filtered, 'DPMAX') if $dp > $dpmax_cutoff; 
	} else {
		push(@filtered, 'DP')
	}

#	if ( $info =~ m/AF1=([\.\d]+);/ ) {
#		my $af = $1;
		
#		push(@filtered, 'AF') if ( $alt eq '.' && $af > 0 ); 
#		push(@filtered, 'AF') if ( $alt ne '.' && $af < $af_cutoff ); 
#	} else {
#		push(@filtered, 'AF');
#	}

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
	#		if ( @filtered == 0 && ( $altf + $altr ) >= $dp_cutoff && (( $altf + $altr ) / $total >= $mixed_cutoff) && $altf >= $strand_min_depth && $altr >= $strand_min_depth ) {
			if ( @filtered == 0 && ( $altf + $altr ) >= $dp_cutoff && $altf >= $strand_min_depth && $altr >= $strand_min_depth ) {
				$het_counter++;
				$het_percent1 = sprintf("%d", ((( $altf + $altr ) / $total) * 100));
				$het_percent2 = sprintf("%d", ((( $reff + $refr ) / $total) * 100));
	#			$het_percent1 = sprintf("%d", ((( $altf + $altr ) / $total) * 100));
#				print $line."\n";

				print "".( $reff + $refr )."\t".( $altf + $altr )."\tref\n" if $verbose;
			}
			
			push(@filtered, 'SBREF') if ( $reff < $strand_min_depth || $refr < $strand_min_depth );
			push(@filtered, 'DP4') if ( ( $reff + $refr ) / $total < ( $dp4_cutoff / 100 ) );
		} else {
			
			# strand bias for variants
#			if ( @filtered == 0 && ( $reff + $refr ) >= $dp_cutoff && (( $reff + $refr ) / $total >= $mixed_cutoff) && $reff >= $strand_min_depth && $refr >= $strand_min_depth ) {
			if ( @filtered == 0 && ( $reff + $refr ) >= $dp_cutoff && $reff >= $strand_min_depth && $refr >= $strand_min_depth ) {
				$het_counter++;
				$het_percent1 = sprintf("%d", ((( $altf + $altr ) / $total) * 100));
				$het_percent2 = sprintf("%d", ((( $reff + $refr ) / $total) * 100));
#	print STDERR "Mixed: ".$pos." : ".$het_percent1." : ".$het_percent2."\n" if $het_percent2 > 70 && $het_percent2 < 85;
#print $pos."\t".$het_percent2."\n";
#				print $line."\n";
				print "".( $reff + $refr )."\t".( $altf + $altr )."\talt\n" if $verbose;
			}
			
			push(@filtered, 'SBALT') if ( $altf < $strand_min_depth || $altr < $strand_min_depth );				
			push(@filtered, 'DP4') if ( ( $altf + $altr ) / $total < ( $dp4_cutoff / 100 ) ); 
		}

#$het_percent1 = sprintf("%d", ((( $altf + $altr ) / $total) * 100));
#$het_percent2 = sprintf("%d", ((( $reff + $refr ) / $total) * 100));
#push(@{$het_percents->{'ref'}}, $het_percent1);
#push(@{$het_percents->{'alt'}}, $het_percent2);
		
		push(@{$het_percents->{'alt'}}, $het_percent1) if $het_percent1 > 0;
		push(@{$het_percents->{'ref'}}, $het_percent2) if $het_percent2 > 0;
		$percent = sprintf("%d", ((( $altf + $altr ) / $total ) * 100 ));

	} else {
		push(@filtered, 'SB');
	}

	# AF1
	if ( @filtered == 0 && $info =~ m/AF1=([\.\d]+);/ ) {
		my $af = $1;
		
		$af1_counter++ if $af > 0 && $af < 1;
	}

	# FQ
	if ( @filtered == 0 && $info =~ m/FQ=([-]*[\.\d]+);/ ) {
		my $fq = $1;
		
		$fq_counter++ if $fq > 0;
	}
	
	if ( @filtered ) {
		$uncertain_total++;
		next;
	}
	
	die if $percent < 0 || $percent > 100;
	$data->[$percent]++;
	$total_percent++;
	
#	last if $pos > 1000;
}

#print STDERR $vcffile."\t".$total_percent."\t".$uncertain_total."\t".$het_counter."\t".$af1_counter."\t".$fq_counter."\n"; #.mean(\@het_percents)."\t".stdev(\@het_percents)."\n";
#print $vcffile."\t".$total_percent."\t".$uncertain_total."\t".$het_counter."\t".$af1_counter."\t".$fq_counter."\n"; #.mean(\@het_percents)."\t".stdev(\@het_percents)."\n";

#for (my $i = 0; $i < @{$het_percents->{'alt'}}; $i++ ) {
#	print $i."\t".$het_percents->{'alt'}->[$i]."\n";
#
#}




if ( $only_total ) {

#	my $R = Statistics::R->new( r_bin => '/homedirs8/share/Tools/x86_64/R-2.15/bin/R' ) ;
	my $R = Statistics::R->new( r_bin => '/usr/bin/R' ) ;
	my $png = $vcffile.".mixed.png";
	my $datafile1 = $vcffile.".hetpercent.ref.txt";
	my $datafile2 = $vcffile.".hetpercent.alt.txt";

#	print STDERR "Starting R\n";
	$R->startR;
	$R->send(q`library(ggplot2)`);
	$R->send(q`library(reshape)`);
	$R->send(q`library(plyr)`);
	
#	print STDERR "Plotting histogram\n";

	my @cmds;
	push(@cmds, "data1 <- data.frame(ref=c(".join(",", @{$het_percents->{'ref'}})."))");
	push(@cmds, "data2 <- data.frame(alt=c(".join(",", @{$het_percents->{'alt'}})."))");
	push(@cmds, "write.table(data1, file='".$datafile1."')") if $verbose;
	push(@cmds, "write.table(data2, file='".$datafile2."')") if $verbose;
	push(@cmds, "data1 <- arrange(data1, ref)");
	push(@cmds, "data2 <- arrange(data2, alt)");
	push(@cmds, "data1 <- melt(data1)");
	push(@cmds, "data2 <- melt(data2)");
#	push(@cmds, "data <- rbind(data1, data2)");
	push(@cmds, "data <- rbind(data2)");
#	push(@cmds, "ggplot(data, aes(x=value, fill=variable)) + geom_histogram(binwidth=1) + xlab('percent') + scale_y_continuous(limits=c(0,100)) + scale_x_continuous(limits=c(0,100), breaks=c(seq(0,100,by=5))) + opts(axis.text.x=theme_text(angle=90, hjust=1, size=15))");
#	push(@cmds, "ggplot(data, aes(x=value, colour=variable)) + geom_freqpoly(binwidth=2) + xlab('percent') + scale_y_continuous(limits=c(0,100)) + scale_x_continuous(limits=c(0,100), breaks=c(seq(0,100,by=5))) + opts(axis.text.x=theme_text(angle=90, hjust=1, size=15))");

	push(@cmds, "ggplot(data, aes(x=value, colour=variable)) + geom_freqpoly(binwidth=2) + xlab('percent') + scale_y_continuous(limits=c(0,100)) + scale_x_continuous(limits=c(0,100), breaks=c(seq(0,100,by=5))) + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))");

#	$R->send('jpeg("'.$jpeg.'", width=1000, height=480)') ;
#	$R->send('pdf("'.$jpeg.'.pdf", width=11, height=8)') ;
	
	foreach ( @cmds ) {
		$R->send($_);
	}
	
	$R->send('ggsave("'.$png.'")') ;

	$R->stopR() ;
	
	exit(0);
}

my @r_data;

foreach ( 0..100 ) {
	print STDERR $_."\t".(defined $data->[$_] ? $data->[$_] : 0 )."\n" if $verbose;
#	print STDERR "".(defined $data->[$_] ? $data->[$_] : 0 ).",";
	
	push(@r_data, (defined $data->[$_] ? $data->[$_] : 0 ));
}

print STDERR "\n\nFails: ".$fails."\n";
#exit(0);
#my $R = Statistics::R->new( r_bin => '/homedirs8/share/Tools/x86_64/R-2.15/bin/R' ) ;
my $R = Statistics::R->new( r_bin => '/usr/bin/R' ) ;
#my $R = Statistics::R->new( r_bin => '/var/local/bin/R' ) ;
my $jpeg = $vcffile.".jpg";
my $plot = $vcffile.".pdf";
my $png = $vcffile.".png";

print STDERR "Starting R\n";
$R->startR;
print STDERR "Loading ggplot2\n";
$R->send(q`library(ggplot2)`);

my $cmd1 = "counts <- matrix(c(".join(",", @r_data)."), nrow=101, ncol=1)";
my $cmd2 = "percents <- matrix(c(".join(",", 0..100)."), nrow=101, ncol=1)";
my $cmd3 = "data <- cbind(percents, counts)";
my $cmd4 = "colnames(data) <- c('percent','counts')";
my $cmd5 = "data <- as.data.frame(data)";
my $cmd6 = "ggplot(data, aes(x=percent, y=counts)) + geom_bar(position='dodge', stat='identity') + scale_y_log10()";

#jpeg("TEST.jpg")

#counts <- matrix(c(1,2,3,4,5,6,7,8,9,10), nrow=10, ncol=1)
#percents <- matrix(c(1,2,3,4,5,6,7,8,9,10), nrow=10, ncol=1)
#data <- cbind(percents, counts)
#colnames(data) <- c('percent','counts')
#data <- as.data.frame(data)
#ggplot(data, aes(x=percent, y=counts)) + geom_bar(position='dodge', stat="identity") + scale_y_log10()

$R->send('jpeg("'.$jpeg.'")') ;
$R->send('pdf("'.$plot.'", width=11, height=8)');

$R->send($cmd1);
$R->send($cmd2);
$R->send($cmd3);
$R->send($cmd4);
$R->send($cmd5);
$R->send($cmd6);

#$R->send(qq`x = 123 \n print(x)`) ;
#my $ret = $R->read ;

$R->stopR() ;

###############################################################################
## Subroutines
###############################################################################

sub mean {
	my($data) = @_;
	if (not @$data) {
		die("Empty array\n");
	}
	
	my $total = 0;
	foreach (@$data) {
		$total += $_;
	}
	
	my $mean = $total / @$data;
	return $mean;
}

sub stdev {
	my($data) = @_;
	if(@$data == 1){
		return 0;
	}
	my $mean = mean($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($mean - $_) ** 2;
	}

	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}

###############################################################################

__END__

=head1 NAME

process_vcf_resistance.pl - Filters VCF file based on various criteria and shows gene locations

=head1 SYNOPSIS

process_vcf_resistance.pl --vcf vcffile --chrom NC_009777

process_vcf_resistance.pl --vcf vcffile --chrom NC_009777 --qual 30 --dp4 75 --dp 10 --dpmax 100 -af 0.75 --showfiltered --noindels --verbose 1

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
  --ignorefile      Specify a text file containing coordinates of regions to ignore e.g. mobile elements
  
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
