#!/usr/bin/perl


## combine information over reps and sort
$in = shift(@ARGV); ## 0 rep file, from fork script
open(OUT, "> mav_$in") or die "failed to write mav $in\n";
$nreps=10;
foreach $rep (0..9){
	$in =~ s/ch\d+/ch$rep/ or die "failed sub for $in to $rep\n";
	print "working on $in\n";
	open(IN, $in) or die "failed to open the infile $in\n";
	<IN>; ## burn header
	while(<IN>){
		chomp;
		@line = split(/\s+/,$_);
		if($rep==0){
			push(@snps,$line[1]);
		}
		$eff{$line[1]} += ($line[4] + $line[5] * $line[6]) * 0.1;## * .1 to divide by 10
	}
	close(IN);
}

foreach $snp (@snps){
	print OUT "$eff{$snp}\n";
}
close(OUT);
