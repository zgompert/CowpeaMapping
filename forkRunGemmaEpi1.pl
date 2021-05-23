#!/usr/bin/perl
#
# core gemma fits for parameter estimation, bslmm and lmm 
#

use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);

$g = shift(@ARGV); ## genotype file
$p = shift(@ARGV); ## phenotype file

foreach $ph (1..4){ ## four traits, starting at 2
	$phn = $ph+1;
	system "sleep 2\n";
	foreach $ch (0..9){
		$pm->start and next;
		$o = "o_cowpea1_epiPc_ph$ph"."_ch$ch";
	#	system "gemma -g $g -p $p -gk 1 -o $o\n";
    		system "gemma -g $g -p $p -k output_stan/o_cowpea1_fit_ph$ph"."_ch0.cXX.txt -bslmm 1 -n $phn -o $o -notsnp -w 200000 -s 1000000\n";
    		#system "gemma -g $g -p $p -k output/$o".".cXX.txt -bslmm 1 -n $phn -o $o -notsnp -w 200000 -s 1000000\n";
		$pm->finish;
	}
}
$pm->wait_all_children;

