#!/usr/bin/perl
#
# core gemma fits for parameter estimation, bslmm and lmm 
#

use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);

$g = shift(@ARGV); ## genotype file
$p = shift(@ARGV); ## phenotype file

foreach $ph (1..2){ ## four traits, starting at 2
	$phn = $ph+1;
	foreach $ch (0..9){
		$pm->start and next;
		$o = "o_cowpea1_larv_fit_ph$ph"."_ch$ch";
		if($ch == 0){ ## lmm only on chain 0
			system "gemma -g $g -p $p -gk 1 -o $o\n";
			system "gemma -g $g -p $p -lmm 4 -n $phn -k output/$o".".cXX.txt -o $o\n";
		}
    		system "gemma -g $g -p $p -bslmm 1 -n $phn -o $o -maf 0 -w 200000 -s 1000000\n";
		$pm->finish;
	}
}
$pm->wait_all_children;

