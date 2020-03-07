#! /usr/bin/perl

use strict;

my %position=();
my %names=();
open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    #my @b=split(/\#/, $a[8]);

    if (defined $position{$a[9]}) {
	my @c=split(/\:/, $position{$a[9]});
	if (($c[0] eq $a[0])&&($a[1] < $c[1])) {
	    $c[1] = $a[1];
	}
	if (($c[0] eq $a[0])&&($a[2] > $c[2])) {
	    $c[2] = $a[2];
	}
	$position{$a[9]}="$c[0]\:$c[1]\:$c[2]";

	my $transposon=$a[3];
	if ($names{$a[9]} !~ /$transposon/) {$names{$a[9]}=$names{$a[9]}.",$transposon";}
    }
    else {
	$position{$a[9]}="$a[0]\:$a[1]\:$a[2]";
	$names{$a[9]}=$a[3];
    }
}
close input;

open (output, ">>temp_for_sort") or die "Can't open temp_for_sort since $!\n";
while ((my $key, my $value) = each (%position)) {
    my @z=split(/\:/, $value);
    print output "$z[0]\t$z[1]\t$z[2]\t$names{$key}\n";
}
close output;

system("sort +0 -1 +1n -2 +2n -3 temp_for_sort > sorted");
system("uniq -c sorted > $ARGV[1]");
system("rm sorted temp_for_sort");
