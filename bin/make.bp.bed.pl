#! /usr/bin/perl

use strict;

my @sample=();
open (in, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
my $line=<in>;
close in;

my %chrs=();
my @a=split(/\t/, $line);
for my $i (0..$#a) {
    if ($a[$i] =~ /_class$/) {
	my $name=$a[$i];
	$name =~ s/_class//;
	my $j=$i+1;
	my $k=$i+2;
	my $l=$i+3;
	system("cut -f7,4,6,$j,$k,$l $ARGV[0] > temp");
	open (input, "<temp") or die "Can't open temp since $!\n";
	open (output, ">$name.insertion.bp.bed") or die "Can't open $name.insertion.bp.bed since $!\n";
	my $header=<input>;
	while (my $line=<input>) {
	    chomp($line);
	    my @b=split(/\t/, $line);
	    if (($b[4] ne "0")||($b[5] ne "0")) {
		my @c=split(/\:/, $b[2]);
		my @d=split(/\./, $c[1]);
		if ($c[0] eq "P") {next;}
		if ($d[0] > $d[1]) {
		    my $temp=$d[0];
		    $d[0]=$d[1];
		    $d[1]=$temp;
		}
		my $lower=$d[0];
		my $upper=$d[1];
		if (($lower >= 0) && ($upper >= 0)) {
		    print output "$c[0]\t$lower\t$upper\t$b[0]\t$b[1]\t$b[3]\t$b[4]\t$b[5]\n";
		}
		$chrs{$c[0]}=1;
	    }
	}
	close input;
	close output;
	system("rm temp");

	if ($ARGV[1] ne "") {
	    open (input, "<$name.insertion.bp.bed") or die "Can't open $name.insertion.bp.bed since $!\n";
	    open (output, ">tmp") or die "Can't tmp since $!\n";
	    while (my $line=<input>) {
		chomp($line);
		my @a=split(/\t/, $line);
		if (($a[0] =~ /^\d{1,2}$/) || ($a[0] eq "X") || ($a[0] eq "Y")) {$a[0]="chr$a[0]";}
		my $strand="+";
		if ($a[4] eq "antisense") {$strand="-";}
		print output "$a[0]\t$a[1]\t$a[2]\t$a[3]\t\.\t$strand\t$a[5]\t$a[6]\t$a[7]\n";
	    }
	    close input;
	    close output;

	    system("bedtools intersect -a tmp -b $ARGV[1] -f 0.1 -wo -s > tmp1");
	    if ($ARGV[2] eq "") {
		system("awk -F \"\\t\" '{OFS=\"\\t\"; if ((\$4==\$13)&&(\$6==\$15)) print \$1,\$2,\$3,\$4,\$5,\$6}' tmp1 > tmp2");
	    }
	    else {
		my %family=();
		open (input, "<$ARGV[2]") or die "Can't open $ARGV[2] since $!\n";
		while (my $line=<input>) {
		    chomp($line);
		    my @a=split(/\t/, $line);
		    $family{$a[0]}=$a[1];
		}
		close input;

		open (input, "<tmp1") or die "Can't open tmp1 since $!\n";
		open (output, ">>tmp2") or die "Can't open tmp2 since $!\n";
		while (my $line=<input>) {
		    chomp($line);
		    my @a=split(/\t/, $line);
		    if (($family{$a[3]} eq $family{$a[12]}) && ($a[5] eq $a[14])) {
			print output "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";
		    }
		}
		close input;
		close output;
	    }
	    
	    if (-s "tmp2") {
		my %to_filter=();
		open (input, "<tmp2") or die "Can't open tmp2 since $!\n";
		while (my $line=<input>) {
		    chomp($line);
		    my @a=split(/\t/, $line);
		    $to_filter{"$a[0]\:$a[1]\:$a[2]\:$a[3]\:$a[5]"}=1;
		}
		close input;
		open (input, "<tmp") or die "Can't open tmp since $!\n";
		open (output, ">$name.insertion.bp.bed") or die "Can't open $name.insertion.bp.bed since $!\n";
		while (my $line=<input>) {
		    chomp($line);
		    my @a=split(/\t/, $line);
		    if (!defined $to_filter{"$a[0]\:$a[1]\:$a[2]\:$a[3]\:$a[5]"}) {
			my $direction="sense";
			if ($a[5] eq "-") {$direction="antisense";}
			my $chr_num=$a[0];
			$chr_num =~ s/chr//;
			if (($chrs{$a[0]} == 1) && (! defined $chrs{$chr_num})) {$chr_num=$a[0];}
			print output "$chr_num\t$a[1]\t$a[2]\t$a[3]\t$direction\t$a[6]\t$a[7]\t$a[8]\n";
		    }
		}
		close input;
		close output;
	    }

	    system("rm tmp*");

	}

    }
}
	
