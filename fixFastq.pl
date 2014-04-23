#!/usr/bin/perl -w
#
use strict;

my ($in, $out, $a, $b, $c, $d, $l, $r, $n);

if (!-d "output") {
    `mkdir output`;
}

foreach $in (`ls *.fastq`) {
    chomp $in;
    $out = "output/new_" . $in;
    open (OUT, ">$out") or die "cannot open $out: $!";
    
    open (IN, "$in") or die "cannot open $in: $!";
    while ($a = <IN>) {
        $b = <IN>;
        $c = <IN>;
        $d = <IN>;
        
        chomp $a;
        ($l, $r) = split( /\s+/, $a );
        $n = substr( $r, 0, 1 );
        print OUT $l, '#0/', "$n \($a\)\n";
        print OUT $b, $c, $d;
    }

    close IN;
    close OUT;
}
