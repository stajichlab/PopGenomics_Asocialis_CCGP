#!/usr/bin/env perl
use strict;
use warnings;

my $prefix = "JNA_AS_CCGP";
my %names;
while(<>) {
    chomp;
    my ($id,$file) = split(',',$_);
    $id =~ s/CCGPMC\d+[a-z]_//;
    if ( $id =~ /JNA_AS_CCGP(\d+[AB]?)(_\d+)?/ ) {
	my $NUM = $1;
	push @{$names{$NUM}}, $file;
    } else {
	print join(",",($id,$file)),"\n";
    }
}
for my $n ( sort { $a <=> $b } keys %names ) {
    print join(",", (sprintf("%s_%s",$prefix,$n), join(";",@{$names{$n}}))),"\n";
}
