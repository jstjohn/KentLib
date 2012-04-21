#!/usr/bin/perl
use strict;

my $file = $ARGV[0];
open(IN,$file) || die "Incorrect file $file. Exiting...\n";

my ($seq, $name)=('','');
while(<IN>){
    chomp;
    my $line = $_;
    $seq.= uc($line) if(eof(IN));
    if (/\>(\S+)/ || eof(IN)){
	if($seq ne ''){
	    my @seqgaps = split(/[N]{100,}/, $seq);
	    if($#seqgaps > 0){
		my $ctgcount=0;
		foreach my $ctgseq (@seqgaps){
		    $ctgcount++;
		    print "$name"."_$ctgcount\n$ctgseq\n";
		    #print "$name contig$ctgcount (size=".length($ctgseq).")\n$ctgseq\n";
		}
	    }else{
		print "$name"."_1\n$seq\n";
	    }
	}
	$seq='';
	my @nparts = split(' ', $_);
	$name = $nparts[0];
    }else{
	$seq.=uc($line);
    }
}
