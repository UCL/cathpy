#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp;
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;

my $dir = "./";

# Tag the Aln file and convert any lowercase AA to uppercase AA (GS does not recognize lowercase AA)

my $ff1 = "1";
my $ff2 = "3";

my $ff1_file = "$dir/$ff1.aln";
my $ff2_file = "$dir/$ff2.aln";

my @ff1_headers = `grep ">" $ff1_file`;

my @ff2_headers = `grep ">" $ff2_file`;

my $alncat = "$dir/$ff1.$ff2.alncat";
my $aln = "$dir/$ff1.$ff2.aln";

system "cat $ff1_file $ff2_file > $alncat";

system "mafft --anysymbol --amino --quiet --maxiterate 2 $alncat > $aln";

my @lines = read_file("$aln");

my $aln_tagged = "$dir/$ff1.$ff2.tagged_for_GS.aln";
open(OUTFILE, ">$aln_tagged") or die "Can't open file $aln_tagged\n";

foreach my $line (@lines){
	chomp($line);
	if($line=~ /\>(.*)/){
        	
		my $header= $1;
		if ( grep( /^>$header$/, @ff1_headers ) ) {
  			print OUTFILE ">$header|$ff1\n";  # FunFam 1 -> ff1
		}
            	elsif ( grep( /^>$header$/, @ff2_headers ) ) {
  			print OUTFILE ">$header|$ff2\n";  # FunFam 2 -> ff2
		}
		else{
			print OUTFILE "SOMETHING WRONG!\n";
		}
	}
	else{
		$line = uc $line;
		print OUTFILE "$line\n";
    	}
}

# Run GS

system("python2 ../group_sim_sdp.py $aln_tagged $ff1 $ff2 > $aln_tagged.gs");

# Get Jalview annotation file 

system("perl ../groupsim2jalview.pl $aln_tagged.gs > $aln_tagged.gs.jlv");