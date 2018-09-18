#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp;
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;

my $dir = "/cath/homes2/ucbtdas/Working/eclipse_workspace/GroupSim";#/3.40.30.10_FunFams

# Tag the Aln file and convert any lowercase AA to uppercase AA (GS does not recognize lowercase AA)

my ($ff1, $ff2) = @ARGV;
chomp($ff1);
chomp($ff2);
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

system("python /cath/homes2/ucbtdas/Working/eclipse_workspace/GroupSim/group_sim_sdp.py $aln_tagged $ff1 $ff2 > $aln_tagged.gs");

# Get Jalview annotation file 

system("perl /cath/homes2/ucbtdas/Working/eclipse_workspace/GroupSim/groupsim2jalview.pl $aln_tagged.gs > $aln_tagged.gs.jlv");

=head
# Get two separate files for each class

my $classA = "$dir/groupsim_analysis/ClassA_forlogo.aln";
open(OUTFILE2, ">$classA") or die "Can't open file $classA\n";

my $classC = "$dir/groupsim_analysis/ClassC_forlogo.aln";
open(OUTFILE3, ">$classC") or die "Can't open file $classC\n";

# create new Bio::SeqIO object    
        my $in = Bio::SeqIO->new( -file   => "<$aln_tagged",
                                  -format => "fasta");
        # loop through each instance in FASTA file
        while(my $seq = $in->next_seq()){
        my $id = $seq->id();
        my $sequence = $seq->seq;
	$id =~ /\w+\/\d+\-\d+\|(\d+)/;
	my $grp = $1;
	#print "$id\t$grp\n";
        if($grp ==1){
            	print OUTFILE2 ">$id\n$sequence\n";    
        }
	elsif($grp ==3){
		print OUTFILE3 ">$id\n$sequence\n";    
        }
    }

mkdir("$dir/groupsim_analysis/SDPs");

# Get SDP sequence logo
	my @sdps; my @cons;
	my @gslines = read_file("$aln_tagged.gs");
	foreach my $line (@gslines){
		chomp($line);
		unless($line =~ /^\#/){
			if($line =~ /(\d+)\t(\d+\.\d+)/){
				my $col =$1;
				$col++;
				my $gs_score = $2;
				#print "$col, $gs_score\n";
				if($gs_score >= 0.7) {
					push(@sdps, $col);
				}
				if($gs_score <= 0.3) {
					push(@cons, $col);
				}
			}
		}
	}
	my $weblogodir = "/cath/homes2/ucbtdas/Documents/PhD/Softwares/weblogo-3.3";
	
	# For Class A
	foreach my $sdp (@sdps){
		system("$weblogodir/weblogo --format png_print -c chemistry < $classA --sequence-type protein --first-index 1 --lower $sdp --upper $sdp --show-xaxis NO --show-yaxis NO --fineprint $sdp > $dir/groupsim_analysis/SDPs/A.sdp.$sdp.png");
	}
	foreach my $sdp (@cons){
		system("$weblogodir/weblogo --format png_print -c chemistry < $classA --sequence-type protein --first-index 1 --lower $sdp --upper $sdp --show-xaxis NO --show-yaxis NO --fineprint $sdp > $dir/groupsim_analysis/SDPs/A.cons.$sdp.png");
	}

	# For Class C
	foreach my $sdp (@sdps){
		system("$weblogodir/weblogo --format png_print -c chemistry < $classC --sequence-type protein --first-index 1 --lower $sdp --upper $sdp --show-xaxis NO --show-yaxis NO --fineprint $sdp > $dir/groupsim_analysis/SDPs/C.sdp.$sdp.png");
	}

	foreach my $sdp (@cons){
		system("$weblogodir/weblogo --format png_print -c chemistry < $classC --sequence-type protein --first-index 1 --lower $sdp --upper $sdp --show-xaxis NO --show-yaxis NO --fineprint $sdp > $dir/groupsim_analysis/SDPs/C.cons.$sdp.png");
	}
=cut
