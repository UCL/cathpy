#!/usr/bin/perl

use strict;

# this script is going to create an alignment features file
# to feed into jalview applet. 

# lets declare some variables
my $CONSfile = shift; # the GroupSim file
my $file = "$CONSfile.processed";
my $bug = "# col_num	score	column";
# now lets get the groupsim score data
my $p=0;
open(CONSFILE, "<$CONSfile")
		or die "Can't open file $CONSfile\n";
open(OUTFILE, ">$file")
	or die "Can't open file $file\n";
while(my $line = <CONSFILE>) {  # reading conservation score file

    if($line=~ /(\d*)\t(\w\w\w\w|\-?\d\.\d*)/){
	my $num = $1 + 1;
	$line =~ s/$1/$num/;
	
	if($p==0){
		$p++;}
	else{
	print OUTFILE "$line";}
}
}
close(CONSFILE);
close(OUTFILE);
# print the header
print STDOUT "JALVIEW_ANNOTATION\n";

# print the graph data
print STDOUT "BAR_GRAPH\tGroupSim\tAlignment SDP Positions based on GroupSim\t";

my $first = 0;
open(INFILE, "<$file")
		or die "Can't open file $file\n";
while(my $line = <INFILE>) {  # reading conservation score file

    if($line=~ /\d*\t(\w\w\w\w|\-?\d\.\d*)/){
	my $score = $1;
   	 if ($first == 0) { # Doesn't prepend on the first field.
       	 $first = 1;
   	 } else {
        	print "|"; # Prepend a pipe as the field separator.
   	 }
	
	if ($score=~ /\-?\d\.\d*/){
    	print STDOUT "$score,$score";
	}
	else {
		print STDOUT "0.000,0.000";
	}
}
}
close(INFILE);
print STDOUT "\n";

