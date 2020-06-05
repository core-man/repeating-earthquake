#!/usr/bin/perl
use strict;
use warnings;
# calculate cc and time shift between sac files

my ($f1, $f2) = ("data/19931201005901500/II.AAK.00.BHZ.SAC",
                 "data/20030906154700205/KN.AAK..BHZ.SAC");
my $phase = "P";    # phase name
my $tmark = "t0";   # reference sac HEADER
my $delta = 0.025;	# data delta
my $tb	  = 5;		# time before tmark
my $te	  = 15;		# time after tmark


my ($ccf, $dt) = split " ", &CalCCF($f1, $f2, $phase, $tmark, $tb, $te, $delta);
print STDERR "$f1 $f2: $ccf (cc) $dt (dt)\n";


# calculate cross-correlation coeffients
# input: file1 file2 phase tmark (t_begin t_end delta)
# output: ccf dt
sub CalCCF {
	my ($file1, $file2, $phase, $tmark) = @_;
	#print STDERR "$file1, $file2\n$phase, $tmark\n";

	my $cc	  = "../ccor_sac";	# cross-correlation

	# default parameters
	my $delta = 0.025;	  # data delta
	my $tb	  = 5;		  # time before tmark
	my $te	  = 15;		  # time after tmark

	if (@_ > 5) {
		$tb = $_[4];
		$te = $_[5];
	}
	$delta = $_[6] if (@_ > 6);
	#print STDERR "$tb	$te	$delta\n";

    # temp filelist with sac files
	my $junk_cc = "junk_cc";
	open(JUNK, "> $junk_cc") || die "Error in opening $junk_cc.\n";
	print JUNK "$file1\n$file2\n";
	close(JUNK);

    # do cross-correlation
    # t : time window with tb second before and te second after tmark
    # p : phase name
    # b : filter flag (1 : bandpass; 0 : no filter)
	`$cc -t $tb,$te -d $delta -p $phase -h $tmark -b 0 $junk_cc`;

    unlink $junk_cc;

    # find cross-correlation coefient and dt
	my $output = "${junk_cc}.${phase}.ccor";
	if (-z $output) {
        unlink $output;
		return "-1000 -1000";
	} else {
		open(JUNK, "< $output") || die "Error in opening $output.\n";
		my @ccfs = split " ", <JUNK>; close(JUNK); chomp @ccfs;
		unlink $output;
		#print STDERR "@ccfs\n";
		return "$ccfs[2] $ccfs[3]";
	}

}

