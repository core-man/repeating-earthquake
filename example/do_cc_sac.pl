#!/usr/bin/perl
use strict;
use warnings;
# calculate cross-correlation coefficient and time shift between sac files

my $cc_sac   = "../bin/cc_sac";	# cross-correlation
my $filelist = "KN_AAK";        # sac file list

my $phase = "P";    # phase name used to name output file only
my $tmark = "t0";   # reference sac HEADER
my $delta = 0.025;	# sampling delta
my $tb	  = 5;		# time before tmark
my $te	  = 15;		# time after tmark


#########################################################
# Cross-correlation of each sac file pair in the filelist
# Notes:
#  1. Only do the calculation for event pair recorded at the same network and station.
#  2. The maximum separation between two events is set to 60 km in the code.
#  3. The depth unit is assumed to be meter in the code.
#  4. The filter parameters are set in the code
#
# Some arguments:
#  t : time window with tb second before and eb second after tmark
#  b : filter (1 - bandpass; 0 - no filter)
#
# Output:
#  KN_AAK.P.ccor : contain file1 file2 ccf dt

`$cc_sac -t $tb,$te -d $delta -p $phase -h $tmark -b 0 $filelist`;


