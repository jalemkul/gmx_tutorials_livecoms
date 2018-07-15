#!/usr/bin/perl

use strict;

# sort_mol2_bonds.pl - a script to reorder the listing in a .mol2 @<TRIPOS>BOND
# section so that the following conventions are preserved:
#   1. Atoms on each line are in increasing order (e.g. 1 2 not 2 1)
#   2. The bonds appear in order of ascending atom number
#   3. For bonds involving the same atom in the first position, the bonds appear
#       in order of ascending second atom
#
# Written by: Justin Lemkul (jalemkul@vt.edu)
#
# Distributed under the GPL-3.0 license

unless (scalar(@ARGV)==2)
{
    die "Usage: perl sort_mol2_bonds.pl input.mol2 output.mol2\n";
}

my $input = $ARGV[0];
my $output = $ARGV[1];

open(IN, "<$input") || die "Cannot open $input: $!\n";
my @in = <IN>;
close(IN);

open(OUT, ">$output") || die "Cannot open $output: $!\n";

# get number of atoms and number of bonds from mol2 file
my @tmp = split(" ", $in[2]);
my $natom = $tmp[0];
my $nbond = $tmp[1];

# print out everything up until the bond section
my $i=0;
while (!($in[$i] =~ /BOND/))
{
    print OUT $in[$i];
    $i++;
}

# print the bond section header line to output
print OUT $in[$i];
$i++;

# read in the bonds and sort them
my $bondfmt = "%6d%6d%6d%5s\n";
my @tmparray; 

# sort the bonds - e.g. the one that has the
# lowest atom number in the first position and then the
# lowest atom number in the second position (swap if necessary)
for (my $j=0; $j<$nbond; $j++)
{
    my @tmp = split(" ", $in[$i+$j]);
    # parse atom numbers
    my $ai = $tmp[1];
    my $aj = $tmp[2];
    # reorder if second atom number < first
    if ($aj < $ai)
    {
        $ai = $tmp[2];
        $aj = $tmp[1];
    }
    # store new lines in a temporary array
    $tmparray[$j] = sprintf($bondfmt, $tmp[0], $ai, $aj, $tmp[3]); 
}

# loop over tmparray to find each atom number
my $nbond = 0;
for (my $x=1; $x<=$natom; $x++)
{
    my @bondarray;
    my $ntmp = scalar(@tmparray);
    for (my $b=0; $b<$ntmp; $b++)
    {
        my @tmp = split(" ", $tmparray[$b]);
        if ($tmp[1] == $x)
        {
            push(@bondarray, $tmparray[$b]);
            splice(@tmparray, $b, 1);
            $ntmp--;
            $b--;
        }
    }

    if (scalar(@bondarray) > 0) # some atoms will only appear in $aj, not $ai
    {
        my $nbondarray = scalar(@bondarray);
        if ($nbondarray > 1)
        {
            # loop over all bonds, find the one with lowest $aj
            # and then print it
            for (my $y=0; $y<$nbondarray; $y++)
            {
                my @tmp2 = split(" ", $bondarray[$y]);
                my $tmpatom = $tmp[2];
                my $lowindex = 0;
                if ($tmp2[2] < $tmpatom)
                {
                    $lowindex = $y; 
                }
                my $keep = splice(@bondarray, $lowindex, 1);
                $y--;
                $nbondarray--;
                my @sorted = split(" ", $keep);
                $nbond++;
                printf OUT $bondfmt, $nbond, $sorted[1], $sorted[2], $sorted[3]; 
            }
        }
        else
        {
            $nbond++;
            my @tmp2 = split(" ", $bondarray[0]);
            printf OUT $bondfmt, $nbond, $tmp2[1], $tmp2[2], $tmp2[3];
        }
    }
}

close(OUT);

exit;
