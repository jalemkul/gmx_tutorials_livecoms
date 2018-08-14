#!/usr/bin/perl

use strict;

# Written by Justin Lemkul, Ph.D.
# Contact: jalemkul@vt.edu
# 6/7/2017
#
# License: GPL-3.0 or later, at user's option

# This script deletes out waters between certain z-coordinate values
# Z-values are specified by a chosen atom on the command line
# The decision on how to split "top" of bilayer from "bottom" - define "middle" atom,
# is taken from the command line.  If z(reference) > z(middle) -> top, and vice versa
# Determine average z-coordinate value for top and bottom, if OW of water falls within
# that range, delete it and its corresponding HW's

# Collect stuff from command line arguments into a hash
unless (@ARGV)
{
    die "Example usage: perl water_deletor.pl -in bilayer.gro -out bilayer_fix.gro -ref O33 -middle C50 -water 3 -v [yes/no]\n";
}

print "\nStarting water deletion process...\n";

my %args = @ARGV;
my $input;
my $output;
my $ref_atom;
my $middle_atom;
my $nwater;
my $verbose;

if (exists($args{"-in"}))
{
    $input = $args{"-in"};
}
else
{
    die "No input specified!\n";
}

# collect all the input lines into an array
open(IN, $input);
my @in = <IN>;
close(IN);

# grab some stuff from the input (non-atom lines)
my $header = shift(@in);
my $orig_natoms = shift(@in);
chomp($orig_natoms);
my $box = pop(@in);

if (exists($args{"-out"}))
{
    $output = $args{"-out"};
}
else
{
    print "Using default output filename of \"bilayer_fix.gro\"\n";
    $output = "bilayer_fix.gro";
}

if (exists($args{"-ref"}))
{
    $ref_atom = $args{"-ref"};
}
else
{
    die "No reference atom specified!\n";
}

if (exists($args{"-middle"}))
{
    $middle_atom = $args{"-middle"};
}
else
{
    die "No middle atom specified!\n";
}

# define number of water atoms for splicing
if (exists($args{"-water"}))
{
    $nwater = $args{"-water"};
}
else
{
    print "Option -nwater not set, assuming 3 atoms in the water model.\n";
    $nwater = 3;
}

# define verbosity
if (exists($args{"-v"}) && (($args{"-v"} eq "yes") || ($args{"-v"} eq "Yes") || ($args{"-v"} eq "y") || ($args{"-v"} eq "Y")))
{
    $verbose = 1;
}
else
{
    $verbose = 0;   # default to non-verbose mode
}

# process input file to strip out all lines corresponding to reference or middle atoms
my @ref_atoms;
my @middle_atoms;

# also count how many water molecules we start with
my $nwater_start = 0;

foreach $_ (@in)
{
    my @data = split(" ", $_);
    my $atom = $data[1];

    if ($atom =~ /$ref_atom/)
    {
        push(@ref_atoms, $_);
    }
    elsif ($atom =~ /$middle_atom/)
    {
        push(@middle_atoms, $_);
    }
    elsif ($atom =~ /OW/)
    {
        $nwater_start++;
    }
}

# determine the middle of the bilayer from the middle atoms
# middle of bilayer is defined as the average z-coordinate of the middle atoms
my $total_z;
my $nmiddle = scalar(@middle_atoms);

foreach $_ (@middle_atoms)
{
    my @data = split(" ", $_);
    my $z = 0;
    # last field on line if no velocities, otherwise n-3
    if (scalar(@data) > 6)
    {
        # there are velocities
        $z = $data[scalar(@data)-4];
    }
    else
    {
        $z = $data[scalar(@data)-1];
    }
    # debug
    #print "$z\n";

    $total_z += $z;
}

my $middle_z = $total_z / $nmiddle;

print "Defining the middle of the bilayer as z = $middle_z.\n";

# now split the reference atoms into top and bottom
my @top_ref;
my @bottom_ref;

foreach $_ (@ref_atoms)
{
    my @data = split(" ", $_);
    my $z = 0;
    # last field on line if no velocities, otherwise n-3
    if (scalar(@data) > 6)
    {
        # there are velocities
        $z = $data[scalar(@data)-4];
    }
    else
    {
        $z = $data[scalar(@data)-1];
    }

    if ($z > $middle_z)
    {
        push(@top_ref, $_);
    }
    elsif ($z < $middle_z)
    {
        push(@bottom_ref, $_);
    }
    else
    {
        print "Weird z-coordinate found!\n";
        print "$_\n";
    }
}

# determine average z-coordinate values for the top and bottom reference atoms
my $top_z;
my $ntop = scalar(@top_ref);

# debug
print "ntop = $ntop\n";

foreach $_ (@top_ref)
{
    my @data = split(" ", $_);
    my $z = 0;
    # last field on line if no velocities, otherwise n-3
    if (scalar(@data) > 6)
    {
        # there are velocities
        $z = $data[scalar(@data)-4];
    }
    else
    {
        $z = $data[scalar(@data)-1];
    }
    # debug
    print "$z\n";

    $top_z += $z;
}

my $top_z_avg = $top_z / $ntop;

my $bottom_z;
my $nbottom = scalar(@bottom_ref);

foreach $_ (@bottom_ref)
{
    my @data = split(" ", $_);
    my $z = 0;
    # last field on line if no velocities, otherwise n-3
    if (scalar(@data) > 6)
    {
        # there are velocities
        $z = $data[scalar(@data)-4];
    }
    else
    {
        $z = $data[scalar(@data)-1];
    }

    $bottom_z += $z;
}

my $bottom_z_avg = $bottom_z / $nbottom;

print "Boundaries for eliminating water are $bottom_z_avg < z < $top_z_avg.\n";

# now actually delete the waters
# go thru the input file and find any OW that has a z-coordinate within the bounds
# defined above.
my @clean_atoms;

my $size = scalar(@in);

my $ndel = 0;   # keep track of how many waters we have deleted

for (my $i=0; $i<$size; $i++)
{
    my @line = split(" ", $in[$i]);
    # account for the fact that at high atom numbers, the string is joined, i.e. OW10000
    my $atom = substr($line[1], 0, 2);
    # get last entry - in the above case, OW10000, there is no distinct atom number
    my $z = $line[(scalar(@line)-1)];

    unless ($atom eq "OW")
    {
        push(@clean_atoms, $in[$i]);
    }
    else
    {
        if (($z > $bottom_z_avg) && ($z < $top_z_avg))
        {
            # assemble the residue numbering for printing/bookkeeping
            if ($verbose == 1)
            {
                my @res = split('', $line[0]);
                pop(@res);  # get rid of "SOL"
                pop(@res);
                pop(@res);
                my $resnr = join('', @res);
                print "Deleting residue $resnr\n";
            }

            # debug
            # my @delete = splice(@in, $i, 3);
            # print "@delete";

            splice(@in, $i, $nwater);
            $i--;
            $size -= $nwater;
            $ndel++;
        }
        else
        {
            push(@clean_atoms, $in[$i]);
        }
    }
}

# Report what happened
print "$ndel water molecules have been deleted.\n";
my $watnew = $nwater_start - $ndel;
print "$watnew water molecules remain. Update your topology!\n";

# print the cleaned output
my $natoms_clean = scalar(@clean_atoms);

open(OUT, ">$output");
print OUT $header;
print OUT "$natoms_clean\n";

foreach $_ (@clean_atoms)
{
    print OUT "$_";
}

print OUT $box;

close(OUT);

print "\nOutput file $output written.  You may want to renumber the file with,\n";
print "for example, genconf in the GROMACS package.\n\n";

exit;
