#!/usr/bin/perl -w



# Rounding procedure of the mighty PerlMonks...
# Halleluja!

sub round {
            my( $num, $prec )= @_;
            return  int( $num/$prec + 0.5 - ($num<0) ) * $prec;
        }





# Procedure to calculate the distance between 2 Vectors


sub abstand {
            my( $a1,$a2,$a3,$b1,$b2,$b3)= @_;
            return  ( (($a1-$b1)**2 + ($a2-$b2)**2 + ($a3-$b3)**2)**0.5);
        }






#Check input

if (@ARGV<7) {die "\n\n
###########################################################################
                               INFLATEGRO 
                  Written by Christian Kandt, (c) 2005-2007

   Kandt C, Ash WL, Tieleman DP (2007): Setting up and running molecular
       dyanmics simulations of membrane proteins. Methods 41:475-488

###########################################################################

INFLATEGRO reads the coordinates of a bilayer and inflates them 
in XY directions using a common scaling factor. To identify
the lipids their actual residue name must be given.

Water coordinates (SOL) will be ingored, everything else will
be centered in the XY plane of the new simulation box.

A distance cutoff in A can be defined: Only lipids with a
P - CA distance exceeding that cutoff will be written. It is currently
assumed that you're actually dealing with phospholipids. However, this can 
be easily extended. Just have a look at the code.

Area per lipid is estimated by caculating the area per protein first.
This is done using a grid-based approach. A grid size of 5 A was found to 
give good results. Output is written as a 3-collumned ASCII file holding
3 area per lipid values: total, upper leaflet & lower leaflet.   

If the additional flag >protein< is set only the protein
coordinates will be written.


USAGE: 
-----

INFLATEGRO  bilayer.gro  scaling_factor  lipid_residue_name  cutoff   inflated_bilayer.gro   gridsize areaperlipid.dat (protein)



...good luck! 

;) Christian Kandt 05/2007


\n\n";}




if (!(open (INPUT, $ARGV[0]))) {print "Eeeeek! No $ARGV[0] at all!\n\n";
                                die;}






# Read lipid coordinates


#$input=$ARGV[0];
$scale=$ARGV[1];
$name=$ARGV[2];
$cutoff=$ARGV[3]*0.1;
$output=$ARGV[4];
$gridsize=$ARGV[5];
$area=$ARGV[6];

$switch=0;

if (@ARGV==8 and $ARGV[7] eq "protein") {
              print "Well, just the protein then....\n";
	      $switch=1;
	      }





$zaehler=1;
$counter=1;
$now=0;
$protein_xmax=-1000;
$protein_ymax=-1000;
$protein_xmin=1000;
$protein_ymin=1000;


print STDOUT "Reading..... \n";
while (<INPUT>) {


      if (/^ +([0-9.-]+) +([0-9.-]+) +([0-9.-]+)/) {
          $box_x=$1;
	  $box_y=$2;
	  $box_z=$3;
      }
	  



unless ($now==0) {

     if (/^ +(\d+)(\S+) +(\S+) +([0-9.-]+) +([0-9.-]+) +([0-9.-]+)/ ) {
      
      
         if ($2 eq $name) {
 	    $resnum_l[$zaehler]=$1;
 	    $resname_l[$zaehler]=$2;
            $atmnamenum=$3;
            $x_l[$zaehler]=$4;
 	    $y_l[$zaehler]=$5;
 	    $z_l[$zaehler]=$6;
 	    
 	    $howlong=length($atmnamenum) - 1;
 	    $atmname_l[$zaehler]=substr($atmnamenum,0,$howlong-4);
 	    $atmnum_l[$zaehler]=substr($atmnamenum,$howlong-4);
 	    
 	    #print "$resnum_l[$zaehler]:    $atmnamenum   $atmname_l[$zaehler] $atmnum_l[$zaehler]\n";   
 	    $zaehler++;
         }
 
 
         if (!($2 eq $name) && !($2 eq "SOL")) {
 	    $resnum_p[$counter]=$1;
 	    $resname_p[$counter]=$2;
            $atmnamenum=$3;
            $x_p[$counter]=$4;
 	    $y_p[$counter]=$5;
 	    $z_p[$counter]=$6;
 	    
 	    $howlong=length($atmnamenum) - 1;
 	    $atmname_p[$counter]=substr($atmnamenum,0,$howlong-4);
 	    $atmnum_p[$counter]=substr($atmnamenum,$howlong-4);
 	    
	    if ($x_p[$counter]>$protein_xmax) {$protein_xmax=$x_p[$counter];}
	    if ($x_p[$counter]<$protein_xmin) {$protein_xmin=$x_p[$counter];}
	    if ($y_p[$counter]>$protein_ymax) {$protein_ymax=$y_p[$counter];}
	    if ($y_p[$counter]<$protein_ymin) {$protein_ymin=$y_p[$counter];}
	    
	    
 	    #print "$resnum_p[$counter]:   $atmnamenum   $atmname_p[$counter]   $atmnum_p[$counter]\n";
 	    $counter++;
         }
 
 
      }


}

      if (/^ +(\d+)(\S+) +(\S+) +(\d+) +([0-9.-]+) +([0-9.-]+) +([0-9.-]+)/ ) {
      
      
         if ($2 eq $name) {
	    $resnum_l[$zaehler]=$1;
	    $resname_l[$zaehler]=$2;
            $atmname_l[$zaehler]=$3;
	    $atmnum_l[$zaehler]=$4;
            $x_l[$zaehler]=$5;
	    $y_l[$zaehler]=$6;
	    $z_l[$zaehler]=$7;
	    
	    if ($atmnum_l[$zaehler]==9999) {$now=1;}
	    #print "$resnum_l[$zaehler]  $now\n";
	    $zaehler++;
	    
         }


         if (!($2 eq $name) && !($2 eq "SOL")) {
	    $resnum_p[$counter]=$1;
	    $resname_p[$counter]=$2;
            $atmname_p[$counter]=$3;
            $atmnum_p[$counter]=$4;
            $x_p[$counter]=$5;
	    $y_p[$counter]=$6;
	    $z_p[$counter]=$7;

	    if ($x_p[$counter]>$protein_xmax) {$protein_xmax=$x_p[$counter];}
	    if ($x_p[$counter]<$protein_xmin) {$protein_xmin=$x_p[$counter];}
	    if ($y_p[$counter]>$protein_ymax) {$protein_ymax=$y_p[$counter];}
	    if ($y_p[$counter]<$protein_ymin) {$protein_ymin=$y_p[$counter];}

	    if ($atmnum_p[$counter]==9999) {$now=1;}
	    #print "$resnum_p[$counter] \n";
	    $counter++;
         }


      }


}


close (INPUT);







$zaehler--;
$counter--;



$totalatmn=$zaehler+$counter;


# Converting nm into A

$protein_xmin=$protein_xmin*10;
$protein_xmax=$protein_xmax*10;
$protein_ymin=$protein_ymin*10;
$protein_ymax=$protein_ymax*10;








# New boxsize

$box_x=$box_x*$scale;
$box_y=$box_y*$scale;










# Scaling P positions & calculating translation vector

print STDOUT "Scaling lipids....\n";

$pcount=1;


for ($k=1; $k<=$zaehler; $k++) {

if (substr($atmname_l[$k],0,1) eq "P") {

    $pxneu=$x_l[$k]*$scale;
    $pyneu=$y_l[$k]*$scale;
    
    $res=$resnum_l[$k];
    $translatex_l[$res]=$pxneu-$x_l[$k];
    $translatey_l[$res]=$pyneu-$y_l[$k];
#    $phos_x[$pcount]=$x_l[$k]+$translatex_l[$res];
#    $phos_y[$pcount]=$y_l[$k]+$translatey_l[$res];
    $phosz[$pcount]=$z_l[$k];
    $pcount++;
    }
}


$pcount--;

print "There are $pcount lipids...\n";

$atomperlipid= $zaehler/$pcount;

print "with $atomperlipid atoms per lipid..\n";










# Determination of upper & lower leaflet

print "\nDetermining upper and lower leaflet...\n";

$middle=0;

for ($p=1;$p<=$pcount; $p++) {

    $middle=$middle + $phosz[$p];
    
    }
    

$middle=$middle / $pcount;



$uppercount=0;
$lowercount=0;
$upper=0;
$lower=0;


for ($p=1;$p<=$pcount; $p++) {

    if ($phosz[$p]>$middle) {$upper=$upper+$phosz[$p];
                             $uppercount++;}
    if ($phosz[$p]<$middle) {$lower=$lower+$phosz[$p];
                             $lowercount++;}
    
    }
$upper=$upper/$uppercount;
$lower=$lower/$lowercount;


print "$uppercount lipids in the upper...\n";
print "$lowercount lipids in the lower leaflet \n\n"; 





#Determining protein XY-COM & calculating translation vector

if ($counter==0) {print "No protein coordinates found...\n";}
if ($counter>0) {

print STDOUT "Centering protein....\n";

$xpsum=0;
$ypsum=0;


for ($k=1; $k<=$counter; $k++) {
     $xpsum=$xpsum+$x_p[$k];
     $ypsum=$ypsum+$y_p[$k];
     }

$xcom=$xpsum/$counter;
$ycom=$ypsum/$counter;


$xcenter=0.5*$box_x;
$ycenter=0.5*$box_y;

$translatex_p=$xcenter - $xcom;
$translatey_p=$ycenter - $ycom;


#print "COM:   $xcom    $ycom\n";
#print "New center: $xcenter   $ycenter\n";
#print "Translation vector: $translatex_p   $translatey_p\n\n";


}










# Checking for protein lipid overlap
$upper_rm=0;
$lower_rm=0;

if ($cutoff>0) {
if ($switch==0) {


print "Checking for overlap....\n";
print "...this might actually take a while....\n";

$overlapcount=0;
for ($k=1; $k<=$zaehler; $k++) {
    $uppercheck=0;
    $lowercheck=0;
    $progress=($k/$zaehler)*100;
    $progress=round ($progress,2);
    print STDOUT "$progress % done...\r";
    if (substr($atmname_l[$k],0,1) eq "P") {
    
             $res=$resnum_l[$k];
	     $overlap[$res]=0;
             
	     for ($i=1; $i<=$counter; $i++) {
	     
	          if ($atmname_p[$i] eq "CA") {
		  
		      $distance=abstand ($x_l[$k]+$translatex_l[$res],$y_l[$k]+$translatey_l[$res],$z_l[$k],$x_p[$i]+$translatex_p,$y_p[$i]+$translatey_p,$z_p[$i]);
		      if ($distance<=$cutoff) {
		           $overlap[$res]=1;
		   	   if ($z_l[$k] > $middle) {$uppercheck=1;}
	                   if ($z_l[$k] < $middle) {$lowercheck=1;}

			   }
		      
		  }

	     }
	     $overlapcount=$overlapcount+$overlap[$res];
	     if ($uppercheck==1) {$upper_rm++;}
	     if ($lowercheck==1) {$lower_rm++;}
	     
	     
    }

}

}
print "\nThere are $overlapcount lipids within cut-off range...\n";
print "$upper_rm will be removed from the upper leaflet...\n";
print "$lower_rm will be removed from the lower leaflet...\n\n";

}



$newlipids=$pcount - $upper_rm - $lower_rm;
$newupper=$uppercount - $upper_rm;
$newlower=$lowercount - $lower_rm;
$totalatmn_new=$totalatmn - ($upper_rm + $lower_rm)*$atomperlipid;







# Writing scaled bilayer & centered protein


print STDOUT "Writing scaled bilayer & centered protein...\n";

open(OUTPUT, ">$output");

print OUTPUT "What you read here has nothing to do with anything. So you don't have to read it. Thank you.\n";
print OUTPUT "$totalatmn_new\n";

for ($k=1; $k<=$counter; $k++) {

    $newx=$x_p[$k]+$translatex_p;
    $newy=$y_p[$k]+$translatey_p;
    
   # print "$newx    $newy\n";
    
     printf OUTPUT "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",$resnum_p[$k],$resname_p[$k],$atmname_p[$k],$atmnum_p[$k],$newx,$newy,$z_p[$k]; 

#     printf OUTPUT "ATOM%7d %4s %-4s%5d %11.3f %7.3f %7.3f %5.2f %5.2f\n",$k,$atmname_p[$k],$resname_p[$k],$resnum_p[$k],$newx,$newy,$z_p[$k],$occupancy_p[$k],$bfactor_p[$k];



}

if ($switch==0) {


for ($k=1; $k<=$zaehler; $k++) {
    $res=$resnum_l[$k];
    $newx=$x_l[$k]+$translatex_l[$res];
    $newy=$y_l[$k]+$translatey_l[$res];


     if ($cutoff>0) {
     if ($overlap[$res]==0) {

         printf OUTPUT "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",$resnum_l[$k],$resname_l[$k],$atmname_l[$k],$atmnum_l[$k],$newx,$newy,$z_l[$k]; 
     
#         printf OUTPUT "ATOM%7d %4s %-4s%5d %11.3f %7.3f %7.3f %5.2f %5.2f\n",$k,$atmname_l[$k],$resname_l[$k],$resnum_l[$k],$newx,$newy,$z_l[$k],$occupancy_l[$k],$bfactor_l[$k];
     }
     }


     if ($cutoff==0) {

         printf OUTPUT "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",$resnum_l[$k],$resname_l[$k],$atmname_l[$k],$atmnum_l[$k],$newx,$newy,$z_l[$k]; 
     
#         printf OUTPUT "ATOM%7d %4s %-4s%5d %11.3f %7.3f %7.3f %5.2f %5.2f\n",$k,$atmname_l[$k],$resname_l[$k],$resnum_l[$k],$newx,$newy,$z_l[$k],$occupancy_l[$k],$bfactor_l[$k];
     }




}

}

printf OUTPUT "%10.5f%10.5f%10.5f\n",$box_x,$box_y,$box_z;
close OUTPUT;






# Translate protein to Xmin & Ymin = 0


print "\n\nCalculating Area per lipid...\n";


$protein_xmax=int ($protein_xmax+1);
$protein_xmin=int ($protein_xmin);
$protein_ymax=int ($protein_ymax+1);
$protein_ymin=int ($protein_ymin);

$xrange=$protein_xmax-$protein_xmin;
$yrange=$protein_ymax-$protein_ymin;

print "Protein X-min/max: $protein_xmin    $protein_xmax\n";
print "Protein Y-min/max: $protein_ymin    $protein_ymax\n";
print "X-range: $xrange A    Y-range: $yrange A\n";



if ($protein_xmin != 0 or $protein_xmax != 0){

for ($k=1; $k<=$counter; $k++) {

    $x_p[$k]=10*$x_p[$k]-$protein_xmin;
    $y_p[$k]=10*$y_p[$k]-$protein_ymin;
    
    }   
}




# Building 2D grid on protein coordinates


print "Building $xrange X $yrange 2D grid on protein coordinates...\n";

#for ($x=0; $x<=$xrange; $x=$x+$gridsize) {
#    for ($y=0; $y<=$yrange; $y=$y+$gridsize) {
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
        $grid[$x][$y]=0;
    }
}



# Calculating area occupied by protein

print "Calculating area occupied by protein..\n";

print "full TMD..\n";

for ($k=1; $k<=$counter; $k++) {
  if ($z_p[$k]>=$lower and $z_p[$k]<=$upper) {
     $x = int( $x_p[$k] / $gridsize);  
     $y = int( $y_p[$k] / $gridsize);
     $grid[$x][$y]=1;
     }
$progress=$k/$counter *100;
$progress=round ($progress,1);
print "$progress % done...\r";
}

$howmany=0;

for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
   
         $howmany=$howmany+$grid[$x][$y];
	 
    }
}



$areaprotein_total=($gridsize)**2 *$howmany *0.01;


$arealipid_total=($box_x * $box_y - $areaprotein_total)/($newlipids*0.5);





print "upper TMD..\n";


for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
        $grid[$x][$y]=0;
    }
}



for ($k=1; $k<=$counter; $k++) {
  if ($z_p[$k]>=$middle and $z_p[$k]<=$upper) {
     $x = int( $x_p[$k] / $gridsize);  
     $y = int( $y_p[$k] / $gridsize);
     $grid[$x][$y]=1;
     }
$progress=$k/$counter *100;
$progress=round ($progress,1);
print "$progress % done...\r";
}

$howmany=0;

for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
   
         $howmany=$howmany+$grid[$x][$y];
	 
    }
}



$areaprotein_upper=($gridsize)**2 *$howmany *0.01;


$arealipid_upper=($box_x * $box_y - $areaprotein_upper)/($newupper);



print "lower TMD..\n";



for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
        $grid[$x][$y]=0;
    }
}



for ($k=1; $k<=$counter; $k++) {
  if ($z_p[$k]>=$lower and $z_p[$k]<=$middle) {
     $x = int( $x_p[$k] / $gridsize);  
     $y = int( $y_p[$k] / $gridsize);
     $grid[$x][$y]=1;
     }
$progress=$k/$counter *100;
$progress=round ($progress,1);
print "$progress % done...\r";
}

$howmany=0;

for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
   
         $howmany=$howmany+$grid[$x][$y];
	 
    }
}



$areaprotein_lower=($gridsize)**2 *$howmany *0.01;


$arealipid_lower=($box_x * $box_y - $areaprotein_lower)/($newlower);


print "Area per protein: $areaprotein_total nm^2\n";
print "Area per lipid: $arealipid_total nm^2\n\n";

print "Area per protein, upper half: $areaprotein_upper nm^2\n";
print "Area per lipid, upper leaflet : $arealipid_upper nm^2\n\n";

print "Area per protein, lower half: $areaprotein_lower nm^2\n";
print "Area per lipid, lower leaflet : $arealipid_lower nm^2\n\n";



print STDOUT "Writing Area per lipid...\n";

open(OUTPUT, ">$area");


# Uncomment this section if you want a plot of the protein area

# for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
#     for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
#        	 if ($grid[$x][$y]==1) {
# 	 $xtmp = $gridsize*$x +0.5*$gridsize;
# 	 $ytmp = $gridsize*$y +0.5*$gridsize;
# 	 print OUTPUT "$xtmp     $ytmp\n";}
 	 #if ($grid[$x][$y]==1) { print OUTPUT "$x     $y\n";}
 	 
#     }
# }

print OUTPUT "$arealipid_total     $arealipid_upper      $arealipid_lower \n";

close OUTPUT;


print "Done!\n\n\n";



