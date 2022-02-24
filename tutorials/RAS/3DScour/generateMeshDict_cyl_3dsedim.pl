#!/usr/bin/perl
# Generate the cylinder blockMeshDict file for the von Karman LES case
#              16          17         18          19
#  y ^        0 - - - - - 1- - - - - - 2 - - - - - 3
#    |        |           |            |           |  
#    --> x    |     0     |    1       |   2       |     ny1
#             |20         |21          |22         |23
#             4 - - - - - 5 - - - - - -6 - - - - - 7
#
#             |           |32/36 33/37 |           |  
#
#             |    3      |34/38 35/39 |   4       |     ny2
#
#             |24         |25          |26         |27
#             8 - - - - - 9 - - - - - 10 - - - - -11
#             |           |            |           |  
#             |    5      |   6        |   7       |     ny1
#             |28         |29          |30         |31
#             12 - - - - 13 - - - - - 14 - - - - - 15
#                                                   
#                 nx1       nx2=ny2       nx3
use strict;

if ($#ARGV != 7)
{
    usage();
    exit;
}

my @input = @ARGV;
my $D = $ARGV[0];
my $H0 = $ARGV[1];
my $H1 = $ARGV[2];
my $H2 = $ARGV[3];
my $H3 = $ARGV[4]; #for the bottom block
my $L1 = $ARGV[5];
my $L2 = $ARGV[6];
my $W = $ARGV[7];

my $w = 0.5*$W;
my $r = 0.5*$D;
my $r0 = 0.5*$D;
my $rinner1 = $r0/4.;
my $rinner = $rinner1/sqrt(2.);
my $xe = $L2-$L1;
my $rext = $r*sqrt(2.);


# aspect ratio for largest/smallest cells for the different blocks
my $g1 = 1.05710412557;#1.0056312014;#1.04286748796;#1.89325432102;
my $g2 = 1.00394312014;#1.04286748796;#1.89325432102;
my $g3 = 1.05710412557;#2.99326472428;#3.66407505143;
my $gC = 10;#7.89;
my $gZ = 100;
my $gZ1 = 0.39576022153092771;
my $gZ1b = 0.29844065970697592;#0.129671569375;#0.0189909804448;#0.13;

my $nx1 = 22;#40;
my $nx2 = 64;#60; 
my $nx3 = 44;#200;

my $ny1 = 11;#40;
my $ny2 = $nx2;

my $nz = 64;
my $nz1 = 100;
my $nz1b = 100;

my $g1i = 1.0/$g1;
my $g2i = 1.0/$g2;
my $g3i = 1.0/$g3;
my $gCi = 1.0/$gC;
my $gZi = 1.0/$gZ;
my $gZi1 = $gZ1;
my $gZi1b = $gZ1b;

printHeader("2.0");
printConvertToMeters("1.0");
printVertices();
printBlocks();
printEdges();
printPatches();
printMerge();

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sub printHeader
{
    print "FoamFile\n";
    print "{\n";
    print "    version\t $_[0];\n";
    print "    format\t ascii;\n";
    print "    class\t dictionary;\n";
    print "    object\t blockMeshDict;\n";
    print "}\n\n";

    print "/*\n";
    print "Automatically generated with $0 @input\n";
    print "*/\n\n";
}

sub printConvertToMeters
{
    print "convertToMeters\t $_[0];\n\n";
}

sub printVertices
{
    print "vertices\n";
    print "(\n";
    print "    ( -$L1  $w $H1 ) //0\n";
    print "    ( -$r   $w $H1 ) //1\n";
    print "    (  $r   $w $H1 ) //2\n";
    print "    (  $xe  $w $H1 ) //3\n";

    print "    ( -$L1  $r $H1 ) //4\n";
    print "    ( -$r   $r $H1 ) //5\n";
    print "    (  $r   $r $H1 ) //6\n";
    print "    (  $xe  $r $H1 ) //7\n";

    print "    ( -$L1  -$r $H1 ) //8\n";
    print "    ( -$r   -$r $H1 ) //9\n";
    print "    (  $r   -$r $H1 ) //10\n";
    print "    (  $xe  -$r $H1 ) //11\n";

    print "    ( -$L1  -$w $H1 ) //12\n";
    print "    ( -$r   -$w $H1 ) //13\n";
    print "    (  $r   -$w $H1 ) //14\n";
    print "    (  $xe  -$w $H1 ) //15\n";

    print "    ( -$L1  $w $H2 ) //16\n";
    print "    ( -$r   $w $H2 ) //17\n";
    print "    (  $r   $w $H2 ) //18\n";
    print "    (  $xe  $w $H2 ) //19\n";

    print "    ( -$L1  $r $H2 ) //20\n";
    print "    ( -$r   $r $H2 ) //21\n";
    print "    (  $r   $r $H2 ) //22\n";
    print "    (  $xe  $r $H2 ) //23\n";

    print "    ( -$L1  -$r $H2 ) //24\n";
    print "    ( -$r   -$r $H2 ) //25\n";
    print "    (  $r   -$r $H2 ) //26\n";
    print "    (  $xe  -$r $H2 ) //27\n";

    print "    ( -$L1  -$w $H2 ) //28\n";
    print "    ( -$r   -$w $H2 ) //29\n";
    print "    (  $r   -$w $H2 ) //30\n";
    print "    (  $xe  -$w $H2 ) //31\n";

# vertices on the cylinder
    print "    ( -$rinner   $rinner $H1 ) //32\n";
    print "    (  $rinner   $rinner $H1 ) //33\n";
    print "    ( -$rinner   -$rinner $H1 ) //34\n";
    print "    (  $rinner   -$rinner $H1 ) //35\n";

    print "    ( -$rinner   $rinner $H2 ) //36\n";
    print "    (  $rinner   $rinner $H2 ) //37\n";
    print "    ( -$rinner   -$rinner $H2 ) //38\n";
    print "    (  $rinner   -$rinner $H2 ) //39\n";

# Lower blocks
    print "    ( -$L1  $w $H0 ) //40\n";
    print "    ( -$r   $w $H0 ) //41\n";
    print "    (  $r   $w $H0 ) //42\n";
    print "    (  $xe  $w $H0 ) //43\n";

    print "    ( -$L1  $r $H0 ) //44\n";
    print "    ( -$r   $r $H0 ) //45\n";
    print "    (  $r   $r $H0 ) //46\n";
    print "    (  $xe  $r $H0 ) //47\n";

    print "    ( -$L1  -$r $H0 ) //48\n";
    print "    ( -$r   -$r $H0 ) //49\n";
    print "    (  $r   -$r $H0 ) //45\n";
    print "    (  $xe  -$r $H0 ) //51\n";

    print "    ( -$L1  -$w $H0 ) //52\n";
    print "    ( -$r   -$w $H0 ) //53\n";
    print "    (  $r   -$w $H0 ) //54\n";
    print "    (  $xe  -$w $H0 ) //55\n";

    print "    ( -$rinner   $rinner $H0 ) //56\n";
    print "    (  $rinner   $rinner $H0 ) //57\n";
    print "    ( -$rinner   -$rinner $H0 ) //58\n";
    print "    (  $rinner   -$rinner $H0 ) //59\n";


# Scour pit blocks

    print "    ( -$r   $r $H3 ) //60\n";
    print "    (  $r   $r $H3 ) //61\n";

    print "    ( -$r   -$r $H3 ) //62\n";
    print "    (  $r   -$r $H3 ) //63\n";

    print "    ( -$rinner   $rinner $H3 ) //64\n";
    print "    (  $rinner   $rinner $H3 ) //65\n";
    print "    ( -$rinner   -$rinner $H3 ) //66\n";
    print "    (  $rinner   -$rinner $H3 ) //67\n";

    print ");\n\n";
}

sub printBlocks
{
    print "blocks\n";
    print "(\n";
    print "    hex (4 5 1 0 20 21 17 16) ($nx1 $ny1 $nz) simpleGrading ($g1i $g2 $gZ) //0\n";
    print "    hex (5 6 2 1 21 22 18 17) ($nx2 $ny1 $nz) simpleGrading (1 $g2 $gZ) //1\n";
    print "    hex (6 7 3 2 22 23 19 18) ($nx3 $ny1 $nz) simpleGrading ($g3 $g2 $gZ) //2\n";

    print "    hex (8 9 5 4 24 25 21 20) ($nx1 $ny2 $nz) simpleGrading ($g1i 1 $gZ) //3\n";
    print "    hex (10 11 7 6 26 27 23 22) ($nx3 $ny2 $nz) simpleGrading ($g3 1 $gZ)//4\n";

    print "    hex (12 13 9 8 28 29 25 24) ($nx1 $ny1 $nz) simpleGrading ($g1i $g2i $gZ)//5\n";
    print "    hex (13 14 10 9 29 30 26 25) ($nx2 $ny1 $nz) simpleGrading (1 $g2i $gZ)//6\n";
    print "    hex (14 15 11 10 30 31 27 26) ($nx3 $ny1 $nz) simpleGrading ($g3 $g2i $gZ)//7\n";

    print "    hex (32 33 6 5 36 37 22 21) ($nx2 $ny2 $nz) simpleGrading (1 $gC $gZ)//8\n";
    print "    hex (35 10 6 33 39 26 22 37) ($nx2 $ny2 $nz) simpleGrading ($gC 1 $gZ)//9\n";
    print "    hex (9 10 35 34 25 26 39 38) ($nx2 $ny2 $nz) simpleGrading (1 $gCi $gZ)//10\n";
    print "    hex (9 34 32 5 25 38 36 21) ($nx2 $ny2 $nz) simpleGrading ($gCi 1 $gZ)//11\n";

    print "    hex (44 45 41 40 4 5 1 0) ($nx1 $ny1 $nz1) simpleGrading ($g1i $g2 $gZ1) //12\n";
    print "    hex (45 46 42 41 5 6 2 1) ($nx2 $ny1 $nz1) simpleGrading (1 $g2 $gZ1) //13\n";
    print "    hex (46 47 43 42 6 7 3 2) ($nx3 $ny1 $nz1) simpleGrading ($g3 $g2 $gZ1) //14\n";

    print "    hex (48 49 45 44 8 9 5 4) ($nx1 $ny2 $nz1) simpleGrading ($g1i 1 $gZ1) //15\n";
    print "    hex (50 51 47 46 10 11 7 6) ($nx3 $ny2 $nz1) simpleGrading ($g3 1 $gZ1)//16\n";

    print "    hex (52 53 49 48 12 13 9 8) ($nx1 $ny1 $nz1) simpleGrading ($g1i $g2i $gZ1)//17\n";
    print "    hex (53 54 50 49 13 14 10 9) ($nx2 $ny1 $nz1) simpleGrading (1 $g2i $gZ1)//18\n";
    print "    hex (54 55 51 50 14 15 11 10) ($nx3 $ny1 $nz1) simpleGrading ($g3 $g2i $gZ1)//19\n";

    print "    hex (56 57 46 45 32 33 6 5) ($nx2 $ny2 $nz1) simpleGrading (1 $gC $gZ1)//20\n";
    print "    hex (59 50 46 57 35 10 6 33) ($nx2 $ny2 $nz1) simpleGrading ($gC 1 $gZ1)//21\n";
    print "    hex (49 50 59 58 9 10 35 34) ($nx2 $ny2 $nz1) simpleGrading (1 $gCi $gZ1)//22\n";
    print "    hex (49 58 56 45 9 34 32 5) ($nx2 $ny2 $nz1) simpleGrading ($gCi 1 $gZ1)//23\n";

    print "    hex (64 65 61 60 56 57 46 45) ($nx2 $ny2 $nz1b) simpleGrading (1 $gC $gZ1b)//24\n";
    print "    hex (67 63 61 65 59 50 46 57) ($nx2 $ny2 $nz1b) simpleGrading ($gC 1 $gZ1b)//25\n";
    print "    hex (62 63 67 66 49 50 59 58) ($nx2 $ny2 $nz1b) simpleGrading (1 $gCi $gZ1b)//26\n";
    print "    hex (62 66 64 60 49 58 56 45) ($nx2 $ny2 $nz1b) simpleGrading ($gCi 1 $gZ1b)//27\n";


    print ");\n\n";
}

sub printEdges
{
    print "edges\n";
    print "(\n";
    print "    arc 32 33 ( 0 $rinner1 $H1 )\n";
    print "    arc 33 35 ( $rinner1 0 $H1 )\n";
    print "    arc 34 35 ( 0 -$rinner1 $H1 )\n";
    print "    arc 34 32 ( -$rinner1 0 $H1 )\n";

    print "    arc 36 37 ( 0 $rinner1 $H2 )\n";
    print "    arc 39 37 ( $rinner1 0 $H2 )\n";
    print "    arc 38 39 ( 0 -$rinner1 $H2 )\n";
    print "    arc 38 36 ( -$rinner1 0 $H2 )\n";

    print "    arc 56 57 ( 0 $rinner1 $H0 )\n";
    print "    arc 59 57 ( $rinner1 0 $H0 )\n";
    print "    arc 58 59 ( 0 -$rinner1 $H0 )\n";
    print "    arc 58 56 ( -$rinner1 0 $H0 )\n";

    print "    arc 64 65 ( 0 $rinner1 $H3 )\n";
    print "    arc 67 65 ( $rinner1 0 $H3 )\n";
    print "    arc 66 67 ( 0 -$rinner1 $H3 )\n";
    print "    arc 66 64 ( -$rinner1 0 $H3 )\n";

    print "    arc 5 6 ( 0 $rext $H1 )\n";
    print "    arc 6 10 ( $rext 0 $H1 )\n";
    print "    arc 9 10 ( 0 -$rext $H1 )\n";
    print "    arc 9 5 ( -$rext 0 $H1 )\n";

    print "    arc 21 22 ( 0 $rext $H2 )\n";
    print "    arc 26 22 ( $rext 0 $H2 )\n";
    print "    arc 25 26 ( 0 -$rext $H2 )\n";
    print "    arc 25 21 ( -$rext 0 $H2 )\n";

    print "    arc 45 46 ( 0 $rext $H0 )\n";
    print "    arc 46 50 ( $rext 0 $H0 )\n";
    print "    arc 49 50 ( 0 -$rext $H0 )\n";
    print "    arc 49 45 ( -$rext 0 $H0 )\n";

    print "    arc 60 61 ( 0 $rext $H3 )\n";
    print "    arc 63 61 ( $rext 0 $H3 )\n";
    print "    arc 62 63 ( 0 -$rext $H3 )\n";
    print "    arc 62 60 ( -$rext 0 $H3 )\n";


    print ");\n\n";
}

sub printPatches
{
    print "patches\n";
    print "(\n";
    print "    wall cylinder\n";
    print "    (\n";
    print "        (32 36 37 33)\n";
    print "        (33 37 39 35)\n";
    print "        (35 39 38 34)\n";
    print "        (34 38 36 32)\n";
    print "        (56 32 33 57)\n";
    print "        (57 33 35 59)\n";
    print "        (59 35 34 58)\n";
    print "        (58 34 32 56)\n";
    print "        (64 56 57 65)\n";
    print "        (65 57 59 67)\n";
    print "        (67 59 58 66)\n";
    print "        (66 58 56 64)\n";
    print "    )\n\n";

    print "    patch inlet\n";
    print "    (\n";
    print "        (4 20 16 0)\n";
    print "        (8 24 20 4)\n";
    print "        (12 28 24 8)\n";
    print "        (44 4 0 40)\n";
    print "        (48 8 4 44)\n";
    print "        (52 12 8 48)\n";
    print "    )\n";

    print "    patch outlet\n";
    print "    (\n";
    print "        (3 19 23 7)\n";
    print "        (7 23 27 11)\n";
    print "        (11 27 31 15)\n";
    #print "        (43 3 7 47)\n";
    #print "        (47 7 11 51)\n";
    #print "        (51 11 15 55)\n";
    print "    )\n";


    print "    wall outletb\n";
    print "    (\n";
    print "        (43 3 7 47)\n";
    print "        (47 7 11 51)\n";
    print "        (51 11 15 55)\n";
    #print "        (19 55 51 23)\n";
    #print "        (23 51 47 27)\n";
    #print "        (27 47 43 31)\n";
    print "    )\n";


    print "    wall pit\n";
    print "    (\n";
    print "        (45 46 61 60)\n";
    print "        (62 63 50 49)\n";
    print "        (61 46 50 63)\n";
    print "        (62 49 45 60)\n";
    print "    )\n";


    print "    patch lateral\n";
    print "    (\n";
    print "        (16 17 1 0)\n";
    print "        (17 18 2 1)\n";
    print "        (18 19 3 2)\n";
    print "        (0 1 41 40)\n";
    print "        (1 2 42 41)\n";
    print "        (2 3 43 42)\n";
    print "    )\n";

    print "    patch symplane\n";
    print "    (\n";
    print "        (12 13 29 28)\n";
    print "        (13 14 30 29)\n";
    print "        (14 15 31 30)\n";
    print "        (52 53 13 12)\n";
    print "        (53 54 14 13)\n";
    print "        (54 55 15 14)\n";
    print "    )\n";

    print "    wall bottom\n";
    print "    (\n";
    print "        (40 41 45 44)\n";
    print "        (41 42 46 45)\n";
    print "        (42 43 47 46)\n";
    print "        (44 45 49 48)\n";
    print "        (46 47 51 50)\n";
    print "        (48 49 53 52)\n";
    print "        (49 50 54 53)\n";
    print "        (50 51 55 54)\n";
#    print "        (45 46 57 56)\n";
#    print "        (46 50 59 57)\n";
#    print "        (58 59 50 49)\n";
#    print "        (45 56 58 49)\n";
    print "        (60 61 65 64)\n";
    print "        (61 63 67 65)\n";
    print "        (66 67 63 62)\n";
    print "        (60 64 66 62)\n";
    print "    )\n";

    print "    patch surface\n";
    print "    (\n";
    print "        (20 21 17 16)\n";
    print "        (21 22 18 17)\n";
    print "        (22 23 19 18)\n";
    print "        (24 25 21 20)\n";
    print "        (26 27 23 22)\n";
    print "        (28 29 25 24)\n";
    print "        (29 30 26 25)\n";
    print "        (30 31 27 26)\n";
    print "        (36 37 22 21)\n";
    print "        (39 26 22 37)\n";
    print "        (25 26 39 38)\n";
    print "        (25 38 36 21)\n";

    print "    )\n";
    print ");\n\n";
}

sub printMerge
{
    print "mergPatchPairs\n";
    print "(\n";
    print ");\n";
}

sub usage
{
    print "usage:\n";
    print "\t $0 <D> <H0> <H1> <H2> <H3> <L1> <L2> <W> <delta>\n";
    print "\t where (recommended value):\n";
    print "\t D  = Width of the refinment region around the cylinder. (4)\n";
    print "\t H0 = Elevation of the fixed sediment bottom. (-1)\n";
    print "\t H1 = Elevation of the sediment bed-fluid interface. (0.)\n";
    print "\t H2 = Elevation of the top boundary. (2.)\n";
    print "\t H3 = Elevation of the pit boundary. (2.)\n";
    print "\t L1 = Length from inlet to cylinder center. (5.)\n";
    print "\t L2 = Length of computational domain. (10.)\n";
    print "\t W = Width of the domain (10)\n";
    print "\t delta = cell size across the sides of the cylinder. (0.001)\n";
}

