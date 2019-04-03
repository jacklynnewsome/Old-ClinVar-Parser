#!/usr/bin/perl -w
use strict;
use warnings;
use IO::Handle;

#get file name
 my $FILE = $ARGV[0];
#open file
	#check file open
	open ( my $readfile, "<", $FILE )
	or die "Could not open file '$FILE' $!";
		#report error

	my @chrom; # chromosome
	my @position; #position
	my @rscombo; #RS###########
	my @consntp; #consensus nucleotide 
	my @varntp; #variation nucleotide
	my @rs; #RS number by itself
	my @rspos; #RSPOS
	my @dbsnpbuild; #dbSNPBuild number
	my @ssr; #SSR number
	my @sao; #SAO number
	my @vp; #VP (complicated set of stuff)
	my @geneinfocombo; #GENEINFO - both the gene symbols and their gene ids
	my @genesym1;  #GENEINFO - gene symbol for the first of the pair
	my @geneid1; #GENEINFO - gene id number for the first of the pair
	my @genesym2; #GENEINFO - gene symbol for the second of the pair
	my @geneid2; #GENEINFO - gene id number for the second of the pair
	my @wgt; #WGT
	my @vc; #
	my @vcdef; #VC definition
		#set 6 of these??
	my @pm; #PM flag
	my @tpa; #TPA flag
	my @pmc; #PMC flag
	my @s3d; #S3D flag
	my @slo; #SLO flag
	my @nsf; #NSF flag
	my @nsm; #NSM flag
	my @nsn; #NSN flag
	my @ref; #REF flag
 	my @syn; #SYN flag
	my @u3; #U3 flag
	my @u5; #U5 flag
	my @ass; #ASS flag
	my @dss; #DSS flag
	my @int; #INT flag
	my @r3; #R3 flag
	my @r5; #R5 flag
	my @oth; #OTH flag
	my @cfl; #CFL flag
	my @asp; #ASP flag
	my @mut; #MUT flag
	my @vld; #VLD flag
	my @g5a; #G5A flag
	my @g5; #G5 flag
	my @hd; #HD flag
	my @gno; #GNO flag
	my @kgval; #KGValidated flag
	my @kgpha; #KGPhase1 flag
	my @kgpil; #KGPilot123 flag
	my @kgprod; #KGPROD flag
	my @otherkg; #OTHERKG flag
	my @ph3; #PH3 flag
	my @cda; #CDA flag
	my @lsd; #LSD flag
	my @mtp; #MTP flag
	my @om; #OM flag
	my @noc; #NOC flag
	my @wtd; #WTD flag
	my @nov; #NOV flag
	my @nc; #NC string
	my @caf; #CAF string
	my @common; #COMMON value
	my @clnhgvs; #CLNHGVS string
	my @clnalle; #CLNSRC String
	my @clnsrc; #CLNSRC String
	my @clndsdbid; #CLNDSDBID string
	my @clndsdb; #CLNDSDB string
	my @clndbn; #CLNDBN string
	my @clnori; #CLNORIGIN string
	my @clnsrcid; #CLNSRCID string	
	my @clnsig; #CLNSIG string
	my @clnacc; #CLNALLE value				
	

	
	# set line array 
	my $line;
	my @lines;


my $i = 0;
#start loop
while ( $line = <$readfile> )
{
	#for each line:

		$lines[$i] = split( /\n/, $line);
	#skip header lines
		#if starts with sharp, skip
	if ( $line =~ /^\s*(\#)/g ) { 
		next; 
	};

#########  p   Preserve the string matched such that ${^PREMATCH}, {$^MATCH}, and ${^POSTMATCH} are available for use after matching.
			#read chrom, position, rscombo, consensus nucleotide(s), variation nucleotides (s)
			if ( $line =~ /^\s*(\w+)\t(\w+)\t(\w+)\t(\w+)\t(\w+)/) { 
					$chrom[$i] = $1;
					$position[$i] = $2;
					$rscombo[$i] = $3;			
					$consntp[$i] = $4;
					$varntp[$i] = $5;
					
				} else {
						print "Could not read line '$line' !";
						last;
			};
			
			
			### Do I need to include "m/" so it searches through the entire line?
			
				
			#read RS number
			if ( $line =~ /.(\bRS=)(\w+)\;/ ) {
					$rs[$i] = $2;
					#if not applicable, set null value?
				} else {
					$rs[$i] = undef;
			};	
			
			#read RSPOS number
			if ( $line =~ /.(\bRSPOS=)(\w+)\;/ ) {
					$rspos[$i] = $2;
					#if not applicable, set null value?
				} else {
					$rspos[$i] = undef;
			};

			
			#read dbSNPBUild number
			if ( $line =~ /.(\bdbSNPBuildID=)(\w+)\;/ ) {
					$dbsnpbuild[$i] = $2;
					#if not applicable, set null value?
				} else {
					$dbsnpbuild[$i] = undef;
			};
				
			#read SSR number
			if ( $line =~ /.(\bSSR=)(\w+)\;/ ) {
					$ssr[$i] = $2;
					#if not applicable, set null value?
				} else {
					$ssr[$i] = undef;
			};
		
			
				
				
				
			#read SAO value
			if ( $line =~ /.(\bSAO=)(\w+)\;/ ) {
					$sao[$i] = $2;
					#if not applicable, set null value?
				} else {
					$sao[$i] = undef;
			};
		
			#read VP                12345678901234567890123456
			if ( $line =~ /.\bVP=(\w{26})\;/ ) {
					$vp[$i] = $1;
					#if not applicable, set null value?
				} else {
					$vp[$i] = undef;
			};
		
			
			
			
			#read GENEINFO set
			if ( $line =~ /(\bGENEINFO=)([a-zA-Z0-9.:|]+)\;/ ) {
					$geneinfocombo[$i] = $2;
					#if not applicable, set null value?
				} else {
					$geneinfocombo[$i] = undef;
			};
			

			
			#read WGT
			if ( $line =~ /.(\bWGT=)(\w+)\;/ ) {
					$wgt[$i] = $2;
					#if not applicable, set null value?
				} else {
					$wgt[$i] = undef;
			};
		
			#read VC value
			if ( $line =~ /.(\bVC=)(\w+)\;/ ) {
					$vc[$i] = $2;
					#if not applicable, set null value?
				#if NA, set null value?
				} else {
					$vc[$i] = undef;
			};
			

			#read PM flag
			if ( $line =~ /.(\bPM)\;/ ) {
					$pm[$i] = "PM";
					#if not applicable, set null value?
				} else {
					$pm[$i] = undef;
			};			

			#read TPA flag
			if ( $line =~ /.(\bTPA)\;/ ) {
					$tpa[$i] = "TPA";
					#if not applicable, set null value?
				} else {
					$tpa[$i] = undef;
			};	
			#read PMC flag
			if ( $line =~ /.(\bPMC)\;/ ) {
					$pmc[$i] = "PMC";
					#if not applicable, set null value?
				} else {
					$pmc[$i] = undef;
			};	
			#read S3D flag
			if ( $line =~ /.(\bS3D)\;/ ) {
					$s3d[$i] = "S3D";
					#if not applicable, set null value?
				} else {
					$s3d[$i] = undef;
			};	
			#read SLO flag
			if ( $line =~ /.(\bSLO)\;/ ) {
					$slo[$i] = "SLO";
				} else {
					$slo[$i] = undef;
			};	
			#read NSF flag
			if ( $line =~ /.(\bNSF)\;/ ) {
					$nsf[$i] = "NSF";
				} else {
					$nsf[$i] = undef;
			};	
			#read NSM flag
			if ( $line =~ /.(\bNSM)\;/ ) {
					$nsm[$i] = "NSM";
				} else {
					$nsm[$i] = undef;
			};
			#read NSN flag
			if ( $line =~ /.(\bNSN)\;/ ) {
					$nsn[$i] = "NSN";
				} else {
					$nsn[$i] = undef;
			};
			#read REF flag
			if ( $line =~ /.(\bREF)\;/ ) {
					$ref[$i] = "REF";
				} else {
					$ref[$i] = undef;
			};
			#read SYN flag
			if ( $line =~ /.(\bSYN)\;/ ) {
					$syn[$i] = "SYN";
				} else {
					$syn[$i] = undef;
			};
			#read U3 flag
			if ( $line =~ /.(\bU3)\;/ ) {
					$u3[$i] = "U3";
				} else {
					$u3[$i] = undef;
			};
			#read U5 flag
			if ( $line =~ /.(\bU5)\;/ ) {
					$u5[$i] = "U5";
				} else {
					$u5[$i] = undef;
			};
			#read ASS flag
			if ( $line =~ /.(\bASS)\;/ ) {
					$ass[$i] = "ASS";
				} else {
					$ass[$i] = undef;
			};
			#read DSS flag
			if ( $line =~ /.(\bDSS)\;/ ) {
					$dss[$i] = "DSS";
				} else {
					$dss[$i] = undef;
			};
			#read INT flag
			if ( $line =~ /.(\bINT)\;/ ) {
					$int[$i] = "INT";
				} else {
					$int[$i] = undef;
			};
			#read R3 flag
			if ( $line =~ /.(\bR3)\;/ ) {
					$r3[$i] = "R3";
				} else {
					$r3[$i] = undef;
			};
			#read R5 flag
			if ( $line =~ /.(\bR5)\;/ ) {
					$r5[$i] = "R5";
				} else {
					$r5[$i] = undef;
			};			
			#read OTH flag
			if ( $line =~ /.(\bOTH)\;/ ) {
					$oth[$i] = $1;
				} else {
					$oth[$i] = undef;
			};			
			#read CFL flag
			if ( $line =~ /.(\CFL)\;/ ) {
					$cfl[$i] = "CFL";
				} else {
					$cfl[$i] = undef;
			};	
			#read ASP flag
			if ( $line =~ /.(\bASP)\;/ ) {
					$asp[$i] = "ASP";
				} else {
					$asp[$i] = undef;
			};	
			#read MUT flag
			if ( $line =~ /.(\bMUT)\;/ ) {
					$mut[$i] = "MUT";
				} else {
					$mut[$i] = undef;
			};	
			#read VLD flag
			if ( $line =~ /.(\bVLD)\;/ ) {
					$vld[$i] = "VLD";
				} else {
					$vld[$i] = undef;
			};	
			#read G5A flag
			if ( $line =~ /.(\bG5A)\;/ ) {
					$g5a[$i] = "G5A";
				} else {
					$g5a[$i] = undef;
			};	
			#read G5 flag
			if ( $line =~ /.(\bG5)\;/ ) {
					$g5[$i] = "G5";
				} else {
					$g5[$i] = undef;
			};	
			#read HD flag
			if ( $line =~ /.(\bHD)\;/ ) {
					$hd[$i] = "HD";
				} else {
					$hd[$i] = undef;
			};
			#read GNO flag
			if ( $line =~ /.(\bGNO)\;/ ) {
					$gno[$i] = "GNO";
				} else {
					$gno[$i] = undef;
			};
			#read KGValidated flag
			if ( $line =~ /.(\bKGValidated)\;/ ) {
					$kgval[$i] = "KGValidated";
				} else {
					$kgval[$i] = undef;
			};
			#read KGPhase1 flag
			if ( $line =~ /.(\bKGPhase1)\;/ ) {
					$kgpha[$i] = "KGPhase1";
				} else {
					$kgpha[$i] = undef;
			};
			#read KGPilot123 flag
			if ( $line =~ /.(\bKGPilot123)\;/ ) {
					$kgpil[$i] = "KGPilot123";
				} else {
					$kgpil[$i] = undef;
			};
			#read KGPROD flag
			if ( $line =~ /.(\bKGPROD)\;/ ) {
					$kgprod[$i] = "KGPROD";
				} else {
					$kgprod[$i] = undef;
			};
			#read OTHERKG flag
			if ( $line =~ /.(\bOTHERKG)\;/ ) {
					$otherkg[$i] = "OTHERKG";
				} else {
					$otherkg[$i] = undef;
			};
			#read PH3 flag
			if ( $line =~ /.(\bPH3)\;/ ) {
					$ph3[$i] = "PH3";
				} else {
					$ph3[$i] = undef;
			};
			#read CDA flag
			if ( $line =~ /.(\bCDA)\;/ ) {
					$cda[$i] = "CDA";
				} else {
					$cda[$i] = undef;
			};
			#read LSD flag
			if ( $line =~ /.(\bLSD)\;/ ) {
					$lsd[$i] = "LSD";
				} else {
					$lsd[$i] = undef;
			};
			#read MTP flag
			if ( $line =~ /.(\bMTP)\;/ ) {
					$mtp[$i] = "MTP";
				} else {
					$mtp[$i] = undef;
			};
			#read OM flag
			if ( $line =~ /.(\bOM)\;/ ) {
					$om[$i] = "OM";
				} else {
					$om[$i] = undef;
			};
			#read NOC flag
			if ( $line =~ /.(\bNOC)\;/ ) {
					$noc[$i] = "NOC";
				} else {
					$noc[$i] = undef;
			};
			#read WTD flag
			if ( $line =~ /.(\bWTD)\;/ ) {
					$wtd[$i] = "WTD";
				} else {
					$wtd[$i] = undef;
			};
			#read NOV flag
			if ( $line =~ /.(\bNOV)\;/ ) {
					$nov[$i] = "NOC";
				} else {
					$nov[$i] = undef;
			};
			
			
			#read NC ??????????
			#It looks like this is a filter that maybe isn't used. There aren't any lines with this field.
			
			#read CAF string
			if ( $line =~ /.\bCAF=(\X+)\;/ ) {
					$caf[$i] = $1;
				} else {
					$caf[$i] = undef;
			};
			#read CLNHGVS string
			if ( $line =~ /.(\CLNHGVS=)(\X{1}|\X{2}|\X{3}|\X{4}|\X{5}|\X{6}|\X{7}|\X{8}|\X{9}|\X{10}|\X{11}|\X{12}|\X{13}|\X{14}|\X{15}|\X{16}|\X{17}|\X{18}|\X{19}|\X{20}|\X{21}|\X{22}|\X{23}|\X{24}|\X{25}|\X{26}|\X{27}|\X{28}\X{29}|\X{30}|\X{31}|\X{31}|\X{33}|\X{34}|\X{35}|\X{36}|\X{37}|\X{38}|\X{39}|\X{40}|\X{41}|\X{42}|\X{43}|\X{44}|\X{45}|\X{46}|\X{47}|\X{48}|\X{49}|\X{50}|\X{51}|\X{52}|\X{53}|\X{54}|\X{55}|\X{56}|\X{57}|\X{58}|\X{59}|\X{60}|\X{61}|\X{62}|\X{63}|\X{64}|\X{65}|\X{66}|\X{68}|\X{69}|\X{70}|\X{71}|\X{72}|\X{73}|\X{74}|\X{75}|\X{76}|\X{77}|\X{78}|\X{79}|\X{80}|\X{81}|\X{82}|\X{83}|\X{84}|\X{85}|\X{86}|\X{87}|\X{88}|\X{89}|\X{90}|\X{91}|\X{92}|\X{93}|\X{94}|\X{95}|\X{96}|\X{97}|\X{98}|\X{99}|\X{100}|\X{101}|\X{102}|\X{103}|\X{104}|\X{104}|\X{106}|\X{107}|\X{108}|\X{109}|\X{110}|\X{111}|\X{112}|\X{113}|\X{114}|\X{115}|\X{116})\;/ ) {
					$clnhgvs[$i] = $2;
				} else {
					$clnhgvs[$i] = undef;
			};
			
			#read COMMON string
			if ( $line =~ /.(\COMMON=)(\d+)/ ) {
					$common[$i] = $2;
				} else {
					$common[$i] = undef;
			};
			#read CLNALLE value
			if ( $line =~ /.(\CLNALLE=)(\d+)\;/ ) {
					$clnalle[$i] = $2;
				} else {
					$clnalle[$i] = undef;
			};
			#read CLNORIGIN string
			if ( $line =~ /.(\bCLNORIGIN=)(\d\,\d|\d)/ ) {
					$clnori[$i] = $2;
				} else {
					$clnori[$i] = undef;
			};
			
			
			#read CLNSRCID string
			if ( $line =~ /.(\CLNSRCID=)(\X{1}|\X{2}|\X{3}|\X{4}|\X{5}|\X{6}|\X{7}|\X{8}|\X{9}|\X{10}|\X{11}|\X{12}|\X{13}|\X{14}|\X{15}|\X{16}|\X{17}|\X{18}|\X{19}|\X{20}|\X{21}|\X{22}|\X{23}|\X{24}|\X{25}|\X{26}|\X{27}|\X{28}\X{29}|\X{30}|\X{31}|\X{31}|\X{33}|\X{34}|\X{35}|\X{36}|\X{37}|\X{38}|\X{39}|\X{40}|\X{41}|\X{42}|\X{43}|\X{44}|\X{45}|\X{46}|\X{47}|\X{48}|\X{49}|\X{50}|\X{51}|\X{52}|\X{53}|\X{54}|\X{55}|\X{56}|\X{57}|\X{58}|\X{59}|\X{60}|\X{61}|\X{62}|\X{63}|\X{64}|\X{65}|\X{66}|\X{68}|\X{69}|\X{70}|\X{71}|\X{72}|\X{73}|\X{74}|\X{75}|\X{76}|\X{77}|\X{78}|\X{79}|\X{80}|\X{81}|\X{82}|\X{83}|\X{84}|\X{85}|\X{86}|\X{87}|\X{88}|\X{89}|\X{90}|\X{91}|\X{92}|\X{93}|\X{94}|\X{95}|\X{96}|\X{97}|\X{98}|\X{99}|\X{100}|\X{101}|\X{102}|\X{103}|\X{104}|\X{104}|\X{106}|\X{107}|\X{108}|\X{109}|\X{110}|\X{111}|\X{112}|\X{113}|\X{114}|\X{115}|\X{116})\;/ ) {
					$clnsrcid[$i] = $2;
				} else {
					$clnsrcid[$i] = undef;
			};		
			
			

		
			
			
			#read CLNSRC string
			if ( $line =~ /.(\CLNSRC=)(\w+)\;/ ) {
					$clnsrc[$i] = $2;
				} else {
					$clnsrc[$i] = undef;
			};
			
			#read CLNSIG string
			if ( $line =~ /.(\bCLNSIG=)(\w+)\;/ ) {
					$clnsig[$i] = $2;
				} else {
					$clnsig[$i] = undef;
			};
			
			#read CLNDSDB string
			if ( $line =~ /.(\bCLNDSDB=)(\w+)\;/ ) {
					$clndsdb[$i] = $2;
				} else {
					$clndsdb[$i] = undef;
			};

			
			
			
			#read CLNDSDBID string
			if ( $line =~ /.(\CLNDSDBID=)(\X{1}|\X{2}|\X{3}|\X{4}|\X{5}|\X{6}|\X{7}|\X{8}|\X{9}|\X{10}|\X{11}|\X{12}|\X{13}|\X{14}|\X{15}|\X{16}|\X{17}|\X{18}|\X{19}|\X{20}|\X{21}|\X{22}|\X{23}|\X{24}|\X{25}|\X{26}|\X{27}|\X{28}\X{29}|\X{30}|\X{31}|\X{31}|\X{33}|\X{34}|\X{35}|\X{36}|\X{37}|\X{38}|\X{39}|\X{40}|\X{41}|\X{42}|\X{43}|\X{44}|\X{45}|\X{46}|\X{47}|\X{48}|\X{49}|\X{50}|\X{51}|\X{52}|\X{53}|\X{54}|\X{55}|\X{56}|\X{57}|\X{58}|\X{59}|\X{60}|\X{61}|\X{62}|\X{63}|\X{64}|\X{65}|\X{66}|\X{68}|\X{69}|\X{70}|\X{71}|\X{72}|\X{73}|\X{74}|\X{75}|\X{76}|\X{77}|\X{78}|\X{79}|\X{80}|\X{81}|\X{82}|\X{83}|\X{84}|\X{85}|\X{86}|\X{87}|\X{88}|\X{89}|\X{90}|\X{91}|\X{92}|\X{93}|\X{94}|\X{95}|\X{96}|\X{97}|\X{98}|\X{99}|\X{100}|\X{101}|\X{102}|\X{103}|\X{104}|\X{104}|\X{106}|\X{107}|\X{108}|\X{109}|\X{110}|\X{111}|\X{112}|\X{113}|\X{114}|\X{115}|\X{116})\;/ ) {
					$clndsdbid[$i] = $2;
				} else {
					$clndsdbid[$i] = undef;
			};
			
			
			#read CLNDBN string
			if ( $line =~ /.(\CLNDBN=)(\X{1}|\X{2}|\X{3}|\X{4}|\X{5}|\X{6}|\X{7}|\X{8}|\X{9}|\X{10}|\X{11}|\X{12}|\X{13}|\X{14}|\X{15}|\X{16}|\X{17}|\X{18}|\X{19}|\X{20}|\X{21}|\X{22}|\X{23}|\X{24}|\X{25}|\X{26}|\X{27}|\X{28}|\X{29}|\X{30}|\X{31}|\X{32}|\X{33}|\X{34}|\X{35}|\X{36}|\X{37}|\X{38}|\X{39}|\X{40}|\X{41}|\X{42}|\X{43}|\X{44}|\X{45}|\X{46}|\X{47}|\X{48}|\X{49}|\X{50}|\X{51}|\X{52}|\X{53}|\X{54}|\X{55}|\X{56}|\X{57}|\X{58}|\X{59}|\X{60}|\X{61}|\X{62}|\X{63}|\X{64}|\X{65}|\X{66}|\X{68}|\X{69}|\X{70}|\X{71}|\X{72}|\X{73}|\X{74}|\X{75}|\X{76}|\X{77}|\X{78}|\X{79}|\X{80}|\X{81}|\X{82}|\X{83}|\X{84}|\X{85}|\X{86}|\X{87}|\X{88}|\X{89}|\X{90}|\X{91}|\X{92}|\X{93}|\X{94}|\X{95}|\X{96}|\X{97}|\X{98}|\X{99}|\X{100}|\X{101}|\X{102}|\X{103}|\X{104}|\X{104}|\X{106}|\X{107}|\X{108}|\X{109}|\X{110}|\X{111}|\X{112}|\X{113}|\X{114}|\X{115}|\X{116}|\X{117}|\X{118}|\X{119}|\X{120}|\X{121}|\X{122}|\X{123}|\X{124}|\X{125}|\X{126}|\X{127}|\X{128}|\X{129}|\X{130}|\X{131}|\X{132}|\X{133}|\X{134}|\X{135}|\X{136}|\X{137}|\X{138}|\X{139}|\X{140}|\X{141}|\X{142}|\X{143}|\X{144}|\X{145}|\X{146}|\X{147}|\X{148}|\X{149}|\X{150}|\X{151}|\X{152}|\X{153}|\X{154}|\X{155}|\X{156}|\X{157}|\X{158}|\X{159}|\X{160}|\X{161}|\X{162}|\X{163}|\X{164}|\X{165}|\X{166}|\X{167}|\X{168}|\X{169}|\X{170}|\X{171}|\X{172}|\X{173}|\X{174}|\X{175}|\X{176}|\X{177}|\X{178}|\X{179}|\X{180}|\X{181}|\X{182}|\X{183}|\X{184}|\X{185}|\X{186}|\X{187}|\X{188}|\X{189}|\X{190}|\X{191}|\X{192}|\X{193}|\X{194}|\X{195}|\X{196}|\X{197}|\X{198}|\X{199}|\X{190}|\X{200})\;/ ) {
					$clndbn[$i] = $2;
				} else {
					$clndbn[$i] = undef;
			};
			
			
			#read CLNACC string

			if ( $line =~ /.(\bCLNACC=)(\X+)/  ) {
					my $tempclnacc = substr($2, 0, index($2, ';'));
					$clnacc[$i] = $tempclnacc;
					chomp( $clnacc[$i] );
				} else {
					$clnacc[$i] = undef;
			};
			
		

			
			
			
			
		#set count to new line
		$i = $i + 1;
	# end loop
	
	
	};
#close input file
	close($readfile)
	#if error, print warning	
	or die "Could not close file '$readfile' $!";
		#if error, print warning


		

		

		
		
	#write new output file
	#read name from user - just list in command line?
	my $outputfile;
	$outputfile = "STDOUT";
	#check file write
	#report error
		open(my $op, '>>', $outputfile)
			or die "Could not write output file '$outputfile' $!";
		#can I append to a file i haven't written yet?


#write header lines to new file

	print "chromosome\tposition\trs_accession_number\tconsensus_nucleotide\tvariation_nucleotide\trs_number\trspos\tdbsnp_build\tssr\tsao\tvp\tgeneinfo\twgt\tvc\tpm\ttpa\tpmc\ts3d\tslo\tnsf\tnsm\tnsn\tref\tsyn\tu3\tu5\tass\tdss\tint\tr3\tr5\toth\tcfl\tasp\tmut\tvld\tg5a\tg5\thd\tgno\tkgvalidated\tkgphase1\tkgpilot123\tkgprod\totherkg\tph3\tcda\tlsd\tmtp\tom\tnoc\twtd\tnov\tnc\tcaf\tcommon\tclnhgvs\tclnalle\tclnsrc\tclnorigin\tclnsrcid\tclnsig\tclndsdbid\tclndbn\tclnacc\n"
		or die "Could not write to output file '$op' $!";

#write hash to output file
$i = 0;
foreach $i (0 .. $#clnacc) {
	print  "$chrom[$i] \t $position[$i] \t $rscombo[$i] \t  $consntp[$i] \t $varntp[$i] \t $rs[$i] \t $rspos[$i] \t $dbsnpbuild[$i] \t $ssr[$i] \t $sao[$i] \t $vp[$i] \t $geneinfocombo[$i] \t $wgt[$i] \t $vc[$i] \t $pm[$i] \t $tpa[$i] \t $pmc[$i] \t $s3d[$i] \t $slo[$i] \t $nsf[$i] \t $nsm[$i] \t $nsn[$i] \t $ref[$i] \t $syn[$i] \t $u3[$i] \t $u5[$i] \t $ass[$i] \t $dss[$i] \t $int[$i] \t $r3[$i] \t $r5[$i] \t $oth[$i] \t $cfl[$i] \t $asp[$i] \t $mut[$i] \t $vld[$i] \t $g5a[$i] \t $g5[$i] \t $hd[$i] \t $gno[$i] \t $kgval[$i] \t $kgpha[$i] \t $kgpil[$i] \t $kgprod[$i] \t $otherkg[$i] \t $ph3[$i] \t $cda[$i] \t $lsd[$i] \t $mtp[$i] \t $om[$i] \t $noc[$i] \t $wtd[$i] \t $nov[$i] \t $nc[$i] \t $caf[$i] \t $common[$i] \t $clnhgvs[$i] \t $clnalle[$i] \t $clnsrc[$i] \t  $clnori[$i] \t $clnsrcid[$i]  \t $clnsig[$i] \t $clndsdbid[$i] \t $clndbn[$i] \t $clnacc[$i]\n";
	$i = $i + 1;
	}
	
	
	
#close write file
		close($op)
			#if error, print warning
			or die "Could not close output file '$outputfile' $!";


