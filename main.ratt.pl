#! /usr/bin/perl -w
#
# File: annotation.correctString.pl
# Time-stamp: <25-Mar-2010 15:52:22 tdo>
# $Id: $
#
# Copyright (C) 2010 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto tdo@sanger.ac.uk and Gary Dilon
#
# Description: Please see http://ratt.sourceforge.net for information
#

use strict;
use Data::Dumper;
use lib $ENV{RATT_HOME};
use ratt_correction;

my $debug=0;

my $SET=1;
my $COLOR_BAD=4;

if (!defined($ENV{RATT_HOME})) {
  print "Please set global variable RATT_HOME in your shell!\n";
  exit 1
}

if (!defined($ARGV[0])) {
  print "Sorry, wrong option.\n";
  print "Tranfer / Correct / Check / EMBLFormatCheck / Mutate / Split / Difference / Embl2Fasta\ncan be used.\n\n";
  exit 1
}
elsif ($ARGV[0] eq "Mutate") {
  if (! defined($ARGV[1])) {
	print "\n\nusage: \$RATT_HOME/main.ratt.pl Mutate <(multi-)fasta-file>\n\n".
	  "Every 250 base pairs a base is changed (mutated). The result is saved as <fastafile>.mutated. This is necessary to recalibrate RATT for similar genomes.\n\n";
	
	exit;
	
  }
  putMutation($ARGV[1]);
  
  exit;
}
elsif ($ARGV[0] eq "Split") {
  if (! defined($ARGV[1])) {
	print "\n\nusage: \$RATT_HOME/main.ratt.pl Split <(multifasta-file>\n\n".
	  "Splits a given multifasta file into individual files containing one sequence. This is necessary as visualization tools (e.g. Artemis) prefer single fasta files.\n\n";
	
	exit (1);
	
  }
  Split($ARGV[1]);
  
  exit;
}
elsif ($ARGV[0] eq "Difference") {
  my $mummerSNP=$ARGV[1];
  my $mummerCoords=$ARGV[2];
  my $resultName=$ARGV[3];
  
  if (!defined($resultName)) {
  	print "\n\nusage: \$RATT_HOME/main.ratt.pl Difference <mummer SNP file> <mummer coord file> <ResultName>\n\n".
	  "Generates files that report the SNP, indels and regions not shared by the reference and query. It also prints a statistic reporting coverage for each replicon.\n\n";
  	exit (1);
  }
  my $stats=getMutations($mummerSNP,$mummerCoords,$resultName);
  
  print $stats;
  exit;
  
}
elsif ($ARGV[0] eq "EMBLFormatCheck") {
  if (scalar(@ARGV) < 3) {
 	print "\n\nusage: \$RATT_HOME/main.ratt.pl EMBLFormatCheck <EMBL file> <ResultName postfix>\n\n".
	  "Some EMBL files have feature positions spanning several lines, this function consolidates these features so they appear on one line. The result name is <EMBL File>.<ResultName postfix>.\n\n";
	
  	exit (1); 
  }
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  my $resultName = shift;
  
  correctEMBL($embl,"tmp.BBA.embl");
}
elsif ($ARGV[0] eq "Correct") {
  if (scalar(@ARGV) < 4) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Correct <EMBL file> <fasta file> <ResultName>\n\n".
	  "Corrects a given annotation, as described previously. The corrections are reported and the new file is saved as <ResultName>.embl.\n\n";	
	exit(1); 
  }
  
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  my $resultName = shift;
  
  correctEMBL($embl,"tmp.BBA.embl");
  
  startAnnotationCorrection( "$embl.tmp.BBA.embl",$fasta,$resultName);
  
  exit;
}
elsif ($ARGV[0] eq "Check") {
  if (scalar(@ARGV) < 4) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Check <EMBL file> <fasta file> <ResultName>\n\n".
	  "Similar to the correct option, but it will only report errors in an EMBL file.\n\n";
 my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  my $resultName = shift;	
	correctEMBL($embl,"tmp.BBA.embl");
	
	startAnnotationCheck( "$embl.tmp.BBA.embl",$fasta,$resultName);
	exit;
  }
}
elsif ($ARGV[0] eq "Embl2Fasta") {
  if (scalar(@ARGV) < 3) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Embl2Fasta <EMBL dir> <fasta file>\n\n".
	  "Extracts the sequence from embl files in the <EMBL directory> and saves it as a <fasta file>.\n\n";
  exit 1;
  }
  
  
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  Embl2Fasta($embl,$fasta);
  exit;
  
}

elsif ($ARGV[0] eq "EMBLFormatCheck") {
  if (scalar(@ARGV) < 3) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl EMBLFormatCheck <EMBL file> <ResultName postfix>\n\n".
	  "Some EMBL files have feature positions spanning several lines, this function consolidates these features so they appear on one line. The result name is <EMBL File>.<ResultName postfix>.\n\n";
	
 	exit(1); 
  }
  
  my $what =shift;
  my $embl=shift;
  my $postfix = shift;
  
  correctEMBL($embl,$postfix);
   
  exit;
}
elsif ($ARGV[0] eq "Transfer") {
  
  if (@ARGV< 5) {
 	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Transfer <embl Directory> <mummer SNP file> <mummer coord file> <ResultName>\n\n".
	  "This functionality uses the mummer output to map the annotation from embl files, which are in the <embl Directory>, to the query. It generates all the new annotation files (ResultName.replicon.embl), as well as files describing which annotations remain untransferred (Replicon_reference.NOTtransfered.embl).\n\n";
	
	
 	exit(1);
  }
  
  
  my $what         = shift;
  my $emblDir     = shift;
  my $mummerSNP   = shift;
  my $mummerCoords= shift;
  my $resultName  = shift;
  
  my $dbg=100;
  

## main hash: %ref_shift{Ref_contig}[pos] [0] query_contig
#										  [1] position
#										  [2] strand
my $ref_shift;


#load the position of the 
my $ref_cdsPosition=loadEmbl($emblDir);

# fill the shift hash with the coords 
$ref_shift = loadCoords($mummerCoords,$ref_shift);

#print Dumper $ref_shift;
# tune the shift hash with the snp file
$ref_shift = loadSNP($mummerSNP,$ref_shift,$ref_cdsPosition);

# clean the space of the annotation position
undef($ref_cdsPosition);
#print Dumper $ref_shift;


# transfer the annotation the annotation
opendir (DIR, $emblDir) or die "Problem to open opendir $emblDir: $!\n";

### will hold the new annotation: $h{queryname}.=annotation as embl
my $ref_results;
my $ref_Counting= {'Partial'         => 0,
				   'ExonNotTransfered' => 0,
				   'Split'         => 0,
				   'NotTransfered' => 0,
				   'Transfered'    => 0,
				   'CDS' => 0,
				   'CDSTransfered' => 0,
				   'CDSNotExons'   => 0,
				   'CDSPartial'   => 0
				  };

map {
	if (/embl$/){
		($ref_results,$ref_Counting)=adaptAnnotationEMBL($emblDir,$_,$ref_shift,$ref_results,$ref_Counting);
	
	}
} readdir(DIR);




  
  ### output results
  print "Overview of transfere of annotation elements:\n$$ref_Counting{Elements}\telements found.\n";
  print "$$ref_Counting{Transfered}\tElements were transfered.\n";
  print "$$ref_Counting{Partial}\tElements could be transfered partially.\n";
  print "$$ref_Counting{Split}\tElements splitted.\n";
  print "$$ref_Counting{ExonNotTransfered}\tParts of elements (i.e.exons tRNA) not transferred.\n";
  print "$$ref_Counting{NotTransfered}\tElements couldn't be transferred.\n";
  
  print "\nCDS:\n$$ref_Counting{CDS}\tGene models to transfer.\n$$ref_Counting{CDSTransfered}\tGene models transferred correctly.\n";
  print "$$ref_Counting{CDSPartial}\tGene models partially transferred.\n";
  print "$$ref_Counting{CDSNotExons}\tExons not transferred from partial CDS matches.\n";
  print ($$ref_Counting{CDS}-$$ref_Counting{CDSTransfered}-$$ref_Counting{CDSPartial});
  print "\tGene models not transferred.\n\n";
  
  
  
  ### then just save it
saveAnnotation($resultName,$ref_results);

}
else {
  print "Sorry, wrong option.\n";
  print "Tranfer / Correct / Check / EMBLFormatCheck / Mutate / Split / Difference / Embl2Fasta\ncan be used.\n\n";
}




###########################################
### subs

############################################
### loadEMBL
############################################
sub loadEmbl{
  my ($emblDir) = @_;
  
  opendir (DIR, $emblDir) or die "Problem to open opendir $emblDir: $!\n";

### will hold the new annotation: $h{queryname}.=annotation as embl
my $ref_annotation;

map {
	if (/embl$/){
		my $embl=$_;
 		open(F, $emblDir."/".$embl) or die "Problems open embl $embl: $!\n";
 
		my ($chr) = $embl =~ /\/{0,1}(\S+)\.embl$/;
		while (<F>){
 			if (/^FT.*CDS.*\d+/){
 				if ($_ =~ /^\W*\(*(\d+)\.\.(\d+)\)*$/) {
					$ref_annotation=
						Mask_CDS($chr,$ref_annotation,$1,$2);
  				}
	  			elsif (/\.\./){
					my @a=split(/,/);
					foreach (@a) {
					  if (/(\d+)\.\.(\d+)/){
						$ref_annotation=
							Mask_CDS($chr,$ref_annotation,$1,$2);
					  }
					  
					}
  				}	
 			}
 		}
	}
	} readdir(DIR);
	return ($ref_annotation)
}

############################################
### Mask_CDS
############################################
sub Mask_CDS{
  my ($chr,$ref,$f,$l,$out) = @_;

  for ($f..$l){
	$$ref{$chr}{$_}=1;
  }
  
  return $ref;
}

############################################
### adaptAnnotationEMBL
############################################
sub adaptAnnotationEMBL{
  my ($DIR,$embl,$ref_shift,$ref_results,$ref_Counting) = @_;
  
  open(F, $DIR."/".$embl) or die "Problems open embl $embl: $!\n";
  
  my $res;
  my $resDeleted='';
  my @ar;
  my $in=1;
  my $maybeWrong=0;

	### get the name of the contig/supercontig/chromosome
  my ($chr) = $embl =~ /\/{0,1}(\S+)\.embl$/;

  
  #  print Dumper %core;
  my $OKCore=1;
  my $queryTarget;
  my $transfer=0;
  
  while (<F>) {
	
	# UTR must be saved
	s/3\'UTR/3TUTR/g;
	s/5\'UTR/5TUTR/g;
	if (/FT   \S+/) {
	  s/<//g;
	  s/>//g;
	}

	my $line=$_;
	# check if entry is over more than one line
	while ($line =~ /^FT   \S+\s{2,}.*\d+,$/) {
	  $_=<F>;
	  chomp($line);
	  	  
	  /^FT   \s{2,}(.*)$/;
	  $line.=$1;
	}

	
	if ($line =~ /^FT   \S+\s{2,}\D+(\d+)\..*\.(\d+)/ ||
		$line =~ /^FT   \S+\s{2,}\D+\d+,(\d+)\..*\.(\d+)/ ||
		$line =~ /^FT   \S+\s{2,}\D+(\d+)/
	   ) {
	  ### This is necessary to not mapped things, which are not covered

	  my $posA=$1;
	  my $posE=$2;
	  if (!defined($posE)){
	  	$posE=$posA	
	  }
	  ### check if CDS
	  if ($line =~ /FT   CDS/) {
		$$ref_Counting{CDS}++
	  }
	  
	  
		  chomp;
	  ($ref_results,$queryTarget,$ref_Counting,$transfer)=doTransfer($ref_shift,$ref_results,$chr,$posA,$line,$ref_Counting);
	  
	  $$ref_Counting{Elements}++;
	  ### case 1, all ok
	  if (defined($$ref_shift{$chr}[$posA][0]) &&
		  defined($$ref_shift{$chr}[$posE][0]) &&
		  ($$ref_shift{$chr}[$posE][0] eq $$ref_shift{$chr}[$posA][0])  # transfer to same query
		 )	
		{

	 	}
	  ### case 2, defined, but gene model is split between two query
	  elsif (defined($$ref_shift{$chr}[$posA][0]) &&
			 defined($$ref_shift{$chr}[$posE][0])
			)
		 {
		   $$ref_Counting{"Split"}++;
		   
#		   print "$chr  position $posA $posE \n $$ref_shift{$chr}[$posA][0]  // ($$ref_shift{$chr}[$posE][0] \n";
		   
		 }  
	  elsif (defined($$ref_shift{$chr}[$posA][0]) ||
			 defined($$ref_shift{$chr}[$posE][0])
			)
		 {
#		   $$ref_Counting{"Partial"}++
			 
		 }

	}
	elsif (/^SQ/) {
	  last;
	} 
	elsif ($transfer==1){
	  $$ref_results[0]{$queryTarget} .= $_;
	}
	elsif ($transfer==0){
	  $$ref_results[1]{$chr} .= $_;
	}
	# this is the case when the annotation could be just mapped partially.
	elsif ($transfer==3){
	  $$ref_results[1]{$chr} .= $_;
	  $$ref_results[0]{$queryTarget} .= $_;
	}
	
  } # end while <F>
  
  close(F);
  return ($ref_results,$ref_Counting);
}

############################################
### doTransfer
############################################
sub doTransfer{
	my $ref_shift=shift;
	my $ref_resultsLocal =shift;
	my $chr =shift;
	my $pos =shift;
	my $line=shift;
    my $ref_Counting = shift;
	
	my $chrqry=0;
	# the zero will hold the no puttable
	my %ResultLine;
	
	### put complement to it, or get rid of it,
	### also will need to reorder the numbers
	my $wasComplement=0;
	if ($line =~ /complement/){
	  $wasComplement=1;
	  $line =~ s/complement\(//g;
	  $line =~ s/\)$//g;
	}
	$line =~ s/join\(//g;
    $line =~ s/\)$//g;

	
	### here to look for missed exons
	
	# 	$line =~ s/(\d+)/($$ref_shift{$chr}[$1][1])/ge;
	my (@parts) = split(/\s+/,$line); 
	my @ar=split(/,/,$parts[2]);
	
	my $mappedOnce=0;
	my $exonMissed=0;
	my $oldQuery;
	my $partialCount=0;
	
	for (my $i=0;$i < scalar(@ar);$i++) {
	  
	  #single base exon
	  if (! ($ar[$i] =~ /\.\./) ) {
		$ar[$i] =~ /(\d+)/;
		my $pos=$1;
#		print " single exon $ar[$i] pos is $pos \n";
#		print Dumper $$ref_shift{$chr};
		
		if (defined($$ref_shift{$chr}[$pos][0])) {
		  $ar[$i] =~ s/(\d+)/($$ref_shift{$chr}[$1][1])/ge;
		  $mappedOnce++;

		  $oldQuery=$$ref_shift{$chr}[$pos][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $pos;
		  
		}		
		else {
		  $ResultLine{0}[0] .= "$ar[$i],";
		  $exonMissed++;
		  
		}
	  }
	  else {
		$ar[$i] =~ /(\d+)\.\.(\d+)/;
		my $posA=$1;
		my $posE=$2;
		if (defined($$ref_shift{$chr}[$posA][0]) &&
			defined($$ref_shift{$chr}[$posE][0]) &&
			$$ref_shift{$chr}[$posE][0] eq $$ref_shift{$chr}[$posA][0] &&
			abs($$ref_shift{$chr}[$posE][1] - $$ref_shift{$chr}[$posA][1])< (2*($posE-$posA))
		   ) {
		  $ar[$i] =~ s/(\d+)/($$ref_shift{$chr}[$1][1])/ge;
		  $mappedOnce++;
		  
		  $oldQuery=$$ref_shift{$chr}[$posA][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $pos;
		  
		}
		elsif (defined($$ref_shift{$chr}[$posA][0]) &&
			defined($$ref_shift{$chr}[($posA+74)][0]) &&
			$$ref_shift{$chr}[($posA+74)][0] eq $$ref_shift{$chr}[$posA][0] &&
			abs($$ref_shift{$chr}[($posA+74)][1] - $$ref_shift{$chr}[$posA][1])< 20000
		   ) {
		  $ar[$i] = $$ref_shift{$chr}[$posA][1]."..".$$ref_shift{$chr}[($posA+74)][1];
		  $mappedOnce++;
		  $partialCount++;
		  
		  $oldQuery=$$ref_shift{$chr}[$posA][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posA;
		  
		}
		elsif (defined($$ref_shift{$chr}[$posE][0]) &&
			defined($$ref_shift{$chr}[($posE-74)][0]) &&
			$$ref_shift{$chr}[($posE-74)][0] eq $$ref_shift{$chr}[$posE][0] &&
			abs($$ref_shift{$chr}[($posE-74)][1] - $$ref_shift{$chr}[$posE][1])< 20000
		   ) {
		  $ar[$i] = $$ref_shift{$chr}[($posE-74)][1]."..".$$ref_shift{$chr}[($posE)][1];
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posE][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posE;
		  
		}
		else {
		  $ResultLine{0}[0] .= "$ar[$i],";
		  $exonMissed++;
		}
	  }
	  
	  
	}

	## default do not transfer
	my $transfer=0;
	#### check the amount
	if ($mappedOnce ==0) {
	  $$ref_Counting{NotTransfered}++
	}
	if (($mappedOnce > 0 && $exonMissed >0)
#		||
#		$partialCount>0
	   ){
	  $$ref_Counting{Partial}++;
	  ### This means, put the annotation to the BB and LB, as it is partial
	  if ($line =~ /FT   CDS/) {
		$$ref_Counting{CDSPartial}++;
		$$ref_Counting{CDSNotExons}+=$exonMissed;
	  }
	  $transfer=3;
	}


	
	if ($exonMissed==0) {
	  $$ref_Counting{Transfered}++;
	  if ($line =~ /FT   CDS/) {
		$$ref_Counting{CDSTransfered}++;
	  }
	  
	  ### Modell fully mapped, put it just to the LB
	  $transfer=1;
	}
	else {
	  $$ref_Counting{ExonNotTransfered}+=$exonMissed;
	  $ResultLine{0}[0]=~ s/,$//g;
	  if ($wasComplement==1){
		$ResultLine{0}[0]="complement(join($ResultLine{0}[0]))";
	  } 
	  else {
	  	$ResultLine{0}[0]="join($ResultLine{0}[0])";
	  }
	  $$ref_resultsLocal[1]{$chr}.=sprintf("%-4s %-15s %s\n",$parts[0],$parts[1],$ResultLine{0}[0]);
	  $ResultLine{0}=undef;
	  undef $ResultLine{0};
	}
	
	
	
	foreach my $trans (keys %ResultLine) {
	
		
	  if ($trans ne '0'){
		
		my ($chrqryLocal,$chrpart) = $trans =~ /^(\S+)::(\S+)$/;
		my $pos =$ResultLine{$trans}[1];
		if (!defined($chrqry)){
		  print $trans."\n";
		  print Dumper %ResultLine;
		  exit;	
		}	  
		
		
		$ResultLine{$trans}[0]=~ s/,$//g;
	  	my @ar=split(/\,/,$ResultLine{$trans}[0]);
	  	my $amountJoined=(scalar(@ar)-1);
	  	if ($amountJoined > 1){
			$ResultLine{$trans}[0]="join($ResultLine{$trans}[0])";
		  }
		  if ($wasComplement==1){
			$ResultLine{$trans}[0]="complement($ResultLine{$trans}[0])";
	  	} 
	  
	### check if the entry must be inversed
 	
	  $chrqry=$chrqryLocal;

	  if (defined($$ref_shift{$chrpart}[$pos][2]) && $$ref_shift{$chrpart}[$pos][2] == -1){
		$line=$ResultLine{$trans}[0];
		
		if ($line =~ /complement/){
		  $wasComplement=1;
		  $line =~ s/complement\(//g;
		  $line =~ s/\)$//g;
		}
		
		my $hadExons=0;
		if ($line =~ /join/){
		  $hadExons=1;
		  $line =~ s/join\(//g;
		}		
		$line =~ s/\)$//g;
		my @numbers = sort {$a <=> $b} split(/,|\.{2}/,$line);
		my $new;
		my $count=0;
		foreach (@numbers){
		  if ($count>0 && ($count%2==1)){
			$new.="..";
		  } elsif ($count>0 && ($count%2==0)){
			$new.=",";
		  }
		  $new.=$_;
		  $count++;
		}
		if ($hadExons){
		  $new="join($new)"	
		}
		if ($wasComplement==0){
		  $new="complement($new)"
		}
		$ResultLine{$trans}[0]=$new;
	  }

	  $$ref_resultsLocal[0]{$chrqry}.=sprintf("%-4s %-15s %s\n",$parts[0],$parts[1],$ResultLine{$trans}[0]);
	}
	}
 	return ($ref_resultsLocal,$chrqry,$ref_Counting,$transfer);
}
############################################
### saveAnnotation
############################################
sub saveAnnotation{
  my ($name,$ref_h) = @_;


  ### map the mapping Stuff
  foreach my $query (sort keys %{$$ref_h[0]}){
	my $nameQry=$query;
	$nameQry =~ s/\|/_/g;
	open (F,"> $name.$nameQry.embl") or die "Couldn't open save file $name: $!\n";
	
	# UTR must be saved
	$$ref_h[0]{$query} =~ s/3TUTR/3\'UTR/g;
	$$ref_h[0]{$query} =~ s/5TUTR/5\'UTR/g;
	print F $$ref_h[0]{$query};
	close(F);
  }
  
  # UTR must be saved
  foreach my $ref (sort keys %{$$ref_h[1]}){
	my $nameRef=$ref;
	$nameRef =~ s/\|/_/g;
	
 	open (F,"> $name.$nameRef.NOTTransfered.embl") or die "Couldn't open save file $name: $!\n";
	
	$$ref_h[1]{$ref} =~ s/3TUTR/3\'UTR/g;
	$$ref_h[1]{$ref} =~ s/5TUTR/5\'UTR/g;
	print F $$ref_h[1]{$ref};
	close(F);
  }

}
############################################
### saveGFF
############################################

sub saveGFF{
  my ($name,$path,$ref_h) = @_;

  if (! -d "$path") {
	!system("mkdir $path") or die "Couldn't create directory $path.\n";
  }

  foreach my $chr (sort keys %$ref_h){
	my $chrName=$chr;
	$chrName =~ s/\|/_/g;
	
	open (F,"> $path/$name.$chrName.Mutations.gff") or die "Couldn't create file $name: $!\n";

  # UTR must be saved
  print F $$ref_h{$chr};
  close(F);
  }
}

############################################
### loadSNP
############################################
sub loadSNP{
  my $fileName        = shift;
  my $ref_shift       = shift;
  my $ref_cdsPosition = shift;
  
  open (F, $fileName) or die "Problem to open SNP File: $fileName \n";
  
  my @File=<F>;
  close(F);

  ## walk through the list. The last 
  for my $pos (0..(scalar(@File)-2)) {

	# get the positions of the snp'indels of before and last
	my @previous;
	if ($pos >0) { @previous = split(/\s+/,$File[($pos-1)]);}
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$refStrand,$queryStrand,$reference,$query) = split(/\s+/,$File[$pos]);
	my @next = split(/\s+/,$File[($pos+1)]);
	
#	my ($refPreviousPos,$refPrevious,$queryPrevious) 
#		= ($previous[0],$previous[3],$previous[8],$previous[9]);
	my ($refNextPos,$queryNextPos,$refNext,$queryNext) 
		= ($next[0],$next[3],$next[8],$next[9]);
	
	# TODO 1: differ is next is not in the same combination	
	# TODO 2: what to do with the last SNP?
	# TODO 3: must be the same contig combination
	

	
	if ($refNext ne $reference) {
	  $ref_shift=walkToEnd($ref_shift,$reference,$refPos,$query,$queryPos,$queryStrand);
	}
	
	#case 1: strand 1 and  update due to annotation
	if (
		($refNextPos - $refPos) < 50000 &&
		!defined($$ref_cdsPosition{$reference}{$refPos})    && # this mutation is not on gene
		defined($$ref_cdsPosition{$refNext}{$refNextPos})   # next mutation is on gene
		
      	){
		for  (my $posLocal=$refNextPos; $posLocal >=($refPos);$posLocal--){
		  if (defined($$ref_shift{$reference}[$posLocal][0])) {
			$$ref_shift{$reference}[$posLocal][1]=$queryNextPos;
		  } 
		  $queryNextPos-=$queryStrand
		}	
	}
	else {
		for my $posLocal ($refPos..($refNextPos-1)){
		  if (defined($$ref_shift{$reference}[$posLocal][0])) {
			$$ref_shift{$reference}[$posLocal][1]=$queryPos;
		  }
		  
		  $queryPos+=$queryStrand
		}	
	}  
  }

  #now have a look a the last line:
  my ($refPos,$refWhat,$queryWhat,$queryPos,
	  $dum1,$dum2,$refStrand,$queryStrand,
	  $reference,$query)
	= split(/\s+/,$File[(scalar(@File)-1)]);

  $ref_shift=walkToEnd($ref_shift,$reference,$refPos,$query,$queryPos,$queryStrand);
  return $ref_shift;	
}

sub  walkToEnd{
  my ($ref_shift,$reference,$refPos,$query,$queryPos,$queryStrand) = @_;
  while (defined($$ref_shift{$reference}[$refPos])) {
	$$ref_shift{$reference}[$refPos][1]=$queryPos;
	$queryPos+=$queryStrand;
	$refPos++
  }
  
  return $ref_shift
}


############################################
### loadCoords
############################################

sub loadCoords{
  my $fileName = shift;
  my $ref_h    = shift;
  
  open (F, $fileName) or die "Problem to open coords File: $fileName \n";
  
  my @File=<F>;
  
  ### position of blocks
 foreach (@File) {
	# 1       115285  837     116121  115285  115285  99.99   5.36    5.31    Neo_chrII       Neo_chrII

	my ($refPos1,$refPos2,$queryPos1,$queryPos2,$overlap1,$overlap2,$identity,$dum1,$dum2,$refLength,$queryLength,$reference,$query) = split(/\s+/);

	### if the alignment is inverted...

	my $maxQuery=$queryPos2; ### as the alignment length might not be
                             ### the same, the querypos cannot be
                             ### bigger than the $queryPos2
	my $minQuery=$queryPos1;

	
	my $strand=1;
	if ($queryPos1 > $queryPos2) {
		$strand=-1;
		$minQuery=$queryPos2;
		$maxQuery=$queryPos1;
	}
	for my $pos ($refPos1..$refPos2){
	  if ($queryPos1< $minQuery ||
		  $queryPos1 > $maxQuery
		 ) {
		$queryPos1-=$strand;
		
	  }
	  
	  @{ $$ref_h{$reference}[$pos]} = ($query,$queryPos1,$strand);
	  $queryPos1+=$strand;
	}
	
  }
  return ($ref_h);
}


############################################
### putMutation
############################################
## Changes every 250 bp the base to G or C
sub putMutation{
  my $file = shift;

  open (F,$file) or die "Please provide a valid query sequence\n";

  my $MUTATION_RATE=250;
  
  my %h;
  my $name;
  while (<F>) {
	chomp;
	if (/^>(\S+)/) {
	  $name=$1;
	}
	else {
	  $h{$name}.=$_;
	}
  }
  close(F);

  # insert the mutation

  ## per chromosome;
  for my $chr (keys %h) {
	my @seq = split "", $h{$chr};	
	my $length=(scalar(@seq)-1);
	
	# mutate the every $MUTATION_RATE base
	
	for (my $i = $MUTATION_RATE; $i < $length ; $i+=$MUTATION_RATE){
	  if ($seq[$i] ne 'G') {
		$seq[$i]='G'
	  }
	  else {
		$seq[$i]='C';
	  }
	}
	### mutate the second and the seconlast
	if ($seq[1] ne 'G') {
	  $seq[1]='G'
	}
	else {
	  $seq[1]='C';
	}
	$length--;
	$length--;
	
	if ($seq[$length] ne 'G') {
	  $seq[$length]='G'
	}
	else {
	  $seq[$length]='C';
	}
	
	# put the sequence back together
	$h{$chr}= join("",@seq)
  }

  
  ### write the result
  open (F, ">$file.mutated") or die "Problems to write the file $file.mutated\n";
   for my $chr (keys %h) {
	 print F ">$chr\n$h{$chr}\n";
   }
  close(F)
}

###########################################
### Functions for gff files
############################################

sub getMutations{
  my $fileNameSNP     = shift;
  my $fileNameCoords  = shift;
  my $resultName      = shift;
  
  open (F, $fileNameSNP) or die "Problem to open SNP File: $fileNameSNP \n";
  
  my @File=<F>;
  close(F);

  my %BB;
  my %LB;

  my (%sizeRef,%sizeQuery,%coveredRef,%coveredQuery,%noCovRef,%noCovQry);
  my %h_sizeRef;
  my %h_sizeQuery;
  
  foreach (@File) {
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$refStrand,$queryStrand,$reference,$query) = split(/\s+/);

	if ($refWhat eq ".") {
	  $BB{$reference} .="unknown\tBBA\tIns\t$refPos\t$refPos\t0\t+\t.\tnote=\"Insertion+in+query:+$queryWhat++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="unknown\tBBA\tDel\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"Deletion+in+reference++++Strand+of+query+is+$queryStrand\"\n";
	}
	elsif($queryWhat eq "."){
	  $BB{$reference} .="unknown\tBBA\tDel\t$refPos\t$refPos\t0\t+\t.\tnote=\"Deletion+in+query++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="unknown\tBBA\tIns\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"Insertion+in+reference:+$refWhat++++Strand+of+query+is+$queryStrand\"\n";
	}
	
	else {
	  $BB{$reference} .="unknown\tBBA\tSNP\t$refPos\t$refPos\t0\t+\t.\tnote=\"SNP+in+query:+$queryWhat++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="unknown\tBBA\tSNP\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"SNP+in+reference:+$refWhat++++Strand+of+query+is+$queryStrand\"\n";
	}
  }

  ## from the coords files we want the regions that are not
  ## covered\n"; 
  
  open (F, $fileNameCoords) or die "Problem to open SNP File: $fileNameCoords \n";
  
  @File=<F>;
  close(F);

  my %coverBB;
  my %coverLB;

  ### get an field, where there is coverage
  foreach (@File) {
	my ($refPos1,$refPos2,$queryPos1,$queryPos2,$overlap1,$overlap2,$identity,$refLength,$queryLength,$dum1,$dum2,$reference,$query) = split(/\s+/);
	$sizeRef{$reference}=$refLength;
	$sizeQuery{$query}=$queryLength;

	$coveredRef{$reference}+=($refPos2-$refPos1+1);
	
	# check on the reference
	$coverBB{$reference}[$refLength]=undef;
	for ($refPos1..$refPos2) {
	  $coverBB{$reference}[$_]=1;
	}

	# check for the query
	$coverLB{$query}[$queryLength]=undef;
	if ($queryPos2<$queryPos1) {  my $tmp=$queryPos2;$queryPos1=$queryPos2;$queryPos2=$tmp	}
	$coveredQuery{$query}+=(abs($queryPos2-$queryPos1)+1);
	for ($queryPos1..$queryPos2) {
	  $coverLB{$query}[$_]=1;
	}
  }
  ### now parse, were there is no coverage
  for my $chr (keys %coverBB) {
	my $start=1;
	my $ok=1;

	for my $pos (1..(scalar(@{$coverBB{$chr}}))) {
	  if ($ok==1 &&
		  !defined($coverBB{$chr}[$pos])) {
		$ok=0;
		$start=$pos
	  }
	  if ($ok==0 &&
		  defined($coverBB{$chr}[$pos])) {
		$ok=1;
		$BB{$chr} .="unknown\tBBA\tSynteny\t$start\t".($pos-1)."\t0\t+\t.\tnote=\"No synteny with query. Possible insert or too divergent\"\n";
		$noCovRef{$chr}+=($pos-1-$start);
		
	  }
	}
	if ($ok==0 && $start <scalar(@{$coverBB{$chr}}) ) {
	  $BB{$chr} .="unknown\tBBA\tSynteny\t$start\t".(scalar(@{$coverBB{$chr}})-1)."\t0\t+\t.\tnote=\"No synteny with query. Possible insert or too divergent\"\n";
	  $noCovRef{$chr}+=($sizeRef{$chr}-1-$start);
	}
  }
  ### now parse, were there is no coverage
  for my $chr (keys %coverLB) {
	my $start=1;
	my $ok=1;

	for my $pos (1..(scalar(@{$coverLB{$chr}})-1)) {
	  if ($ok==1 &&
		  !defined($coverLB{$chr}[$pos])) {
		$ok=0;
		$start=$pos
	  }
	  if ($ok==0 &&
		  defined($coverLB{$chr}[$pos])) {
		$ok=1;
		$LB{$chr} .="unknown\tBBA\tSynteny\t$start\t".($pos-1)."\t0\t+\t.\tnote=\"No synteny with reference. Possible insert or too divergent\"\n";
		$noCovQry{$chr}+=($pos-1-$start);
	  }
	}
	if ($ok==0 && $start < scalar(@{$coverLB{$chr}})) {
	  $LB{$chr} .="unknown\tBBA\tSynteny\t$start\t".(scalar(@{$coverLB{$chr}})-1)."\t0\t+\t.\tnote=\"No synteny with reference. Possible insert or too divergent\"\n";
	  $noCovQry{$chr}+=($sizeQuery{$chr}-1-$start);
	}
  }
  saveGFF($resultName,"Reference",\%BB);
  saveGFF($resultName,"Query",\%LB);

  my $res;
  
  foreach my $chr (sort keys %sizeRef ) {
	$res.=sprintf("Of the reference chromosome $chr\t%.2f per cent\thas no synteny with the query\n",(( $sizeRef{$chr} - $coveredRef{$chr})*100/$sizeRef{$chr}) );
  }
  foreach my $chr (sort keys %sizeQuery ) {
	$res.=sprintf("Of the query chromosome $chr\t %.2f per cent\thas no synteny with the reference\n",((( $sizeQuery{$chr}-$coveredQuery{$chr} ))*100/$sizeQuery{$chr}) );
  }
  return $res
}



###########################################
### Seperate Replicon
############################################

sub Split{
  my $FileName = shift;
  my $Amount = shift;
  
  if (!defined($Amount)) {
	$Amount = 99999999999;
  }
  open (FILE, $FileName) or die "Couldn't open file to Seqparate $FileName $!\n";
  
  my $Counter = 0;
  my $Out;
  my @NameFiles;
  while (<FILE>) {
	if ($Counter < $Amount) {
	  
	  if (/^>(\S+)/) {
		close(FILEOUT);
		my $Name = $1;#."_$Counter.fasta";
		open (FILEOUT, "> $Name") or die "Problem do $Name file  $!\n";
		push @NameFiles, $Name;
		$Counter++;
		print FILEOUT $_;
	  }
	  else {
		print FILEOUT $_;
	  }
	}
  }
  print "$Counter Sequences where generated\nDone.\n";
  return \@NameFiles;
}

###########################################
### Seperate Replicon
############################################

sub Embl2Fasta{
  my $emblDir = shift;
  my $fastaResult = shift;

  
  opendir (DIR, $emblDir) or die "Problem to open opendir $emblDir: $!\n";
  
  ### will hold the new annotation: $h{queryname}.=annotation as embl

  my $fasta;
  
  map {
	if (/embl$/){
	  open F, "$emblDir/$_ " or die "Problem open $emblDir/$_ in Embl2Fasta:$! \n";
	  my ($name)=$_ =~ /^(\S+)\.embl/;
	  
	  while (<F>) {
		if (/^SQ/) {
		  $fasta.=">$name\n";
		  while (<F>) {
			if (/\/\//) {
			  
			  last;
			}
			## get away space and number of seq
			s/\d+//g;
			s/\s+//g;
			$fasta.=$_."\n";
						
		  }
		}
	  }
	  
	}
  } readdir(DIR);
  open F, "> $fastaResult" or die "Couldn't write file $fastaResult in Embl2Fasta: $!\n";
  print F $fasta;
  close(F);
  
}
