package ratt_correction;

# tdo 01.11.09: include that when new contig, the contig with max
# overhang left, is take

#To DO lkist
# inlcude %lefthandPosContig

use 5.008008;
use strict;
use warnings;
use Data::Dumper;
no warnings "recursion";
require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration       use Mapping::MappingTools ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = (
    'all' => [
        qw(
		   correctEMBL
          )
    ]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( startAnnotationCorrection
  startAnnotationCheck
  correctEMBL

);

our $VERSION = '0.01';
my $DEBUG      = 0;
my %STOP_CODON = (
    'TGA' => 1,
    'TAG' => 1,
    'TAA' => 1
);
my $CORRECT_SPLICESITE=0;
my %START_CODON     = ( 'ATG' => 1 );
my %SPLICE_DONOR    = ( 'GT'  => 1 );
my %SPLICE_ACCEPTOR = ( 'AG'  => 1 );

# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Protocols::Assembly::TillingRaph - Perl extension for performing
assemblies guided with a reference.

=head1 SYNOPSIS

use ratt_Correction;

=head1 DESCRIPTION



=head1 EXPORT

startTiling

=head1 SEE ALSO

http://scratchy.internal.sanger.ac.uk/wiki/index.php/Team_133

# Preloaded methods go here.

=head1 FUNCTIONS

=head2 startTilling - This function will call all programs

        Usage           :

        Arg [1]         : LSF queue for jobs

        Arg [2]         : file of mapped lanes to be summarised

        Arg [3]         : faidx index file for reference

        Example         :

        Description     : Summarises the coverage of the mapped lanes specified in the lanes.fofn file

        Returntype      : none

        Author          : Jacqueline McQuillan  E<lt>tdo@sanger.ac.uk<gt>

=cut

sub startAnnotationCheck {
    my $embl         = shift;
    my $sequenceFile = shift;
    my $ResultName   = shift;

    doit( $embl, $sequenceFile, $ResultName, "Check" );
}

sub startAnnotationCorrection {
    my $embl         = shift;
    my $sequenceFile = shift;
    my $ResultName   = shift;

    doit( $embl, $sequenceFile, $ResultName, "Correct" );

}

sub doit {
    my $embl         = shift;
    my $sequenceFile = shift;
    my $ResultName   = shift;
    my $method       = shift;

    loadConfig();
    my $sequence = loadSequence($sequenceFile);

    ### will have the new annotation
    my $newAnnotation;

    ### will have the stats per genes, as hash, see description
    my $ref_stats;

    ### gff file including position of errors and done changes
    my $GFFfile;

    my ( $revAnnotation, $ref_revStats, $revGFFfile, $newrevEMBL, $newEMBL );
    if ( $method eq "Check" ) {
        my $revSequence = reverseEMBL( $embl, $sequence, "reverse.embl" );
        ( $newAnnotation, $ref_stats, $GFFfile ) =
          checkEMBL( $embl, $sequence );
        ( $revAnnotation, $ref_revStats, $revGFFfile ) =
          checkEMBL( $embl . ".reverse.embl", $revSequence );
        $revGFFfile = reverseGFF( $revGFFfile, $sequence );
    }
    elsif ( $method eq "Correct" ) {
        ( $newAnnotation, $ref_stats, $GFFfile, $newEMBL ) =
          correctModel( $embl, $sequence );
        $revGFFfile   = "";
        $ref_revStats = "";
    }

    #  foreach my $gene (keys %$ref_stats) {
    #	print Dumper $$ref_stats{$gene}
    #  }
    # foreach my $gene (keys %$ref_revStats) {
    #	print Dumper $$ref_revStats{$gene}
    #  }

    open( F, "> $ResultName.final.embl" )
      or die "Couldn't write the gff file...\n";
    print F $newEMBL;
    close(F);

    # write results:
    open( F, "> $ResultName.Report.gff" )
      or die "Couldn't write the gff file...\n";
    print F $GFFfile . $revGFFfile;
    close(F);

    my $res .= printErrorStats($ref_stats);
    if ( $method eq "Check" ) {
        $res .= printErrorStats($ref_revStats);
    }
    print $res;

    open( F, "> $ResultName.Report.txt" )
      or die "Couldn't write the gff file...\n";
    print F $res;
    close(F);
    print
"done.\nPlease see the file $ResultName.Report.txt for reporting the errors and the file $ResultName.Report.gff for reporting the error in Artemis or other sequence viewer.\n";

}

####################
### correctModel
####################
sub correctModel {
    my $emblFile         = shift;
    my $originalSequence = shift;

    my $revSequence  = revcomp($originalSequence);
    my $sequence     = $originalSequence;
    my $isComplement = 0;
    my $ref_annotation;
    my $ref_stats;
    my $GFFfile;

    my $ref_structure;

    open( F, $emblFile )
      or die "Sorry, couldn't open annotation file $emblFile: $_\n";

    my @EMBL = <F>;

    my $pos = 0;
    my $res;
    foreach (@EMBL) {
	  chomp;
	  if (/^FT   CDS\s{2,}(\S+)$/) {
		
		my $line = $1;
		if (/complement/) {
		  $isComplement = 1;
		  $sequence     = $revSequence;
		  $line = reverseCDS( $line, length($sequence) );
		}
		my $ref_structure = getStructure($line);
		my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
		my ( $id, $product ) = getID( \@EMBL, $pos );
		if ( !defined($id) ) {
		  $id = $$ref_structure{start};
		}
		
		$$ref_stats{$id}{product} = $product;
		
		debug(1,"\nWorking on $id ".$$ref_structure{pos}[0]);
		
		
					if ($DEBUG > 500) {
		  print Dumper  \@{ $$ref_structure{pos} };
		  
		}	
		if ( !isStartOK($cds) ) {
		  debug(50,"Start is wrong");
		  #		  $GFFfile.=doGFF($$ref_structure{start},"BadStart","Start wrong",$isComplement,length($sequence));
		  $$ref_stats{$id}{StartBad} = 1;
		  $$ref_stats{$id}{error}++;
		  $ref_structure = correctStart(
										$ref_structure,         $sequence,
										\%{ $$ref_stats{$id} }, \$GFFfile
									   );
		  my $ok = 0;
		  
		  if ($DEBUG>50) {
			print Dumper @{$$ref_structure{pos}}
		  }
		  
		  ( $ref_structure, $ok ) = extentModelUpstreamStart(
															 $ref_structure,         $sequence,
															 \%{ $$ref_stats{$id} }, \$GFFfile
																);
		  if ($DEBUG>50) {
			print Dumper @{$$ref_structure{pos}}
		  }
		  if ( !$ok ) {
			( $ref_structure, $ok ) = extentModelDownstreamStart(
																 $ref_structure,         $sequence,
																 \%{ $$ref_stats{$id} }, \$GFFfile
																);
		  }
		  
		  if ( !$ok ) {
			$GFFfile .=
			  doGFF( $$ref_structure{start}, "BadStart", "Start wrong",
					 $isComplement, length($sequence) );
			$$ref_stats{$id}{StartStillBad} = 1;
			$$ref_stats{$id}{errorStill}++;
		  }
		  if ($ok) {
			$$ref_stats{$id}{CorrectionLog} .= " // Corrected Start";
			$GFFfile .= doGFF(
							  $$ref_structure{start}, "CorrectStart",
							  "Corrected Start",      $isComplement,
							  length($sequence)
							 );
		  }
		  
		  if ($DEBUG>50) {
			print Dumper @{$$ref_structure{pos}}
		  }
		}
		
		
		### check to undo old frameshifts
		$ref_structure =
		  checkIntrons( $ref_structure, $sequence, \%{ $$ref_stats{$id} },
						\$GFFfile, $isComplement );
		
		my ( $amount, $lastPos ) =
		  amountWrongspliceDonors( $ref_structure, $sequence );
		  
		if ($amount) {
		  debug( 1, "Splice donor wrong" );
		  
		  $GFFfile .=
			doGFF( $lastPos, "Wrong_Splice",
				   "Model has $amount wrong splice sites.",
				   $isComplement, length($sequence) );
		  $$ref_stats{$id}{splicesites} = $amount;
		  $$ref_stats{$id}{error}++;
		}
		
		### try to correct one splice site.
		if ($amount ==1 && $CORRECT_SPLICESITE){
			my $result = correctSplicesite(
										   $ref_structure,         $sequence,
										   \%{ $$ref_stats{$id} }, \$GFFfile, $isComplement 
										  );
			
		  }
		
		
		##update the cds
		$cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
		debug(10,"Check Frame shifts");
			if ($DEBUG > 500) {
		  print Dumper  \@{ $$ref_structure{pos} };
		  
		}	
		if ( ( my $amount = getAmountFrameshifts($cds) ) ) {

		  $$ref_stats{$id}{frameshifts} = $amount;
		  $$ref_stats{$id}{error}++;
		  my ($ok) = 0;
		  $ref_structure =
			checkFrameShits( $ref_structure, \%{ $$ref_stats{$id} },
							 $sequence );

		  ### check if model still have framesfhits
		  $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
		  my $amount2 = getAmountFrameshifts($cds);
		  
		  if ( $amount2 == 0 ) {
			$GFFfile .= doGFF(
							  $$ref_structure{start},
							  "FrameshiftCorrected",
							  ( $amount - $amount2 ) . " frameshifts corrected.",
							  $isComplement,
							  length($sequence)
							 );
		  }
		}
		
		if ($DEBUG > 500) {
		  print Dumper  \@{ $$ref_structure{pos} };
		  
		}
		
		$cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
		debug(10,"Check Length");
			
            if ( !isMod3Length($cds) ) {
			  $GFFfile .= doGFF(
								$$ref_structure{start},  "Error",
								"Model length is wrong - shiftet it to correct length", $isComplement,
								length($sequence)
							   );
			  $$ref_stats{$id}{length} = 1;
			  $$ref_stats{$id}{error}++;
			  my $count_=0;
			  
			  while (! isMod3Length($cds) && $count_ <6 ){
				$count_++;
				$ref_structure=shiftStructureLast($ref_structure,1);
				$cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );	
			  }
			  
            }
			
			debug(10,"Check Stop");
			
			
           $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );


			my $amountFrame = getAmountFrameshifts($cds);
			
            if ( !isStopOK($cds) || ($amountFrame > 0 && $amountFrame < 5)) {
                ### two case, no stop codon, so look at the end.
                ### more than one, find an earlier stop...
                debug( 50, "Stop is wrong" );

                $$ref_stats{$id}{StopBad} = 1;
                $$ref_stats{$id}{error}++;
                $cds=getLastExon($ref_structure,$sequence);
                my $amount = getAmountFrameshifts($cds);
                my $ok     = 0;
                if ( $amount == 0 ) {
                    ( $ref_structure, $ok ) = extentModelDownstreamStop(
                        $ref_structure,         $sequence,
                        \%{ $$ref_stats{$id} }, \$GFFfile
                    );
                }
                else {
				  
                    ( $ref_structure, $ok ) = extentModelUpstreamStop(
                        $ref_structure,         $sequence,
                        \%{ $$ref_stats{$id} }, \$GFFfile
                    );
                }
                ### check if still stop
                $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
                if ( !isStopOK($cds) ) {
				  $$ref_stats{$id}{StopStillBad} = 1;
				  $$ref_stats{$id}{errorStill}++;
				  $GFFfile .= doGFF(
									$$ref_structure{end}, "BadStop",
									"Stop wrong",         $isComplement,
									length($sequence)
								   );
                }
                else {
				  $$ref_stats{$id}{CorrectionLog} .= " // Corrected Stop";
				  $GFFfile .= doGFF(
									$$ref_structure{end}, "CorrectStop",
									"Corrected Stop",     $isComplement,
									length($sequence)
								   );
                }

				### check if model still have framesfhits
				$cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
				my $amount2 = getAmountFrameshifts($cds);

				if ( $amount2 > 0 ) {
				  $GFFfile .= doGFF(
									$$ref_structure{start},              "Frameshift",
									"The model has $amount frameshifts", $isComplement,
									length($sequence)
								   );
				  $$ref_stats{$id}{frameshiftsStill} = $amount2;
				}
				
			  }
		
            if ($isComplement) {
                  my $tmp = printStructurePos($ref_structure);
                $tmp =~ /^(FT   \S+\s+)(\S+)$/;
                my $CDS = reverseCDS( $2, length($sequence) );
                $res .= $1 . $CDS . "\n";
                $sequence     = $originalSequence;
                $isComplement = 0;
            }
            else {
                $res .= printStructurePos($ref_structure) . "\n";
            }

        }    # enf if
        else {
            $res .= $_ . "\n";
        }    # else if CDS
        $pos++;
    }    # foreach
    return ( $ref_annotation, $ref_stats, $GFFfile, $res );

}

####################
### checkIntrons
####################
sub checkFrameShits {
    my ( $ref_structure, $ref_statsGene, $sequence ) = @_;

    # 1. find the exon that has the frameshift
    # 2. get the real sequence position of it. 
    # 3. find the needed amount of frames, as for the start
    # 4. rebuild the model and call checkFrameshifts
    
 	### go through the exons
   	my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );

    my $pos=0;
    while (! isStopOK(uc(substr($cds,$pos,3))) && $pos <= length($cds) ){
    	$pos+=3;
    }

	# as the function is just called, if there are frameshifts, $pos holds the start of the stop
	### never the less check, if a stop was found
	if ($pos >= length($cds)){
		warn "checkFrameShifts was called, without a cds in the coding sequence\n";
		return $ref_structure	
	}
	else {
		   my $ref_structPos = getStructure2ar( $ref_structure );
 
			### exon must stop $posFirstStop-1-1
			my ($posStop,$numExon) = $$ref_structPos[$pos] =~ /^(\d+)\s(\d+)$/; 		
			
			### the next "exon" will start +3 bases. All three frames should be tested.
			my $foundNewExon=0;
			   # holds the shift position 1-3. If it is 0, means that two or more are equal
			my ($start,$stop) = $$ref_structure{pos}[($numExon-1)] =~ /^(\d+)\.\.(\d+)$/;
			# ( cannot be one base exon
			
			my $midStop  = ($posStop - 3); # position (absolute in str) of last ok codon
			my $midStart = ($posStop + 3) ; # first possible start of new exon			
   			my $minAllpos = 0;
			my $minAll    = 9999999;
    		my $amountF   = 0;

    		my $newExon;
    		
    		if (($stop-$midStart+1) > 150){
	    	for my $shift ( 0 .. 2 ) {
    		   $newExon = substr($sequence,($midStart-1+$shift),($stop-$midStart+1));
    		   $amountF = getAmountFrameshifts($newExon);
        		if ($amountF < $minAll) {
	        		( $minAll, $minAllpos ) = min( $amountF, $shift, $minAll, $minAllpos );
        		}
    		}
    		
    		### change the $ref_structure, include this new exon
			$ref_structure = updateStructure($ref_structure,$numExon,($midStop+2),($midStart+$minAllpos));
			
			### check if the model has still frameshift. 
			### as it can be a dogy end, we have the variable $foundNewExon
			### Todo This is a possible problem,  
			$cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
            my $amountFrameShift = getAmountFrameshifts($cds);

			if ($foundNewExon && $amountFrameShift > 0){
				$ref_structure = checkFrameShits($ref_structure, $ref_statsGene, $sequence);
			}
    		}
	}    
    return $ref_structure;
}



####################
### getStructure2ar
####################
sub getStructure2ar{
	my $ref_structure = shift;
	
	my $amountExons = scalar(@{$$ref_structure{pos} });
	my @ar;
	my $pos;
	foreach my $exons (1..$amountExons){
		if ($$ref_structure{pos}[($exons-1)] =~ /^(\d+)\.\.(\d+)$/) {
            ### if a gene model got shiftet, the last exon might be gone, so it should be rebuild
    		for ($1..$2){
    			push @ar, "$_ $exons"	
    		}
        }
        elsif ($$ref_structure{pos}[($exons-1)] =~ /^(\d+)$/) {
    			push @ar, "$_ $exons"	
        }
        else {
			die "Problem with getStructure2ar\n";
        }
	}
	return \@ar
}

####################
### updateStructure
####################
sub updateStructure{
	my $ref_structure = shift;
	my $exon = shift;
	my $midStop = shift;
	my $midStart = shift;
	
	my @new;
	my $amountExons = scalar(@{$$ref_structure{pos}});
	foreach my $num (1..$amountExons){
		if ($num == $exon){
			my ($s,$e) = $$ref_structure{pos}[($num-1)] =~ /^(\d+)\.\.(\d+)$/;
			push @new, "$s..$midStop";
			push @new, "$midStart..$e";
			$$ref_structure{exons}++
		}
		else {
			push @new, $$ref_structure{pos}[($num-1)];	
		}
	}
	
	### put the new coordantes onto the structure
	$$ref_structure{pos}=\@new;
		
	return $ref_structure
}

####################
### checkIntrons
####################
sub checkIntrons {
    my ( $ref_structure, $sequence, $ref_statsGene, $ref_GFFfile,
        $isComplement ) = @_;

    ### if an intron in %mod3=0 long, and has no stop at all, it can be taken away.
    my $exons = ( scalar( @{ $$ref_structure{pos} } ) );
    my @ar;
    my $found = 0;
    for my $number ( 0 .. ( $exons - 2 ) ) {
        my ($start) = $$ref_structure{pos}[$number] =~ /(\d+)$/;
        my ($end) = $$ref_structure{pos}[ ( $number + 1 ) ] =~ /^(\d+)/;

        my $intron = substr( $sequence, ($start), ( $end - $start - 1 ) );
        if (
               isMod3Length($intron)
            && !isStopOK($intron)
            && getAmountFrameshifts($intron) == 0
            && !defined( $SPLICE_DONOR{ uc( substr( $intron, 0, 2 ) ) } )
            && !defined(
                $SPLICE_ACCEPTOR{ uc( substr( $sequence, ( $end - 3 ), 2 ) ) }
            )
          )
        {
            my ($oldstart) = $$ref_structure{pos}[$number] =~ /^(\d+)/;
            my ($oldend) = $$ref_structure{pos}[ ( $number + 1 ) ] =~ /(\d+)$/;
            $$ref_structure{pos}[ ( $number + 1 ) ] =
              $oldstart . ".." . $oldend;
            $$ref_structure{pos}[ ($number) ] = '';
            $found++;
        }
    }
    if ($found) {
        $$ref_statsGene{JoinExons} = $found;
        $$ref_statsGene{error}--;
        $$ref_GFFfile .=
          doGFF( $$ref_structure{start}, "JoinedExons",
            "$found introns/frameshitf were eliminated",
            $isComplement, length($sequence) );
        $$ref_statsGene{CorrectionLog} .= " // Joined $found introns";
        my @new;
        foreach ( @{ $$ref_structure{pos} } ) {
            if ($_) {
                push @new, $_;
            }
        }
        $$ref_structure{pos} = \@new;
    }
    return $ref_structure;
}
####################
### correctStart
####################
sub printStructurePos {
    my $ref_structure = shift;
    my $Tag           = "CDS";
    if ( !defined($Tag) ) {
        $Tag = "CDS";
    }
    my $pos = join( ',', @{ $$ref_structure{pos} } );
    if ( scalar( @{ $$ref_structure{pos} } ) > 1 ) {
        $pos = "join($pos)";
    }
    return "FT   $Tag             " . $pos;
}

####################
### correctStart
####################
sub correctStart {
    my $ref_structure = shift;
    my $sequence      = shift;
    my $ref_statsGene = shift;
    my $ref_GFF       = shift;

    ### we want to have less possible frameshift in the first 50aa and in total.
    my $minAll = 99999;

    # holds the shift position 1-3. If it is 0, means that two or more are equal
    my $minAllpos = 0;
    my $min150    = 99999;
    my $min150pos = 0;
    my $amountF   = 0;

    my $ref_structureN = copyStructure($ref_structure);

    for ( 1 .. 3 ) {
        $ref_structureN = shiftStructure( $ref_structureN, 1 );
        my $cds = buildGene( \@{ $$ref_structureN{pos} }, $sequence );
        $amountF = getAmountFrameshifts($cds);
        ( $minAll, $minAllpos ) = min( $amountF, $_, $minAll, $minAllpos );
        $amountF = getAmountFrameshifts( substr( $cds, 0, 150 ) );
        ( $min150, $min150pos ) = min( $amountF, $_, $min150, $min150pos );
    }

    if ( $min150pos == "0" ) {
        if ( $minAllpos != "0" ) {
            $ref_structure =
              correctModelShift( $ref_structure, $minAllpos, $ref_statsGene,
                $ref_GFF );
        }
        else {
            $$ref_statsGene{'CorrectionLog'} .=
              " // No unique shift for model found... MANUAL check";

        }
    }
    ### one shift has a minimum
    else {
        $ref_structure =
          correctModelShift( $ref_structure, $min150pos, $ref_statsGene,
            $ref_GFF );
    }

    return $ref_structure;
}

####################
### copyStructure
####################
sub copyStructure {
    my $ref_structure = shift;

    my %new;

    foreach my $key ( keys %$ref_structure ) {
        if ( defined( @{ $$ref_structure{$key} } ) ) {
            foreach my $val ( @{ $$ref_structure{$key} } ) {
                push @{ $new{$key} }, $val;
            }
        }
        else {
            $new{$key} = $$ref_structure{$key};
        }
    }

    return \%new;
}

####################
### findFirstStop
####################
sub findFirstStop {
    my ( $ref_structure, $sequence, $ref_statsGene, $refGFF ) = @_;
	die "not implemented\n";
	
    my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );

	my $pos=0;
    while ( !isMod3Length($cds) && $pos <=5 ) {
        debug( 5, "shift in extentModelDownstreamStop" );

        $ref_structure = shiftStructureLast( $ref_structure, -1 );
        $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
        $pos++
    }
    if ($pos > 3){
    	debug(1,"Couldn't find a good shift in extentModelDownstreamStop")	
    }
    my $AmountExons = scalar( @{ $$ref_structure{pos} } );
    my ($stopPos) = $$ref_structure{pos}[ ( $AmountExons - 1 ) ] =~ /(\d+)$/;
    my $ok = walkFirstStop( ( $stopPos + 1 ), \$sequence, 1, 9991, 3 );
    if ( $ok != 9991 ) {
        ### the +2 is necessary as it is the first position of the stopcodon...
        $ref_structure = shiftStructureLast( $ref_structure, ( $ok + 2 ) );
        $ok = 1;
    }
    else { $ok = 0 }
    return ( $ref_structure, $ok );
}


####################
### extentModelDownstreamStop
####################
sub extentModelDownstreamStop {
    my ( $ref_structure, $sequence, $ref_statsGene, $refGFF ) = @_;

    my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );

	my $pos=0;
    while ( !isMod3Length($cds) && $pos <=5 ) {
        debug( 5, "shift in extentModelDownstreamStop" );

        $ref_structure = shiftStructureLast( $ref_structure, -1 );
        $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
        $pos++
    }
    if ($pos > 3){
	  debug(1, "Couldn't find a good shift in extentModelDownstreamStop")	
    }
    my $AmountExons = scalar( @{ $$ref_structure{pos} } );
    my ($stopPos) = $$ref_structure{pos}[ ( $AmountExons - 1 ) ] =~ /(\d+)$/;
    my $ok = walkFirstStop( ( $stopPos + 1 ), \$sequence, 1, 9991, 3 );
    if ( $ok != 9991 ) {
        ### the +2 is necessary as it is the first position of the stopcodon...
        $ref_structure = shiftStructureLast( $ref_structure, ( $ok + 2 ) );
        $ok = 1;
    }
    else { $ok = 0 }
    return ( $ref_structure, $ok );
}

####################
### extentModelUpstreamStop
####################
sub extentModelUpstreamStop {
    my ( $ref_structure, $sequence, $ref_statsGene, $refGFF ) = @_;

    my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
	
    while ( !isMod3Length($cds) ) {
        $ref_structure = shiftStructureLast( $ref_structure, -1 );
        $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
    }
	
    my $AmountExons = scalar( @{ $$ref_structure{pos} } );
    my ($stopPos) = $$ref_structure{pos}[ ( $AmountExons - 1 ) ] =~ /(\d+)$/;
    my $ok = walkFirstStop( ( $stopPos - 2 ), \$sequence, -2, 9991, -3 );
    if ( $ok != 9991 ) {
        ### the +2 is necessary as it is the first position of the stopcodon...
        $ref_structure = shiftStructureLast( $ref_structure, ( $ok + 2 ) );
        $ok = 1;
    }
    else {
	  $ok = 0
	}
	
    $cds = getLastExon( $ref_structure, $sequence );
	my $peace=substr( $cds, (length($cds)-150));
	
    my $amount = getAmountFrameshifts( substr( $cds, (length($cds)-150)) );
	
	# check this
    if ($amount > 0 && $amount < 5 && length($cds) > 50 ) {
        $ref_structure = shiftStructureLast( $ref_structure, -1 );
        ( $ref_structure, $ok ) =
          extentModelUpstreamStop( $ref_structure, $sequence, $ref_statsGene,
            $refGFF );
    }


    return ( $ref_structure, $ok );
}
####################
### extentModelUpstreamStart
####################
sub extentModelUpstreamStart {
    my ( $ref_structure, $sequence, $ref_statsGene, $refGFF ) = @_;

    my ($startPos) = $$ref_structure{pos}[0] =~ /^(\d+)/;

    my $ok = walkLastStart( $startPos, \$sequence, 0, 9991, -3 );
    if ( $ok != 9991 ) {
        $ref_structure = shiftStructure( $ref_structure, $ok );
        $ok = 1;
    }
    else { $ok = 0 }

    return ( $ref_structure, $ok );
}

####################
### extentModelDownstreamStart
####################
sub extentModelDownstreamStart {
    my ( $ref_structure, $sequence, $ref_statsGene, $refGFF ) = @_;

    ###
	my $ok=0;
    my ($startPos,$endPos) = $$ref_structure{pos}[0] =~ /^(\d+)\.\.(\d+)/;
	if (defined($startPos)) { # one exon genes

	  # if maximal distance a start codon can be looked for - half the
	  # gene model
	  my  $maxSearch=int(($endPos-$startPos)/2);
	  
	  $ok = walkFirstStart( $startPos, \$sequence, 0, 9991, 3, $maxSearch );
	  if ( $ok != 9991 ) {
        $ref_structure = shiftStructure( $ref_structure, $ok );
        $ok = 1;
	  }
	  else { $ok = 0 }
	}
	
    return ( $ref_structure, $ok );
}

####################
### walkStart
####################
sub walkLastStart {
    my ( $check, $ref_sequence, $dist, $ok, $step ) = @_;

    ## recursive procedure, stops, if
    ### has no sequence start end sequence
    ### finds a stop in frame
    ### return the difference to the possible start to the first starting point
    ###

    if ( $check > 3 && $check < ( length($$ref_sequence) - 2 ) ) {

        my $codon = uc( substr( $$ref_sequence, ( $check - 1 ), 3 ) );

        if ( isStopOK($codon) ) {
            return $ok;
        }
        else {
            if ( defined( $START_CODON{$codon} ) ) {
                $ok = $dist;
            }
            $dist += $step;
            $ok = walkLastStart( ( $check + $step ),
                $ref_sequence, $dist, $ok, $step );
        }
    }
    return $ok;
}
####################
### walkStart
####################
sub walkFirstStart {
    my ( $check, $ref_sequence, $dist, $ok, $step, $maxRun ) = @_;

    ## recursive procedure, stops, if
    ### has no sequence start end sequence
    ### finds a stop in frame
    ### return the difference to the possible start to the first starting point
    ###
	
    if ( $maxRun > 1 && $check > 3 && $check < ( length($$ref_sequence) - 2 ) ) {
        my $codon = uc( substr( $$ref_sequence, ( $check - 1 ), 3 ) );
		$maxRun-=$step;

		
        #	if (isStopOK($codon)){#
        #		return $ok;
        #	}
        if ( defined( $START_CODON{$codon} ) ) {
            $ok = $dist;

            return $ok;
        }
        else {
            $dist += $step;
            $ok = walkFirstStart( ( $check + $step ),
                $ref_sequence, $dist, $ok, $step, $maxRun )

        }
    }
    return $ok;
}
####################
### walkStart
####################
sub walkFirstStop {
    my ( $check, $ref_sequence, $dist, $ok, $step ) = @_;

    ## recursive procedure, stops, if
    ### has no sequence start end sequence
    ### finds a stop in frame
    ### return the difference to the possible start to the first starting point
    ###
	
    if ( $check < 3 || $check > ( length($$ref_sequence) - 2 ) ) {

        return $ok;
    }
    else {
        my $codon = uc( substr( $$ref_sequence, ( $check - 1 ), 3 ) );

        if ( isStopOK($codon) ) {

            return $dist;
        }
        else {
            $dist += $step;

            $ok = walkFirstStop( ( $check + $step ),
                $ref_sequence, $dist, $ok, $step )

        }
    }
    return $ok;
}

####################
### min
####################
sub min {
    my ( $new, $newPos, $min, $minPos ) = @_;

    if ( $new == $min ) {
        return ( $new, 0 );
    }
    elsif ( $new > $min ) {
        return ( $min, $minPos );
    }
    else {
        return ( $new, $newPos );
    }
}

####################
### correctModelShift
####################
sub correctModelShift {
    my ( $ref_structure, $shift, $ref_statsGene, $ref_GFF ) = @_;
    if ( ( $shift % 3 ) > 0 ) {
        $$ref_statsGene{'CorrectionLog'} .= " // Modell shifted $shift base(s)";

#		$$ref_GFF.=doGFF($$ref_structure{start},"Modell_Shifted","Modell shifted $shift bases");
    }
    return shiftStructure( $ref_structure, ( $shift % 3 ) );
}

####################
### shiftStructure
####################
sub shiftStructure {
    my ( $ref_structure, $shift ) = @_;
    ### important: We don't handle complements, as reversed before!

    ## change the first digit of the first exon
    $$ref_structure{pos}[0] =~ s/^(\d+)/($1+$shift)/e;
    ### one exon bases ok...

    return $ref_structure;
}

####################
### shiftStructureLast
####################
sub shiftStructureLast {
    my ( $ref_structure, $shift ) = @_;
    ### important: We don't handle complements, as reversed before!

    my $amountExon = ( scalar( @{ $$ref_structure{pos} } ) );
    ## change the last digit of the last exon
    ## change the last digit of the last exon
    if ( $$ref_structure{pos}[ ( $amountExon - 1 ) ] =~ /(\d+)\.\.(\d+)/ ) {
        if ( $1 < ( $2 + $shift ) ) {
            $$ref_structure{pos}[ ( $amountExon - 1 ) ] =~
              s/(\d+)$/($1+$shift)/e;
        }
    }
    elsif ( $$ref_structure{pos}[ ( $amountExon - 1 ) ] =~ /(\d+)/ ) {
        if ( abs($shift) < 1000 ) {
            $$ref_structure{pos}[ ( $amountExon - 1 ) ] =~
              s/(\d+)$/($1+$shift)/e;
        }
    }
    else {
        print Dumper $ref_structure;

        die "couldn't transfer at all... \n";

    }

    ###	$$ref_structure{pos}[($amountExon-1)] =~  s/(\d+)$/($1+$shift)/e;
    ### one exon bases ok...

    return $ref_structure;
}

####################
### printErrorStats
####################
sub printErrorStats {
    my $ref_stats = shift;

    my @flags = (
        'error',            'StartBad',
        'StopBad',          'frameshifts',
        'splicesites',      'length',
        'product',          'errorStill',
        'StartStillBad',    'StopStillBad',
        'frameshiftsStill', 'JoinExons',
        'PossiblePseudo',
        'CorrectionLog'
    );
    my $res = "\nGene_ID\t";
    foreach (@flags) {
        $res .= $_ . "\t";
    }
    $res .= "\n";
    foreach my $gene ( sort keys %$ref_stats ) {
        if ( defined( $$ref_stats{$gene}{error} ) ) {
            $res .= "$gene\t";

            foreach (@flags) {
                if ( defined( $$ref_stats{$gene}{$_} ) ) {
                    $res .= "$$ref_stats{$gene}{$_}\t";
                }
                else {
                    $res .= "0\t";

                }
            }
            $res .= "\n";
        }

    }
    return $res;
}

####################
### checkEMBL
####################
sub checkEMBL {
    my $emblFile = shift;
    my $sequence = shift;

    my $ref_annotation;
    my $ref_stats;
    my $GFFfile;

    my $ref_structure;

    open( F, $emblFile )
      or die "Sorry, couldn't open annotation file $emblFile: $_\n";

    my @EMBL = <F>;

    my $pos = 0;
    foreach (@EMBL) {
        chomp;
        if (/^FT   CDS\s{2,}(\S+)$/) {

            my $line = $1;

            if ( !/complement/ ) {

                $ref_structure = getStructure($line);
                my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
                my ( $id, $product ) = getID( \@EMBL, $pos );
                if ( !defined($id) ) {
                    $id = $$ref_structure{start};
                }

                $$ref_stats{$id}{product} = $product;

                if ( !isStopOK($cds) ) {
                    $GFFfile .=
                      doGFF( $$ref_structure{end}, "BadStop", "Stop wrong" );
                    $$ref_stats{$id}{StopBad} = 1;
                    $$ref_stats{$id}{error}++;
                }
                if ( !isStartOK($cds) ) {
                    $GFFfile .=
                      doGFF( $$ref_structure{start}, "BadStart",
                        "Start wrong" );
                    $$ref_stats{$id}{StartBad} = 1;
                    $$ref_stats{$id}{error}++;
                }
                if ( my $amount = getAmountFrameshifts($cds) ) {
                    $GFFfile .=
                      doGFF( $$ref_structure{start}, "Frameshift",
                        "The model has $amount frameshifts" );
                    $$ref_stats{$id}{frameshifts} = $amount;
                    $$ref_stats{$id}{error}++;
                }
                if ( !isMod3Length($cds) ) {
                    $GFFfile .=
                      doGFF( $$ref_structure{start}, "Error",
                        "Model length is wrong" );
                    $$ref_stats{$id}{length} = 1;
                    $$ref_stats{$id}{error}++;
                }
                my ( $amount, $lastPos ) =
                  amountWrongspliceDonors( $ref_structure, $sequence );
                if ($amount) {
                    $GFFfile .= doGFF( $lastPos, "Wrong_Splice",
                        "Model has $amount wrong splice sites." );
                    $$ref_stats{$id}{splicesites} = $amount;
                    $$ref_stats{$id}{error}++;
                }
            }

        }
        $pos++;
    }
    return ( $ref_annotation, $ref_stats, $GFFfile );

}

####################
### doGFF
####################
sub doGFF {
    my $pos       = shift;
    my $tag       = shift;
    my $str       = shift;
    my $isReverse = shift;
    my $length    = shift;

    $str =~ s/\s/+/g;

    if ( defined($isReverse) and $isReverse ) {
        return
            "unknown\tBBA\t$tag\t"
          . ( $length - $pos + 1 ) . "\t"
          . ( $length - $pos + 1 )
          . "\t0\t-\t.\tnote=\"$str\"\n";
    }
    else {
        return "unknown\tBBA\t$tag\t$pos\t$pos\t0\t+\t.\tnote=\"$str\"\n";
    }
}

####################
### getID
####################

sub getID {
    my $ref_ar = shift;
    my $pos    = shift;

    my $notFound = 1;
    my $id;

    my $product = 'Not Set';
    $pos++;

    while (
        defined( $$ref_ar[$pos] )
        && !( $$ref_ar[$pos] =~ /^FT   \S+\s+/ )

      )
    {
        chomp( $$ref_ar[$pos] );

        if ( $$ref_ar[$pos] =~ /FT\s+.*.locus_tag="(\S+)\"/ ) {
            $id = $1;
        }
        if ( $$ref_ar[$pos] =~ /FT\s+.*c_id=\"(\S+)\"/ ) {
            $id = $1;
        }
        elsif ( $$ref_ar[$pos] =~ /FT\s{4,}\/product=\"(.*\S+)\"{0,1}$/ ) {
            $product = $1;
        }
        $pos++;
    }
    return ( $id, $product );

}

####################
### getStructure
####################
sub getStructure {
    my $str = shift;

    my %struct;

    $str =~ s/join\(//g;
    $str =~ s/\)//g;
    $str =~ s/<//g;
    $str =~ s/</>/g;
    my $count = 0;
    my @part = split( /,/, $str );

    ### get the start of the gff
    $part[0] =~ /^(\d+)/;
    $struct{start} = $1;
    $part[ ( scalar(@part) - 1 ) ] =~ /(\d+)$/;
    $struct{end} = $1;

    foreach (@part) {
        if (/(\d+..\d+)/) {
            $count++;
            push @{ $struct{pos} }, $1;
        }
        elsif (/^(\d+)$/) {
            $count++;
            push @{ $struct{pos} }, $1;
        }
        elsif (/^\s*/) {
            print "empty line... bad annotation way... $_ $str \n";
        }
        else {
            die "Error in getStructure $_ \n";
        }
    }
    $struct{exons} = $count;

    return \%struct;
}


####################
### getLastExon
####################
sub getLastExon {
    my $ref_ar   = shift;
    my $sequence = shift;

    my $geneSeq;

	my $exons = (scalar(@{ $$ref_ar{pos} } ));
	
	if ($$ref_ar{pos}[($exons-1)] =~ /^(\d+)\.\.(\d+)$/) {
            ### if a gene model got shiftet, the last exon might be gone, so it should be rebuild
            if ( $2 > $1 ) {
                $geneSeq .= substr( $sequence, ( $1 - 1 ), ( $2 - $1 + 1 ) );
            }
        }
        elsif ($$ref_ar{pos}[($exons-1)] =~ /^(\d+)$/) {
            $geneSeq .= substr( $sequence, ( $1 - 1 ), 1 );
        }
        else {
			print "Amount of Exons $exons \n";
        	print Dumper $ref_ar;
            die "Problems $_ in getLastExon\n";
        }
    
    return $geneSeq;
}
####################
### buildGene
####################
sub buildGene {
    my $ref_ar   = shift;
    my $sequence = shift;

    my $geneSeq;

    foreach (@$ref_ar) {
        if (/^(\d+)\.\.(\d+)$/) {
            ### if a gene model got shiftet, the last exon might be gone, so it should be rebuild
            if ( $2 > $1 ) {
                $geneSeq .= substr( $sequence, ( $1 - 1 ), ( $2 - $1 + 1 ) );
            }
        }
        elsif (/^(\d+)$/) {
            $geneSeq .= substr( $sequence, ( $1 - 1 ), 1 );
        }
        else {
			print Dumper $ref_ar;
            die "Problems $_ in buildGene\n";
 
        }
    }
    return $geneSeq;
}

####################
### isStartOK
####################
sub isStartOK {
    my $seq = shift;

	if (!defined($seq)) {
	  return 0
	}
	
    if ( defined( $START_CODON{ uc( substr( $seq, 0, 3 ) ) } ) ) {
        return 1;
    }
    return 0;
}

####################
### isStopOK
####################
sub isStopOK {
    my $seq = shift;

    my $length = length($seq);

    if ( defined( $STOP_CODON{ uc( substr( $seq, ( $length - 3 ), 3 ) ) } ) ) {
        return 1;
    }
    return 0;
}


sub findSpliceSiteDonor{
	my $ref_structure = shift;
	my $sequence      = shift;
	my $newPos        = shift;
	my $stepSize      = shift;
	my $exonNumber    = shift;
	my $exonStart     = shift;
	my $maxStep       = shift;
	
	my $MIN_SIZE=10;
	my ($posNextExon) = $$ref_structure{pos}[($exonNumber)] =~ /^(\d+)/;
	
	### change the model dowanstream
  if ($posNextExon> ($newPos + $MIN_SIZE ) &&
  	defined( $SPLICE_DONOR{ uc( substr( $sequence, $newPos, 2 ) ) } )){
  	### possible splice site found
  	#update structure
  	$$ref_structure{pos}[($exonNumber-1)] = "$exonStart..$newPos";
  	# check of zero stop codons
	my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
	my $amount = getAmountFrameshifts($cds);
	if ($amount==0){
		return ($newPos,$ref_structure);
	}
  }
  elsif(($exonStart < ($newPos-$stepSize-$MIN_SIZE)) && defined( $SPLICE_DONOR{ uc( substr( $sequence, ($newPos-$stepSize), 2 ) ) } ))
  {
  		### possible splice site found
  	#update structure
  	$$ref_structure{pos}[($exonNumber-1)] = "$exonStart..".($newPos-$stepSize);
  	# check of zero stop codons
	my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
	my $amount = getAmountFrameshifts($cds);
	if ($amount==0){
		return (($newPos-$stepSize),$ref_structure);
	}
  }
	my $ok;
  if (($stepSize <= $maxStep) && 
  ($posNextExon > ($newPos + $MIN_SIZE) || ($exonStart < ($newPos-$stepSize-$MIN_SIZE)))){
  		($ok, $ref_structure)=findSpliceSiteDonor($ref_structure,$sequence,($newPos+1),($stepSize+2),$exonNumber,$exonStart,$maxStep); 
	 }	
	 else {
	 	return (0,$ref_structure)
	 }
	 return ($ok,$ref_structure);
}


sub findSpliceSiteAcceptor{
	my $ref_structure = shift;
	my $sequence      = shift;
	my $newPos        = shift;
	my $stepSize      = shift;
	my $exonNumber    = shift;
	my $exonEnd       = shift;
	my $maxStep       = shift;
	
	my $MIN_SIZE=10;
	my ($posNextExon) = $$ref_structure{pos}[($exonNumber-2)] =~ /^(\d+)/;
	
	### change the model dowanstream
  if ($exonEnd > ($newPos + $MIN_SIZE ) &&
  	defined( $SPLICE_ACCEPTOR{ uc( substr( $sequence, ($newPos-3), 2 ) ) } )){
  	### possible splice site found

  	#update structure
  	$$ref_structure{pos}[($exonNumber-1)] = "$newPos..$exonEnd";
  	# check of zero stop codons
	my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
	my $amount = getAmountFrameshifts($cds);
	if ($amount==0){
		return ($newPos,$ref_structure);
	}
  }
  elsif(($posNextExon < ($newPos-$stepSize-$MIN_SIZE)) && defined( $SPLICE_ACCEPTOR{ uc( substr( $sequence, ($newPos-$stepSize-3), 2 ) ) } ))
  {
  		### possible splice site found
  	#update structure
  	$$ref_structure{pos}[($exonNumber-1)] = ($newPos-$stepSize)."..$exonEnd";
  	# check of zero stop codons
	my $cds = buildGene( \@{ $$ref_structure{pos} }, $sequence );
	my $amount = getAmountFrameshifts($cds);
	if ($amount==0){
		return (($newPos-$stepSize),$ref_structure);
	}
  }
my $ok;
  if (($stepSize <= $maxStep) && 
  ($exonEnd > ($newPos + $MIN_SIZE ) || ($posNextExon < ($newPos-$stepSize-$MIN_SIZE)))){
  		($ok, $ref_structure)=findSpliceSiteAcceptor($ref_structure,$sequence,($newPos+1),($stepSize+2),$exonNumber,$exonEnd,$maxStep); 
	 }	
	 else {
#	 	print "Nothing found $stepSize <= $maxStep  $posNextExon > ($newPos + $MIN_SIZE)     ($exonStart < ($newPos-$stepSize-$MIN_SIZE)\n";
	 	return (0,$ref_structure)
	 }
	 return ($ok,$ref_structure);
}

####################
### correctSplicesite
####################
sub correctSplicesite {
    my $ref_structure = shift;
    my $sequence      = shift;
    my $ref_statsGene = shift;
    my $ref_GFF       = shift;
	my $isComplement  = shift;
	
    my $count     = 1;
    my $wrong     = 0;
    my $lastwrong = 0;

    foreach ( @{ $$ref_structure{pos} } ) {
        if ( $count != scalar(@{ $$ref_structure{pos} } )) {
            /(\d+)$/;
            # splice_donor is wrong GT
            if ( !defined( $SPLICE_DONOR{ uc( substr( $sequence, $1, 2 ) ) } ) )
            {
            	my $old=$_;
				 /(\d+)\.\.(\d+)$/;

            	### look on the downstream site for splice site
	            my ($ok, $ref_structure)=findSpliceSiteDonor($ref_structure,$sequence,$2,1,$count,$1,200); 
				if ($ok != 0){
					  $$ref_statsGene{'CorrectionLog'} .= " // Correct Splice Site (Donor)";

					  $$ref_GFF.=doGFF($ok,"SpliceSiteCorrect","Splice site correct.",$isComplement, length($sequence) );
					  
					}
				else {
				$$ref_structure{pos}[($count-1)	]=$old;
				}
            }
        }
        if ( $count != 1 ) {
            /^(\d+)/;
            if (
                !defined(
                    $SPLICE_ACCEPTOR{ uc( substr( $sequence, ( $1 - 3 ), 2 ) ) }
                )
              )
            {
            	my $old=$_;
				 /(\d+)\.\.(\d+)$/;

            	### look on the downstream site for splice site
	            my ($ok, $ref_structure)=findSpliceSiteAcceptor($ref_structure,$sequence,$1,1,$count,$2,200); 
				if ($ok != 0){
					  $$ref_statsGene{'CorrectionLog'} .= " // Correct Splice Site (Acceptor)";

				$$ref_GFF.=doGFF($ok,"SpliceSiteCorrect","Splice site correct.",$isComplement, length($sequence));

				}
				else {
		   			$$ref_structure{pos}[($count-1)	]=$old;
				}
		    }
        }
        $count++;
    }

    return ( $wrong, $lastwrong );
}

####################
### amountWrongspliceDonors
####################
sub amountWrongspliceDonors {
    my $ref_structure = shift;
    my $sequence      = shift;
	my $ref_stats     = shift;
	my $refGFF        = shift;

    # we know it is one wrong
    my $count     = 1;
    my $wrong     = 0;
    my $lastwrong = 0;

    foreach ( @{ $$ref_structure{pos} } ) {
        if ( $count != scalar(@{ $$ref_structure{pos} } )) {
            /(\d+)$/;
            if ( !defined( $SPLICE_DONOR{ uc( substr( $sequence, $1, 2 ) ) } ) )
            {
            	
                $wrong++;
                $lastwrong = $1;

            }
        }
        if ( $count != 1 ) {
            /^(\d+)/;
            if (
                !defined(
                    $SPLICE_ACCEPTOR{ uc( substr( $sequence, ( $1 - 3 ), 2 ) ) }
                )
              )
            {
                $wrong++;
                $lastwrong = $1;
            }
        }
        $count++;
    }

    return ( $wrong, $lastwrong );
}



####################
### isMod3Length
####################
sub isMod3Length {
    my $seq = shift;
	
	if (!defined($seq)) {
	  return 0
	}
    my $length = length($seq);

    return !( $length % 3 );
}
####################
### getAmountFrameshifts
####################
sub getAmountFrameshifts {
    my $seq = shift;

	if (!defined($seq)) {
	  return 0;
	}
    my $length = length($seq);

    my $amountFrameshifts = 0;

    # we don't want to count the stop codon
    for ( my $i = 0 ; $i < ( $length - 3 ) ; $i += 3 ) {
        if ( defined( $STOP_CODON{ uc( substr( $seq, $i, 3 ) ) } ) ) {
            $amountFrameshifts++;
        }
    }

    return $amountFrameshifts;
}

####################
### loadConfig
####################
sub loadConfig {

    if ( defined( $ENV{RATT_CONFIG} )
        && -f "$ENV{RATT_CONFIG}" )
    {
        open( F, "$ENV{RATT_CONFIG}" )
          or die
          "Couldn't open Config file $ENV{RATT_CONFIG}\n";

		print "Using the $ENV{RATT_CONFIG} file for specifications.\n";
		
		
        my $count = 0;
        while (<F>) {
            chomp;
            if (/#START/) {
                $count = 1;
            }
            elsif (/#STOP/) {
                $count = 2;
            }
            elsif (/#SPLICE/) {
                $count = 3;
            }
            elsif (/#CORRECTSPLICESITE/){
            	$count = 4
            }
            elsif ( $count == 1 ) {
                $START_CODON{$_} = 1;
            }
            elsif ( $count == 2 ) {
                $STOP_CODON{$_} = 1;
            }
            elsif ( $count == 3 ) {
                /(\d+)\.\.(\d+)/;

                $SPLICE_DONOR{$_}    = 1;
                $SPLICE_ACCEPTOR{$_} = 1;
            }
            elsif ( $count == 4 ) {
            	chomp;
            	$CORRECT_SPLICESITE = $_ 	
            }
            
        }
    }
	else {
	  print "Using the default specifications for start/codons and splice sites.\n";

	}
}

####################
### loadSequence
####################
sub loadSequence {
    my $file = shift;

    my $seq;
    my $amount;

    open( F, $file )
      or die "Couldn't open Sequence file $file: $!\n";

    while (<F>) {
        if (/^>/) {
            $amount++;
        }
        else {
            chomp;
            $seq .= $_;
        }
    }
    if ( $amount > 1 ) {
        die
		  "Please do not enter a multifasta file, in $file --  but just the sequence file for the used replicon.\n";
    }
    return $seq;
}

####################
### revSeq
####################
sub revcomp {
    my $str = uc(shift);
    $str =~ tr/ATGC/TACG/;

    return reverse($str);
}

####################
### revSeq
####################
sub reverseEMBL {
    my ( $embl, $sequence, $postfix ) = @_;

    my $length = length($sequence);
    my $res;
    open( F, $embl ) or die "Couldn't open EMBL file $embl.\n";

    while (<F>) {
        if (/^FT   CDS\s+(\S+)/) {
            $res .= "FT   CDS             " . reverseCDS( $1, $length ) . "\n";
        }
        else {
            $res .= $_;
        }
    }
    open( F, "> $embl.$postfix" )
      or die "Couldn't write EMBL file $embl.$postfix\n";
    print F $res;
    close(F);

    return ( revcomp($sequence) );
}

sub reverseCDS {
    my $line   = shift;
    my $length = shift;

    chomp($line);
    $line =~ s/>//g;
    $line =~ s/<//g;

    my $wasComplement = 0;
    if ( $line =~ /complement/ ) {
        $wasComplement = 1;
        $line =~ s/complement\(//g;
        $line =~ s/\)$//g;
    }

    my $hadExons = 0;
    if ( $line =~ /join/ ) {
        $hadExons = 1;
        $line =~ s/join\(//g;
        $line =~ s/\)$//g;
    }
    $line =~ s/\.\./_/g;

    my @parts = split( /([,|_])/, $line );

    my $res;
    for ( my $i = ( scalar(@parts) - 1 ) ; $i >= 0 ; $i-- ) {
        if ( $parts[$i] =~ /\d+/ ) {
            $res .= ( $length - $parts[$i] + 1 );
        }
        else {
            $res .= $parts[$i];
        }
    }
    $res =~ s/_/../g;
    if ($hadExons) {
        $res = "join($res)";
    }
    if ( !$wasComplement ) {
        $res = "complement($res)";
    }
    return $res;
}

sub reverseGFF {
    my $gff      = shift;
    my $sequence = shift;

    my $length = length($sequence);
    my $res;

    my @lines = split( /\n/, $gff );
    foreach (@lines) {
        my @ar = split(/\s+/);
        $ar[3] = ( $length - $ar[3] + 1 );
        $ar[4] = ( $length - $ar[4] + 1 );
        $ar[6] = '-';

        $res .= join( "\t", @ar ) . "\n";
    }
    return $res;

}

sub correctEMBL {
    my $embl    = shift;
    my $postfix = shift;

    my $res;
    open( F, $embl ) or die "Couldn't open EMBL file $embl.\n";

    #  my @F=<F>;

    my $getNext = 0;

    while (<F>) {
        chomp;
        if (/^FT   \S+\s+(\S+,)$/) {
            $getNext = 1;
            $res .= $_;
        }
        elsif ($getNext) {
            if (/^FT\s+(.*,)$/) {
                $getNext = 1;
                $res .= $1;
            }
            else {
                /^FT\s+(.*)$/;
                $getNext = 0;
                $res .= $1 . "\n";
            }
        }
        else {
            $res .= $_ . "\n";
        }
    }
    open( F, "> $embl.$postfix" )
      or die "Couldn't write EMBL file $embl.$postfix\n";
    print F $res;
    close(F);
}

sub debug {
    my $val = shift;

    my $str = shift;

    if ( $DEBUG >= $val ) {
        print $str. "\n";

    }

}
