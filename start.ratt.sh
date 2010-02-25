#!/bin/bash


refembl=$1
query=$2
result=$3
parameterSet=$4
ref=$5

## check the entrance
if [ -z "$parameterSet" ]; then 
	echo "Please use RATT with the following options:
  $RATT_HOME/start.ratt.sh <Directory with embl-files> <Query-fasta sequence> <Resultname> <Transfer type> <optional: reference (multi) Fasta>

  
  Directory name with 
  embl-annotation files  - This directory contains all the embl files that should be transfered to the query.
  Query.fasta            - A multifasta file to, which the annotation will be mapped.
  ResultName             - The prefix you wish to give to each result file.
  Transfer type          - Following parameters can be used (see below for the different used sets)
       (i)   Assembly:             Transfer between different assemblies. 
       (ii)  Assembly.Repetitive:  As before, but the genome is extremely repetitive. 
                   This should be run, only if the parameter Assembly doesn't return good results (misses too many annotation tags).  
       (iii) Strain:              Transfer between strains. Similarity is between 95-99%.
       (iv)  Strain.Repetitive:   As before, but the genome is extremely repetitive. 
                   This should be run, only if the parameter Strain doesn't return good results (misses too many annotation tags).
       (v)   Species:              Transfer between species. Similarity is between 50-94%.
       (vi)  Species.Repetitive:   As before, but the genome is extremely repetitive. 
                   This should be run, only if the parameter Species doesn't return good results (misses too many annotation tags).
       (vii) Multiple:             When many annotated strains are used as a reference, and you assume the newly sequenced genome has many insertions
                   compared to the strains in the query (reference?). This parameter will use the best regions of each reference strain to transfer tags.    
       (viii)Free:                 The user sets all parameter individually.

  reference fasta        - Name of multi-fasta. VERY I M P O R T A N T The name of each sequence in the fasta description, 
                   MUST be the same name as its corresponding embl file. So if your embl file is call Tuberculosis.embl, in your reference.fasta file, 
                   the description has to be 
                           >Tuberculsosis
                           ATTGCGTACG
                           ..."

   exit
 fi


### nucmer call as function
function doNucmer {
	if [ -f "$verbose" ] 
		then 
		echo "nucmer  $other_nucmer -g $g -p $name -c $c -l $l $ref $query"
		echo "delta-filter $rearrange -i $minInd $name.delta > $name.filter.delta"
	fi

	nucmer $other_nucmer -g $g -p $name -c $c -l $l $ref $query &> /dev/null
	delta-filter $rearrange -i $minInd $name.delta > $name.filter.delta
	show-snps -CHTr $name.filter.delta > $name.snp
	show-coords -clHT $name.delta > $name.coords
	show-coords -clHT $name.filter.delta > $name.filter.coords
}

### check path of nucmer and the perl transferprogram
NUCMER_EXE=${NUCMER_EXE-`which nucmer 2>/dev/null`}
EMBL_files=${EMBL_files-`ls $refembl* 2>/dev/null`}
#PERL_SCRIPT=${NUCMER_EXE-`which samtools 2>/dev/null`}
if (test "$NUCMER_EXE" == "");
then
   echo "nucmer is not on the PATH"
   exit; 
elif(test "$EMBL_files" == "");
then
   echo "Cannot find any EMBL files in $refembl is not on the PATH."
   echo "So far just embl formatted files can be used. They must end with .embl."
   exit; 
fi

# check for the query
if [ ! -f "$query" ]                # be sure the directory /mnt exists
  then
    echo "Query $query doesn't exist"
	exit
fi


# check for the reference
if [ -f "$ref" ] 
	then 
	# if the reference exist, use this
	echo "I am using the reference $ref. Please make sure that the description line of each fasta entry is the same than in the embl file name!";
	head -n 1 $ref
	echo "should be the same name as the embl file in embl"
	ls $refembl
	echo
	echo
else
	perl $RATT_HOME/main.ratt.pl Embl2Fasta $refembl Reference.fasta
	ref="Reference.fasta"
fi
	
### check if files ok
if [ ! -f "$ref" ]                # be sure the directory /mnt exists
  then
    echo "Sorry the reference file wasn't generated corretly."
	exit
fi;

name=nucmer.$result

### paramters to set for the nucmer
# l - "k-mer size" 10 short 20 assembly set
# c - cluster size
# g - extending of cluster

# rearrange = set for delta-filter -1    1-to-1 alignment allowing for rearrangements
#                                  -g    1-to-1 global alignment not allowing rearrangements



doneDifference=0;

# minInd - minimal indentity
# minLength - minlength of a hit   -- not used so far...

if [ "$parameterSet" == "Assembly" ] || [ "$parameterSet" == "Assembly.Repetitive" ] ;
	then 
	c=400;
	l=30;
	g=1000;
	
	if [ "$parameterSet" == "Assembly.Repetitive" ] ;
		then
		other_nucmer=" --maxmatch "
	else
		other_nucmer="  "
	fi
	rearrange=" -r -o 10 ";
	minInd=99;
	
	### get real SNP before mutate
	doNucmer

	perl $RATT_HOME/main.ratt.pl Difference  $name.snp $name.filter.coords $result
	doneDifference=1;

	### insert of mutation to have better anchors
	perl $RATT_HOME/main.ratt.pl Mutate $query
	return=$?
	if [ "$return" != "0" ] ;
		then 
		echo "See Error in BBA.main script, to mutate the query for Assembly to Assembly annotation transfer."
		exit 1;
	fi;
	### update name of query
	query=$query."mutated"

elif [ "$parameterSet" == "Strain" ] ||  [ "$parameterSet" == "Strain.Repetitive" ] ;
	then 
	c=400;
	l=20;
	g=500;

	if [ "$parameterSet" == "Strain.Repetitive" ] ;
		then
		other_nucmer=" --maxmatch "
	else
		other_nucmer="  "
	fi
	
	rearrange=" -r -o 10 ";
	minInd=95;
	
	### get real SNP before mutate
	doNucmer

	perl $RATT_HOME/main.ratt.pl Difference  $name.snp $name.filter.coords $result
	doneDifference=1;

	### insert of mutation to have better anchors
	perl $RATT_HOME/main.ratt.pl Mutate $query
	return=$?
	if [ "$return" != "0" ] ;
		then 
		echo "See Error in BBA.main script, to mutate the query for Assembly to Assembly annotation transfer."
		exit 1;
	fi;
	### update name of query
	query=$query."mutated"


elif [ "$parameterSet" == "Species" ]  || [ "$parameterSet" == "Species.Repetitive" ] ;
	then 
	if [ "$parameterSet" == "Species.Repetitive" ] ;
		then
		other_nucmer=" --maxmatch "
	else
		other_nucmer="  "
	fi

	c=400;
	l=10;
	g=500;
    other_nucmer="--maxmatch -o 1"
	rearrange="-q";
	minInd=40;

elif [ "$parameterSet" == "Multiple" ] ;
	then 
	c=400;
	l=25;
	g=1000;
	
	other_nucmer=" --maxmatch "
	rearrange=" -q -o 1";
	minInd=98;
elif [ "$parameterSet" == "Free" ] ;
	then 
	c=$RATT_c;
	l=$RATT_l;
	g=$RATT_g;
	rearrange=$RATT_rearrange;
	minInd=$RATT_minInd;
	other_nucmer=$RATT_anchor;
else
	echo "Plese set: Transfer type: Assembly / Strain / Free "
	exit
fi

### do the comparison using nucmer
doNucmer

# print the synteny stats and do the difference files.
if [ "$doneDifference" == "0" ]
	then 
	perl $RATT_HOME/main.ratt.pl Difference  $name.snp $name.filter.coords $result
fi

### do the tranfer
if [ -f "$verbose" ] 
	then
	echo "Nucmer is done. Now transfer the annotation."
	echo "perl $RATT_HOME/main.ratt.pl $refembl $name.snp $name.filter.coords $result"
fi

perl $RATT_HOME/main.ratt.pl Transfer $refembl $name.snp $name.filter.coords $result

### do ther correction
if [ -f "$verbose" ] 
	then
	echo "Nucmer is done. Now transfer the annotation."
	echo "perl $RATT_HOME/main.ratt.pl $refembl $name.snp $name.filter.coords $result"
fi

for name in `grep '>' $query | perl -nle 's/\|/_/g;/>(\S+)/; print $1'` ; do
	perl $RATT_HOME/main.ratt.pl Correct $result.$name.embl $query $name
done


echo "If you want to start art (assume just one replicon):"
x=$(ls $result*final.embl);
file=$(echo $x | sed 's/\.embl$//g';)
echo "art $query + $x + Query/$file.Mutations.gff"

### clean the files
if [ ! -f "$verbose" ]
	then 
	rm $result*embl.tmp.BBA.embl
	rm $name*
fi