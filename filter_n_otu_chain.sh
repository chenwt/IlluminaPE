#!/bin/bash

# ###########################################################################
# Usage:
#      bash filter_n_otu_chain.sh <sample_name>.aligned.composite.gz [STEP]
#
# Requires Qiime (1.3.0)
#
# STEP 0: run all steps
# STEP 1: filter out low qual seqs (has 1 or more base with phred < 20)
# STEP 2: pick_otus' prefix_suffix
# STEP 3: pick_otus' cd-hit
# STEP 4: pick rep set from otu output
# STEP 5: combine the chained otus to .chained_otu.[list|txt]
# STEP 6: remove intermediate files
# ###########################################################################


SRC=$EBB/IlluminaPE
QIIME=~/software_download/Qiime-1.3.0/scripts
INPUT=$1
STEP=$2
INDIR=`dirname $1`
INBASE=`basename $1`
LOG=$1.filter_n_chain_otu.log

if [ $STEP -eq 1 -o $STEP -eq 0 ]; then
	rm $LOG
	cmd="python $SRC/filter_low_count_low_qual_seqs.py $1"
	echo `date` >> $LOG
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG
fi

if [ $STEP -eq 2 -o $STEP -eq 0 ]; then
	cmd="python $QIIME/pick_otus.py -i $INPUT.phred20_passed.unique.fasta -o $INDIR/picked_otus1/ -m prefix_suffix -p 50 -u 50"
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG

	cmd="python $QIIME/pick_rep_set.py -m longest -i $INDIR/picked_otus1/$INBASE.phred20_passed.unique_otus.txt -f $INPUT.phred20_passed.unique.fasta -o $INPUT.phred20_passed.unique.prefix50suffix50otu.fasta"
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG
fi

if [ $STEP -eq 3 -o $STEP -eq 0 ]; then
	cmd="python $QIIME/pick_otus.py -m cdhit -s 0.97 -n 100 -i $INPUT.phred20_passed.unique.prefix50suffix50otu.fasta -o $INDIR/picked_otus2/"
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG
	
	cmd="python $QIIME/pick_rep_set.py -m longest -i $INDIR/picked_otus2/$INBASE.phred20_passed.unique.prefix50suffix50otu_otus.txt -f $INPUT.phred20_passed.unique.prefix50suffix50otu.fasta -o $INPUT.phred20_passed.unique.prefix50suffix50otu.cdhit.fasta"
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG
fi

if [ $STEP -eq 4 -o $STEP -eq 0 ]; then
	cmd="python $SRC/combine_chained_otus.py $INDIR/picked_otus1/$INBASE.phred20_passed.unique_otus.txt $INDIR/picked_otus2/$INBASE.phred20_passed.unique.prefix50suffix50otu_otus.txt $INPUT.chained.tmp"
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG

	cmd="python $SRC/combine_chained_otus.py $INPUT.phred20_passed.unique.count $INPUT.chained.tmp $INPUT.chained_otu.txt"
	echo "Running command: $cmd" >> $LOG
	$cmd
	rm $INPUT.chained.tmp
	echo `date` >> $LOG
fi

if [ $STEP -eq 5 -o $STEP -eq 0 ]; then
	cmd="python $SRC/make_fn_list_based_on_uniquified_otus.py $INPUT.chained_otu.txt $INPUT.chained_otu.list"
	echo "Running command: $cmd" >> $LOG
	$cmd
	echo `date` >> $LOG
fi

if [ $STEP -eq 6 -o $STEP -eq 0 ]; then
	gzip -f $INPUT.phred20_passed.unique.fasta
	gzip -f $INPUT.phred20_passed.unique.prefix50suffix50otu.fasta
	gzip -f $INPUT.phred20_passed.unique.prefix50suffix50otu.cdhit.fasta
	gzip -f $INPUT.phred20_passed.unique.count
	gzip -f $INDIR/picked_otus1/$INBASE.phred20_passed.unique_otus.txt
	gzip -f $INDIR/picked_otus2/$INBASE.phred20_passed.unique.prefix50suffix50otu_otus.txt
	gzip -f $INPUT.chained_otu.txt
fi

