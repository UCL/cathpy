#!/bin/bash

# This is a simple wrapper around `cath-sc-merge-alignment` that merges
# all the structural cluster files for a given CATH version (eg v4.2) and 
# cluster_type (eg FF_SSG5) 
#
# You probably want to check the variables below, then run this with 
# something like:
#
#   nohup scripts/run_all_sc.sh sfam.list > run_all_sc.out >& run_all_sc.err &
#
# (where 'sfam.list' is a list of all sfams you want to include)


if [[ $# -eq 0 ]]; then
	echo "usage: $0 sfam.list"
	exit 1
fi

SFAM_LIST="$1"

# make sure these are okay
CATH_VERSION=v4_2_0
FF_DIR=/cath/data/${CATH_VERSION}/funfam/families/
CLUSTER_TYPE=FF_SSG5
PROC=4
SC_SUFFIX='.aln_reps.cora.fa'
OVERWRITE=0

# shouldn't need to change anything else
JOB_DIR=merge.$$
LOG_DIR=${JOB_DIR}
ALN_DIR=sc-alignments
LOG_TAG=$JOB_DIR
LOGSTASH_PORT=5514

CODE_BASE=/cath/homes2/ucbcisi/git/cathpy
MERGE_SCRIPT=$CODE_BASE/scripts/cath-sc-merge-alignment

export MERGE_SCRIPT
export CATH_VERSION
export JOB_DIR
export OVERWRITE

log(){
	echo $(date) $1
}
export -f log

logevent() {
	printf "%s %-20s %s\n" "$(date)" "$1" "$2"
	#logger -t $JOB_DIR -P 9876 "$1 - $2"
}
export -f logevent

domerge() {
	task_id="$1"
	sc_file="$2"
	out_file="$3"
	log_stem="$4"

	stderr_file="$log_stem.stderr"
	stdout_file="$log_stem.stdout"

	logevent 'START' $task_id
	if [ -f $out_file ] && [ "$OVERWRITE" = "0" ]; then
		logevent 'ALREADY_EXISTS' $out_file
	else
		logevent 'RUN' $task_id
		( "$MERGE_SCRIPT" --cv "$CATH_VERSION" --in $sc_file --out $out_file > $stdout_file >& $stderr_file ) \
			|| logevent 'ERROR' $stderr_file
	fi
	logevent 'END' $task_id
}
export -f domerge

log "SFAM_LIST:     $SFAM_LIST"
log "CATH_VERSION:  $CATH_VERSION"
log "CLUSTER_TYPE:  $CLUSTER_TYPE"
log "FF_DIR:        $FF_DIR"
log "OVERWRITE:     $OVERWRITE"
log "PROCESSES:     $PROC"
log "LOG_DIR:       $LOG_DIR"
log "ALN_DIR:       $ALN_DIR"

mkdir -p $JOB_DIR $LOG_DIR $ALN_DIR

#set -x
find $FF_DIR -name "*${CLUSTER_TYPE}*${SC_SUFFIX}" -print | \
    grep -F -f $SFAM_LIST | \
#	head -n 1 | \
	sed "s#$FF_DIR##" | \
	parallel --line-buffer -j$PROC domerge \
		"{/.}" "$FF_DIR{}" "$ALN_DIR/{/.}.merged.sto" "$LOG_DIR/{/.}"
		# name sc_path out_file log_stub

