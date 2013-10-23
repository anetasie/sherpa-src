#! /bin/sh

# July 2010

# This is the official template for tool regression test scripts.
# In addition to supporting the "SHORTTEST" option, this script also
# allows the user to run individual subtests from the command line.
# The script will accept a series of test identifiers, generally of
# the form "test1" "test2" ... which are to be run.

# Portions of the script which must be customized are marked with "!!",
# below.

# The complete list of tests must be placed in "alltests", !!3, below.
# The test[s] for the SHORTTEST must be placed in "shortlist", !!4 below.


######################################################################
# subroutine
# error_exit <message>
# Fatal error exit

error_exit()
{
  echo "$1" | tee -a $LOGFILE
  echo "${toolname} : FAIL" | tee -a $LOGFILE
  exit 1
}

######################################################################
# subroutine
# keyfilter infile outfile
# filters out CHECKSUM, Dataset, CREATOR, HISTORY, DATASUM,
#             ASCDSVER, HISTNUM, and DATE
# To filter additional keywords, add s/KEYWORD/Dataset/g; for each.

keyfilter()
{
  cat $1 | sed -e 's/CHECKSUM/Dataset/g;s/COMMENT/Dataset/g;
  s/DATE/Dataset/g;s/CREATOR/Dataset/g;s/HISTORY/Dataset/g;
  s/DATASUM/Dataset/g;s/ASCDSVER/Dataset/g;s/HISTNUM/Dataset/g' | \
  grep -v Dataset > $2
  zerotest $2
}

######################################################################
# subroutine
# find_tool <toolname>
# checks that tool exists and is runnable

find_tool()
{
  s1=`type $1`
  s2=`echo $s1 | awk -F" " '{ print $3}'`
  if test -x $s2 ; then
    :
  else
    error_exit "tool $1 not found"
  fi
  punlearn $1
}

######################################################################
# subroutine
# find_tool_py <toolname>.py
# checks that tool exists and is runnable

find_tool_py()
{
  s1=`type $1.py`
  s2=`echo $s1 | awk -F" " '{ print $3}'`
  if test -x $s2 ; then
    :
  else
    error_exit "tool $1 not found"
  fi
  punlearn $1
}

######################################################################
# subroutine
# zerotest <file>
# Makes sure that file is not 0 length.
# Use this to protect yourself against empty files  (which will
# 'diff' without error).  This can happen when the input file to
# cat $infile | do_something >> $outfile
# is missing.  This is used by keyfilter(), above.

zerotest()
{
 if test -s $1 ;
 then
   :
 else
   echo "ERROR: file $1 is of zero length" >> $LOGFILE
   #  Indicate failure, but do not exit.
   mismatch=0
 fi
}


######################################################################
# Initialization

# !!3
toolname="lev3_iss"

# set up list of tests
# !!4
alltests="lev3_iss1"

# "short" test to run
# !!5
shortlist=""


# compute date string for log file
DT=`date +'%d%b%Y_%T'`


# convenience definitions
OUTDIR=$TESTOUT/$toolname
SAVDIR=$TESTSAV/$toolname
INDIR=$TESTIN/$toolname
LOGDIR=$TESTLOG/$toolname

# set up log file name
LOGFILE=$LOGDIR/${toolname}_log.$DT

#get rid of old logs
rm -f $LOGDIR/${toolname}_log.*

# Any tests specified on command line?
if test $# -gt 0; then
  # yes, do those tests
  testlist=$*
else
  # No, see if we are to do "short" test
  if test "x$SHORTTEST" = "x" ; then
    # No, do everything
    testlist=$alltests
  else
    # yes, do short test
    testlist=$shortlist
  fi
fi


# Make sure we have a log directory
if test -d $LOGDIR ; then
 :
else
  mkdir -p $LOGDIR
  if test $? -ne 0 ; then
    error_exit ""
  fi
fi


# Make sure we have an output directory
if test -d $OUTDIR ; then
 :
else
  mkdir -p $OUTDIR >> $LOGFILE 2>&1
  if test $? -ne 0 ; then
    error_exit "can't create output directory $OUTDIR"
  fi
fi

# check for directory environment variables
if test "x${TESTIN}" = "x" -o "x${TESTOUT}" = "x" -o "x${TESTSAV}" = "x" \
   -o "x${TESTLOG}" = "x" ; then
  error_exit "one or more of TESTIN/TESTOUT/TESTSAV/TESTLOG not defined"
fi


# check for tools
# if a utility is used in the form "utility <args> > outfile", and 'utility'
# cannot be run, 'outfile' will still be created.  If utility is used on
# both the output and reference files of a tool the resultant utility output
# files will both exist and be empty, and will pass a diff.


# announce ourselves
echo ""
echo "${toolname} regression" | tee $LOGFILE
echo ""

# All parameters except verbose should be set anyway, but clear them
# to be safe.
#bose=`pget $toolname verbose`
punlearn $toolname
#pset $toolname verbose=$bose

script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do

  # delete old outputs
  rm -f $OUTDIR/${toolname}_${testid}.param.txt

  # Set up file names
  srcfile="$INDIR/acisf00635_000N001_r0001b_imho3.par"
  psffile="$INDIR/acisf00635_000N001_r0001b_pmho3.par"
  #outfile="$OUTDIR/acisf00635_000N001_r0001b_iss3.par"
  outfile="$OUTDIR/${toolname}_${testid}.param.txt"

  echo "" | tee -a $LOGFILE
  echo "Running test:$testid" | tee -a $LOGFILE
  echo "" | tee -a $LOGFILE
  ####################################################################
  # run the tool
  case ${testid} in
    # !!1
      lev3_iss1 ) test_string="lev3_iss.py srcfile=$srcfile psffile=$psffile \
              outfile=$outfile"
      	  ;;

  esac

  echo $test_string | tee -a $LOGFILE
  eval $test_string | tee -a  $LOGFILE  2>&1                #Original
  #eval pdump $toolname >> $parfile
  ####################################################################
  # check the outputs

  # Init per-test error flag
  mismatch=1

#    diff $OUTDIR/${toolname}_${testid}.txt $SAVDIR/${toolname}_${testid}.ref.txt > \
#          /dev/null 2>> ${LOGFILE}

#    if  test $? -ne 0 ; then
#      echo "ERROR: TEXT MISMATCH in $OUTDIR/${toolname}_${testid}.txt" >> $LOGFILE
#      mismatch=0
#    fi

   diff $OUTDIR/${toolname}_${testid}.param.txt $SAVDIR/${toolname}_${testid}.param.ref.txt > \
       /dev/null 2>> ${LOGFILE}

   if  test $? -ne 0 ; then
     echo "ERROR: TEXT MISMATCH in $OUTDIR/${toolname}_${testid}.par" >> $LOGFILE
     mismatch=0
   fi


  ####################################################################
  # Did we get an error?
  echo ""
  if test $mismatch -eq 0 ; then
    # Yes
    echo "${testid} NOT-OK"
    script_succeeded=1
  else
    # No
    echo "${testid} OK"
  fi

done
# end per-test loop
######################################################################


######################################################################
# report results

# blank line
echo ""

if test $script_succeeded -eq 0; then
    echo "${toolname} : PASS" | tee -a $LOGFILE
else
    echo "${toolname} : FAIL" | tee -a $LOGFILE
fi

echo "log file in ${LOGFILE}"


exit $script_succeeded
