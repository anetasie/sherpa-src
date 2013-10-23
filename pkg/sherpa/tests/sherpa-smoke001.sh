#! /bin/sh

TESTID=sherpa-smoke001
export TESTID

INFILE=''
OUTFILE=''
OUTLIST='${OUTFILE}'
TOLFILE=''
#
# Synopsis: Sherpa Self Test
TESTCMD='python -c "import sherpa; sherpa.test()"'

# Setup CIAO environment
#   (replace with variable via configure)
#
ASCDS_QUIET=1; export ASCDS_QUIET
ASCDS_OVERRIDE=1; export ASCDS_OVERRIDE
. setup.sh
  
#
# Setup test env
#
TESTDIR=$ASCDS_INSTALL/test/smoke/data; export TESTDIR
OUTDIR=$ASCDS_WORK_PATH/smoke.${LOGNAME}/$TESTID; export OUTDIR
SAVEDIR=$TESTDIR; export SAVEDIR

INDIR=$TESTDIR

#
# --!!  check $ASCDS_WORK_PATH exist, writeable, has enough space
#
# Setup local parameter location
PFILES="$OUTDIR/param;$ASCDS_INSTALL/param"; export PFILES

if test -d $OUTDIR
then
  \rm -rf $OUTDIR
fi

mkdir -p $OUTDIR/param

#
# Run command
#
printf "Running test $TESTID"
printf " ."

fail=1
(eval "$TESTCMD")  > $OUTDIR/diff.log 2>&1
if test $? != 0
then
  printf "."
  fail=1
else
  failure_found=0
  # Hack -- check output log for FAIL string
  /bin/grep -i "FAIL" $OUTDIR/diff.log > /dev/null 2>&1
  if test $? = 0
  then
    failure_found=1
  fi
  printf "."
  # Hack -- check output log for Error string
  /bin/grep -i "Error" $OUTDIR/diff.log > /dev/null 2>&1
  if test $? = 0
  then
    failure_found=1
  fi
  if test $failure_found = 1
  then
    fail=1
  else
    fail=0
  fi  
fi

#
# Exit condition
#
printf ". "
if test $fail = 0
then
  echo "PASS"
  \rm -rf $OUTDIR
  xx=`\ls $ASCDS_WORK_PATH/smoke.${LOGNAME}`
  if test "x$xx" = "x"
  then
    rmdir $ASCDS_WORK_PATH/smoke.${LOGNAME}
  fi
else
  echo " FAIL"
fi

exit $fail
