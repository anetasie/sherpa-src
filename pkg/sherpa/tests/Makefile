# Makefile for Sherpa smoke test suite
MK_TOP = ../../..
include $(MK_TOP)/include/Makefile.cxc

TEST_SCRIPTS = sherpa-smoke001.sh
TEST_FILES   = $(TEST_SCRIPTS)
TEST_TARGS   = Makefile.test.sherpa

# Making use of install macros from Makefile.cxc to 
# install the smoke test Makefile and scripts to their
# respective homes.
INSTALL_TEST_SCRIPTS = $(TEST_FILES:%=$(TST_DIR)/smoke/bin/%)
INSTALL_TEST_BASELINES=$(TEST_TARGS:%=$(TST_DIR)/smoke/%)

include $(MK_TOP)/Makefile.all

# DO NOT DELETE THIS LINE -- make depend depends on it.

announce1:
