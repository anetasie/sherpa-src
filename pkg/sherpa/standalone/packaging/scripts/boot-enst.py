# Copyright (c) 2009-2011 by Enthought, Inc.
# All rights reserved.
"""\
This code is used to bootstrap Enstaller (into sys.prefix of the Python
running this script).  The functionality is implemented in Enstaller itself.

Install (bootstrap) Enstaller:
$ python custom_tools/boot-enst.py path/to/Enstaller...egg

Uninstall (unbootstrap) Enstaller:
$ python custom_tools/boot-enst.py --remove
"""
import sys


def install():
    sys.path.insert(0, sys.argv[1])
    from egginst.bootstrap import main
    sys.exit(main())

def remove():
    import egginst
    ei = egginst.EggInst('Enstaller', sys.prefix)
    ei.remove()
    sys.exit(0)

def help():
    print __doc__
    sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        help()

    if sys.argv[1] == '--remove':
        remove()

    if sys.argv[1].endswith('.egg'):
        install()

    help()
