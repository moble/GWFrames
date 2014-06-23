# encoding: utf-8
"""
This file's functions are mostly stolen from IPython.
"""

from __future__ import print_function

import errno
import os
import sys
import subprocess

from distutils.command.build_py import build_py
from distutils.command.build_scripts import build_scripts
from distutils.command.install import install
from distutils.command.install_scripts import install_scripts
from distutils.cmd import Command
from fnmatch import fnmatch
from glob import glob
from subprocess import call

#-------------------------------------------------------------------------------
# Useful globals and utility functions
#-------------------------------------------------------------------------------

# A few handy globals
isfile = os.path.isfile
pjoin = os.path.join
repo_root = os.path.dirname(os.path.abspath(__file__))

# Py3 compatibility hacks, without assuming IPython itself is installed with
# the full py3compat machinery.

try:
    execfile
except NameError:
    def execfile(fname, globs, locs=None):
        locs = locs or globs
        exec(compile(open(fname).read(), fname, "exec"), globs, locs)

def SubmoduleUpdater(submodules):
    class UpdateSubmodules(Command):
        """Update git submodules

        IPython's external javascript dependencies live in a separate repo.
        """
        description = "Update git submodules"
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            failure = False
            try:
                self.spawn('git submodule init'.split())
                self.spawn('git submodule update --recursive'.split())
            except Exception as e:
                failure = e
                print(e)

            if not check_submodule_status(repo_root, submodules) == 'clean':
                print("submodules could not be checked out")
                sys.exit(1)

    return UpdateSubmodules

#-------------------------------------------------------------------------------
# Stolen from IPython.utils.submodule
#-------------------------------------------------------------------------------

def check_submodule_status(root, submodules):
    """check submodule status

    Has three return values:

    'missing' - submodules are absent
    'unclean' - submodules have unstaged changes
    'clean' - all submodules are up to date
    """

    if hasattr(sys, "frozen"):
        # frozen via py2exe or similar, don't bother
        return 'clean'

    if not os.path.exists(pjoin(root, '.git')):
        # not in git, assume clean
        return 'clean'

    for submodule in submodules:
        if not os.path.exists(submodule):
            return 'missing'

    # Popen can't handle unicode cwd on Windows Python 2
    if sys.platform == 'win32' and sys.version_info[0] < 3 \
        and not isinstance(root, bytes):
        root = root.encode(sys.getfilesystemencoding() or 'ascii')
    # check with git submodule status
    proc = subprocess.Popen('git submodule status',
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True,
                            cwd=root,
    )
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")

    for line in status.splitlines():
        if line.startswith('-'):
            return 'missing'
        elif line.startswith('+'):
            return 'unclean'

    return 'clean'

def update_submodules(repo_dir):
    """update submodules in a repo"""
    subprocess.check_call("git submodule init", cwd=repo_dir, shell=True)
    subprocess.check_call("git submodule update --recursive", cwd=repo_dir, shell=True)



#-------------------------------------------------------------------------------
# Make sure we aren't trying to run without submodules
#-------------------------------------------------------------------------------
def require_clean_submodules(repo_root, submodules):
    """Check on git submodules before distutils can do anything

    Since distutils cannot be trusted to update the tree
    after everything has been set in motion,
    this is not a distutils command.
    """
    # PACKAGERS: Add a return here to skip checks for git submodules

    # don't do anything if nothing is actually supposed to happen
    for do_nothing in ('-h', '--help', '--help-commands', 'clean', 'submodule'):
        if do_nothing in sys.argv:
            return

    status = check_submodule_status(repo_root, submodules)

    if status == "missing":
        print("checking out submodules for the first time")
        update_submodules(repo_root)
    elif status == "unclean":
        print('\n'.join([
            "Cannot build / install with unclean submodules.",
            "Please update submodules with",
            "    git submodule init",
            "    git submodule update",
            "in the top-level directory, or commit any",
            "submodule changes you have made."
        ]))
        sys.exit(1)
