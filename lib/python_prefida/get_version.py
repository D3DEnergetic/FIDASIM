#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess


def get_version(fidasim_dir):
    #+#get_version
    #+ Gets FIDASIM version number from git.
    #+ Falls back to reading VERSION file when git is not available
    #+***
    #+##Input Arguments
    #+    **fidasim_dir**: FIDASIM install directory
    #+
    #+##Output Arguments
    #+     **version**: FIDAIM version number.
    #+
    #+##Example Usage
    #+```python
    #+>>> version = get_version(get_fidasim_dir())
    #+```
    version = ''
    alt = False

    # Location of .git folder
    git_dir = r'{}{}.git'.format(fidasim_dir, os.path.sep)

    # Check that git is installed

    # Check that .git folder is present
    if os.path.isdir(git_dir):
        try:
            subprocess.check_output(['git', '--git-dir={}'.format(git_dir), 'describe', '--tags', '--always', '--dirty'])
        except:
            alt = True
    else:
        alt = True

    if alt:
        # Git 'version' filepath
        ver_file = '{}{}VERSION'.format(fidasim_dir, os.path.sep)

        if os.path.isfile(ver_file):
            with open(ver_file) as f:
                version = f.read()

    return version
