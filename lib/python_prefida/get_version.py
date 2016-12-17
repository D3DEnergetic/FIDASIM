#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import platform


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

    if platform.system() == 'Windows':
        alt = True
    else:
        # Location of .git folder
        git_dir = r'{}{}.git'.format(fidasim_dir, os.path.sep)

        # git is installed if git_file is a file
        proc = subprocess.Popen('command -v git', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        git_file = proc.communicate()[0]
        git_file = git_file.replace('\n', '')

        # Check that .git folder is present and git is installed
        if os.path.isfile(git_file) and os.path.isdir(git_dir):
            try:
                version = subprocess.check_output(['git', '--git-dir={}'.format(git_dir), 'describe', '--tags', '--always', '--dirty'])
                version = version.replace('\n', '')
            except:
                alt = True
        else:
            alt = True

    # If above didn't work, read version file
    if alt:
        # Git 'version' filepath
        ver_file = '{}{}VERSION'.format(fidasim_dir, os.path.sep)

        if os.path.isfile(ver_file):
            with open(ver_file) as f:
                version = f.read()

    return version
