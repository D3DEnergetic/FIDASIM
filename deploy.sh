#!/bin/bash
# This script lets Travis CI deploy the FORD generated documentation website

set -e

if [ ! "$TRAVIS" ]; then
    echo "Documentation can only be deployed by Travis CI"
    exit 0
fi

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    echo "Skipping documentation deployment"
    exit 0
fi

REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
echo $REPO
echo $SSH_REPO

openssl aes-256-cbc -K $encrypted_60fa3eb19c83_key -iv $encrypted_60fa3eb19c83_iv -in .travis_key.enc -out travis_key -d

if [ ! -f $TRAVIS_BUILD_DIR/travis_key ]; then
    echo "Missing Travis Deploy Key"
    exit 1
fi

chmod 600 travis_key
eval `ssh-agent -s`
ssh-add travis_key

git clone --branch=gh-pages --depth 1 $REPO gh-pages

if [[ $TRAVIS_BRANCH = release-* ]]; then
    cd gh-pages
    rm -rf css favicon.png fonts index.html interface \
       js lists media module page proc program search.html \
       sourcefile src tipuesearch type
    cp -r "$TRAVIS_BUILD_DIR"/docs/html/* .
    git add -A
    git commit -m "Development Documentation updated by travis job $TRAVIS_JOB_NUMBER for commits $TRAVIS_COMMIT_RANGE"
    git push $SSH_REPO gh-pages
fi

if [[ $TRAVIS_BRANCH == "master" ]]; then
    cd gh-pages
    rm -rf $TRAVIS_BRANCH
    mkdir $TRAVIS_BRANCH
    cp -r "$TRAVIS_BUILD_DIR"/docs/html/* $TRAVIS_BRANCH
    git add -A
    git commit -m "$TRAVIS_BRANCH documentation updated by travis job $TRAVIS_JOB_NUMBER for commits $TRAVIS_COMMIT_RANGE"
    git push $SSH_REPO gh-pages
fi
