#!/bin/bash
# This script lets Travis CI deploy the FORD generated documentation website

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
if [ ! -f $TRAVIS_BUILD_DIR/travis_key ]; then
    echo "Missing Travis Key"
    exit 0
fi

chmod 600 travis_key
eval `ssh-agent -s`
ssh-add travis_key

if [ "$TRAVIS_BRANCH" = "master" ] && \
   [ "(ls -A $TRAVIS_BUILD_DIR/docs/html)" ]; then
    git clone --branch=gh-pages $REPO gh-pages
    cd gh-pages
    rm -rf css favicon.png fonts index.html interface \
       js lists media module page proc program search.html \
       sourcefile src tipuesearch type
    cp -r "$TRAVIS_BUILD_DIR"/docs/html/* .
    git add -A
    git commit -m "Development Documentation updated by travis job $TRAVIS_JOB_NUMBER for commits $TRAVIS_COMMIT_RANGE"
    git push $SSH_REPO gh-pages
fi

if [[ $TRAVIS_BRANCH == release-* ]]; then
    git clone --branch=gh-pages $REPO gh-pages
    cd gh-pages
    rm -rf $TRAVIS_BRANCH
    mkdir $TRAVIS_BRANCH
    cp -r "$TRAVIS_BUILD_DIR"/docs/html/* $TRAVIS_BRANCH
    git add -A
    git commit -m "$TRAVIS_BRANCH documentation updated by travis job $TRAVIS_JOB_NUMBER for commits $TRAVIS_COMMIT_RANGE"
    git push $SSH_REPO gh-pages
fi
