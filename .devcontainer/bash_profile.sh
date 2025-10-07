#!/bin/bash

DISCLAIMER="""
By running commands in this container, you automatically agree to the Terms of Services (ToS)
of various anaconda channels. To review them, run the command : conda tos view. To suppress
this message, either :
- set the environment variable NFNEURO_ACCEPT_CONDA_TOS to true
- switch to a different VSCode terminal profile
- configure the settings under Terminal › Integrated › Default Profile › Linux
"""

if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi

if [ ! $NFNEURO_ACCEPT_CONDA_TOS ];
then

echo "$DISCLAIMER"

fi
