#!/bin/sh

cd $HOME

BASH_RC="$HOME"/.bashrc
BASH_FILE="$HOME"/.profile
BASH_PROFIEL="$HOME"/.bash_profile

if [ -f "$BASH_RC" ]; then
    printf "\\n"
    printf "Initializing PyFrag in %s\\n" "$BASH_RC"
    printf "\\n"
    BASHSET="$BASH_RC"
elif [ -f "$BASH_FILE" ]; then
    printf "\\n"
    printf "Initializing PyFrag in %s\\n" "$BASH_FILE"
    printf "\\n"
    BASHSET="$BASH_FILE"
else
    printf "\\n"
    printf "Initializing PyFrag in %s\\n" "$BASH_PROFIEL"
    BASHSET="$BASH_PROFIEL"
fi

cat <<EOF >> "$BASH_SET"
# added by conda installer
export HOSTPYFRAG="$HOME/pyfrag"
export PATH="\$HOSTPYFRAG/bin:\$PATH"
# added by conda installer
EOF

source "$HOME"/.bashrc
source "$HOME"/.profile
source "$HOME"/.bash_profile


curl -L -o pyfrag.zip  https://github.com/sunxb05/PyFrag/zipball/master/
PREZIP="$HOME/pyfrag.zip"
unzip "$PREZIP"
rm -rf $PREZIP
mv "$HOME"/sunxb05-PyFrag*/host $HOME/pyfrag
rm -r "$HOME"/sunxb05-PyFrag*
