#!/bin/sh

cd $HOME

BASH_RC="$HOME"/.bashrc
BASH_FILE="$HOME"/.profile
BASH_PROFILE="$HOME"/.bash_profile

if [ -f "$BASH_RC" ]; then
    printf "\\n"
    printf "Initializing PyFrag in %s\\n" "$BASH_RC"
    printf "\\n"

cat <<EOF >> "$BASH_RC"
# added by PyFrag installer
export PATH="$HOME/miniconda/bin:\$PATH"
export HOSTPYFRAG="$HOME/pyfrag"
export PATH="\$HOSTPYFRAG/bin:\$PATH"
export QMWORKS="$HOME/miniconda/envs/qmworks"
export USERNAME="$USER"
# added by PyFrag installer
EOF
source "$HOME"/.bashrc

elif [ -f "$BASH_FILE" ]; then
    printf "\\n"
    printf "Initializing PyFrag in %s\\n" "$BASH_FILE"
    printf "\\n"
cat <<EOF >> "BASH_FILE"
# added by PyFrag installer
export PATH="$HOME/miniconda/bin:\$PATH"
export HOSTPYFRAG="$HOME/pyfrag"
export PATH="\$HOSTPYFRAG/bin:\$PATH"
export QMWORKS="$HOME/miniconda/envs/qmworks"
export USERNAME="$USER"
# added by PyFrag installer
EOF
source "$HOME"/.profile

else
    printf "\\n"
    printf "Initializing PyFrag in %s\\n" "$BASH_PROFILE"
    BASHSET="$BASH_PROFIEL"
cat <<EOF >> "$BASH_PROFILE"
# added by PyFrag installer
export PATH="$HOME/miniconda/bin:\$PATH"
export HOSTPYFRAG="$HOME/pyfrag"
export PATH="\$HOSTPYFRAG/bin:\$PATH"
export QMWORKS="$HOME/miniconda/envs/qmworks"
export USERNAME="$USER"
# added by PyFrag installer
EOF

source "$HOME"/.bash_profile
fi


curl -L -o pyfrag.zip  https://github.com/TheoChem-VU/PyFrag/zipball/master/
PREZIP="$HOME/pyfrag.zip"
unzip "$PREZIP"
rm -rf $PREZIP
mv "$HOME"/TheoChem-VU-PyFrag*/host $HOME/pyfrag
chmod +x $HOME/pyfrag/bin/*
rm -r "$HOME"/TheoChem-VU-PyFrag*
