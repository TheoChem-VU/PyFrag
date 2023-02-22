#!/bin/sh

cd $HOME

if  ! which conda >/dev/null 2>&1; then
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh || exit 1
  echo "download finished"
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $HOME/miniconda || exit 1
  echo "conda figure finished"
fi


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
export HOST

PYFRAG="$HOME/pyfrag"
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


conda create -n qmworks python=3.5        || exit 1
source activate qmworks

conda install h5py==2.6.0 || exit 1
conda install -c rdkit rdkit==2018.03.4.0  || exit 1


pip install --no-cache-dir git+https://github.com/TheoChem-VU/PyFrag@master#egg=qmworks-0.0.1    || exit 1


curl -L -o pyfrag.zip  https://github.com/TheoChem-VU/PyFrag/zipball/master/
PREZIP="$HOME/pyfrag.zip"
unzip "$PREZIP"
rm -rf $PREZIP
mv "$HOME"/TheoChem-VU-PyFrag*/host $HOME/pyfrag
rm -r "$HOME"/TheoChem-VU-PyFrag*
