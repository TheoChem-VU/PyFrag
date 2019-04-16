#!/bin/sh

cd $HOME

if  [ -e "conda -h" ] >/dev/null 2>&1; then
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh || exit 1
  echo "download finished"
  ./miniconda.sh -b -p $HOME/miniconda || exit 1
  echo "conda figure finished"
fi


BASH_RC="$HOME"/.bashrc
cat <<EOF >> "$BASH_RC"
# added by conda installer
export PATH="$HOME/miniconda/bin:\$PATH"
export HOSTPYFRAG="$HOME/bin/host"
export QMWORKS="$HOME/miniconda/envs/qmworks"
export USERNAME="$USER"
# added by conda installer
EOF

source "$HOME"/.bashrc

conda create -n qmworks python=3.5        || exit 1

source activate qmworks

conda install h5py==2.6.0 || exit 1
conda install -c rdkit rdkit==2018.03.4.0  || exit 1


pip install --no-cache-dir git+https://github.com/sunxb05/PyFrag@master#egg=qmworks-0.0.1    || exit 1


curl -L -o pyfrag.zip  https://github.com/sunxb05/PyFrag/zipball/master/


PREZIP="$HOME/pyfrag.zip"

unzip "$PREZIP"

rm -f $PREZIP

mv "$HOME"/sunxb05-PyFrag*/host $HOME/bin

rm -r "$HOME"/sunxb05-PyFrag*
