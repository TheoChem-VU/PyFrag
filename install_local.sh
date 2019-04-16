#!/bin/sh

THIS_DIR=$(DIRNAME=$(dirname "$0"); cd "$DIRNAME"; pwd)
THIS_FILE=$(basename "$0")
THIS_PATH="$THIS_DIR/$THIS_FILE"
PREFIX=$HOME/bin


USAGE_ADF="
set up should look like this in you bashrc file
########ADF

# export ADFHOME="/Applications/ADF2017.app/Contents/Resources/adfhome"         # here ADFHOME is the directory where you unpacked the downloaded files
# export ADFBIN="/Applications/ADF2017.app/Contents/Resources/adfhome/bin"             # the directory for run scripts and program binaries
# export ADFRESOURCES="/Applications/ADF2017.app/Contents/Resources/adfhome/atomicdata" # the location of the ADF database
# export SCMLICENSE="/Library/Application_Support/SCM/license.txt"      # the location of your SCM license file (when you get it)
# export PATH=\$ADFBIN:\$PATH

########ADF
For more information, please check https://www.scm.com/doc/Installation/Installation.html?highlight=installation#set-up-the-environment
"

if  ! [ -e "$ADFBIN" ] >/dev/null 2>&1; then
    printf "\\nWARNING: ADF does not appear to be installed or properly set up, you need resolve this problem to continue" >&2
    printf "%s \\n"  "$USAGE_ADF"  >&2
    exit 1
fi



USAGE_Python3="
The easy way to install python3 is to download python3 package from
https://www.python.org/downloads/release/python-368/
"

if  ! [ -e "python3 -h" ] >/dev/null 2>&1; then
    printf "\\nWARNING: python3 does not appear to be installed or properly set up, you need resolve this problem to continue" >&2
    printf "%s \\n"  "$USAGE_Python3"  >&2
    exit 1
fi


if  ! [ -e "pip3 -h" ] >/dev/null 2>&1; then
    printf "\\nWARNING: Pip3 does not appear to be installed or properly set up, you need resolve this problem to continue" >&2
    exit 1
fi


xcode-select --install


pip3 install bokeh
pip3 install opencv-python


export PATH="/System/Library/Frameworks/Python.framework/Versions/2.7/bin:$PATH"
sudo easy_install pip &>/dev/null
sudo pip install Fabric==1.12.2 || exit 1



sudo apachectl start &>/dev/null
sudo cp /etc/apache2/httpd.conf /etc/apache2/httpd.conf_original
script=/etc/apache2/httpd.conf
coordinates=$HOME/Sites
coordinate="LoadModule php7_module libexec/apache2/libphp7.so"
sudo awk -v r="$coordinates" '{sub("/Library/WebServer/Documents",r,$0); print}' $script > tmpfile
sudo awk -v r="$coordinate" '{sub("#LoadModule php7_module libexec/apache2/libphp7.so",r,$0); print}' tmpfile > inputfile

sudo mv inputfile /etc/apache2/httpd.conf
sudo rm tmpfile
sudo rm inputfile
sudo apachectl restart


if [ -e "$HOME/Sites" ]; then
    cd "$HOME/Sites"
    mkdir "$HOME/Sites/video"
else
    mkdir "$HOME/Sites"
    cd "$HOME/Sites"
    mkdir "$HOME/Sites/video"
fi


printf "PyFrag will now be installed into this location:\\n"
printf "%s\\n" "$PREFIX"
printf "\\n"
printf "  - Press ENTER to confirm the location\\n"
printf "  - Press CTRL-C to abort the installation\\n"
printf "  - Or specify a different location below\\n"
printf "\\n"
printf "[%s] >>> " "$PREFIX"
read -r user_prefix
if [ "$user_prefix" != "" ]; then
    case "$user_prefix" in
        *\ * )
            printf "ERROR: Cannot install into directories with spaces\\n" >&2
            exit 1
            ;;
        *)
            eval PREFIX="$user_prefix"
            ;;
    esac
fi


if [ -e "$PREFIX" ]; then
    printf "ERROR: File or directory already exists: '%s'\\n" "$PREFIX" >&2
    printf "If you want to complish pyfrag installation, you need to remove this existing directory\\n" >&2
    exit 1
fi


if ! mkdir -p "$PREFIX"; then
    printf "ERROR: Could not create directory: '%s'\\n" "$PREFIX" >&2
    exit 1
fi

PREFIX=$(cd "$PREFIX"; pwd)
export PREFIX

cd "$PREFIX"

if ! curl --help >/dev/null 2>&1; then
    printf "WARNING: curl does not appear to be installed, you need install it to continue\\n" >&2
    exit 1
fi

curl -L -o pyfrag.zip  https://github.com/sunxb05/PyFrag/zipball/master/


PREZIP="$PREFIX/pyfrag.zip"

unzip "$PREZIP"

rm -f $PREZIP

mv "$PREFIX"/sunxb05-PyFrag* "$PREFIX"/pyfrag


PYFRAG_RC="$PREFIX"/pyfrag/.pyfragrc
CONFIG_FILE="$PREFIX"/pyfrag/utils/configure.py

printf "Please give your remote host machine host name, like cartesius.surfsara.nl\\n"
read -r hostname

printf "Please give your remote host machine user name, like xnt453\\n"
read -r username

printf "Please give your remote host machine home directory, like /home/xnt453\\n"
read -r hostdirectory

cat <<EOF >> "$PYFRAG_RC"
# added by PyFrag installer
export REMOTEBASE="$hostdirectory"
# added by PyFrag installer
EOF

cat <<EOF >> "$CONFIG_FILE"
# added by PyFrag installer
USERNAME = "$username"
HOSTNAME = "$hostname"
# added by PyFrag installer
EOF


if True; then
    BASH_RC="$HOME"/.bashrc
    BASH_FILE="$HOME"/.profile
    DEFAULT=no
    printf "Do you wish the installer to initialize PyFrag\\n"
    printf "in your %s ? [yes|no]\\n" "$BASH_RC"
    printf "or in your %s ? [yes|no]\\n" "$BASH_FILE"
    read -r ans
    if [ "$ans" = "" ]; then
        ans=$DEFAULT
    fi
    if [ "$ans" != "yes" ] && [ "$ans" != "Yes" ] && [ "$ans" != "YES" ] && \
       [ "$ans" != "y" ]   && [ "$ans" != "Y" ]
    then
        printf "\\n"
        printf "You may wish to edit your $BASH_RC to setup PyFrag:\\n"
        printf "\\n"
        printf "export PATH=\"$PREFIX/bin:\$PATH\"\\n"
        printf "\\n"
    else
        if [ -f "$BASH_RC" ]; then
            printf "\\n"
            printf "Initializing PyFrag in %s\\n" "$BASH_RC"
            printf "\\n"
        elif [ -f "$BASH_FILE" ]; then
            printf "\\n"
            printf "Initializing PyFrag in %s\\n" "$BASH_FILE"
            printf "\\n"
        else
            printf "\\n"
            printf "Initializing PyFrag in newly created %s\\n" "$BASH_RC"
        fi
        cat <<EOF >> "$BASH_RC"
# added by PyFrag installer
export PYFRAGHOME="$PREFIX/pyfrag"
export PATH="\$PYFRAGHOME/bin:\$PATH"
# added by PyFrag installer
EOF
        printf "\\n"
        printf "For this change to become active, you have to open a new terminal.\\n"
        printf "\\n"
    fi

    printf "Thank you for installing PyFrag!\\n"
fi # !BATCH


exit 0
