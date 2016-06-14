#!/bin/bash
case $HOSTNAME in
  karel-work | karel-macbook | karel-super)
    ARCH=archlinux
    ;;
  shannon)
    ARCH=shannon
    ;;
  rs*)
    ARCH=rsrijnhuizen
    ;;
  edison*)
    ARCH=edison
    ;;
esac

echo Detected machine \'$ARCH\' with root \'$PWD\'. Is this correct? \[Y\/n]
read yno
if [ "$yno" == "" ]
then
  yno=Y
fi

case $yno in
  Y | y)
    echo Creating and adjusting src/Makefile.inc
    cp $PWD/src/make.inc/Makefile.$ARCH $PWD/src/Makefile.inc
    sed -i '/QUALIKIZ*/c\QUALIKIZ='$PWD $PWD/src/Makefile.inc
    ;;
  N | n)
    echo "Please manually copy src/make.inc/Makefile.* for your arch to src/Makefile.inc and edit QUALIKIZ and compiler flags"
    exit 1
    ;;
esac
