#!/bin/bash
case $HOSTNAME in
  karel-work | karel-macbook | karel-super | karel-differ)
    ARCH=archlinux
    ;;
  rs*)
    ARCH=rsrijnhuizen
    ;;
  edison*)
    ARCH=edison
    ;;
  fusionsim)
    ARCH=fusionsim
    ;;
  freia*)
    ARCH=freia
    ;;
  andromede*)
    ARCH=andromede
    ;;
esac

if [ "$ARCH" == "" ]
then
    echo "Machine not found"
    yno=N
else
    echo Detected machine \'$ARCH\' with root \'$PWD\'. Is this correct? \[Y\/n]
    read yno
    if [ "$yno" == "" ]
    then
      yno=Y
    fi
fi

case $yno in
  Y | y)
    echo Creating and adjusting src/Makefile.inc
    cp $PWD/src/make.inc/Makefile.$ARCH $PWD/src/Makefile.inc
    sed -i '/QUALIKIZ=*/c\QUALIKIZ='$PWD $PWD/src/Makefile.inc
    ;;
  N | n)
    echo "Please manually copy src/make.inc/Makefile.* for your arch to src/Makefile.inc and edit QUALIKIZ and compiler flags"
    exit 1
    ;;
esac
