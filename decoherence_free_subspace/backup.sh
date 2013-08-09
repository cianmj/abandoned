#!/bin/sh
rm *.tar
name=`date | awk '{print $2_$3}'`
tar -cf "$name".tar *.f90 Makefile *.sh r_matrix.*
echo "Files:" `tar -tvf "$name".tar | awk '{print $6}'`
echo "Saved to "$name".tar"
cp "$name".tar temp/
