#!/bin/sh
echo 'Character :'
read origString
echo 'to be replaced with :'
read newString

for i in `grep -rlI -e "${origString}" *`
do
  sed "s/${origString}/${newString}/g" $i ">" $i.__tmp
  echo mv -f $i.__tmp $i
done
