#!/bin/bash
rm -rf *.dat_edit
ls *.dat0 > files.txt
for i in $(cat files.txt); do
    echo $i
    sed 's/,/ /g' $i > $i"_1"
    sed -e 's/"//g' $i"_1" > $i"_2"
    sed -e 's/\\//g' $i"_2" > $i"_end"
    wc -l $i"_end" > num.txt
    nl=$(awk '{print $1}' num.txt)
    nl=$(($nl+1))
    rm -rf num.txt
    for (( j=1;j<$nl;j++ )); do
	val1=$(awk 'NR=='$j'{print $1}' $i"_end")
	val2=$(awk 'NR=='$j'{print $2}' $i"_end")
	val3=$(awk 'NR=='$j'{print $3}' $i"_end")
	val4=$(awk 'NR=='$j'{print $4}' $i"_end")
	j1=$(($j+1))
	val1_1=$(awk 'NR=='$j1'{print $1}' $i"_end")
	val1_2=$(awk 'NR=='$j1'{print $2}' $i"_end")
	if [ -z $val2 ]; then
	    continue;
	else
	    if [ -z $val1_2 ]; then
		echo $val1 $val2 $val3 $val4$val1_1 >> $i"_out"
	    else
		echo $val1 $val2 $val3 $val4 >> $i"_out"
	    fi
	fi
    done
    sed -e 's/## ##/ /g' $i"_out" > $i"_3"
    sed -e 's/##/ \n 0. 0. 0. 1. /g' $i"_3" > $i"_edit"
    #
    rm -rf $i"_1" $i"_2" $i"_3" $i"_end" $i"_out"
    rm -rf files.txt
done

exit