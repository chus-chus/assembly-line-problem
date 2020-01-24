#!/bin/sh
rm -f sols.txt

END=10
if [ "$1" = "hard" ]; then
    END=20
fi

for index in $(seq 1 $END)
do
    ./greedy.exe ./public_benchs/$1-$index.txt out.txt
    ./check ./public_benchs/$1-$index.txt out.txt
    echo $index >> prova_hard1.txt
    head -2 out.txt >> prova_hard1.txt
    echo "" >> prova_hard1.txt
    echo "" >> prova_hard1.txt
done
