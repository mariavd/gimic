#!/bin/bash

gimicdir="$1"

if [ -z $gimicdir ]
then
    echo "Gimic directory not specified as a command line argument"
    exit
fi

testname=benzene

printf "\nPerforming test on $testname bond integral\n"

cp ./$testname/XDENS ./$testname/int/XDENS
cp ./$testname/MOL ./$testname/int/MOL

(cd ./$testname/int && $gimicdir/gimic gimic.inp > gimic.test.out )

diatropic=$(grep -A 2 "Induced current" ./$testname/int/gimic.test.out | awk '{ if (NR == 2) printf("% f\n", $5); }')
paratropic=$(grep -A 2 "Induced current" ./$testname/int/gimic.test.out | awk '{ if (NR == 3) printf("% f\n", $5); }')
total=$(grep "Induced current (nA/T)" ./$testname/int/gimic.test.out | awk '{printf("% f\n", $5); }')

# Debugging
echo $diatropic $paratropic $total

dia_ref=$(grep -A 2 "Induced current" ./$testname/int/gimic.out | awk '{ if (NR == 2) printf("% f\n", $5); }')
para_ref=$(grep -A 2 "Induced current" ./$testname/int/gimic.out | awk '{ if (NR == 3) printf("% f\n", $5); }')
total_ref=$(grep "Induced current (nA/T)" ./$testname/int/gimic.out | awk '{printf("% f\n", $5); }')

# Debugging
echo $dia_ref $para_ref $total_ref

# avoid returning success value by default
test1=2
test2=2
test3=2

test1=$( awk -v diatropic="$diatropic" -v dia_ref="$dia_ref" 'BEGIN{ if ((dia_ref - diatropic) < 1e-5) { print 0; } else {print 1; } }'  )
test2=$( awk -v paratropic="$paratropic" -v para_ref="$para_ref" 'BEGIN{ if ((para_ref - paratropic) < 1e-5) { print 0; } else {print 1; } }' )
test3=$( awk -v total="$total" -v total_ref="$total_ref" 'BEGIN{ if ((total_ref - total) < 1e-5) { print 0; } else {print 1; } }' )

#Debugging
echo $test1 $test2 $test3

if [ $test1 -eq 0 ] && [ $test2 -eq 0 ] && [ $test3 -eq 0 ]
then
    echo "0" # success
else
    echo "1" # failure
fi

rm -rf ./$testname/int/XDENS ./$testname/int/MOL
