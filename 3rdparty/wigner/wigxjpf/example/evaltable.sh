#!/bin/sh

if [ -z "$STATS" ]
then
    STATS=--stats
fi

EVAL=bin/wigxjpf

# Small memory use and fast evaluation
if [ "x$MEMUSE" = "xsmall" ]
then
    $EVAL $STATS --3j=1.5,1.5,1,1.5,-0.5,-1
    $EVAL $STATS --3j=1,2,4,1,-2,1
    $EVAL $STATS --3j=1,2,1,1,1,-2
    $EVAL $STATS --3j=1,2,3,0.5,-0.5,0

    $EVAL $STATS --6j=2,2,1,2,1,1
    $EVAL $STATS --6j=1,1,0,2,2,1
    $EVAL $STATS --6j=1,1,1,2,2,1
    $EVAL $STATS --6j=1,1,2,2,2,1
    $EVAL $STATS --6j=1,1,0,2,2,2
    $EVAL $STATS --6j=1,1,0,2,2,1
    $EVAL $STATS --6j=1,1,1,2,2,2
    $EVAL $STATS --6j=1,1,1,2,2,1
    $EVAL $STATS --6j=1,1,2,2,2,2
    $EVAL $STATS --6j=1,1,2,2,2,1

    $EVAL $STATS --9j=1,2,3,1,2,3,0,2,2
    $EVAL $STATS --9j=0.5,1.5,2,1.5,0.5,2,1,1,0
    $EVAL $STATS --9j=1,2,3,3,2,4,2,3,5
    $EVAL $STATS --9j=1,2,3,3,2,4,2,4,5
    $EVAL $STATS --9j=1,2,3,3,2,4,3,2,5
    $EVAL $STATS --9j=1,2,3,3,2,4,3,3,5
    $EVAL $STATS --9j=1,2,3,3,2,4,4,1,5
    $EVAL $STATS --9j=1,2,3,3,2,4,4,2,5
    $EVAL $STATS --9j=1,2,3,3,2,4,4,3,5
    $EVAL $STATS --9j=1,2,3,3,2,4,4,4,5
    $EVAL $STATS --9j=1,2,3,2,3,4,3,4,5

    $EVAL $STATS --3j=1,1,1,1,0,-1
    $EVAL $STATS --3j=3,3,3,-1,-1,2
    $EVAL $STATS --3j=3.5,2.5,2,3.5,-1.5,-2
    $EVAL $STATS --3j=7.5,7.5,0,1.5,-1.5,0
    $EVAL $STATS --3j=8,5.5,4.5,2,-3.5,1.5
    $EVAL $STATS --3j=8,7,6,3,0,-3

    $EVAL $STATS --3j=15,30,40,2,2,-4
    $EVAL $STATS --3j=30,30,30,0,15,-15
    $EVAL $STATS --3j=143,100,60,-10,60,-50
    $EVAL $STATS --3j=160,100,60,-10,60,-50
    $EVAL $STATS --3j=200,200,200,-10,60,-50

    $EVAL $STATS --3j=70,75,80,20,-40,20

    $EVAL $STATS --6j=2,2,1,2,1,2
    $EVAL $STATS --6j=3.5,3,1.5,1,1.5,3
    $EVAL $STATS --6j=4,4,1,4,3,1
    $EVAL $STATS --6j=6,6,4,4.5,3.5,5.5
    $EVAL $STATS --6j=6.5,6,1.5,3,3.5,6
    $EVAL $STATS --6j=8,8,8,8,7,7

    $EVAL $STATS --6j=8,8,8,8,8,8
    $EVAL $STATS --6j=20,20,20,20,20,20
    $EVAL $STATS --6j=128,120,72,112,48,80
    $EVAL $STATS --6j=230,80,150,190,230,120
    j=200 ; $EVAL $STATS --6j=$j,$j,$j,$j,$j,$j

    $EVAL $STATS --6j=80,70,60,60,70,80

    $EVAL $STATS --9j=1,3,2,1.5,2.5,2,0.5,0.5,0
    $EVAL $STATS --9j=3,4,2,3.5,3.5,2,0.5,0.5,1
    $EVAL $STATS --9j=1,3,4,0.5,3.5,3,0.5,0.5,1
    $EVAL $STATS --9j=3,3.5,3.5,2.5,3,3.5,0.5,0.5,1

    $EVAL $STATS --9j=8.5,9.5,7.0,12.5,8.0,8.5,8.0,10.5,9.5
    $EVAL $STATS --9j=15,15,15,15,3,15,15,18,10
    $EVAL $STATS --9j=45,30,20,45,15,35,90,45,45
    $EVAL $STATS --9j=60,70,130,50,70,120,60,50,40
    $EVAL $STATS --9j=100,80,50,50,100,70,60,50,100

    $EVAL $STATS --9j=5.0,0.5,4.5,5.0,0.5,5.5,9.0,1.0,10.0
    $EVAL $STATS --9j=15.0,15.0,30.0,15.0,3.0,15.0,15.0,18.0,30.0
    $EVAL $STATS --9j=20.0,10.0,30.0,30.0,30.0,60.0,10.0,20.0,30.0
    $EVAL $STATS --9j=30.0,20.0,10.0,30.0,10.0,20.0,60.0,30.0,30.0
    $EVAL $STATS --9j=45.0,30.0,20.0,45.0,15.0,35.0,90.0,45.0,45.0
    $EVAL $STATS --9j=50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0
    $EVAL $STATS --9j=74.0,67.0,60.0,82.0,60.0,70.0,10.0,20.0,30.0
    $EVAL $STATS --9j=84.0,90.0,61.0,60.0,64.0,40.0,52.0,40.0,30.0
    $EVAL $STATS --9j=94.0,67.0,86.0,61.0,73.0,52.0,70.0,60.0,80.0
    $EVAL $STATS --9j=60.0,70.0,115.0,50.0,70.0,110.0,60.0,50.0,40.0
    $EVAL $STATS --9j=100.0,80.0,50.0,50.0,100.0,70.0,60.0,50.0,100.0

    $EVAL $STATS --9j=20,20,40,20,20,40,20,20,40
fi

# Reasonable memory use (< 100 MB)
if [ "x$MEMUSE" = "xlarge" ]
then
    j=600 ; $EVAL $STATS --6j=$j,$j,$j,$j,$j,$j

    j=200 ; $EVAL $STATS --9j=$j,$j,$j,$j,$j,$j,$j,$j,$j
    j=1000 ; $EVAL $STATS --9j=$j,$j,$j,$j,$j,$j,$j,$j,$j

    $EVAL $STATS --9j=250,250,500,250,250,500,500,500,1000
fi

# Huge memory use (< 32 GB)
if [ "x$MEMUSE" = "xhuge" ]
then
    $EVAL $STATS --3j=50000,50000,50000,1000,-6000,5000

    j=10000 ; $EVAL $STATS --6j=$j,$j,$j,$j,$j,$j
    j=50000 ; $EVAL $STATS --6j=$j,$j,$j,$j,$j,$j

    j=2000 ; $EVAL $STATS --9j=$j,$j,$j,$j,$j,$j,$j,$j,$j
fi
