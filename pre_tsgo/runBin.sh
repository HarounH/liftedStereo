#!/bin/bash

# Setup variables

	# GM variables
NLABELS=40
LAMBDA=50
UNARYPOTENTIALRADIUS=3
	# Alpha Expansion Variables
NSTEPS=100

	# Counting BP Variables
NMAXITER=0
L1THRESHOLD=7 # its multiplied by nLabels within the function.
NBINS=48

	# Output variables
NOSYMOUTPUT=nosyms
L1SYMOUTPUT=l1syms
BINSYMOUTPUT=binsyms

# Compile
# echo "Compiling nosyms.out"
# g++ stereoInfer.cpp -std=c++11 -lpng -Wfatal-errors -o nosyms.out
# echo "Compiling l1syms.out"
# g++ stereoInfer2.cpp -std=c++11 -lpng -Wfatal-errors -o l1syms.out
# echo "Compiling binsyms.out"
# g++ stereoInfer3.cpp -std=c++11 -lpng -Wfatal-errors -o binsyms.out

# Run all code on one photo at a time.



# Tsubaka - average photo
LOCATION="all2/repeatedTsu"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	


# Tsubaka - average photo
LOCATION="all2/Tsukuba"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	echo "Running binsyms"
	# Run bins alpha expansion
	./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	

# Art 
LOCATION="all2/Art"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Books
LOCATION="all2/Books"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Computer
LOCATION="all2/Computer"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Dolls
LOCATION="all2/Dolls"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Drumsticks
LOCATION="all2/Drumsticks"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Dwarves
LOCATION="all2/Dwarves"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Laundry
LOCATION="all2/Laundry"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Moebius
LOCATION="all2/Moebius"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

# Reindeer
LOCATION="all2/Reindeer"
	# echo "Running normal"
	# Run normal alpha expansion
	# ./nosyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS  ../../data/${LOCATION}/disp1.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NOSYMOUTPUT}.txt
	# echo "Running l1syms"
	# Run counting bp alpha expansion
	# ./l1syms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $L1THRESHOLD ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${L1THRESHOLD}_${L1SYMOUTPUT}.txt	
	# echo "Running binsyms"
	# Run bins alpha expansion
	# ./binsyms.out single ../../data/${LOCATION}/view1.png ../../data/${LOCATION}/view5.png ../../data/${LOCATION}/output_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png $NLABELS $LAMBDA $UNARYPOTENTIALRADIUS  $NSTEPS $NMAXITER $NBINS ../../data/${LOCATION}/disp1.png ../../data/${LOCATION}/syms_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.png > ../../data/${LOCATION}/log_${NLABELS}_${LAMBDA}_${UNARYPOTENTIALRADIUS}_${NSTEPS}_${NMAXITER}_${NBINS}_${BINSYMOUTPUT}.txt	
	

