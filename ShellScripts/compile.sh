# This script goes into the proper folder, compiles, executes, and moves on
# This is my first script that I've written, and I have a for loop written in 
# two ways, one uses a while and one is just a for loop
# The for loop 


for (( i=1; i<=6; i++ ))
do 
	cd /nobackup/leftynm/Nov2011/Mix

	cd $i

# 	if [ $i -eq 1 ]; then
# 	  cd Exponential/0
# 	elif [ $i -eq 2 ]; then
# 	  cd Quadratic/0
# 	elif [ $i -eq 3 ]; then
# 	  cd SinLog/0
# 	fi
	
	cd 0
	
	for (( j=0; j<=4; j++ ))
	do
		cd ../$j		
		echo 'compling folder' $i '/' $j
		make
		qsub submit.sh
	done
done



# i=1
# while [ $i -le 5 ]
# do
# 	cd ../$i
# 	gfortran S1* ../S2* ../S3*
# 	./a.out
# 	(( i++ ))
# done

# for ((  i = 0 ;  i <= 5;  i++  ))
# do
# 	cd ../$i
# 	gfortran S1* ../S2* ../S3*
# 	./a.out
# done