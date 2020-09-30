counter=000
while [ $counter -le 22300 ]
	do
	if [[ $counter -lt 1000 ]]
	then
        append=000$counter
	elif [[ $counter -lt 10000 ]]
	then
        append=00$counter
	elif [[ $counter -lt 100000 ]]
	then
        append=0$counter
	elif [[ $counter -lt 1000000 ]]
	then
        append=$counter
	fi
        inputs="inputs"
        dir="dir"
        inputs+=$append
        dir+=$append

	srun -n 216 ./clvsdeg.ex $inputs
	mkdir $dir
	mv CL.* $dir
	let counter=counter+100
done
