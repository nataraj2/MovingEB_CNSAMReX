for counter in {0..22300..100}
	do
	echo $counter
		if [[ $counter -lt 100 ]]
		then
			append=0000$counter
		elif [[ $counter -lt 1000 ]]
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
		dir="dir"
		dir+=$append
		CL="CL"
		CL+=$append
		CL+=".txt"

		cat $dir/CL.*.txt > $CL
     	let counter=counter+100
done
