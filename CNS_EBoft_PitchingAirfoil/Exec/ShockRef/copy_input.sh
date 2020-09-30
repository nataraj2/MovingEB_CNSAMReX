counter=000
counterplus1=1
while [ $counter -le 22300 ]
do
echo $counter

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
	chk="chk"
	inputs+=$append
	chk+=$append

	cp inputs $inputs
        sed -i "s/50000/$counterplus1/g" $inputs
        sed -i "s/chk147000/$chk/g" $inputs

	let counter=counter+100
	let counterplus1=counterplus1+100
done

echo All done
