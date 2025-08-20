mkdir -p applied
for f in *.py; do
	[ -f "$f" ] || break
	echo $f
	grep -o '^[^#]*' $f | expand -t 2 -i > applied/$f 	
done
