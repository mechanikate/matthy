mkdir -p $2 
grep -o '^[^#]*' $1.py > $2/$1.py  
