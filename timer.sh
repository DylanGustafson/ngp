#!/bin/bash
x="./$1"
if ! [ -f $x ]; then
	echo Rebuilding executable...
	make $x
	echo Done! Running now...
fi

echo "Running $x..."
a=$(date +%s.%N)
$x $2 $3 $4
b=$(date +%s.%N)

echo Total time\: $(echo "$b-$a" | bc) seconds
