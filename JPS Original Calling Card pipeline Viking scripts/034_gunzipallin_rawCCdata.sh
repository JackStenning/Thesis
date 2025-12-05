
#!/bin/sh
for file in jack/rawCCdata/*/; do
   	echo “$file”
   	cd "$file"
	gunzip *.gz
	echo "$file done"
	cd ../..
	done
