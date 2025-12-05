
#!/bin/sh
for file in rawChipdata/*/; do
   	echo “$file”
   	cd "$file"
	gunzip *.gz
	echo "$file done"
	cd ..
	done
	cd ..
