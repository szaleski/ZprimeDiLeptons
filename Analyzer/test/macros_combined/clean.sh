#!/bin/bash

rm -f att.txt

for file in `ls -trd1 jobs/*.sh.* | head -500`; do
 echo $file
 echo "rm -f $file" >> att.txt
done

bash att.txt
