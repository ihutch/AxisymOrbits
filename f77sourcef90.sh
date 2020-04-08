#!/bin/bash
# Convert f77 source form to f90.
#Convert all line-initial c,  C or * character to !.
sed -i -e 's/^[cC*]/!/' $1
#Convert f77 continuations to f90 continuations.
# Multi-line sed analysis, by paragraphs terminated by blanks or comments.
#For non-empty lines not starting !, c, or C, just append pattern to hold space
#then delete the pattern space and restart. 
#For others, swap hold/pattern (x), and analyse with subsequent s command.
sed -i '1{h;d} ; /^[^!cC]/{H;$!d} ; x ; s/\n     [^ ]/\&\n     \&/g' $1
# Then we are in f90 form but it is better get rid of tabs and
#to align the continuation in column 74, at the line end, like this:
sed -i -e 's/\t/   /' $1
sed -i -e :a -e 's/\(^.\{1,72\}\)\&$/\1 \&/;ta' $1
# This is supposed to tidy up anything that has generated multiple & 
# by multiple conversion. A line with & in column 74 is truncated there.
sed -i 's/\(^.\{74\}\)\(\&\+\)$/\1/' $1
