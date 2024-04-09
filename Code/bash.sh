# replace all occurences of '__' with '_' in all files in the current directory

for file in *; do
    sed -i 's/__/_/g' $file
done

# As a one-liner

for file in *; do sed -i 's/__/_/g' $file; done

# Above does not work, try using mv instead

for file in *; do mv $file $(echo $file | sed 's/__/_/g'); done