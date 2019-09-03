for i in *.cso; do
./maxcso --decompress "$i";
files=$(echo $i | sed 's/\.cso$//g');
7z a "$files.iso.gz" "$files.iso";
rm "$i"
rm "$files.iso"
done
