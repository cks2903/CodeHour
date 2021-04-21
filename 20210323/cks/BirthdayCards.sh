# automatic generation of invitations



for i; do 
cp "BirthdayDraft.txt" "BirthdayDraft$i.txt"
sed "s/Name/$i/g" "BirthdayDraft$i.txt"  > "BirthdayDraft_$i.txt"
rm "BirthdayDraft$i.txt"

done

