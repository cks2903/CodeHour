# automatic generation of workflows

for i in `seq 1 99`; do 
cp "workflow1.py" "workflow$[i+1].py"
sed "s/Round1/Round$[i+1]/g" workflow$[i+1].py > workflow$[i+1]_.py
sed "s/grouping1/grouping$[i+1]/g" workflow$[i+1]_.py > workflow$[i+1]__.py
mv workflow$[i+1]__.py workflow$[i+1].py
rm workflow$[i+1]_.py;
done

