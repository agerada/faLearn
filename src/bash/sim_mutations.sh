infile="/volatile/agerada/molecularMIC/sims/esco_25922/genome/sequence.fasta"
outpath="/volatile/agerada/molecularMIC/sims/esco_25922/mutations/"
n=2
option="args -sn 0.01 -in 0.005 -inmin 1 -inmax 5"

for i in $(seq 1 $n); do
    mutation-simulator -o ${outpath}mut$i ${infile} ${option}
    mv "${outpath}mut${i}_ms.fasta" "${outpath}mut${i}_ms.fna"
done
