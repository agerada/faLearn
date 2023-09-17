infile="/volatile/agerada/molecularMIC/sims/esco_25922/genome/sequence.fasta"
outpath="/volatile/agerada/molecularMIC/sims/esco_25922/mutations/"
n=10
option="args -sn 0.01 -in 0.005 -inmin 1 -inmax 5"

for i in $(seq -f "%04g" 1 $n); do
    mutation-simulator -o ${outpath}${i}mut ${infile} ${option}
    mv "${outpath}${i}mut_ms.fasta" "${outpath}${i}mut_ms.fna"
done
