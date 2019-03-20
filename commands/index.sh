for k in {7..31..4}
do
for w in $(seq 2 4 $k)
do
echo "K:$k W:$w"
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/indexes/k$k.w$w
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/indexes/k$k.w$w
cp /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/short_prg.fa .
bsub.py 4 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg short_prg.fa --num_prg 13397 --w $w --k $k --chunk_size 500 --max_forks 20
done
done
