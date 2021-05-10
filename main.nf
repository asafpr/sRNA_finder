#!/usr/bin.env nextflow

params.outdir="."
params.split=1
params.colors="#aab4bb,#003e6a,#f5b400,#f78a88"
params.shapes="1,1,1,1"
params.label="sRNAs"
params.apikey="ebb6c6d7dd3133cc27642320fb9c039fcd08"
params.source="PATRIC"
Channel.value(params.apikey).set{apikey}
Channel.fromPath(params.genomes).set{genomes_lst}
Channel.fromPath(params.genomes).splitCsv(header:false, sep:"\t").map{row -> [row[0], row[1]]}.set{genomes}
Channel.value(file(params.sequences)).set{seqs}
Channel.fromPath(params.tree).set{tree_ch}
// genomes.splitText(by:params.split).set{genomes_spl}

if (params.source == "PATRIC"){
process download_genome_patric{
  publishDir path:params.outdir, mode:'copy'
  label 'single_cpu'
  input:
    tuple gen, name from genomes

  output:
    tuple gen, name, file("${gen}.fna") into ffiles

  script:
"""
wget ftp://ftp.patricbrc.org/genomes/${gen}/*.fna
"""
}

}else{
process download_genomes{
  publishDir path:params.outdir, mode:'copy'

  label 'eutils'
  label 'long_run'
  input:
    tuple gen, name from genomes
    env NCBI_API_KEY from apikey
  output:
    tuple gen, name, file ("library/gtdb/*.fna") into ffiles

  script:
  """
  gtdbToTaxonomy.pl --genome $gen 
  """
}
}

process build_db{
  publishDir path:params.outdir, mode:'copy'
  
  label 'blast'
  label 'high_mem'

  input:
    tuple fn, name, file(f) from ffiles

  output:
    tuple fn, name, file(f), file("$f.*") into bdb

  script:
  """
  makeblastdb -in $f -dbtype nucl
  """
}

process run_blast{
  publishDir path:params.outdir, mode:'copy'

  label 'blast'
  label 'high_mem'

  input:
  tuple fn, name, file(fname), file(fall ) from bdb
  file seq from seqs

  output:
    tuple fn, name, file("${fname}.blout"), file(fname) into blout

  script:
  """
  blastn -db $fname -query $seq -out ${fname}.blout -outfmt '6 std qlen slen'  -template_type coding_and_optimal -template_length 16 -evalue 1000 -num_threads ${task.cpus} -task dc-megablast 
  """
}

process filter_blast{
  publishDir path:params.outdir, mode:'copy'
  label 'biopython'
  tag "$fn"
  input:
    tuple fn, name, file(blo), file(orig) from blout
    file sq from seqs
  output:
    file "${fn}.txt" into digest
    file "${fn}*.fa" optional true into outseqs
  script:
"""
#!/usr/bin/env python

from Bio import SeqIO
# Read the sequences file and get the names in order
ord = []
with open("${sq}") as fin:
  for line in fin:
    if (line.startswith(">")):
      ord.append(line.strip()[1:])

# Read the blast input
hasb = {o:-1 for o in ord}
outseq = {}
records = SeqIO.to_dict(SeqIO.parse("$orig", "fasta"))
maxid = {o : 80 for o in ord}
maxcov = {o : 0.8 for o in ord}
with open("${blo}") as bin:
  for line in bin:
    spl = line.strip().split("\\t")
    if ((float(spl[2]) > maxid[spl[0]]) & (float(spl[3]) > float(spl[12])*maxcov[spl[0]])):
      maxid[spl[0]] = float(spl[2])
      maxcov[spl[0]] = float(spl[3])/float(spl[12])
      hasb[spl[0]] = 1
      mfrom = min(int(spl[8]), int(spl[9]))-1
      mto = max(int(spl[8]), int(spl[9]))
      rev = int(spl[9]) < int(spl[8])
      outseq[spl[0]] = str(records[spl[1]].seq[mfrom:mto])
      if rev:
        outseq[spl[0]] = str(records[spl[1]].seq[mfrom:mto].reverse_complement())
outstr = "${fn}," + ",".join(str(hasb[o]) for o in ord)
with open("${fn}.txt", 'w') as fout:
  fout.write(outstr + "\\n")
for k in  outseq:
  with open("${fn}_" + k + ".fa", "w") as sout:
    sout.write(">${fn} ${name}\\n" + outseq[k] + "\\n") 
"""
}

process msa{
  publishDir path:params.outdir, mode:'copy'

  label 'mafft'
  input:
    file fa from outseqs.collect()
  output:
    file "*.mafft" into msa

  script:
"""
for m in \$(ls *.fa | awk -F"_" '{print \$NF}'| sort |uniq) ; do cat *"\$m" > all"\$m".fa &&  mafft all"\$m".fa > MSA_"\$m".mafft; done
"""
}

process combine_out{
  publishDir path:params.outdir, mode:'copy'
  input:
    file allout from digest.collect()
    file sq from seqs
  
  output:
    file "all_srnas_status.txt" into allout

  script:
"""
cat <<EOF > all_srnas_status.txt
DATASET_BINARY
SEPARATOR COMMA

FIELD_COLORS,${params.colors}
FIELD_SHAPES,${params.shapes}
DATASET_LABEL,${params.label}
EOF

printf "FIELD_LABELS" >> all_srnas_status.txt
grep ">" $sq | tr  ">"  "," | tr -d "\\n" >> all_srnas_status.txt
echo "" >> all_srnas_status.txt
echo "DATA" >> all_srnas_status.txt

cat $allout >> all_srnas_status.txt
"""
}

process prune_tree{
  publishDir path:params.outdir, mode:'copy'
  label 'ete3'

  input:
    file tree from tree_ch
    file leaves from genomes_lst

  output:
    file "${tree.baseName}.pruned.newick" into prntree

  script:
"""
#!/usr/bin/env python
from ete3 import Tree
from ete3.coretype.tree import TreeError
# Read the tree
t = Tree("$tree", format=1, quoted_node_names=True)

# Read the genomes (first column)
lv = []
with open("$leaves") as fin:
    for line in fin:
        lv.append(line.strip().split("\\t")[0])
tlv = t.get_leaf_names()
try:
  t.prune(set(lv) & set(tlv))
except TreeError:
  raise
t.write(format=1, outfile="${tree.baseName}.pruned.newick")
"""
} 


