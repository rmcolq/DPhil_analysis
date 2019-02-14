from Bio import SeqIO
import argparse

def sim_genomes(random_paths):
    record_dict = SeqIO.to_dict(SeqIO.parse(random_paths, "fasta"))
    g1 = []
    g2 = []
    g3 = []
    g4 = []
    g4_p = []
    t = [] 
    count = 0
    for record in record_dict:
        if record.startswith('GC0') and record.endswith('_0'):
            s1 = str(record_dict[record].seq)
            try:
                pair = record[:-1] + '1'
                name = record[:-2]
                s2 = str(record_dict[pair].seq)
                if 0 <= count < 200:
                    g1.append(s1)
                    g2.append(s2)
                    if 0 <= count < 50:
                        g3.append(s1)
                    elif 50 <= count < 100:
                        g4_p.append(s2)
                elif 200 <= count < 250:
                    g3.append(s1)
                elif 250 <= count < 300:
                    g4.append(s1)
                    if count == 299:
                        g4.extend(g4_p)
                elif 300 <= count < 400:
                    g3.append(s1)
                    g4.append(s2)
                elif count >= 400:
                    break
                t.append(name)
                count += 1
            except:
                continue
                
    genome1 = "".join(g1)
    genome2 = "".join(g2)
    genome3 = "".join(g3)
    genome4 = "".join(g4)
    
    with open("genomes.fa", 'w') as f:
        f.write(">genome1\n%s\n" %genome1)
        f.write(">genome2\n%s\n" %genome2)
        f.write(">genome3\n%s\n" %genome3)
        f.write(">genome4\n%s\n" %genome4)
    with open("genomes_truth.txt", 'w') as f:
        f.write(" ".join(t))
    
parser = argparse.ArgumentParser(description='Constructs 4 genomes from random paths')
parser.add_argument('--in_fa', type=str,
                    help='Input FASTA of random paths through PanRG')
args = parser.parse_args()

sim_genomes(args.in_fa)
