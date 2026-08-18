```python
from Bio import SeqIO

records = list(SeqIO.parse("/Users/.../dna2.fasta", "fasta"))
print(len(records))
```

    18



```python

lengths = [len(r.seq) for r in records]

print("Longest:", max(lengths))
print("Shortest:", min(lengths))
```

    Longest: 4894
    Shortest: 115



```python
def find_orfs(seq, frame):
    seq = str(seq)
    start_pos = frame - 1
    orfs = []

    for i in range(start_pos, len(seq) - 2, 3):
        if seq[i:i+3] == "ATG":
            for j in range(i, len(seq) - 2, 3):
                codon = seq[j:j+3]
                if codon in ["TAA", "TAG", "TGA"]:
                    length = j + 3 - i
                    orfs.append((length, i + 1))
                    break
    return orfs
```


```python


# ---------- Q3: ORFs ----------
all_orfs_f1 = []
all_orfs_f2 = []
all_orfs_f3 = []

for r in records:
    all_orfs_f1.extend(find_orfs(r.seq, 1))
    all_orfs_f2.extend(find_orfs(r.seq, 2))
    all_orfs_f3.extend(find_orfs(r.seq, 3))

# longest ORF in frame 2
print("4. Longest ORF (frame 2):", max(all_orfs_f2)[0])

# longest ORF in frame 3 + start position
longest_f3 = max(all_orfs_f3, key=lambda x: x[0])
print("5. Start position longest ORF (frame 3):", longest_f3[1])
```

    4. Longest ORF (frame 2): 1458
    5. Start position longest ORF (frame 3): 636



```python


# ORF finder (forward frames only)
def find_orfs(seq, frame):
    seq = str(seq)
    orfs = []

    for i in range(frame - 1, len(seq) - 2, 3):
        if seq[i:i+3] == "ATG":
            for j in range(i, len(seq) - 2, 3):
                if seq[j:j+3] in ["TAA", "TAG", "TGA"]:
                    length = j + 3 - i
                    orfs.append(length)
                    break
    return orfs


# ---------------------------
# Q1: longest ORF in ALL sequences + ALL frames
# ---------------------------
all_orfs = []

for r in records:
    for frame in [1, 2, 3]:
        all_orfs.extend(find_orfs(r.seq, frame))

print("Longest ORF overall:", max(all_orfs))


# ---------------------------
# Q2: longest ORF in specific sequence
# ---------------------------
target_id = "gi|142022655|gb|EQ086233.1|16"

target_seq = None
for r in records:
    if r.id == target_id:
        target_seq = r.seq
        break

target_orfs = []
for frame in [1, 2, 3]:
    target_orfs.extend(find_orfs(target_seq, frame))

print("Longest ORF in target sequence:", max(target_orfs))
```

    Longest ORF overall: 2394
    Longest ORF in target sequence: 1644



```python
from collections import defaultdict

# ---------- helper: count repeats ----------
def count_repeats(seq, n):
    counts = defaultdict(int)

    for i in range(len(seq) - n + 1):
        subseq = seq[i:i+n]
        counts[subseq] += 1

    return counts


# =========================
# Q1: repeat length 6
# =========================
total_6_counts = defaultdict(int)

for r in records:
    c = count_repeats(str(r.seq), 6)
    for k, v in c.items():
        total_6_counts[k] += v

max_repeat_6 = max(total_6_counts, key=total_6_counts.get)
print("Most frequent 6-mer total occurrences:", total_6_counts[max_repeat_6])


# =========================
# Q2: repeat length 12
# =========================
total_12_counts = defaultdict(int)

for r in records:
    c = count_repeats(str(r.seq), 12)
    for k, v in c.items():
        total_12_counts[k] += v

max_freq_12 = max(total_12_counts.values())

# how many different 12-mers occur max times
num_max_12 = sum(1 for v in total_12_counts.values() if v == max_freq_12)

print("Number of 12-mers occurring max times:", num_max_12)
```

    Most frequent 6-mer total occurrences: 153
    Number of 12-mers occurring max times: 4



```python

counts = defaultdict(int)

for r in records:
    seq = str(r.seq)
    for i in range(len(seq) - 6):
        kmer = seq[i:i+7]
        counts[kmer] += 1

# check only the given options
options = ["CGCGCCG", "AATGGCA", "TGCGCGC", "CATCGCC"]

for opt in options:
    print(opt, counts[opt])
```

    CGCGCCG 63
    AATGGCA 2
    TGCGCGC 36
    CATCGCC 13


New file


```python
from Bio import SeqIO
records = list(SeqIO.parse("/Users/manalali/Downloads/lambda_virus.fa", "fasta"))
print(len(records))
```

    1



```python
genome = open("/Users/manalali/Downloads/lambda_virus.fa").read().splitlines()
genome = "".join(line for line in genome if not line.startswith(">"))

```


```python
def reverseComplement(s):
    complement = {'A':'T','T':'A','C':'G','G':'C', 'N':'N'}
    return "".join(complement[b] for b in s[::-1])

P = "AGGT"
P_rc = reverseComplement(P)

count = 0
for i in range(len(genome) - len(P) + 1):
    if genome[i:i+4] == P:
        count += 1
    if P_rc != P and genome[i:i+4] == P_rc:
        count += 1

print(count)
```

    306



```python
reverseComplement("TTAA")
P = "TTAA"

count = 0
for i in range(len(genome) - len(P) + 1):
    if genome[i:i+4] == P:
        count += 1

print(count)
```

    195



```python
def reverseComplement(s):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return "".join(complement[b] for b in s[::-1])

P = "ACTAAGT"
P_rc = reverseComplement(P)
print(P_rc)
```

    ACTTAGT



```python
min_pos = None

for pattern in [P, P_rc]:
    for i in range(len(genome) - len(pattern) + 1):
        if genome[i:i+len(pattern)] == pattern:
            if min_pos is None or i < min_pos:
                min_pos = i

print(min_pos)
```

    26028



```python
def reverseComplement(s):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return "".join(complement[b] for b in s[::-1])

P = "AGTCGA"
P_rc = reverseComplement(P)

print(P_rc)

```

    TCGACT



```python
min_pos = None

for pattern in [P, P_rc]:
    for i in range(len(genome) - len(pattern) + 1):
        if genome[i:i+len(pattern)] == pattern:
            if min_pos is None or i < min_pos:
                min_pos = i

print(min_pos)
```

    450



```python
def naive_2mm(P, T):
    occurrences = []
    m = len(P)

    for i in range(len(T) - m + 1):
        mismatches = 0

        for j in range(m):
            if P[j] != T[i + j]:
                mismatches += 1
                if mismatches > 2:
                    break

        if mismatches <= 2:
            occurrences.append(i)

    return occurrences
```


```python
P = "TTCAAGCC"

matches = naive_2mm(P, genome)

print(len(matches))
```

    191



```python
def naive_2mm_leftmost(P, T):
    m = len(P)

    for i in range(len(T) - m + 1):
        mismatches = 0

        for j in range(m):
            if P[j] != T[i + j]:
                mismatches += 1
                if mismatches > 2:
                    break

        if mismatches <= 2:
            return i   # first (leftmost) match

    return None
```


```python

P = "AGGAGGTT"

print(naive_2mm_leftmost(P, genome))
```

    49



```python
import numpy as np

def readFastq(filename):
    sequences = []
    qualities = []

    with open(filename) as fh:
        while True:
            fh.readline()  # name
            seq = fh.readline().rstrip()
            fh.readline()  # +
            qual = fh.readline().rstrip()

            if len(seq) == 0:
                break

            sequences.append(seq)
            qualities.append(qual)

    return sequences, qualities


def phred33_to_q(qual):
    return np.array([ord(c) - 33 for c in qual])


# Use your full file path here
filename = "/Users/manalali/Downloads/ERR037900_1.first1000.fastq"

sequences, qualities = readFastq(filename)

read_len = len(qualities[0])
profile = np.zeros(read_len)

for qual in qualities:
    profile += phred33_to_q(qual)

profile /= len(qualities)

bad_cycle = np.argmin(profile)

print(bad_cycle)
```

    66


A new practace


```python
from Bio import SeqIO
dff = list(SeqIO.parse('/Users/manalali/Downloads/dna2.fasta', 'fasta'))
print('The Length of the file:' ,len(dff))
```

    The Length of the file: 18



```python
count = 0
with open('/Users/manalali/Downloads/dna2.fasta') as f:
    for line in f:
        if line.startswith('>'):
            count += 1
print('Number of Sequences:', count)
```

    Number of Sequences: 18



```python
from collections import defaultdict

def count_repeats(seq,n):
    counts = defaultdict(int)
    for i in range(len(seq) - n + 1):
        subseq= seq[i:i+n]
        counts[subseq] += 1
    return counts


total_5_count = defaultdict(int)
for r in dff:
    c = count_repeats(str(r.seq), 5)
    for k, v in c.items():
        total_5_count[k] += v

max_freq_5 = max(total_5_count, key=total_5_count.get)
print('Most frequent 5-mer total occurrences:', total_5_count[max_freq_5])

total_6_count = defaultdict(int)
for r in dff:
    c = count_repeats(str(r.seq), 6)
    for k, v in c.items():
        total_6_count[k] += v
max_freq_6 = max(total_6_count, key=total_6_count.get)
print('Most frequent 6-mer total occurrences:', total_6_count[max_freq_6])

```

    Most frequent 5-mer total occurrences: 418
    Most frequent 6-mer total occurrences: 153



```python
from collections import defaultdict

def count_repeats(seq,n):
    counts = defaultdict(int)
    for i in range(len(seq) - n + 1):
        subseq= seq[i:i+n]
        counts[subseq] += 1
    return counts
for r in dff:
    c = count_repeats(str(r.seq), 5)
    print(max(c, key=c.get), max(c.values()))


r = dff[0]
c = count_repeats(str(r.seq), 5)
max_kmer = max(c, key=c.get)
print('The max 5-mer',c[max_kmer])

```

    CGGCG 54
    GCGCG 14
    CGCGC 47
    GCGCG 33
    CGCGC 62
    GCGCG 37
    CGCGC 10
    CGACG 12
    CGCCG 11
    GCGCG 29
    GCCGC 16
    GCGCG 51
    GCGCG 49
    TCGTC 6
    CGTCG 16
    GCGCG 18
    CCGGC 3
    CGCCG 27
    The max 5-mer 54


Nivee exact matching


```python
def longestcommonprefix(s1,s2):
    i=0
    while i< len(s1) and i<len(s2) and s1[i]==s2[i]:
        i+=1
    return(s1[:i])
longestcommonprefix('ACCATGT', 'ACCAGTA')
```




    'ACCA'




```python

```
