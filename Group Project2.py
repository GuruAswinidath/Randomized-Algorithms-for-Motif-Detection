###############
def count_motif(motif):
    count = {}  # initializing the count dictionary
    k = len(motif[0])  # the length of the first row in motifs is actually the column number of motifs matrix
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(motif)  # t is equal to the row number of motifs matrix
    for i in range(t):
        for j in range(k):
            symbol = motif[i][j]
            count[symbol][j] += 1
    return count


motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
counts_dict = count_motif(motifs)
counts_dict


def profile_motif(motif):
    t = len(motif)
    k = len(motif[0])
    profile = {}
    count = count_motif(motif)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)
    return profile


motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
profile_dict = profile_motif(motifs)
profile_dict
#############################



def consensus_motif(motif):
    consensus_ = ""
    count = count_motif(motif)
    k = len(motif[0])
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus_ += frequent_symbol
    return consensus_


motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
consensus = consensus_motif(motifs)
consensus


def score_motif(motifs):
    consensus = consensus_motif(motifs)
    count = count_motif(motifs)
    t = len(motifs)
    k = len(motifs[0])
    max_com = 0
    for j in range(k):
        max_com += count[consensus[j]][j]
    score = t * k - max_com
    # another way that using a former function hamming_distance in replication.py
    # score = 0
    #     for i in range(len(motifs)):
    #         score += hamming_distance(motifs[i], consensus_motif(motifs))
    #     return score
    return score


motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
score = score_motif(motifs)
score


text = "ACGGGGATTACC"
motif_profile = {'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
                 'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
                 'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
                 'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}

def probability_profile(text, motif_profile):
    pe = 1
    for i in range(len(text)):
        pe *= motif_profile[text[i]][i]
    return pe

pr = probability_profile(text, motif_profile)


def profile_most_probable_kmer(text, k, motif_profile):
    probability_kmer = {}
    most_probable_kmer = []
    for i in range(len(text) - k + 1):
        probability_kmer[text[i:i + k]] = probability_profile(text[i:i + k], motif_profile)
    prob_max = max(probability_kmer.values())
    for key, value in probability_kmer.items():
        if probability_kmer[key] == prob_max:
            most_probable_kmer.append(key)
    return most_probable_kmer[0]


text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
k = 5
motif_profile = {'A': [0.2, 0.2, 0.3, 0.2, 0.3], 'C': [0.4, 0.3, 0.1, 0.5, 0.1], 'G': [0.3, 0.3, 0.5, 0.2, 0.4],
                 'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
kmer = profile_most_probable_kmer(text, k, motif_profile)
kmer

####################
def greedy_motif_search(dna, k, t):
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
    for m in range(len(dna[0]) - k + 1):  # loop through the first dna string, motifs[0] = the first k-mer in dna[0]
        motifs = [dna[0][m:m + k]]  # motifs is the m-th k-mer got through looping the first dna string
        for j in range(1, t):
            profile = profile_motif(motifs[0:j])
            motifs.append(profile_most_probable_kmer(dna[j], k, profile))
        if score_motif(motifs) < score_motif(best_motifs):
            best_motifs = motifs
    return best_motifs


dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
k = 3
t = 5
t_kmers = greedy_motif_search(dna, k, t)
t_kmers


def count_with_pseudocount(motifs):
    t = len(motifs)  # t is equal to the row number of motifs matrix
    k = len(motifs[0])  # the length of the first row in motifs is actually the column number of motifs matrix
    count = {}  # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count
#######################

motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
pseudocounts_dict = count_with_pseudocount(motifs)
pseudocounts_dict


def profile_with_pseudocount(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = {}
    count = count_with_pseudocount(motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / (t+4))
    return profile

motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
profile_pseudo_dict = profile_with_pseudocount(motifs)
profile_pseudo_dict


def consensus_with_pseudocount(motifs):
    consensus = ""
    count = count_with_pseudocount(motifs)
    k = len(motifs[0])
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus



def score_with_pseudocount(motifs):
    consensus = consensus_with_pseudocount(motifs)
    count = count_with_pseudocount(motifs)
    t = len(motifs)
    k = len(motifs[0])
    max_com = 0
    for j in range(k):
        max_com += count[consensus[j]][j]
    score = t * k - max_com
    # another way that using a former function hamming_distance in replication.py
    # score = 0
    #     for i in range(len(motifs)):
    #         score += hamming_distance(motifs[i], consensus_motif(motifs))
    #     return score
    return score


def probability_with_pseudocount(text, motif_profile):
    probability = 1
    for i in range(len(text)):
        probability *= motif_profile[text[i]][i]
    return probability



def profile_most_probable_kmer_with_pseudocount(text, k, motif_profile):
    probability_kmer = {}
    most_probable_kmer = []
    for i in range(len(text)-k+1):
        probability_kmer[text[i:i+k]] = probability_with_pseudocount(text[i:i+k], motif_profile)
    prob_max = max(probability_kmer.values())
    for key, value in probability_kmer.items():
        if probability_kmer[key] == prob_max:
            most_probable_kmer.append(key)
    return most_probable_kmer[0]


def greedy_motif_search_with_pseudocount(dna, k, t):
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
        # best_motifs initialized with all the first k-mer (total number is t) from each dna string

    for m in range(len(dna[0]) - k + 1):  # loop through the first dna string, motifs[0] = the first k-mer in dna[0]
        motifs = [dna[0][m:m + k]]  # motifs is the m-th k-mer got through looping the first dna string
        for j in range(1, t):
            profile = profile_with_pseudocount(motifs[0:j])
            motifs.append(profile_most_probable_kmer_with_pseudocount(dna[j], k, profile))
        if score_with_pseudocount(motifs) < score_with_pseudocount(best_motifs):
            best_motifs = motifs
    return best_motifs


dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
k = 3
t_kmers_pseudo= greedy_motif_search_with_pseudocount(dna, k, t)
t_kmers_pseudo



def profile_probable_motifs(motif_profile, dna):
    k = len(motif_profile["A"])  # get the column number from the profile matrix
    motifs = []
    for i in range(len(dna)):
        motifs.append(profile_most_probable_kmer(dna[i], k, motif_profile))
    return motifs


motif_profile = {'A': [0.8, 0.0, 0.0, 0.2], 'C': [0.0, 0.6, 0.2, 0.0], 'G': [0.2, 0.2, 0.8, 0.0],
                 'T': [0.2, 0.2, 0.0, 0.8]}
dna = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
kmers = profile_probable_motifs(motif_profile, dna)
kmers


import random
def random_motifs(dna, k, t):
    motifs = []
    for i in range(t):
        m = random.randint(0, len(dna[0]) - k)
        motifs.append(dna[i][m:m + k])
    return motifs


dna = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
k = 3
t = 5
random_t_kmers = random_motifs(dna, k, t)
random_t_kmers
