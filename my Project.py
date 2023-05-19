# motif Discovery : Randomized Motif Search Algorithm
# BA2F

import random              #creating an list of motifs of length k
def randomkmers(dna,k):                              
    motifs = []
    for one in dna:
        pos = random.randrange(0,len(dna[0])-k+1)
        motifs.append(one[pos:pos+k])
    return motifs

# creating an best matif and comparing with the other motifs to find the best one 

def randomizedmotifsearch(dna,k,t):
    motifs = randomkmers(dna,k)
    best = motifs
    while True:
        matrix = motifsToProfile(motifs)
        motifs = []
        for m in range(t):
            motifs.append(profileMostProbablekmer(dna[m],k,matrix))
        if score(motifs) < score(best):
            best = motifs
        else:
            return best
        
        
#counting the number of times each kmers appears in the dna  and returning the most frequent one

def score(motifs):
    z = zip(*motifs)
    thescore = 0
    for string in z:
        score = len(string) - max([string.count('A'), string.count('C'), string.count('G'), string.count('T')])
        thescore += score
    return thescore

#creating an profile matrix for the count of motifs

def motifsToProfile(motifs):
    d = {}
    n = float(len(motifs))
    z = list(zip(*motifs))
    for i in range(len(z)):
        d.setdefault('A', []).append((z[i].count('A')+1)/n/2)
        d.setdefault('C', []).append((z[i].count('C')+1)/n/2)
        d.setdefault('G', []).append((z[i].count('G')+1)/n/2)
        d.setdefault('T', []).append((z[i].count('T')+1)/n/2)
    return d



def profileMostProbablekmer(text, k , matrix):
    maxp = None
    probablekmer = None
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
    
        pt = 1
        for j in range(k):
            p = matrix[kmer[j]][j]
            pt *=p
        if maxp == None or pt > maxp:
            maxp = pt
            probablekmer = kmer
    return probablekmer



def runrandomtimes(dna,k,t,times):  
    bestmotifs = []
    highscore = None
    for i in range(int(times)):
        tempmotifs = randomizedmotifsearch(dna, k, t)
        tempscore = score(tempmotifs)
        if highscore == None or highscore > tempscore:
            highscore = tempscore
            bestmotifs = tempmotifs
    return bestmotifs


#inputing the file and slipting the dna and the numbers in the file
with open('C:\\Users\\C Gangadhara\\OneDrive\\Documents\\Desktop\\Amrita college\\Biology\\programs\\Test.txt') as f:
    k,t = map(int,f.readline().rstrip().split(' '))
    strings = [st.rstrip() for st in f.readlines()]
print('\n'.join(runrandomtimes(strings,k,t,1000)))