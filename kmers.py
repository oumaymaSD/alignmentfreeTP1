
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def encode(letter):
    #letters={'A':0, 'C':1, 'T': 2, 'G':3}
    letters = ['A', 'C', 'T', 'G']
    
    return letters.index(letter)

  

def mask(k):
    return( (1<< 2*(k-1))-1)


def str2kmer(nucl,scores):
    if nucl in scores:
        return scores[nucl]
    else:
        return 0

def stream_kmers(text, k):
    list_kmer=[]
    kmer=0
    revkmer=0
    mask=(1<<((k-1)*2))-1
    revmask=(1<<(k*2))-1-3
    scores={"A":0,"C":1,"T":2,"G":3}
    rscores={"A":2,"C":3,"T":0,"G":1}
    for i in range(k-1):
        kmer=kmer<<2
        revkmer=revkmer>>2
        kmer+=str2kmer(text[i],scores)
        revkmer+=str2kmer(text[i],rscores)<<(k-1)*2
    for nucl in text[k-1:]:
        kmer=kmer&mask
        #revkmer=revkmer&revmask
        kmer=kmer<<2
        revkmer=revkmer>>2
        kmer+=str2kmer(nucl,scores)
        revkmer+=str2kmer(nucl,rscores)<<(k-1)*2
        list_kmer.append(min(kmer, revkmer))
    return list_kmer



k=3
a=stream_kmers('TGTA', k)
print (a)
for j in a:
    print(kmer2str(j, k))
