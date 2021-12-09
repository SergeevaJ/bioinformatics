def read_fasta_file(fasta_file):
    reads=[]
    file=open(fasta_file, 'r')
    for line in file:
        if line[0]!='>' and line[0]!=';':
            reads.append(line[:-1])
    return reads

def make_fasta_file(contigi):
    file=open('out.fa','w')
    for i in range(len(contigi)):
        file.write('>seq'+str(i)+'\n')
        file.write(contigi[i]+'\n')
    return

def kmers_maker(reads, k): #на случай, если данные заданы в виде ридов, разбиваем их н к-меры
    return [read[i:k + i] for read in reads for i in range(0, len(read) - k + 1)]

# строим граф де брёйна
def de_bruijn_graph_maker(k_mers):
    vertexes=set()
    not_starts=set()
    ss=dict()
    st =set() #множество вершин, которые не могут быть промежуточными в пути
    for kmer in k_mers:
        p=kmer[:-1]
        s=kmer[1:]
        vertexes.add(p)
        vertexes.add(s)
        if s in not_starts:
            st.add(s)
        not_starts.add(s)
        if p not in ss.keys():
            ss[p]=[s]
        else:
            st.add(p)
            ss[p].append(s)
    for v in vertexes:
        if v not in ss.keys():
            ss[v]=[]
    return (ss, vertexes - not_starts, st)

def find_contig(ss, st, starts, start):
    next=ss[start].pop()
    way=[start, next]
    while len(ss[next])==1 and (next not in st):
        way.append(ss[next][0])
        new_next=ss[next][0]
        next=new_next
    if len(ss[next])!=0:
        starts.add(next)
    if len(ss[start])!=0:
        starts.add(start)
    return way

reads=read_fasta_file("proba.fa")
k_mers=kmers_maker(reads,4) # если данные заданы в виде  ридов, если нет то убираем эту строчку
ss, starts, st= de_bruijn_graph_maker(k_mers)
trail=[]

if starts==set():
    starts.add(st.pop())
while len(starts)!=0:
    v=starts.pop()
    contig=find_contig(ss, st, starts, v)
    trail.append(contig)

def assembler(trail):
    contigi=[]
    for way in trail:
        t=way[0]
        for i in range(1,len(way)):
            t=t+way[i][-1]
        contigi.append(t)
    return contigi


contigi=assembler(trail)
make_fasta_file(contigi)
