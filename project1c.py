def import_reads(read_file, read_length):
    reads = []
    data = open(read_file).readlines()
    for line in data:
        if line[0] == '>':
            pass
        else:
            i = read_length - len(line)
            if i < 0:
                reads.append(line[:i])
            else:
                reads.append(line)
            #reads.append(line[:-1])
    return reads

def read_genome(genome_file):
    genome = ''
    data = open(genome_file).readlines()
    for line in data[1:]:
        genome+=line[:-1]
    return genome

def generate_index(reference_genome, kmer_length):
    index = {}
    for i in range(len(reference_genome) - kmer_length + 1):
        window = reference_genome[i:i+kmer_length]
        if window in index:
            index[window].append(i)
        else:
            index[window] = [i]
    return index

def get_kmers(read): #Split read into 3 kmers
    kmers = []
    kmers.append(read[:16])
    kmers.append(read[16:32])
    kmers.append(read[32:])
    return kmers

def get_positions(kmers, index, kmer_length): #use the index to get positions of a kmer within reference genome
    possible_positions = []
    for i in range(len(kmers)):
        kmer = kmers[i]
        if kmer in index:
            positions = index[kmer]
            for pos in positions:
                if pos-kmer_length*i >=0:
                    possible_positions.append(pos-kmer_length*i)
                else:
                    pass
        else:
            pass
    return possible_positions

def hamming_distance(str1, str2):
    distance = 0
    if len(str1) != len(str2):
        print('unequal lengths, cannot calculate hamming distance')
        print(len(str1))
        print(len(str2))
        return 
    else:
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                distance+=1
            else:
                pass
    return distance

def output(output_file, dict):
    with open(output_file, 'w') as file:
        for read_num in dict:
            genome_num = dict[read_num][0]
            file.write(f'>read_{read_num} Genome_Number_{genome_num}\n')
    
    return


kmer_length = 16
read_length = 48
read_file = 'project1c_reads.fasta'
reads = import_reads(read_file, read_length)
read_matches = {}
for i in range(800000):
    read_matches[i] = []

scores = []
for i in range(100):
    counter = 0
    '''
    for read in read_matches:
        if len(read_matches[read]) != 0:
            counter+=1
    print(counter)
    '''
    genome_file = 'project1c_genome_' + str(i) +'.fasta'
    genome = read_genome(genome_file)
    index = generate_index(genome, kmer_length)

    for j in range(len(reads)): #read in reads:
        read = reads[j]
        if len(read) != read_length:
            continue
        else:
            pass
        kmers = get_kmers(read)
        possible_positions = get_positions(kmers, index, kmer_length)
        if len(possible_positions) == 0:
            pass
        else:
            boolean = False
            for pos in possible_positions:
                window = genome[pos:pos+read_length]
                if hamming_distance(window, read) <= 2:
                    read_matches[j].append(i)
                    boolean = True
                    counter+=1
                else:
                    window_ins_1 = window[:-1]
                    window_ins_2 = window[1:]
                    window_del_1 = genome[pos:pos+read_length+1]
                    if pos > 0:
                        window_del_2 = genome[pos-1:pos+read_length]
                    else:
                        window_del_2 = genome[pos:pos+read_length+1]
                    for k in range(len(read)):
                        read_ins = read[:k]+read[k+1:]
                        read_del = read[:k]+'-'+read[k:]    
                        if hamming_distance(window_ins_1, read_ins) <= 2 or hamming_distance(window_ins_2, read_ins) <= 2:
                            read_matches[j].append(i)
                            boolean = True
                            counter+=1
                            break
                        elif hamming_distance(window_del_1, read_del) <= 2 or hamming_distance(window_del_2, read_del) <= 2:
                            read_matches[j].append(i)
                            boolean = True
                            counter+=1
                            break
                        else:
                            pass
                if boolean == True:
                    break
                else:
                    pass
    scores.append(counter)
    print(i)
    #print(counter)


possible_genome_list = []
max_genome_matched = 0
max_matches = 0
for i in range(len(scores)):
    score = scores[i]
    if score>= 50000:
        possible_genome_list.append(i)
    else:
        pass
    if score > max_matches:
        max_matches = score
        max_genome_matched = i
    else:
        pass

for read in read_matches:
    if len(read_matches[read]) == 0:
        read_matches[read] = [max_genome_matched]
    read_matches[read] = [x for x in read_matches[read] if x in possible_genome_list]
    if len(read_matches[read]) > 1:
        read_matches[read] = [read_matches[read][0]]
    if len(read_matches[read]) == 0:
        read_matches[read] = [max_genome_matched]
    if len(read_matches[read]) != 1:
        print('something wrong')
    else:
        pass

output_file = 'predictions.txt'
output(output_file, read_matches)

print('finished')


#around 50k matches for genomes that aren't present in sample