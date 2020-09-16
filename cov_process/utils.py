# Determine if the mutation is a transition mutation.
def is_transition(ori_allele, alt_allele):
    p = False
    transition = [("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")]
    transversion = [("C", "A"), ("C", "G"), ("G", "T"), ("T", "G"), ("G", "C"), ("A", "C"), ("A", "T"), ("T", "A")]
    if (ori_allele, alt_allele) in transition:
        p = True
    elif (ori_allele, alt_allele) in transversion:
        p = False
    else:
        print ori_allele, alt_allele
        raise CategError("NEITHER TRANSITION NOR TRANSVERSION MUTATION")

    return p


# if base in gc dinucleotide
def is_cpg(base, base_before, base_after):
    flag = False
    if (base + base_after == gc_di_nuc) or (base_before + base == gc_di_nuc):
        flag = True

    return flag


# if base is A:T pair
def is_at(base):
    flag = False
    if base == "A" or base == "T":
        flag = True

    return flag


# if a mutation cause silent codon
def is_silent(before_context, after_context):
    before_aa = table[before_context]
    after_aa = table[after_context]
    flag = False
    if before_aa == after_aa:
        flag = True

    return flag


# just add count to whatever dictionary passed
def add_count(base, base_before, base_after, mut_bases, dict_cov):
    dict_cov[6] += float(1) / 3
    if is_at(base):
        if is_transition(base, mut_bases):
            dict_cov[4] += float(1) / 3
        else:
            dict_cov[5] += float(1) / 3
    elif is_cpg(base, base_before, base_after):  # This is a CG dinucleotide
        if is_transition(base, mut_bases):
            dict_cov[0] += float(1) / 3
        else:
            dict_cov[1] += float(1) / 3
    else:
        if is_transition(base, mut_bases):
            dict_cov[2] += float(1) / 3
        else:
            dict_cov[3] += float(1) / 3

    return dict_cov


# Transform all the elements in the input list into numeric.
def list_tonumeric(test_list):
    test_list = list(map(int, test_list))
    return test_list


# Group the positions in a list into position start and end pairs in a gene.
def list_to_pairs(list_old):
    new_list = []
    for i in range(0, len(list_old), 2):
        new_list.append(list_old[i:i + 2])
    return new_list


# This code is for process sequences on negative strand For the sequences in list, this code will get the reverse
# complement for each sequence and then reverse the sequence order.
def negative_process(list_sequence):
    new_list = []
    for sequences in list_sequence:
        sequences = sequences.reverse_complement()
        new_list.append(sequences)
    new_list = new_list[::-1]

    return new_list


# Find out if the base fall into non-coding region or coding region.
def get_region(region_list, base_pos, region_f):
    r_flag = "nc"
    j_flag = 0
    for regions in region_list:
        if regions[1] - 1 >= base_pos >= regions[0]:
            r_flag = region_f
            if r_flag == "coding":
                if base_pos == regions[1] - 1 or base_pos == regions[1] - 2 or base_pos == regions[0] or base_pos == \
                        regions[0] + 1:
                    j_flag = 1
            break

    return r_flag, j_flag


def get_annotation(annotation_file_path):
    name_dict = {}  # {transcript:[gene_name]}
    transcript_dict = {}  # {'gene_name':['ENST00000']}
    transcript_info_dict = {}  # {transcriptID:'exon'[]; 'CDS'[];'gene'[];'UTR

    with open(annotation_file_path, 'r') as annotation_f:
        for lines in annotation_f:
            if lines.startswith("chr"):
                line_split = lines.split('\t')
                chr_n, source, feature, start_pos, end_pos = line_split[0:5]
                strand, frame, attribute = line_split[6:9]
                position_list = [int(start_pos), int(end_pos)]

                # feature = line_split[2] # 'exon' / 'CDS'/ 'UTR'/ 'start_codon'/ 'stop_codon'/ 'transcript'/ 'gene'

                # if this is a protein coding gene and is the primary transcript:
                if 'protein_coding' in attribute and "KNOWN" in attribute and 'appris_principal' in attribute:
                    # parse the attributes
                    attribute_split = attribute.split(';')
                    gene_id_col, transcript_id_col, gene_type_col, gene_status_col, gene_name_col, transcript_type_col, \
                    transcript_status_col, transcript_name_col = attribute_split[0:8]
                    # not include the version number
                    gene_id = re.findall(r'(ENSG\d+|ENSGR\d+)', gene_id_col)[0]
                    transcript_id = re.findall(r'(ENST\d+|ENSTR\d+)', transcript_id_col)[0]
                    gene_status = re.findall(r'\"(.*?)\"', gene_status_col)[0]  # KNOWN or not
                    gene_name = re.findall(r'\"(.*?)\"', gene_name_col)[0]
                    transcript_status = re.findall(r'\"(.*?)\"', transcript_status_col)[0]  # KNOWN or not

                    if gene_name not in transcript_dict:
                        transcript_dict[gene_name] = []
                    if transcript_id not in transcript_dict[gene_name]:
                        transcript_dict[gene_name].append(transcript_id)

                    if transcript_id not in transcript_info_dict:
                        transcript_info_dict[transcript_id] = {}
                        transcript_info_dict[transcript_id]['strand'] = strand
                        transcript_info_dict[transcript_id]['chr'] = chr_n
                        transcript_info_dict[transcript_id]['exon'] = []
                        transcript_info_dict[transcript_id]['CDS'] = []
                        transcript_info_dict[transcript_id]['UTR'] = []
                        transcript_info_dict[transcript_id]['transcript'] = []

                        # if the transcript and genes are known, parse according to feature type
                    if feature == 'exon':
                        # exon_n = re.findall(r'\d+', attribute_split[8])
                        transcript_info_dict[transcript_id]['exon'].append(position_list)
                    elif feature == 'CDS':
                        transcript_info_dict[transcript_id]['CDS'].append(position_list)
                    elif feature == 'UTR':
                        transcript_info_dict[transcript_id]['UTR'].append(position_list)
                    elif feature == 'transcript':
                        transcript_info_dict[transcript_id]['transcript'] = position_list

                    if transcript_id not in name_dict:
                        name_dict[transcript_id] = gene_name

    # delete the transcript record that have 2 or more principle transcripts, only keep the longest transcript
    for names in transcript_dict:
        if len(transcript_dict[names]) > 1:
            store_max_t = 0
            for t in transcript_dict[names]:
                len_t = abs(transcript_info_dict[t]['transcript'][1] - transcript_info_dict[t]['transcript'][0])
                if len_t >= store_max_t:
                    store_max_t = len_t
                else:
                    del transcript_info_dict[t]

    return name_dict, transcript_info_dict


# get mrna positions in gene, then delete the utr positions
# the returned result is the cds positions
def get_mrna_position(transcript, info_dict, strand_gene):
    list_all_pos = info_dict[transcript]['UTR'] + info_dict[transcript]['CDS']
    list_all_pos.sort()

    # get the index of utr pair
    utr_index_list = []
    for pairs in info_dict[transcript]['UTR']:
        utr_idx = list_all_pos.index(pairs)
        utr_index_list.append(utr_idx)

    # if negative strand, the start position is the biggest position
    list_gene_position = []
    if strand_gene == '-':
        transcript_start = list_all_pos[-1][1]
        for pairs in list_all_pos:
            for positions in pairs:
                gene_position = -(positions - transcript_start)
                list_gene_position.append(gene_position)
    else:
        transcript_start = list_all_pos[0][0]
        for pairs in list_all_pos:
            for positions in pairs:
                gene_position = positions - transcript_start
                list_gene_position.append(gene_position)

    # put the mrna positions into pairs
    list_gene_position.sort()
    mrna_list = list_to_pairs(list_gene_position)
    # negative strand need to turn the utr index list around
    if strand_gene == '-':
        exon_n = len(mrna_list) - 1
        utr_index_list = [abs(i - exon_n) for i in utr_index_list]

    # remove utr positions according to the utr index
    for utr_index in sorted(utr_index_list,
                            reverse=True):  # reverse the list so that index in list won't change after deletion
        del mrna_list[utr_index]

    return mrna_list


# get the transcript sequence
def get_transcript_sequence(transcript, info_dict, fasta_dict, strand_gene):
    chromosome = info_dict[transcript]['chr']
    transcript_start_pos = info_dict[transcript]['transcript'][0]
    transcript_end_pos = info_dict[transcript]['transcript'][1]
    if strand_gene == '-':
        transcript_seq = fasta_dict[chromosome][transcript_start_pos - 1:transcript_end_pos].reverse_complement()
    else:
        transcript_seq = fasta_dict[chromosome][transcript_start_pos - 1:transcript_end_pos]

    return transcript_seq


# get the cdna sequence from the transcript sequence, inputs are cdna positions within the gene
def get_cdna_sequence(cds_list, strand_gene, sequence_transcript):
    list_cds_seq = []
    for cds_pairs in cds_list:
        cds_start = cds_pairs[0]
        cds_end = cds_pairs[1]
        if strand_gene == '-':
            list_cds_seq.append(sequence_transcript[cds_start:cds_end + 1])
        else:
            list_cds_seq.append(sequence_transcript[cds_start:cds_end + 1])

    concatenated_cds = Seq("", generic_dna)
    for s in list_cds_seq:
        concatenated_cds += s

    return concatenated_cds.seq