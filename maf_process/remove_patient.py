list_p = []
with open('/Users/stella/Downloads/MutSigrelated_local_input/Biliary-AdenoCA_merged.coverage.txt', 'r') as f:
    list_p = next(f).split('\t')[3:]


with open('/Users/stella/Downloads/MutSigrelated_local_input/merged_consensus_passonly_snv_mnv_indel_Biliary-AdenoCA_05252020.maf.output.txt', 'r') as m:
    with open('/Users/stella/Downloads/MutSigrelated_local_input/merged_Biliary-AdenoCA.maf',
            'wb') as output:
        output.write(next(m))
        for lines in m:
            patient = lines.split('\t')[12]
            # patient = patient.replace('-','')
            if patient in list_p:
                print 'T'
                output.write(lines)

