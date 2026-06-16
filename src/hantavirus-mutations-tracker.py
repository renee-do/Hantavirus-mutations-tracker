import Bio
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Input user email. This is necessary to search the NCBI database using Bio.Entrez.
user_email = input('Enter user email.')

Entrez.email = user_email

# Search NCBI database for virus and genome.
# Search: Andes virus [Organism] AND complete genome
search_term = input('Define search term.')

handle = Entrez.esearch(db = 'nucleotide', term = search_term, retmax = 15)

search_results = Entrez.read(handle)
handle.close()

id_list = search_results['IdList']
print('Found {} sequences matching the query.'.format(len(id_list)))
print('The NCBI IDs are: {}'.format(id_list))

# Fetch IDs of sequences and save to new .fasta file.
if id_list:
    ids_to_fetch = ','.join(id_list)

    fetch_handle = Entrez.efetch(db = 'nucleotide', id = ids_to_fetch, rettype = 'fasta', retmode = 'text')

    raw_fasta_data = fetch_handle.read()
    fetch_handle.close()

    output_file = input('Name of output file:')

    fasta_sequences = open(output_file, 'w')
    fasta_sequences.write(raw_fasta_data)
    fasta_sequences.close

    print('Success. Files saved to {}.'.format(output_file))
else:
    print('No IDs found to fetch.')

# Parse through sequences .fasta file. 
for record in SeqIO.parse(output_file, 'fasta'):
    print('Sequence ID: {}'.format(record.id))
    print('Description: {}'.format(record.description))
    print('Length of Sequence: {} base pairs'.format(len(record.seq)))
    print(repr(record.seq))
    print('-' * 50)

isolated_s_segments = []

for record in SeqIO.parse(output_file, 'fasta'):
    description = record.description.lower()

    if ('segment s' or 's segment' or 'small segment' or 'segment small') in description:
        isolated_s_segments.append(record)

SeqIO.write(isolated_s_segments, 'isolated_s_segments.fasta', 'fasta')
print('There are {} S segments for analysis'.format(len(isolated_s_segments)))

# Andes Hantavirus Segment S = NC_003466(.1)
# Compare segments with reference sequence using remote BLAST. 
print('Submitting sequences to NCBI BLAST servers.')

sequence_file = input('Name of file containing isolated sequences:')
isolated_sequences = open(sequence_file, 'r')
query_sequences = isolated_sequences.read()

result_handle = NCBIWWW.qblast(program = 'blastn', database = 'nt', sequence = query_sequences, entrez_query = 'NC_003466 [Accession]')

output_file = input('Name of output file:')

alignment_results = open(output_file, 'w')
alignment_results.write(result_handle.read())
alignment_results.close
print('Success. BLAST results saved to {}.'.format(output_file))

# Parse through BLAST results to find alignments and view mutations.
result_handle = open(output_file)
blast_records = NCBIXML.parse(result_handle)

mutation_data = []

for record in blast_records:
    query_id = record.query

    for alignment in record.alignments:
        for hsp in alignment.hsps:
            query_seq = hsp.query
            match_seq = hsp.match
            sbjct_seq = hsp.sbjct

            ref_start_pos = hsp.sbjct_start

            for index, match_char in enumerate(match_seq):
                if match_char == ' ':
                    actual_genome_pos = ref_start_pos + index
                    
                    ref_base = sbjct_seq[index]
                    variant_base = query_seq[index]

                    if ref_base != '-' and variant_base != '-':
                        mutation_data.append({'Variant_ID': query_id, 'Genome_Position': actual_genome_pos, 'Ref_Base': ref_base, 'Variant_Base': variant_base, 'Mutation_Type': '{} -> {}'.format(ref_base, variant_base)})

result_handle.close()

df_mutations = pd.DataFrame(mutation_data)

df_mutations.to_csv('Hantavirus_mutations.csv', index = False)
print('Found {} mutations.'.format(len(df_mutations)))
print(df_mutations.head())

df = pd.read_csv('Hantavirus_mutations.csv')

window_size = 50
df['Genome_Window'] = (df['Genome_Position'] // window_size) * window_size

hotspots = df.groupby('Genome_Window').size().reset_index(name = 'Mutation_Count')

# Visualize data.
sns.set_theme(style = 'whitegrid')
plt.figure(figsize = (12, 6))

plt.scatter(hotspots['Genome_Window'], hotspots['Mutation_Count'], color = '#b2182b', label = 'Mutation_ Count')
plt.plot(hotspots['Genome_Window'], hotspots['Mutation_Count'], color = '#b2182b', linewidth = 2, label = 'Mutation_ Count')
plt.fill_between(hotspots['Genome_Window'], hotspots['Mutation_Count'], color = '#b1282b', alpha = 0.15)

plt.title('Hantavirus S-Segment Mutation Hotspots', fontsize = 16, fontweight = 'bold',  pad = 15)
plt.xlabel('Genome Position (Base Pairs)', fontsize = 12, labelpad = 10)
plt.ylabel('Number of Mutations Detected', fontsize = 12, labelpad = 10)

plt.tight_layout()
plt.savefig('Hantavirus_mutations_hotspots.png', dpi = 300)

print('Top 5 Mutation Types')
type_counts = df['Mutation_Type'].value_counts().head(5)
print(type_counts)
