query_lines = []
with open('blast_queries.fasta', 'rt') as queries:
    for line in queries:
        query_lines.append(line.strip())

with open('blast_queries2.fasta', 'wt') as new_queries:
    for line in query_lines:
        if line.startswith('>'):
            new_queries.write(line)
        else:
            new_queries.write(line[1000:-1000])
        new_queries.write('\n')
