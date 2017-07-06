import sys

query_lines = []

total_length = 0
total_error_count = 0

with open(sys.argv[1], 'rt') as blast_hits:
    for line in blast_hits:
        line_parts = line.split('\t')
        percent_identity = float(line_parts[2])
        alignment_length = int(line_parts[3])

        if percent_identity < 95.0:
            continue
        if alignment_length < 9000:
            continue

        error_count = round(alignment_length * (100.0 - percent_identity) / 100.0)

        total_length += alignment_length
        total_error_count += error_count

overall_percent_identity = 100.0 * (1.0 - (total_error_count / total_length))
overall_error_rate = 100.0 - overall_percent_identity

try:
    mean_distance_between_errors = 1.0 / (1.0 - (overall_percent_identity / 100.0))
    mean_distance_between_errors_str = str(round(mean_distance_between_errors))
except ZeroDivisionError:
    mean_distance_between_errors_str = 'inf'

print('%.4f' % overall_percent_identity + '\t' + '%.4f' % overall_error_rate + '\t' + mean_distance_between_errors_str, flush=True, end='')
