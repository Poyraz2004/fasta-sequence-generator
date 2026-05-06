# ============================================================
# Album number : s29979
# Date         : 2026-05-06
# Description  : Random nucleotide sequence generator in FASTA format.
#                The program creates a DNA sequence, calculates statistics,
#                and provides additional bioinformatics analyses.
# ============================================================



import random
import csv




def generate_sequence(length: int) -> str:

    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(nucleotides) for _ in range(length))


def calculate_stats(sequence: str) -> dict:


    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0
    for char in sequence:
        if char in counts:
            counts[char] += 1
            total += 1

    if total == 0:
        return {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'GC': 0.0}

    stats = {nuc: round((counts[nuc] / total) * 100, 2) for nuc in counts}
    stats['GC'] = round(((counts['G'] + counts['C']) / total) * 100, 2)
    return stats


def insert_name(sequence: str, name: str) -> str:

    name_lower = name.lower()
    # Choose a random position between 0 and len(sequence)
    position = random.randint(0, len(sequence))
    return sequence[:position] + name_lower + sequence[position:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:


    if description:
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"

    # Split sequence into lines of 80 characters
    lines = []
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    return header + '\n' + '\n'.join(lines) + '\n'


def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    while True:
        raw = input(prompt)
        try:
            value = int(raw)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


def validate_id(prompt: str) -> str:

    while True:
        seq_id = input(prompt).strip()
        if not seq_id:
            print("Error: ID cannot be empty.")
        elif any(c.isspace() for c in seq_id):
            print("Error: ID cannot contain whitespace.")
        else:
            return seq_id



def find_motif(sequence: str, motif: str) -> list:

    clean_seq = ''.join(c for c in sequence if c.isupper())
    positions = []
    start = 0
    motif_upper = motif.upper()
    while True:
        idx = clean_seq.find(motif_upper, start)
        if idx == -1:
            break
        positions.append(idx + 1)  # 1-based index
        start = idx + 1
    return positions


def run_motif_search(sequence: str) -> None:

    motif = input("Enter the motif to search for (e.g. ATG): ").strip().upper()
    if not motif:
        print("Motif cannot be empty, search skipped.")
        return
    positions = find_motif(sequence, motif)
    if positions:
        print(f"Motif '{motif}' was found {len(positions)} time(s). Positions (1-based):")
        print(', '.join(map(str, positions)))
    else:
        print(f"Motif '{motif}' was not found in the sequence.")



COMPLEMENT_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')


def complement(sequence: str) -> str:

    return sequence.translate(COMPLEMENT_TABLE)


def reverse_complement(sequence: str) -> str:

    return complement(sequence)[::-1]


def add_complement_records(seq_id: str, description: str, sequence: str) -> str:

    clean_seq = ''.join(c for c in sequence if c.isupper())
    comp_seq = complement(clean_seq)
    rev_comp_seq = reverse_complement(clean_seq)

    comp_record = format_fasta(
        f"{seq_id}_complement",
        f"Complement of {description}" if description else f"Complement of {seq_id}",
        comp_seq
    )
    rev_comp_record = format_fasta(
        f"{seq_id}_revcomp",
        f"Reverse complement of {description}" if description else f"Reverse complement of {seq_id}",
        rev_comp_seq
    )
    return comp_record + rev_comp_record



def transcribe_to_mrna(sequence: str) -> str:

    clean_seq = ''.join(c for c in sequence if c.isupper())
    return clean_seq.replace('T', 'U')


def add_mrna_record(seq_id: str, description: str, sequence: str) -> str:

    mrna = transcribe_to_mrna(sequence)
    return format_fasta(
        f"{seq_id}_mRNA",
        f"mRNA transcript of {description}" if description else f"mRNA transcript of {seq_id}",
        mrna
    )


def sliding_window_gc(sequence: str, window_size: int, step: int = 1) -> list:

    clean_seq = ''.join(c for c in sequence if c.isupper())
    results = []
    for start in range(0, len(clean_seq) - window_size + 1, step):
        window = clean_seq[start:start + window_size]
        gc_count = window.count('G') + window.count('C')
        gc_pct = round((gc_count / window_size) * 100, 2)
        results.append((start + 1, gc_pct))  # 1-based start position
    return results


def save_sliding_window_csv(seq_id: str, results: list) -> str:

    filename = f"{seq_id}_gc_window.csv"
    with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['start_position', 'gc_content'])
        writer.writerows(results)
    return filename


def run_sliding_window(seq_id: str, sequence: str) -> None:

    clean_len = sum(1 for c in sequence if c.isupper())
    if clean_len < 2:
        print("The sequence is too short for sliding window analysis.")
        return

    window_size = validate_positive_int(
        f"Enter sliding window width [1, {clean_len}]: ",
        min_val=1,
        max_val=clean_len
    )
    results = sliding_window_gc(sequence, window_size)
    csv_file = save_sliding_window_csv(seq_id, results)
    print(f"Sliding window GC analysis completed. Results saved to: {csv_file}")
    # Preview the first few rows
    preview_count = min(5, len(results))
    print(f"  First {preview_count} window(s):")
    for pos, gc in results[:preview_count]:
        print(f"    Position {pos}: {gc:.2f}% GC")


# Codon table (standard genetic code)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def find_orfs(sequence: str, min_length: int = 100) -> list:

    clean_seq = ''.join(c for c in sequence if c.isupper())
    orfs = []


    for frame in range(3):
        i = frame
        while i < len(clean_seq) - 2:
            codon = clean_seq[i:i + 3]
            if codon == 'ATG':

                start = i
                j = i + 3
                found_stop = False
                while j < len(clean_seq) - 2:
                    stop_codon = clean_seq[j:j + 3]
                    if stop_codon in STOP_CODONS:
                        end = j + 3
                        orf_length = end - start
                        if orf_length >= min_length:
                            orfs.append((
                                start + 1,
                                end,
                                orf_length,
                                frame + 1
                            ))
                        found_stop = True
                        i = j + 3
                        break
                    j += 3
                if not found_stop:
                    i += 3
            else:
                i += 3

    # Sort by start position
    orfs.sort(key=lambda x: x[0])
    return orfs


def run_orf_finder(sequence: str) -> None:

    clean_len = sum(1 for c in sequence if c.isupper())
    min_len = validate_positive_int(
        "Enter minimum ORF length in nucleotides [1, 10000]: ",
        min_val=1,
        max_val=min(10000, clean_len)
    )
    orfs = find_orfs(sequence, min_length=min_len)
    if orfs:
        print(f"\n{len(orfs)} ORF(s) found (minimum length: {min_len} nt):")
        for start, end, length, frame in orfs:
            print(f"  Frame +{frame}: position {start}–{end}, length {length} nt")
    else:
        print(f"No ORFs with minimum length {min_len} nt were found.")


def print_stats(stats: dict, length: int) -> None:

    print(f"\nSequence statistics (n={length}):")
    for nuc in ['A', 'C', 'G', 'T']:
        print(f"  {nuc}: {stats[nuc]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")


def main():

    print("=" * 50)
    print("  FASTA Sequence Generator")
    print("=" * 50)

    # ── Get input parameters ──────────────────────────
    length = validate_positive_int("Enter sequence length: ")
    seq_id = validate_id("Enter sequence ID: ")
    description = input("Enter sequence description (can be empty): ").strip()
    name = input("Enter your name: ").strip()


    sequence = generate_sequence(length)


    sequence_with_name = insert_name(sequence, name)

    stats = calculate_stats(sequence_with_name)
    print_stats(stats, length)


    fasta_content = format_fasta(seq_id, description, sequence_with_name)


    fasta_content += add_complement_records(seq_id, description, sequence_with_name)


    fasta_content += add_mrna_record(seq_id, description, sequence_with_name)


    output_file = f"{seq_id}.fasta"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(fasta_content)

    print(f"\nSequence saved to file: {output_file}")


    print(f"\nFirst lines of the '{output_file}' file:")
    lines = fasta_content.splitlines()
    for line in lines[:4]:
        print(f"  {line}")
    if len(lines) > 4:
        print("  ...")


    print("\n--- Motif Search ---")
    run_motif_search(sequence_with_name)


    print("\n--- Sliding Window GC Analysis ---")
    clean_len = sum(1 for c in sequence_with_name if c.isupper())
    if clean_len >= 2:
        run_sliding_window(seq_id, sequence_with_name)
    else:
        print("The sequence is too short for GC analysis.")


    print("\n--- ORF Detection ---")
    run_orf_finder(sequence_with_name)

    print("\nProgram completed.")


if __name__ == "__main__":
    main()