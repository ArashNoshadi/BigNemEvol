import pandas as pd
import os
import csv
import sys
import itertools
from Bio.Seq import Seq
from collections import defaultdict

# ==================== 1. CONFIGURATION & PARAMETERS ====================
INPUT_FILES_LIST = [
    r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\Silva + NCBI\mearg\counsensus\ITS1.xlsx",
    r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\Silva + NCBI\mearg\counsensus\ITS2.xlsx",
    r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\Silva + NCBI\mearg\counsensus\18S.xlsx",
    r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\Silva + NCBI\mearg\counsensus\28S.xlsx",
    r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\Silva + NCBI\mearg\counsensus\5.8S.xlsx",
]

OUTPUT_DIR = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\Silva + NCBI\mearg\counsensus\Final_counsensus"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TEMP_CSV_PATH = os.path.join(OUTPUT_DIR, 'live_results_temp.csv')
FINAL_EXCEL_PATH = os.path.join(OUTPUT_DIR, 'final_concatenated_consensus.xlsx')
FINAL_FASTA_PATH = os.path.join(OUTPUT_DIR, 'final_sequences.fasta')
STATS_EXCEL_PATH = os.path.join(OUTPUT_DIR, 'Overlap_vs_Concat_Statistics.xlsx')

MARKER_MAP = {
    'ITS1.xlsx': 'ITS1', 'ITS2.xlsx': 'ITS2', '18S.xlsx': '18S',
    '28S.xlsx': '28S', '5.8S.xlsx': '5.8S',
}

ORDERED_MARKERS = ['18S', 'ITS1', '5.8S', 'ITS2', '28S']
MARKER_RANK = {m: i for i, m in enumerate(ORDERED_MARKERS)}

ESTIMATED_LENGTHS = {
    '5.8S': 160, 'ITS2': 400, 'ITS1': 730, '18S': 1700, '28S': 3500
}

SEARCH_WINDOW = 800
MIN_OVERLAP_LEN = 40
MAX_MISMATCH_PCT = 0.10

IUPAC_CODES = {
    frozenset(['A']): 'A', frozenset(['C']): 'C', frozenset(['G']): 'G', frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
    frozenset(['A', 'C']): 'M', frozenset(['G', 'T']): 'K',
    frozenset(['G', 'C']): 'S', frozenset(['A', 'A', 'T']): 'W',
    frozenset(['A', 'C', 'T']): 'H', frozenset(['G', 'C', 'T']): 'B',
    frozenset(['A', 'G', 'C']): 'V', frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'G', 'C', 'T']): 'N'
}

# ==================== 2. CORE CLASSES ====================
class Contig:
    def __init__(self, sequence, markers_contained, history_log, connections=None):
        self.seq = sequence
        self.markers = set(markers_contained)
        self.log = history_log
        self.connections = connections if connections is not None else {}

    def min_rank(self):
        ranks = [MARKER_RANK[m] for m in self.markers if m in MARKER_RANK]
        return min(ranks) if ranks else 999
    
    def max_rank(self):
        ranks = [MARKER_RANK[m] for m in self.markers if m in MARKER_RANK]
        return max(ranks) if ranks else -1

# ==================== 3. HELPER FUNCTIONS ====================
def get_iupac_consensus(b1, b2):
    b1, b2 = b1.upper(), b2.upper()
    if b1 == b2: return b1
    if b1 in ['-', 'N', '']: return b2
    if b2 in ['-', 'N', '']: return b1
    valid = set('ACGT')
    if b1 not in valid or b2 not in valid: return 'N'
    return IUPAC_CODES.get(frozenset([b1, b2]), 'N')

def clean_dna(val):
    if pd.isna(val): return ""
    val = str(val).upper().replace('\n', '').replace(' ', '').replace('\r', '').replace('\xa0', '').replace('\u200b', '')
    seq = "".join([c for c in val if c in "ATGCNRYMKSWBDHV-"]).strip('N')
    return seq

def safe_str(val):
    return str(val).strip().replace(' ', '_').replace('|', '-').replace(':', '').replace('/', '-').replace('\\', '-')

def get_normalized_key(row):
    fields = ['Main_Organism', 'Family', 'ID_Name', 'Level']
    return tuple((str(row.get(f, '')).strip().lower() if str(row.get(f, '')).strip().lower() != 'nan' else '') for f in fields)

def is_biologically_sound(contig_a, contig_b):
    ranks_a = set([MARKER_RANK[m] for m in contig_a.markers if m in MARKER_RANK])
    ranks_b = set([MARKER_RANK[m] for m in contig_b.markers if m in MARKER_RANK])
    if not ranks_a or not ranks_b: return True 
    for r_a in ranks_a:
        for r_b in ranks_b:
            if abs(r_a - r_b) <= 1: return True
    return False

def consolidate_marker_sequences(marker_name, sequences_list):
    if not sequences_list: return ""
    if len(sequences_list) == 1: return sequences_list[0]
    pool = [Contig(seq, [marker_name], [f"Part_{i}"]) for i, seq in enumerate(sequences_list)]
    final_pool = run_assembly_core(pool, check_biology=False) 
    return max(final_pool, key=lambda c: len(c.seq)).seq if final_pool else max(sequences_list, key=len)

# ==================== 4. OVERLAP LOGIC ====================
def check_overlap_directional(seq1, seq2):
    def _one_way(s1, s2):
        max_check = min(len(s1), len(s2), SEARCH_WINDOW)
        for ov in range(max_check, MIN_OVERLAP_LEN - 1, -1):
            suffix, prefix = s1[-ov:], s2[:ov]
            if all(c in 'N-' for c in suffix) or all(c in 'N-' for c in prefix): continue
            mismatches, consensus_part, possible = 0, [], True
            for c1, c2 in zip(suffix, prefix):
                consensus_part.append(get_iupac_consensus(c1, c2))
                if c1 != c2 and c1 not in 'N-' and c2 not in 'N-':
                    mismatches += 1
                    if mismatches > (ov * MAX_MISMATCH_PCT):
                        possible = False
                        break
            if possible:
                return s1[:-ov] + "".join(consensus_part) + s2[ov:], ov, f"{ov}bp_Mis:{mismatches}"
        return None, 0, ""

    m, l, log = _one_way(seq1, seq2)
    if m: return m, l, f"Fwd_Join({log})"
    m_rc, l_rc, log_rc = _one_way(seq1, str(Seq(seq2).reverse_complement()))
    if m_rc: return m_rc, l_rc, f"Rev_Join({log_rc})"
    return None, 0, ""

# ==================== 5. ASSEMBLY PHASE ====================
def run_assembly_core(pool, check_biology=True):
    while True:
        best_merge_seq, best_merge_info, best_merge_len, best_indices = None, "", -1, None
        n = len(pool)
        if n < 2: break

        for idx_a, idx_b in itertools.permutations(range(n), 2):
            c_a, c_b = pool[idx_a], pool[idx_b]
            if check_biology and not is_biologically_sound(c_a, c_b): continue

            merged, ov_len, info = check_overlap_directional(c_a.seq, c_b.seq)
            if merged and ov_len > best_merge_len:
                best_merge_len, best_merge_seq, best_merge_info, best_indices = ov_len, merged, info, (idx_a, idx_b)

        if best_indices:
            idx_a, idx_b = best_indices
            contig_a, contig_b = pool[idx_a], pool[idx_b]
            
            max_a, min_a = contig_a.max_rank(), contig_a.min_rank()
            max_b, min_b = contig_b.max_rank(), contig_b.min_rank()
            pair_name = f"{ORDERED_MARKERS[max_a]}-{ORDERED_MARKERS[min_b]}" if max_a < min_b and max_a != -1 and min_b != 999 else (f"{ORDERED_MARKERS[max_b]}-{ORDERED_MARKERS[min_a]}" if max_b < min_a and max_b != -1 and min_a != 999 else None)

            new_conn = {**contig_a.connections, **contig_b.connections}
            if pair_name: new_conn[pair_name] = 'Overlap'

            new_contig = Contig(best_merge_seq, contig_a.markers.union(contig_b.markers), contig_a.log + contig_b.log + [f"Merge[{best_merge_info}]"], new_conn)
            pool.pop(max(idx_a, idx_b))
            pool.pop(min(idx_a, idx_b))
            pool.append(new_contig)
        else: break
    return pool

def run_greedy_assembly(markers_dict):
    pool = [Contig(seq, [m_name], [f"Origin({m_name})"]) for m_name, seq in markers_dict.items() if seq]
    return run_assembly_core(pool, check_biology=True) if pool else []

# ==================== 6. SCAFFOLDING PHASE ====================
def scaffold_contigs(contigs):
    if not contigs: return "", "", "", {}
    sorted_contigs = sorted(contigs, key=lambda c: c.min_rank())
    
    curr_c = sorted_contigs[0]
    final_seq, final_log, final_markers = curr_c.seq, list(curr_c.log), [f"[{'+'.join(curr_c.markers)}]"]
    covered_markers, final_connections = set(curr_c.markers), dict(curr_c.connections)

    first_rank = curr_c.min_rank()
    if first_rank > 0:
        prefix_gap = "".join("N" * ESTIMATED_LENGTHS.get(ORDERED_MARKERS[r], 100) for r in range(0, first_rank))
        final_seq = prefix_gap + final_seq
        final_log = [f"GapStart_{ORDERED_MARKERS[r]}" for r in range(0, first_rank)] + final_log

    for i in range(1, len(sorted_contigs)):
        next_contig = sorted_contigs[i]
        last_max_rank = max([MARKER_RANK[m] for m in covered_markers if m in MARKER_RANK], default=-1)
        next_min_rank = next_contig.min_rank()
        is_neighbor = (next_min_rank - last_max_rank) <= 1
        pair_name = f"{ORDERED_MARKERS[last_max_rank]}-{ORDERED_MARKERS[next_min_rank]}" if last_max_rank != -1 and next_min_rank != 999 else "Unknown_Pair"

        merged_res, ov_len, info = check_overlap_directional(final_seq, next_contig.seq) if is_neighbor else (None, 0, "")

        if merged_res:
            final_seq = merged_res
            final_log.extend([f"Scaffold_Merge[{info}]"] + next_contig.log)
            final_markers.append(f"[{'+'.join(next_contig.markers)}]")
            final_connections[pair_name] = 'Overlap'
        else:
            overlapping_markers = covered_markers.intersection(next_contig.markers)
            merged_successfully = False
            if overlapping_markers and is_neighbor:
                original_seq, MIN_KEEP_LEN = final_seq, 450
                for cut_len in range(20, 800, 20):
                    if len(original_seq) - cut_len < MIN_KEEP_LEN: break
                    m_res, ov_len_new, info_new = check_overlap_directional(original_seq[:-cut_len], next_contig.seq)
                    if m_res:
                        final_seq = m_res
                        final_log.extend([f"Fixed_Duplicate_Trim_{cut_len}bp[{info_new}]"] + next_contig.log)
                        final_markers.append(f"[{'+'.join(next_contig.markers)}]")
                        final_connections[pair_name] = 'Overlap'
                        merged_successfully = True
                        break
                if not merged_successfully:
                    safe_trim = max(0, min(300, len(original_seq) - MIN_KEEP_LEN))
                    final_seq = original_seq[:-safe_trim] + "NNNNN" + next_contig.seq
                    final_log.extend([f"Forced_Join_Safe[{list(overlapping_markers)}]"] + next_contig.log)
                    final_markers.append(f"[{'+'.join(next_contig.markers)}]")
                    final_connections[pair_name] = 'Concat'
            else:
                last_max = max([MARKER_RANK[m] for m in covered_markers if m in MARKER_RANK], default=-1)
                gap_seq = ""
                if next_min_rank > last_max + 1:
                    gap_seq = "".join("N" * ESTIMATED_LENGTHS.get(ORDERED_MARKERS[r], 100) for r in range(last_max + 1, next_min_rank))
                    final_log.extend(f"Gap_{ORDERED_MARKERS[r]}" for r in range(last_max + 1, next_min_rank))
                elif not gap_seq:
                    gap_seq = "N" * 50
                    final_log.append("Gap_Spacer_50N")
                final_seq += gap_seq + next_contig.seq
                final_log.extend(next_contig.log)
                final_markers.append(f"[{'+'.join(next_contig.markers)}]")
                final_connections[pair_name] = 'Concat'
        covered_markers.update(next_contig.markers)
        final_connections.update(next_contig.connections)

    if covered_markers:
        last_max = max([MARKER_RANK[m] for m in covered_markers if m in MARKER_RANK], default=-1)
        if last_max < len(ORDERED_MARKERS) - 1:
            final_seq += "".join("N" * ESTIMATED_LENGTHS.get(ORDERED_MARKERS[r], 100) for r in range(last_max + 1, len(ORDERED_MARKERS)))
            final_log.extend(f"GapEnd_{ORDERED_MARKERS[r]}" for r in range(last_max + 1, len(ORDERED_MARKERS)))

    return final_seq, "|".join(final_markers), "; ".join(final_log), final_connections

# ==================== 7. MAIN EXECUTION & STATS ====================
def main():
    print("Step 1: Loading Data...")
    marker_lengths = defaultdict(list)
    tax_data = {}

    for fpath in INPUT_FILES_LIST:
        fname = os.path.basename(fpath)
        if fname not in MARKER_MAP: continue
        m_name = MARKER_MAP[fname]
        try:
            df = pd.read_excel(fpath)
            seq_col = next((c for c in df.columns if 'seq' in c.lower() or 'consensus' in c.lower()), None)
            if not seq_col: continue
            for _, row in df.iterrows():
                key = get_normalized_key(row)
                if not any(key): continue
                clean_seq = clean_dna(row[seq_col])
                if len(clean_seq) > 20:
                    if key not in tax_data: tax_data[key] = {'info': row.to_dict(), 'markers': defaultdict(list)}
                    tax_data[key]['markers'][m_name].append(clean_seq)
                    if len(clean_seq) > len(''.join(tax_data[key]['markers'][m_name])[:-1] or ""):
                        tax_data[key]['info'] = row.to_dict()
        except Exception as e:
            print(f" ❌ Error reading {fname}: {e}")

    print(f"Loaded {len(tax_data)} organisms.")
    
    # متغیرهای آماری برای تولید گزارش دقیق برای استیو
    global_pair_stats = defaultdict(lambda: {'Overlap': 0, 'Concat': 0})
    seq_stats_rows = []
    marker_presence_counts = defaultdict(int)
    num_markers_counts = defaultdict(int)

    csv_file = open(TEMP_CSV_PATH, 'w', newline='', encoding='utf-8')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['Main_Organism', 'Family', 'ID_Name', 'Level', 'Consensus_Sequence', 'Length', 'Markers_Structure', 'Provenance', 'Long_Sequence_Flag'])
    fasta_file = open(FINAL_FASTA_PATH, 'w', encoding='utf-8')

    count = 0
    for key, data in tax_data.items():
        raw_markers_dict = data['markers']
        info = data['info']
        if not raw_markers_dict: continue

        processed_markers = {}
        for m_name, seq_list in raw_markers_dict.items():
            if seq_list:
                final_marker_seq = consolidate_marker_sequences(m_name, seq_list)
                if final_marker_seq and len(final_marker_seq) > 20:
                    processed_markers[m_name] = final_marker_seq
                    marker_lengths[m_name].append(len(final_marker_seq))

        if not processed_markers: continue

        org = str(info.get('Main_Organism', 'Unknown')).strip() or 'Unknown'
        family = str(info.get('Family', 'Unknown')).strip() or 'Unknown'
        id_name = str(info.get('ID_Name', 'NA')).strip() or 'NA'
        level = str(info.get('Level', 'NA')).strip() or 'NA'

        contigs = run_greedy_assembly(processed_markers)
        final_seq, m_str, l_str, organism_connections = scaffold_contigs(contigs)
        
        real_bases = len(final_seq.replace('N', '').replace('-', ''))
        if real_bases < 150: continue

        # جمع‌آوری آمارهای دقیق در سطح توالی
        total_len = len(final_seq)
        n_count = final_seq.count('N')
        n_pct = (n_count / total_len) * 100 if total_len > 0 else 0
        completeness = 100 - n_pct
        num_markers = len(processed_markers)
        markers_present = ", ".join(sorted(processed_markers.keys()))
        
        base_header = f"{safe_str(org)}_{safe_str(id_name)}_{safe_str(level)}"
        safe_header = f"{safe_str(family)}_{base_header}" if family != 'Unknown' else base_header

        seq_stats_rows.append({
            'Taxon_Name': safe_header,
            'Total_Length(bp)': total_len,
            'Real_Bases(bp)': real_bases,
            'N_Count': n_count,
            'N_Percentage(%)': round(n_pct, 2),
            'Completeness_Occupancy(%)': round(completeness, 2),
            'Num_Markers_Present': num_markers,
            'Markers_List': markers_present
        })

        for m in processed_markers.keys():
            marker_presence_counts[m] += 1
        num_markers_counts[num_markers] += 1

        for pair, j_type in organism_connections.items():
            if pair != "Unknown_Pair":
                global_pair_stats[pair][j_type] += 1

        long_flag = 1 if len(final_seq) > 13000 else 0
        csv_writer.writerow([org, family, id_name, level, final_seq, len(final_seq), m_str, l_str, long_flag])
        fasta_file.write(f">{safe_header}\n{final_seq}\n")

        count += 1
        if count % 50 == 0: print(f" ⏳ Processed {count}...")

    csv_file.close()
    fasta_file.close()

    # ==================== GENERATING REVIEWER STATS REPORT ====================
    print("\nGenerating Multi-sheet Reviewer Statistics Report...")
    
    # Dataframe 1: Supermatrix Completeness
    total_taxa = len(seq_stats_rows)
    total_bases_dataset = sum(r['Total_Length(bp)'] for r in seq_stats_rows)
    total_n_dataset = sum(r['N_Count'] for r in seq_stats_rows)
    overall_n_pct = (total_n_dataset / total_bases_dataset * 100) if total_bases_dataset > 0 else 0
    overall_matrix_occupancy = 100 - overall_n_pct

    df_supermatrix = pd.DataFrame([{
        'Total_Taxa_In_Matrix': total_taxa,
        'Total_Dataset_Length(bp)': total_bases_dataset,
        'Total_Missing_Data_Ns(bp)': total_n_dataset,
        'Overall_N_Percentage(%)': round(overall_n_pct, 2),
        'Overall_Matrix_Occupancy(%)': round(overall_matrix_occupancy, 2)
    }])

    # Dataframe 2: Connection Stats
    stats_rows = []
    for pair, counts in global_pair_stats.items():
        total = counts['Overlap'] + counts['Concat']
        if total > 0:
            stats_rows.append({
                'Marker_Pair': pair,
                'Total_Connections': total,
                'Real_Overlap_Count': counts['Overlap'],
                'Concatenated_Gap_Count': counts['Concat'],
                'Real_Overlap(%)': round((counts['Overlap'] / total) * 100, 2),
                'Concatenated(%)': round((counts['Concat'] / total) * 100, 2)
            })
    df_connections = pd.DataFrame(stats_rows).sort_values(by='Total_Connections', ascending=False) if stats_rows else pd.DataFrame()

    # Dataframe 3: Sequence Level Stats
    df_seq_stats = pd.DataFrame(seq_stats_rows).sort_values(by='Completeness_Occupancy(%)', ascending=False)

    # Dataframe 4: Marker Presence Summary
    marker_summary_rows = [
        {'Metric': f"Taxa with exactly {k} markers", 'Count': v, 'Percentage(%)': round((v/total_taxa)*100, 2)}
        for k, v in sorted(num_markers_counts.items(), reverse=True)
    ]
    marker_summary_rows.append({'Metric': '---', 'Count': '---', 'Percentage(%)': '---'})
    for m in ORDERED_MARKERS:
        if m in marker_presence_counts:
            count_m = marker_presence_counts[m]
            marker_summary_rows.append({
                'Metric': f"Total Taxa containing {m}", 
                'Count': count_m, 
                'Percentage(%)': round((count_m/total_taxa)*100, 2)
            })
    df_marker_summary = pd.DataFrame(marker_summary_rows)

    # Writing all to a single Excel file with multiple sheets
    try:
        with pd.ExcelWriter(STATS_EXCEL_PATH, engine='openpyxl') as writer:
            df_supermatrix.to_excel(writer, sheet_name='1_Supermatrix_Completeness', index=False)
            df_connections.to_excel(writer, sheet_name='2_Connections_Stats', index=False)
            df_seq_stats.to_excel(writer, sheet_name='3_Sequence_Level_Stats', index=False)
            df_marker_summary.to_excel(writer, sheet_name='4_Marker_Presence_Summary', index=False)
        print(f" 📊 Detailed Validation Report for Reviewer saved to:\n   {STATS_EXCEL_PATH}")
    except Exception as e:
        print(f" ❌ Could not save stats excel: {e}")
    # =========================================================================

    avg_lengths = {m: int(round(sum(marker_lengths[m]) / len(marker_lengths[m]))) if marker_lengths[m] else ESTIMATED_LENGTHS.get(m, 500) for m in ORDERED_MARKERS}
    
    with open(os.path.join(OUTPUT_DIR, 'charset_partitions.txt'), 'w', encoding='utf-8') as f:
        pos = 1
        for m in ORDERED_MARKERS:
            f.write(f"charset {m} = {pos}-{pos + avg_lengths[m] - 1}; ")
            pos += avg_lengths[m]

    try:
        pd.read_csv(TEMP_CSV_PATH).to_excel(FINAL_EXCEL_PATH, index=False)
        print(f"\n🎉 Done! Final files:\n 📄 {FINAL_EXCEL_PATH}\n 🧬 {FINAL_FASTA_PATH}")
    except Exception as e:
        print(f"Could not convert final data to Excel: {e}")

if __name__ == "__main__":
    main()