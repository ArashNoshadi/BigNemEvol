import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ================= CONFIGURATION =================

# Input file paths (Update these paths as needed)
INPUT_ANNOTATION = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\whole genome\Caenorhabditis elegans\genomic.gff"
FASTA_FILE = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\whole genome\Caenorhabditis elegans\GCA_000002985.3_WBcel235_genomic.fna"

# Output file paths
OUTPUT_PARTS = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\whole genome\Caenorhabditis elegans\Caenorhabditis elegans_parts.fasta"
OUTPUT_FULL = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\whole genome\Caenorhabditis elegans\Caenorhabditis elegans_clusters.fasta"
# --- NEW: Output path for Padded sequences ---
OUTPUT_PADDED = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\whole genome\Caenorhabditis elegans\Caenorhabditis elegans_padded_structure.fasta"


# =================================================

# --- MANUAL GENE LENGTHS FOR PADDING ---
# مقادیر زیر را بر اساس طول دقیق یا میانگین ژن‌های گونه خود تنظیم کنید
# این اعداد تعیین می‌کنند که چند تا N برای هم‌تراز کردن ژن‌ها اضافه شود
MANUAL_GENE_LENGTHS = {
    '18S': 1600,  # مثال: طول تقریبی 18S
    'ITS1': 500,  # مثال
    '5.8S': 130,  # مثال
    'ITS2': 300,  # مثال
    '28S': 3300   # مثال
}

# Maximum gap allowed between rRNA genes to consider them part of the same cluster
MAX_CLUSTER_GAP = 7000 

# Valid feature types to parse from GFF/GTF
VALID_TYPES = {
    'gene', 'exon', 'rRNA', 'rRNA_gene', 'transcript', 'ncRNA', 'misc_RNA', 
    'sequence_feature', 'mRNA', 'CDS', 'pseudogene', 'region'
}
# =================================================

def detect_format(file_path):
    """
    Detects if the annotation file is in GFF or GTF format based on the attribute column style.
    """
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            attributes = parts[8]
            if "=" in attributes and ";" in attributes:
                return "GFF"
            else:
                return "GTF"
    return "GTF" 

def detect_rrna_from_attributes(attr_string):
    """
    GREEDY SEARCH: Scans the ENTIRE attribute string for keywords.
    """
    attr_lower = attr_string.lower()
    
    # --- 18S (Small Subunit) ---
    if re.search(r'\b18s\b', attr_lower) or \
       'small subunit' in attr_lower or \
       'ssu' in attr_lower or \
       'rrn-1.1' in attr_lower or \
       'rn18s' in attr_lower:
        return '18S'
    
    # --- 5.8S ---
    if '5.8s' in attr_lower or \
       'rrn-2.1' in attr_lower or \
       'rn5-8s' in attr_lower or \
       'rn5.8s' in attr_lower:
        return '5.8S'
    
    # --- 28S (Large Subunit) ---
    if re.search(r'\b28s\b', attr_lower) or \
       re.search(r'\b26s\b', attr_lower) or \
       re.search(r'\b25s\b', attr_lower) or \
       'large subunit' in attr_lower or \
       'lsu' in attr_lower or \
       'rrn-3.1' in attr_lower or \
       'rn28s' in attr_lower:
        return '28S'
        
    return None

def resolve_chrom_name(genome_dict, chrom_name):
    """
    Attempts to match the chromosome name from annotation to the FASTA file keys.
    """
    if chrom_name in genome_dict: return chrom_name
    
    if chrom_name.startswith("chr"):
        alt_name = chrom_name[3:]
        if alt_name in genome_dict: return alt_name
    else:
        alt_name = f"chr{chrom_name}"
        if alt_name in genome_dict: return alt_name
    
    for k in genome_dict.keys():
        if chrom_name in k or k in chrom_name:
            return k
            
    return None

def safe_get_seq(genome_dict, chrom_resolved, start, end, strand):
    """
    Safely retrieves sequence slicing. Handles strand orientation (Reverse Complement).
    """
    if start >= end: return Seq("") 
    if chrom_resolved not in genome_dict: return Seq("")
    
    raw_seq = genome_dict[chrom_resolved][start:end].seq
    if strand == '-':
        return raw_seq.reverse_complement()
    return raw_seq

def debug_gff_content(file_path, file_fmt):
    """
    DIAGNOSTIC FUNCTION
    """
    print("\n--- DIAGNOSTIC MODE: Analyzing GFF file content ---")
    print(f"File Format Detected: {file_fmt}")
    
    keywords = ["18s", "28s", "5.8s", "rrna", "ribosom", "ssu", "lsu", "rrn-"]
    scan_hits = 0
    with open(file_path, 'r') as f:
        for line in f:
            line_lower = line.lower()
            if any(k in line_lower for k in keywords):
                if scan_hits < 5: 
                    print(f"Found keyword match: {line.strip()[:100]}...")
                scan_hits += 1
    
    if scan_hits == 0:
        print("RESULT: No ribosomal keywords found.")
    else:
        print(f"\nTotal lines containing keywords: {scan_hits}")
    print("---------------------------------------------------\n")

def parse_annotation(file_path):
    """
    Parses the annotation file using Greedy Attribute Search.
    """
    features = []
    file_fmt = detect_format(file_path)
    print(f"Reading {file_fmt} file: {file_path}...")

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            if parts[2] not in VALID_TYPES: continue
            
            attributes = parts[8]
            found_common_name = detect_rrna_from_attributes(attributes)
            
            if found_common_name:
                features.append({
                    'name': found_common_name,
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6]
                })
                
    if len(features) == 0:
        debug_gff_content(file_path, file_fmt)
        
    unique_features = {(f['chrom'], f['start'], f['end'], f['name']): f for f in features}
    return sorted(list(unique_features.values()), key=lambda x: (x['chrom'], x['strand'], x['start']))

def cluster_features(features):
    """
    Groups features. 
    Accepts partial clusters (e.g., just 18S, or 18S+28S).
    """
    clusters = []
    current_cluster = {'chrom': None, 'strand': None, 'features': [], 'last_end': -1}
    
    for feat in features:
        start_new = True
        
        if current_cluster['features']:
            if (feat['chrom'] == current_cluster['chrom'] and 
                feat['strand'] == current_cluster['strand']):
                
                dist = feat['start'] - current_cluster['last_end']
                if dist < MAX_CLUSTER_GAP:
                    start_new = False
        
        if start_new:
            if current_cluster['features']: clusters.append(current_cluster)
            current_cluster = {
                'chrom': feat['chrom'], 'strand': feat['strand'],
                'features': [feat], 'last_end': feat['end']
            }
        else:
            current_cluster['features'].append(feat)
            current_cluster['last_end'] = max(current_cluster['last_end'], feat['end'])
            
    if current_cluster['features']: clusters.append(current_cluster)
    
    valid_clusters = []
    for c in clusters:
        names = {f['name'] for f in c['features']}
        if len(names) > 0:
            valid_clusters.append(c)
            
    print(f"Total clusters identified (Complete or Partial): {len(valid_clusters)}")
    return valid_clusters

def process_and_save(clusters, fasta_path):
    """
    Extracts sequences. 
    Checks which genes are present.
    ADDS LENGTH TO FASTA HEADER.
    CREATES PADDED SEQUENCES BASED ON **MANUAL** GENE LENGTHS.
    """
    print(f"Indexing Genome: {fasta_path}...")
    try:
        genome_dict = SeqIO.index(fasta_path, "fasta")
    except Exception as e:
        print(f"Error indexing FASTA: {e}")
        return

    parts_records = []
    full_records = []
    padded_records = [] # List for the new requested functionality

    for idx, cluster in enumerate(clusters):
        cluster_id = idx + 1
        chrom_raw = cluster['chrom']
        strand = cluster['strand']
        
        chrom = resolve_chrom_name(genome_dict, chrom_raw)
        if not chrom:
            print(f"Warning: Chromosome {chrom_raw} not found in FASTA file.")
            continue
            
        gene_map = {f['name']: f for f in cluster['features']}
        
        has_18s = '18S' in gene_map
        has_58s = '5.8S' in gene_map
        has_28s = '28S' in gene_map
        
        # --- Temp Storage for Padded Logic ---
        current_seqs = {
            '18S': "",
            'ITS1': "",
            '5.8S': "",
            'ITS2': "",
            '28S': ""
        }

        # --- Helper to create record with LENGTH in header ---
        def make_rec(seq, name, sub_id):
            seq_len = len(seq)
            # Format: >ID Description (including len=...)
            return SeqRecord(
                seq, 
                id=f"{name}_Copy{cluster_id}", 
                description=f"len={seq_len} Cluster{cluster_id} {chrom}:{strand}"
            )

        # --- 1. EXTRACT INDIVIDUAL GENES (PARTS) ---
        if has_18s:
            f18s = gene_map['18S']
            s18s = safe_get_seq(genome_dict, chrom, f18s['start']-1, f18s['end'], strand)
            parts_records.append(make_rec(s18s, "18S", cluster_id))
            current_seqs['18S'] = str(s18s) 
            
        if has_58s:
            f58s = gene_map['5.8S']
            s58s = safe_get_seq(genome_dict, chrom, f58s['start']-1, f58s['end'], strand)
            parts_records.append(make_rec(s58s, "5.8S", cluster_id))
            current_seqs['5.8S'] = str(s58s) 
            
        if has_28s:
            f28s = gene_map['28S']
            s28s = safe_get_seq(genome_dict, chrom, f28s['start']-1, f28s['end'], strand)
            parts_records.append(make_rec(s28s, "28S", cluster_id))
            current_seqs['28S'] = str(s28s) 
            
        # --- 2. EXTRACT ITS REGIONS ---
        all_feats = sorted(cluster['features'], key=lambda x: x['start'])
        
        if has_18s and has_58s:
            f_a = gene_map['18S']
            f_b = gene_map['5.8S']
            pair = sorted([f_a, f_b], key=lambda x: x['start'])
            gap_seq = safe_get_seq(genome_dict, chrom, pair[0]['end'], pair[1]['start']-1, strand)
            parts_records.append(make_rec(gap_seq, "ITS1", cluster_id))
            current_seqs['ITS1'] = str(gap_seq) 
            
        if has_58s and has_28s:
            f_a = gene_map['5.8S']
            f_b = gene_map['28S']
            pair = sorted([f_a, f_b], key=lambda x: x['start'])
            gap_seq = safe_get_seq(genome_dict, chrom, pair[0]['end'], pair[1]['start']-1, strand)
            parts_records.append(make_rec(gap_seq, "ITS2", cluster_id))
            current_seqs['ITS2'] = str(gap_seq) 
            
        # --- 3. EXTRACT FULL CLUSTER ---
        if all_feats:
            full_start = all_feats[0]['start'] - 1 
            full_end = all_feats[-1]['end']
            
            full_seq = safe_get_seq(genome_dict, chrom, full_start, full_end, strand)
            
            content_label = []
            if has_18s: content_label.append("18S")
            if has_58s: content_label.append("5.8S")
            if has_28s: content_label.append("28S")
            desc_str = "+".join(content_label)
            
            full_records.append(SeqRecord(
                full_seq, 
                id=f"Full_Cluster_Copy{cluster_id}", 
                description=f"len={len(full_seq)} {chrom}:{strand} Contains:{desc_str}"
            ))

        # --- 4. GENERATE PADDED SEQUENCES (USING MANUAL LENGTHS) ---
        
        # Retrieve Manual Lengths (default to 0 if key missing)
        l_18s = MANUAL_GENE_LENGTHS.get('18S', 0)
        l_its1 = MANUAL_GENE_LENGTHS.get('ITS1', 0)
        l_58s = MANUAL_GENE_LENGTHS.get('5.8S', 0)
        l_its2 = MANUAL_GENE_LENGTHS.get('ITS2', 0)
        # l_28s = MANUAL_GENE_LENGTHS.get('28S', 0)
        
        # 18S (No padding needed, it is first)
        if current_seqs['18S']:
            padded_records.append(make_rec(Seq(current_seqs['18S']), "18S_Padded", cluster_id))
            
        # ITS1 (Padding = Manual len 18S)
        if current_seqs['ITS1']:
            pad = "N" * l_18s
            new_seq = Seq(pad + current_seqs['ITS1'])
            padded_records.append(make_rec(new_seq, "ITS1_Padded", cluster_id))
            
        # 5.8S (Padding = Manual len 18S + Manual len ITS1)
        if current_seqs['5.8S']:
            pad = "N" * (l_18s + l_its1)
            new_seq = Seq(pad + current_seqs['5.8S'])
            padded_records.append(make_rec(new_seq, "5.8S_Padded", cluster_id))
            
        # ITS2 (Padding = Manual len 18S + Manual len ITS1 + Manual len 5.8S)
        if current_seqs['ITS2']:
            pad = "N" * (l_18s + l_its1 + l_58s)
            new_seq = Seq(pad + current_seqs['ITS2'])
            padded_records.append(make_rec(new_seq, "ITS2_Padded", cluster_id))
            
        # 28S (Padding = Manual len 18S + Manual len ITS1 + Manual len 5.8S + Manual len ITS2)
        if current_seqs['28S']:
            pad = "N" * (l_18s + l_its1 + l_58s + l_its2)
            new_seq = Seq(pad + current_seqs['28S'])
            padded_records.append(make_rec(new_seq, "28S_Padded", cluster_id))


    # --- SAVE FILES ---
    if parts_records:
        with open(OUTPUT_PARTS, "w") as f: SeqIO.write(parts_records, f, "fasta")
        with open(OUTPUT_FULL, "w") as f: SeqIO.write(full_records, f, "fasta")
        # Save the new padded file
        with open(OUTPUT_PADDED, "w") as f: SeqIO.write(padded_records, f, "fasta")
        
        print(f"\nSuccess! Files saved:\n{OUTPUT_PARTS}\n{OUTPUT_FULL}\n{OUTPUT_PADDED}")
        
        counts = {'18S':0, '5.8S':0, '28S':0, 'ITS1':0, 'ITS2':0}
        for r in parts_records:
            for k in counts:
                if r.id.startswith(k): counts[k] += 1
        print("\nSummary of extracted parts:")
        print(f"  18S sequences: {counts['18S']}")
        print(f"  5.8S sequences: {counts['5.8S']}")
        print(f"  28S sequences: {counts['28S']}")
        print(f"  ITS1 regions:  {counts['ITS1']}")
        print(f"  ITS2 regions:  {counts['ITS2']}")
    else:
        print("No rRNA features found to save.")

# ================= MAIN =================
if __name__ == "__main__":
    if os.path.exists(INPUT_ANNOTATION) and os.path.exists(FASTA_FILE):
        features = parse_annotation(INPUT_ANNOTATION)
        print(f"Found {len(features)} rRNA related features.")
        clusters = cluster_features(features)
        process_and_save(clusters, FASTA_FILE)
    else:
        print("Error: Input files not found. Please check the paths in CONFIGURATION.")