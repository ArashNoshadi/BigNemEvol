import os
import time
import re
import html
import numpy as np
import pandas as pd
from Bio import Entrez

# ================= تنظیمات کاربر =================
# ایمیل خود را برای رعایت قوانین NCBI وارد کنید
Entrez.email = "your_email@example.com" 

# مسیر ذخیره‌سازی فایل خروجی
output_dir = r"G:\Paper\nema-Nanopore-Sequencing\Nature Genetics"
output_filename = "Nematoda_Genome_Unique_Stats.xlsx"

# ایجاد پوشه در صورت عدم وجود
os.makedirs(output_dir, exist_ok=True)
full_output_path = os.path.join(output_dir, output_filename)

def run_analysis():
    print("--- 1. Connecting to NCBI Assembly Database ---")
    
    search_term = "Nematoda[Organism]"
    
    try:
        # جستجو برای یافتن تمام IDها
        search_handle = Entrez.esearch(db="assembly", term=search_term, retmax=10000)
        search_results = Entrez.read(search_handle)
        search_handle.close()
    except Exception as e:
        print(f"Error connecting to NCBI: {e}")
        return

    id_list = search_results["IdList"]
    total_found = len(id_list)
    print(f"Found {total_found} assemblies. Fetching metadata in batches...")

    if total_found == 0:
        print("No assemblies found. Please check your internet connection.")
        return

    # لیست برای ذخیره داده‌های خام
    raw_data = []
    
    # پردازش دسته‌ای (Batch Processing) برای جلوگیری از قطع شدن سرور
    batch_size = 200
    
    for start in range(0, total_found, batch_size):
        end = min(total_found, start + batch_size)
        batch_ids = id_list[start:end]
        print(f"Processing batch {start} to {end}...")
        
        try:
            # دریافت اطلاعات خلاصه (Summary)
            summary_handle = Entrez.esummary(db="assembly", id=",".join(batch_ids))
            summaries = Entrez.read(summary_handle)
            summary_handle.close()
            
            # مدیریت تفاوت ساختار خروجی در نسخه‌های مختلف Biopython
            if 'DocumentSummarySet' in summaries:
                records = summaries['DocumentSummarySet']['DocumentSummary']
            else:
                records = summaries

            for rec in records:
                try:
                    # --- 1. استخراج نام و اطلاعات پایه ---
                    org_name = rec.get("Organism", "Unknown")
                    tax_id = rec.get("Taxid", "")
                    
                    # --- 2. تمیزکاری متادیتا (Decode HTML Entities) ---
                    # داده‌های NCBI اغلب به صورت &lt;Stat... ذخیره می‌شوند
                    raw_meta_xml = rec.get("Meta", "")
                    clean_meta = html.unescape(raw_meta_xml)
                    
                    # --- 3. استخراج اعداد با Regex (غیر حساس به بزرگی/کوچکی حروف) ---
                    
                    # الف) استخراج طول ژنوم (Genome Size)
                    size_match = re.search(r'category="total_length".*?>(\d+)<', clean_meta, re.IGNORECASE)
                    if not size_match: # الگوی جایگزین
                        size_match = re.search(r'<Stat category="total_length".*?>(\d+)</Stat>', clean_meta, re.IGNORECASE)
                    
                    genome_size = int(size_match.group(1)) if size_match else 0

                    # ب) استخراج درصد GC
                    gc_match = re.search(r'category="gc_percent".*?>([\d\.]+)<', clean_meta, re.IGNORECASE)
                    gc_pct = float(gc_match.group(1)) if gc_match else np.nan

                    # ج) استخراج تعداد ژن‌های کدکننده (Coding Genes)
                    # نام این فیلد متغیر است (total_coding_genes یا gene_count)
                    gene_match = re.search(r'category="total_coding_genes".*?>(\d+)<', clean_meta, re.IGNORECASE)
                    if not gene_match:
                        gene_match = re.search(r'category="gene_count".*?>(\d+)<', clean_meta, re.IGNORECASE)
                    
                    coding_genes = int(gene_match.group(1)) if gene_match else 0

                    # فقط رکوردهایی که حداقل سایز ژنوم را دارند ذخیره کن
                    if genome_size > 0:
                        raw_data.append({
                            "Organism": org_name,
                            "TaxID": tax_id,
                            "Genome_Size_bp": genome_size,
                            "GC_Percent": gc_pct,
                            "Coding_Genes": coding_genes
                        })
                        
                except Exception as inner_e:
                    continue # رد کردن رکورد خراب

        except Exception as e:
            print(f"Batch Error: {e}")
            time.sleep(2)

    # =========================================================
    # پردازش نهایی، میانگین‌گیری و ذخیره
    # =========================================================
    if raw_data:
        print("\n--- 2. Aggregating Data (Removing Duplicates) ---")
        df = pd.DataFrame(raw_data)
        
        # تبدیل 0 به NaN در ستون ژن‌ها
        # (تا ژنوم‌های ناقص که تعداد ژن ندارند، میانگین ژنوم‌های کامل را خراب نکنند)
        df['Coding_Genes'] = df['Coding_Genes'].replace(0, np.nan)

        # *** میانگین‌گیری بر اساس نام گونه (Organism) ***
        df_unique = df.groupby('Organism', as_index=False).agg({
            'Genome_Size_bp': 'mean',    # میانگین سایز
            'GC_Percent': 'mean',        # میانگین GC
            'Coding_Genes': 'mean',      # میانگین تعداد ژن (NaNها نادیده گرفته می‌شوند)
            'TaxID': 'first'             # نگه داشتن شناسه تاکسونومی
        })
        
        # اضافه کردن ستون تعداد اسمبلی‌های یافت شده برای هر گونه
        count_series = df.groupby('Organism').size().reset_index(name='Assembly_Count')
        df_unique = df_unique.merge(count_series, on='Organism')

        # --- محاسبات فرمولی (روی داده‌های میانگین‌گیری شده) ---
        print("--- 3. Calculating Derived Metrics (Log & Density) ---")
        
        # 1. Log Genome Size
        df_unique['Log_Genome_Size'] = np.log10(df_unique['Genome_Size_bp'])
        
        # 2. Log Coding Genes (با شرط اینکه مقدار وجود داشته باشد)
        df_unique['Log_Coding_Genes'] = df_unique['Coding_Genes'].apply(
            lambda x: np.log10(x) if pd.notnull(x) and x > 0 else 0
        )
        
        # 3. Gene Density (تعداد ژن تقسیم بر مگابایت)
        df_unique['Gene_Density_per_Mb'] = df_unique.apply(
            lambda row: row['Coding_Genes'] / (row['Genome_Size_bp'] / 1_000_000) 
            if pd.notnull(row['Coding_Genes']) and row['Genome_Size_bp'] > 0 else 0,
            axis=1
        )

        # گرد کردن اعداد برای زیبایی فایل اکسل
        df_unique = df_unique.round({
            'Genome_Size_bp': 0,
            'GC_Percent': 2,
            'Coding_Genes': 0,
            'Log_Genome_Size': 4,
            'Log_Coding_Genes': 4,
            'Gene_Density_per_Mb': 2
        })

        # ذخیره در اکسل
        df_unique.to_excel(full_output_path, index=False)
        
        print(f"\n✅ SUCCESS! Process Completed.")
        print(f"   - Total Raw Assemblies Processed: {len(df)}")
        print(f"   - Unique Species (Rows in Excel): {len(df_unique)}")
        print(f"   - File Saved to: {full_output_path}")
        
    else:
        print("\n❌ No valid data parsed from NCBI.")

if __name__ == "__main__":
    run_analysis()