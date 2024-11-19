# sort_mags.py

import os
import csv
from pathlib import Path

def process_sample(sample, checkm_path, metabat2_bins_path, sorted_mags_path, summary_writer):
    quality_file = checkm_path / f"{sample}/{sample}_quality.txt"
    
    # Check if the quality file exists before proceeding
    if not quality_file.exists():
        print(f"Quality file not found for sample {sample}, skipping...")
        return
    
    bins_dir = next(metabat2_bins_path.glob(f"{sample}_*.metabat-bins-*"), None)

    if not bins_dir:
        print(f"Bins directory not found for sample {sample}, skipping...")
        return

    with quality_file.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            bin_id = row['Bin Id']
            new_bin_id = f"{sample}__{bin_id.replace('bin.', 'bin_')}"
            completeness = float(row['Completeness'])
            contamination = float(row['Contamination'])

            if completeness >= 90.0 and contamination <= 10.0:
                quality = "high_quality"
            elif 50.0 <= completeness < 90.0 and contamination <= 10.0:
                quality = "medium_quality"
            else:
                quality = "low_quality"

            fa_file_path = bins_dir / f"{bin_id}.fa"
            if fa_file_path.exists():
                quality_dir = sorted_mags_path / quality
                quality_dir.mkdir(parents=True, exist_ok=True)
                symlink_filename = f"{new_bin_id}.fa"
                symlink_path = quality_dir / symlink_filename
                if not symlink_path.exists():
                    symlink_path.symlink_to(fa_file_path.resolve())

                summary_writer.writerow([bin_id, sample, new_bin_id, completeness, contamination, quality])

def main(checkm_path_str, metabat2_bins_path_str, sorted_mags_path_str):
    checkm_path = Path(checkm_path_str)
    metabat2_bins_path = Path(metabat2_bins_path_str)
    sorted_mags_path = Path(sorted_mags_path_str)
    
    # Create sorted_mags directory and its subdirectories
    sorted_mags_path.mkdir(parents=True, exist_ok=True)

    summary_file_path = sorted_mags_path / "checkm_summary.txt"
    with summary_file_path.open('w', newline='') as summary_file:
        writer = csv.writer(summary_file, delimiter='\t')
        writer.writerow(["BinID", "Sample", "NewBinID", "Completeness", "Contamination", "Quality"])

        for sample_dir in checkm_path.iterdir():
            if sample_dir.is_dir():
                sample = sample_dir.name
                process_sample(sample, checkm_path, metabat2_bins_path, sorted_mags_path, writer)

if __name__ == "__main__":
    checkm_path = "checkm"
    metabat2_bins_path = "metabat2_bins"
    sorted_mags_path = "sorted_mags"
    main(checkm_path, metabat2_bins_path, sorted_mags_path)
