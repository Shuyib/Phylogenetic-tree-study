#!/usr/bin/env python3
"""
Fix errors in pathogenic E. coli CSV data by comparing with original file
"""

import csv
import re
import pandas as pd
from pathlib import Path

def fix_pathogenic_ecoli_csv(original_file, error_file, output_file):
    """Fix errors in pathogenic E. coli CSV file by referencing original data"""
    
    # Read the original file
    original_df = pd.read_csv(original_file)
    
    # Read the error file
    error_df = pd.read_csv(error_file)
    
    # Create a copy for fixing
    fixed_df = error_df.copy()
    
    print(f"Original file has {len(original_df)} rows")
    print(f"Error file has {len(error_df)} rows")
    print("\nIdentified errors and fixes:")
    
    # Fix specific data errors by row
    for idx, row in fixed_df.iterrows():
        original_row = original_df.iloc[idx] if idx < len(original_df) else None
        
        if original_row is None:
            continue
            
        # Fix missing pathogenesis_mechanisms (row 6 in original had this issue)
        if pd.isna(row['pathogenesis_mechanisms']) or row['pathogenesis_mechanisms'] == '':
            if pd.notna(original_row['pathogenesis_mechanisms']):
                fixed_df.at[idx, 'pathogenesis_mechanisms'] = original_row['pathogenesis_mechanisms']
                print(f"✓ Row {idx+2}: Fixed missing pathogenesis_mechanisms")
        
        # Fix missing reference_or_ncbi_id (row 24 in original had this issue)
        if pd.isna(row['reference_or_ncbi_id']) or row['reference_or_ncbi_id'] == '':
            if pd.notna(original_row['reference_or_ncbi_id']):
                fixed_df.at[idx, 'reference_or_ncbi_id'] = original_row['reference_or_ncbi_id']
                print(f"✓ Row {idx+2}: Fixed missing reference_or_ncbi_id")
        
        # Fix missing strain_name (row 31 in original had this issue)
        if pd.isna(row['strain_name']) or row['strain_name'] == '':
            if pd.notna(original_row['strain_name']):
                fixed_df.at[idx, 'strain_name'] = original_row['strain_name']
                print(f"✓ Row {idx+2}: Fixed missing strain_name")
        
        # Fix invalid pathogenic_type
        if pd.notna(row['pathogenic_type']) and str(row['pathogenic_type']).upper() == 'INVALID_TYPE':
            if pd.notna(original_row['pathogenic_type']):
                fixed_df.at[idx, 'pathogenic_type'] = original_row['pathogenic_type']
                print(f"✓ Row {idx+2}: Fixed invalid pathogenic_type")
        
        # Fix invalid reference_or_ncbi_id
        if pd.notna(row['reference_or_ncbi_id']) and str(row['reference_or_ncbi_id']).upper() == 'INVALID_ID':
            if pd.notna(original_row['reference_or_ncbi_id']):
                fixed_df.at[idx, 'reference_or_ncbi_id'] = original_row['reference_or_ncbi_id']
                print(f"✓ Row {idx+2}: Fixed invalid reference_or_ncbi_id")
        
        # Fix inconsistent quoting - remove unnecessary quotes
        for col in ['strain_name', 'pathogenic_type', 'infection_type', 'pathogenesis_mechanisms']:
            if pd.notna(row[col]) and isinstance(row[col], str):
                # Remove quotes if they're inconsistent with original
                if row[col].startswith('"') and row[col].endswith('"'):
                    if not (original_row[col].startswith('"') and original_row[col].endswith('"')):
                        fixed_df.at[idx, col] = row[col].strip('"')
                        print(f"✓ Row {idx+2}: Fixed inconsistent quoting in {col}")
    
    # Remove any extra rows that don't exist in original
    if len(fixed_df) > len(original_df):
        rows_to_remove = len(fixed_df) - len(original_df)
        fixed_df = fixed_df.iloc[:-rows_to_remove]
        print(f"✓ Removed {rows_to_remove} extra rows")
    
    # Save the fixed file
    fixed_df.to_csv(output_file, index=False)
    print(f"\n✓ Fixed file saved as: {output_file}")
    
    # Validation summary
    print(f"\nValidation Summary:")
    print(f"- Total rows processed: {len(fixed_df)}")
    print(f"- Missing values in strain_name: {fixed_df['strain_name'].isna().sum()}")
    print(f"- Missing values in pathogenic_type: {fixed_df['pathogenic_type'].isna().sum()}")
    print(f"- Missing values in pathogenesis_mechanisms: {fixed_df['pathogenesis_mechanisms'].isna().sum()}")
    print(f"- Missing values in reference_or_ncbi_id: {fixed_df['reference_or_ncbi_id'].isna().sum()}")
    
    return fixed_df

def main():
    # File paths
    original_file = "/Users/maibee/Downloads/open-data-editor/pathogenic_e_coli_30strains.csv"
    error_file = "/Users/maibee/Downloads/open-data-editor/pathogenic_e_coli_30strains_with_errors.csv"
    output_file = "/Users/maibee/Downloads/open-data-editor/pathogenic_e_coli_30strains_errors_fixed.csv"
    
    # Check if files exist
    if not Path(original_file).exists():
        print(f"Error: Original file not found: {original_file}")
        return
    
    if not Path(error_file).exists():
        print(f"Error: Error file not found: {error_file}")
        return
    
    print("Pathogenic E. coli CSV Error Fixer")
    print("=" * 40)
    
    # Fix the errors
    fixed_df = fix_pathogenic_ecoli_csv(original_file, error_file, output_file)
    
    print("\n" + "=" * 40)
    print("Error fixing completed successfully!")

if __name__ == "__main__":
    main()