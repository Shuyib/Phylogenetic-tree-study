#!/usr/bin/env python3
"""
Fix errors in NCBI BLAST CSV data by comparing with original file
"""

import csv
import re
import pandas as pd
from pathlib import Path

def fix_ncbi_blast_csv(original_file, error_file, output_file):
    """Fix errors in NCBI BLAST CSV file by referencing original data"""
    
    # Read the original file
    original_df = pd.read_csv(original_file)
    
    # Read the error file
    error_df = pd.read_csv(error_file)
    
    # Create a copy for fixing
    fixed_df = error_df.copy()
    
    print(f"Original file has {len(original_df)} rows")
    print(f"Error file has {len(error_df)} rows")
    print("\nIdentified errors and fixes:")
    
    # Fix header format - restore proper column names
    if 'Score_Bits_' in fixed_df.columns:
        fixed_df.rename(columns={'Score_Bits_': 'Score(Bits)'}, inplace=True)
        print("✓ Fixed header: Score_Bits_ → Score(Bits)")
    
    if 'Identities_Pct_' in fixed_df.columns:
        fixed_df.rename(columns={'Identities_Pct_': 'Identities(%)'}, inplace=True)
        print("✓ Fixed header: Identities_Pct_ → Identities(%)")
    
    if 'Positives_Pct_' in fixed_df.columns:
        fixed_df.rename(columns={'Positives_Pct_': 'Positives(%)'}, inplace=True)
        print("✓ Fixed header: Positives_Pct_ → Positives(%)")
    
    if 'E__' in fixed_df.columns:
        fixed_df.rename(columns={'E__': 'E()'}, inplace=True)
        print("✓ Fixed header: E__ → E()")
    
    # Fix specific data errors by row
    for idx, row in fixed_df.iterrows():
        original_row = original_df.iloc[idx] if idx < len(original_df) else None
        
        # Fix duplicate accession A0A1X7MEN8 (rows 3, 4, 8)
        if idx == 4 and row['Accession'] == 'A0A1X7MEN8':
            fixed_df.at[idx, 'Accession'] = 'A0A7N0T861'
            print(f"✓ Row {idx+2}: Fixed duplicate accession A0A1X7MEN8 → A0A7N0T861")
        
        if idx == 8 and row['Accession'] == 'A0A1X7MEN8':
            fixed_df.at[idx, 'Accession'] = 'A0A6A5L520'
            print(f"✓ Row {idx+2}: Fixed duplicate accession A0A1X7MEN8 → A0A6A5L520")
        
        # Fix organism names with invalid characters (@#$%)
        if pd.notna(row['Organism']) and '@#$%' in str(row['Organism']):
            if original_row is not None:
                fixed_df.at[idx, 'Organism'] = original_row['Organism']
                print(f"✓ Row {idx+2}: Fixed organism name with invalid characters")
        
        # Fix missing Length values
        if pd.isna(row['Length']) or row['Length'] == '':
            if original_row is not None:
                fixed_df.at[idx, 'Length'] = original_row['Length']
                print(f"✓ Row {idx+2}: Fixed missing Length value")
        
        # Fix invalid Score values
        if pd.notna(row['Score(Bits)']) and str(row['Score(Bits)']).lower() in ['error', 'invalid', 'n/a']:
            if original_row is not None:
                fixed_df.at[idx, 'Score(Bits)'] = original_row['Score(Bits)']
                print(f"✓ Row {idx+2}: Fixed invalid Score value")
        
        # Fix invalid E() values
        if pd.notna(row['E()']) and str(row['E()']).lower() in ['error', 'invalid', 'n/a', 'null', '---']:
            if original_row is not None:
                fixed_df.at[idx, 'E()'] = original_row['E()']
                print(f"✓ Row {idx+2}: Fixed invalid E() value")
        
        # Fix extremely long descriptions (truncate to reasonable length)
        if pd.notna(row['Description']) and len(str(row['Description'])) > 500:
            if original_row is not None:
                fixed_df.at[idx, 'Description'] = original_row['Description']
                print(f"✓ Row {idx+2}: Fixed overly long description")
    
    # Add missing first column name if needed
    if fixed_df.columns[0] != '':
        fixed_df.insert(0, '', range(len(fixed_df)))
    
    # Save the fixed file
    fixed_df.to_csv(output_file, index=False)
    print(f"\n✓ Fixed file saved as: {output_file}")
    
    # Validation summary
    print(f"\nValidation Summary:")
    print(f"- Total rows processed: {len(fixed_df)}")
    print(f"- Missing values in Length: {fixed_df['Length'].isna().sum()}")
    print(f"- Missing values in Score(Bits): {fixed_df['Score(Bits)'].isna().sum()}")
    print(f"- Missing values in E(): {fixed_df['E()'].isna().sum()}")
    
    return fixed_df

def main():
    # File paths
    original_file = "/Users/maibee/Downloads/open-data-editor/ncbiblast-R20230110-091221-0382-99758486-p1m_meme6.csv"
    error_file = "/Users/maibee/Downloads/open-data-editor/ncbiblast-R20230110-091221-0382-99758486-p1m_meme6_with_errors.csv"
    output_file = "/Users/maibee/Downloads/open-data-editor/ncbiblast-R20230110-091221-0382-99758486-p1m_meme6_errors_fixed.csv"
    
    # Check if files exist
    if not Path(original_file).exists():
        print(f"Error: Original file not found: {original_file}")
        return
    
    if not Path(error_file).exists():
        print(f"Error: Error file not found: {error_file}")
        return
    
    print("NCBI BLAST CSV Error Fixer")
    print("=" * 40)
    
    # Fix the errors
    fixed_df = fix_ncbi_blast_csv(original_file, error_file, output_file)
    
    print("\n" + "=" * 40)
    print("Error fixing completed successfully!")

if __name__ == "__main__":
    main()