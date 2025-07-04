#!/usr/bin/env python3
"""
Enhanced NCBI BLAST CSV Data Cleaning and Issue Detection Script

This script identifies and fixes common issues in NCBI BLAST CSV files:
- Missing values
- Inconsistent data types
- Malformed columns
- Duplicate entries
- Invalid E-values and scores
- Creates a file with intentional errors for testing
"""

import pandas as pd
import numpy as np
import re
import sys
from pathlib import Path


class NCBIBlastCleaner:
    def __init__(self, csv_path):
        self.csv_path = Path(csv_path)
        self.df = None
        self.issues = []
        
    def load_data(self):
        """Load CSV data with error handling"""
        try:
            self.df = pd.read_csv(self.csv_path, index_col=0)
            print(f"Loaded {len(self.df)} rows from {self.csv_path}")
        except Exception as e:
            print(f"Error loading CSV: {e}")
            sys.exit(1)
    
    def detect_issues(self):
        """Detect various data quality issues"""
        print("\n=== DETECTING ISSUES ===")
        
        # Check for missing values
        missing_counts = self.df.isnull().sum()
        if missing_counts.sum() > 0:
            self.issues.append("Missing values detected")
            print(f"Missing values per column:\n{missing_counts[missing_counts > 0]}")
        
        # Check column names
        expected_cols = ['Hit', 'DB', 'Accession', 'Description', 'Organism', 'Length', 
                        'Score(Bits)', 'Identities(%)', 'Positives(%)', 'E()']
        if list(self.df.columns) != expected_cols:
            self.issues.append("Column names don't match expected format")
            print(f"Expected: {expected_cols}")
            print(f"Actual: {list(self.df.columns)}")
        
        # Check for duplicate accessions
        if 'Accession' in self.df.columns:
            duplicates = self.df['Accession'].duplicated().sum()
            if duplicates > 0:
                self.issues.append(f"Found {duplicates} duplicate accessions")
                print(f"Duplicate accessions: {duplicates}")
        
        # Check numeric columns
        numeric_cols = ['Length', 'Score(Bits)', 'Identities(%)', 'Positives(%)']
        for col in numeric_cols:
            if col in self.df.columns:
                non_numeric = pd.to_numeric(self.df[col], errors='coerce').isnull().sum()
                if non_numeric > 0:
                    self.issues.append(f"Non-numeric values in {col}: {non_numeric}")
                    print(f"Non-numeric values in {col}: {non_numeric}")
        
        # Check E-values format
        if 'E()' in self.df.columns:
            invalid_evals = 0
            for val in self.df['E()']:
                if pd.isna(val):
                    continue
                try:
                    float(val)
                except (ValueError, TypeError):
                    invalid_evals += 1
            if invalid_evals > 0:
                self.issues.append(f"Invalid E-values: {invalid_evals}")
                print(f"Invalid E-values: {invalid_evals}")
        
        # Check for extremely long descriptions
        if 'Description' in self.df.columns:
            long_desc = (self.df['Description'].str.len() > 500).sum()
            if long_desc > 0:
                self.issues.append(f"Extremely long descriptions: {long_desc}")
                print(f"Extremely long descriptions (>500 chars): {long_desc}")
        
        # Check for unusual characters in organism names
        if 'Organism' in self.df.columns:
            special_chars = self.df['Organism'].str.contains(r'[^\w\s\(\)\-\.]', na=False).sum()
            if special_chars > 0:
                self.issues.append(f"Organism names with unusual characters: {special_chars}")
                print(f"Organism names with unusual characters: {special_chars}")
        
        print(f"\nTotal issues found: {len(self.issues)}")
        return self.issues
    
    def fix_issues(self):
        """Fix detected issues"""
        print("\n=== FIXING ISSUES ===")
        
        if self.df is None:
            print("No data loaded")
            return
        
        # Fix column names
        expected_cols = ['Hit', 'DB', 'Accession', 'Description', 'Organism', 'Length', 
                        'Score(Bits)', 'Identities(%)', 'Positives(%)', 'E()']
        if list(self.df.columns) != expected_cols:
            print("Fixing column names...")
            self.df.columns = expected_cols
        
        # Remove duplicate accessions (keep first occurrence)
        if 'Accession' in self.df.columns:
            before_count = len(self.df)
            self.df = self.df.drop_duplicates(subset='Accession', keep='first')
            after_count = len(self.df)
            if before_count != after_count:
                print(f"Removed {before_count - after_count} duplicate accessions")
        
        # Fix numeric columns
        numeric_cols = ['Length', 'Score(Bits)', 'Identities(%)', 'Positives(%)']
        for col in numeric_cols:
            if col in self.df.columns:
                # Convert to numeric, replace invalid with NaN
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
                # Fill NaN with median for that column
                if self.df[col].isnull().sum() > 0:
                    median_val = self.df[col].median()
                    self.df[col].fillna(median_val, inplace=True)
                    print(f"Fixed non-numeric values in {col} using median: {median_val}")
        
        # Fix E-values
        if 'E()' in self.df.columns:
            def fix_evalue(val):
                if pd.isna(val):
                    return 1.0
                try:
                    return float(val)
                except (ValueError, TypeError):
                    return 1.0
            
            self.df['E()'] = self.df['E()'].apply(fix_evalue)
            print("Fixed invalid E-values")
        
        # Truncate extremely long descriptions
        if 'Description' in self.df.columns:
            long_mask = self.df['Description'].str.len() > 500
            if long_mask.sum() > 0:
                self.df.loc[long_mask, 'Description'] = self.df.loc[long_mask, 'Description'].str[:500] + '...'
                print(f"Truncated {long_mask.sum()} long descriptions")
        
        # Clean organism names
        if 'Organism' in self.df.columns:
            # Remove unusual characters, keep only alphanumeric, spaces, parentheses, hyphens, dots
            self.df['Organism'] = self.df['Organism'].str.replace(r'[^\w\s\(\)\-\.]', '', regex=True)
            print("Cleaned organism names")
        
        # Fill any remaining NaN values
        self.df.fillna('', inplace=True)
        
        print("Data cleaning completed")
    
    def save_cleaned_data(self, output_path=None):
        """Save cleaned data to new CSV file"""
        if output_path is None:
            output_path = self.csv_path.stem + '_cleaned.csv'
        
        self.df.to_csv(output_path, index=False)
        print(f"Cleaned data saved to: {output_path}")
        return output_path
    
    def create_error_file(self, error_output_path=None):
        """Create a version of the file with intentional errors for testing"""
        if error_output_path is None:
            error_output_path = self.csv_path.stem + '_with_errors.csv'
        
        # Create a copy of the original data
        error_df = self.df.copy()
        
        # Introduce various errors
        print(f"\n=== CREATING ERROR FILE: {error_output_path} ===")
        
        # Add missing values (10% of rows)
        sample_size = int(len(error_df) * 0.1)
        missing_indices = np.random.choice(error_df.index, sample_size, replace=False)
        error_df.loc[missing_indices, 'Length'] = np.nan
        print(f"Added {sample_size} missing values in Length column")
        
        # Add invalid E-values
        eval_indices = np.random.choice(error_df.index, 5, replace=False)
        error_df['E()'] = error_df['E()'].astype(str)
        error_df.loc[eval_indices, 'E()'] = ['invalid', 'N/A', 'error', '---', 'null']
        print("Added 5 invalid E-values")
        
        # Add duplicate accessions
        if len(error_df) > 10:
            dup_indices = np.random.choice(error_df.index[:10], 3, replace=False)
            error_df.loc[dup_indices, 'Accession'] = error_df.loc[dup_indices[0], 'Accession']
            print("Added 3 duplicate accessions")
        
        # Add non-numeric values to numeric columns
        score_indices = np.random.choice(error_df.index, 3, replace=False)
        error_df['Score(Bits)'] = error_df['Score(Bits)'].astype(str)
        error_df.loc[score_indices, 'Score(Bits)'] = ['N/A', 'error', 'invalid']
        print("Added 3 non-numeric values in Score(Bits) column")
        
        # Add extremely long descriptions
        long_desc_indices = np.random.choice(error_df.index, 2, replace=False)
        long_desc = "This is an extremely long description that exceeds normal limits and should be truncated. " * 10
        error_df.loc[long_desc_indices, 'Description'] = long_desc
        print("Added 2 extremely long descriptions")
        
        # Add unusual characters to organism names
        org_indices = np.random.choice(error_df.index, 3, replace=False)
        error_df.loc[org_indices, 'Organism'] = error_df.loc[org_indices, 'Organism'] + ' @#$%'
        print("Added unusual characters to 3 organism names")
        
        # Mess up column names - use same number of columns as original
        new_cols = []
        for col in error_df.columns:
            new_cols.append(col.replace('(', '_').replace(')', '_').replace('%', 'Pct'))
        error_df.columns = new_cols
        print("Modified column names to non-standard format")
        
        # Save error file
        error_df.to_csv(error_output_path, index=False)
        print(f"Error file created: {error_output_path}")
        return error_output_path


def main():
    if len(sys.argv) < 2:
        print("Usage: python tidy_ncbiblast_enhanced.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    
    # Create cleaner instance
    cleaner = NCBIBlastCleaner(csv_file)
    
    # Load data
    cleaner.load_data()
    
    # Create error file first (before cleaning)
    error_file = cleaner.create_error_file()
    
    # Detect issues
    issues = cleaner.detect_issues()
    
    # Fix issues
    cleaner.fix_issues()
    
    # Save cleaned data
    cleaned_file = cleaner.save_cleaned_data()
    
    print(f"\n=== SUMMARY ===")
    print(f"Original file: {csv_file}")
    print(f"Cleaned file: {cleaned_file}")
    print(f"Error file: {error_file}")
    print(f"Issues found and fixed: {len(issues)}")
    
    return cleaned_file, error_file


if __name__ == "__main__":
    main()