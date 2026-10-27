import os
import pandas as pd
from pathlib import Path

def compile_fragle_data(fragle_dir):
    # Initialize an empty list to store all data
    all_data = []
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(fragle_dir):
        # Skip if not a leaf directory (has subdirectories)
        if dirs:
            continue
            
        # Get the sample ID from the directory name
        sample_id = os.path.basename(root)
        
        # Initialize a dictionary for this sample
        sample_data = {'Sample_ID': sample_id}
        
        # Check for required CSV files
        fragle_csv = os.path.join(root, 'Fragle.csv')
        ht_csv = os.path.join(root, 'HT.csv')
        lt_csv = os.path.join(root, 'LT.csv')
        
        if os.path.exists(fragle_csv) and os.path.exists(ht_csv) and os.path.exists(lt_csv):
            # Read each CSV file
            fragle_df = pd.read_csv(fragle_csv)
            ht_df = pd.read_csv(ht_csv)
            lt_df = pd.read_csv(lt_csv)
            
            # Extract ctDNA burden values
            sample_data['fragle_ctDNAburden'] = fragle_df['ctDNA_Burden'].iloc[0]
            sample_data['HT_ctDNAburden'] = ht_df['ctDNA_Burden'].iloc[0]
            sample_data['LT_ctDNAburden'] = lt_df['ctDNA_Burden'].iloc[0]
            
            # Add to the list
            all_data.append(sample_data)
    
    # Convert to DataFrame
    result_df = pd.DataFrame(all_data)
    
    # Save to CSV
    output_path = os.path.join(os.path.dirname(fragle_dir), 'compiled_fragle_data.csv')
    result_df.to_csv(output_path, index=False)
    print(f"Compiled data saved to: {output_path}")
    return result_df

if __name__ == "__main__":
    # Get the FRAGLE directory path from the command line argument
    import sys
    if len(sys.argv) != 2:
        print("Usage: python compile_fragle_csv.py <path_to_FRAGLE_directory>")
        sys.exit(1)
    
    fragle_dir = sys.argv[1]
    if not os.path.exists(fragle_dir):
        print(f"Error: Directory {fragle_dir} does not exist")
        sys.exit(1)
    
    result = compile_fragle_data(fragle_dir)
    print("\nCompiled data preview:")
    print(result.head()) 