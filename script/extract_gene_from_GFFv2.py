#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import re
from collections import defaultdict

# Script information
SCRIPT_NAME = "gff_to_bed"
VERSION = "2.0.0"
DESCRIPTION = "Extract gene ranges from CDS features in GFF file and convert to BED format"

def show_help():
    """Display help information"""
    help_text = f"""
    {SCRIPT_NAME} (v{VERSION}) - {DESCRIPTION}

    Usage: {SCRIPT_NAME} <input.gff>

    Arguments:
    <input.gff>     Input GFF/GFF3 file containing CDS annotations

    Options:
    --idItem        Key(s) to extract gene ID from attributes (comma-separated)
    -o, --output    Output prefix (default: "output")

    Output:
    - output_genefromgff.bed  BED format file containing gene coordinates and IDs

    Description:
    This script extracts CDS features from a GFF file, groups them by parent gene,
    calculates the min/max coordinates, and outputs gene ranges in BED format.

    Example:
    {SCRIPT_NAME} annotations.gff --idItem ID,Name
    """
    print(help_text)

def process_gff(input_gff, id_items, output):
    """Process GFF file and extract gene ranges from CDS features"""
    try:
        # Dictionary to store CDS ranges by gene ID
        gene_dict = defaultdict(lambda: {
            'chrom': None,
            'starts': [],
            'ends': [],
            'strand': None,
            'attributes': ""
        })

        with open(input_gff, 'r') as infile:
            for line in infile:
                if line.startswith('#') or line.strip() == '':
                    continue  # Skip comments and empty lines
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue  # Skip invalid lines
                
                feature_type = fields[2]
                if feature_type.lower() == "cds":                
                    # Extract basic fields
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    attributes = fields[8]
                    
                    # Extract parent ID
                    parent_id = None
                    for item in id_items:
                        pattern = rf"{item}[=\s][\"']?([A-Za-z0-9_.\-:|]+)[\"']?"
                        match = re.search(pattern, attributes)
                        if match:
                            parent_id = match.group(1)
                            break
                    
                    if not parent_id:
                        continue  # Skip CDS without parent
                    
                    # Extract gene ID from parent's attributes (if available)
                    gene_id = parent_id  # Default to parent ID if no gene ID found
                    # Update gene dictionary
                    gene_dict[gene_id]['chrom'] = chrom
                    gene_dict[gene_id]['starts'].append(start)
                    gene_dict[gene_id]['ends'].append(end)
                    gene_dict[gene_id]['strand'] = strand
                    gene_dict[gene_id]['attributes'] = attributes
            
                if feature_type.lower() == "transcript":                
                    # Extract basic fields
                    attributes = fields[8]
                    
                    # Extract parent ID
                    parent_id = None
                    for item in ["gene_id","Parent"]:
                        pattern = rf"{item}[=\s][\"']?([A-Za-z0-9_.\-:|]+)[\"']?"
                        match = re.search(pattern, attributes)
                        if match:
                            parent_id = match.group(1)
                            break
                    if not parent_id:
                        continue  # Skip gene without parent
                        
                    transcript_id = None
                    for item in ["transcript_id","ID"]:
                        pattern = rf"{item}[=\s][\"']?([A-Za-z0-9_.\-:|]+)[\"']?"
                        match = re.search(pattern, attributes)
                        if match:
                            transcript_id = match.group(1)
                            break

                    gene_dict[transcript_id]['gene'] = parent_id
                    
        # Write to BED file
        count = 0
        with open(f'{output}_genefromgff.bed', 'w') as outfile:
            for gene_id, data in gene_dict.items():
                if not data['starts']:
                    continue
                
                min_start = min(data['starts'])
                max_end = max(data['ends'])
                gene=data['gene']
                
                outfile.write(f"{data['chrom']}\t{min_start}\t{max_end}\t{gene_id}\t.\t{data['strand']}\t{gene}\n")
                count += 1
        
        return count
    
    except Exception as e:
        print(f"Error: Failed to process files - {str(e)}", file=sys.stderr)
        return None

def main(args_list=None):
    """Main function to handle command-line arguments and execution"""
    parser = argparse.ArgumentParser(
        prog=SCRIPT_NAME,
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-g", "--gxf", required=True, help='Input GFF/GFF3 file')
    parser.add_argument("--idItem", 
                       help='Comma-separated keys to extract gene ID from attributes')
    parser.add_argument("-o", "--output", default="output", 
                       help='Output prefix (default: "output")')
    
    #args = parser.parse_args(args_list) if args_list else parser.parse_args()
    args, _ = parser.parse_known_args(args_list) if args_list else parser.parse_args()
    # Check input file exists
    if not os.path.isfile(args.gxf):
        print(f"Error: Input GFF file '{args.gxf}' not found", file=sys.stderr)
        show_help()    
        sys.exit(1)

    print(f"Processing GFF file: {args.gxf}")
    print("Extracting gene ranges from CDS features...")
    
    # Process ID items
    id_items = ["transcript_id","ID", "gene_id", "Name"]  # Default search terms
    if args.idItem:
        id_items = [item.strip() for item in args.idItem.split(",")]
    
    # Process the file
    count = process_gff(args.gxf, id_items, args.output)
    
    if count is None:
        sys.exit(1)  # Error already printed by process_gff

    if count == 0:
        print("Warning: No valid CDS features with parent genes found", file=sys.stderr)
        if os.path.exists(f'{args.output}_genefromgff.bed'):
            os.remove(f'{args.output}_genefromgff.bed')
        sys.exit(1)
    
    print(f"Conversion completed successfully (extracted {count} gene ranges)")
    print(f"Output saved to: {args.output}_genefromgff.bed")
    
    # Show sample output
    print("\nFirst 5 lines of output:")
    try:
        with open(f"{args.output}_genefromgff.bed", 'r') as f:
            for i, line in enumerate(f):
                if i >= 5:
                    break
                print(line.strip().expandtabs(10))
    except IOError:
        print("Could not display sample output")

if __name__ == "__main__":
    main()
