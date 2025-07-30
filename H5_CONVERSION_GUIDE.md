# PLINK to HDF5 Conversion Guide

This guide explains how to use the PLINK to HDF5 conversion utilities in the gensim package.

## Overview

The `convert_bed_to_h5.py` script converts PLINK binary files (.bed, .bim, .fam) to HDF5 format for more efficient storage and faster random access. This is particularly useful for large genomic datasets where you need to frequently access subsets of data.

## Installation Requirements

```bash
pip install h5py pandas-plink
```

Or install all requirements:
```bash
pip install -r requirements.txt
```

## Basic Usage

### Command Line Usage

**Auto mode (recommended) - automatically determines processing method:**
```bash
python convert_bed_to_h5.py --input-dir data/simulations --auto
```

Convert all .bed files in a directory:
```bash
python convert_bed_to_h5.py --input-dir data/simulations
```

Convert files and save to a different directory:
```bash
python convert_bed_to_h5.py --input-dir data/simulations --output-dir data/h5_files
```

**Use chunked processing for large files:**
```bash
python convert_bed_to_h5.py --input-dir data --chunked --memory-usage 4.0
```

Use chunked processing with specific chunk size:
```bash
python convert_bed_to_h5.py --input-dir data --chunked --chunk-size 1000
```

Enable verbose logging:
```bash
python convert_bed_to_h5.py --input-dir data --verbose
```

Verify existing H5 files without converting:
```bash
python convert_bed_to_h5.py --input-dir data --verify-only
```

Save logs to a file:
```bash
python convert_bed_to_h5.py --input-dir data --log-file conversion.log
```

### Programmatic Usage

**Auto mode (recommended):**
```python
from convert_bed_to_h5 import BedToH5Converter

# Initialize converter
converter = BedToH5Converter(
    input_dir="data/simulations",
    output_dir="data/h5_files"
)

# Convert all files using automatic processing mode selection
converted_files = converter.convert_all_auto()

# Verify converted files
for h5_file in converted_files:
    is_valid = converter.verify_h5_file(h5_file)
    print(f"{h5_file}: {'Valid' if is_valid else 'Invalid'}")
```

**Standard mode:**
```python
from convert_bed_to_h5 import BedToH5Converter

# Initialize converter
converter = BedToH5Converter(
    input_dir="data/simulations",
    output_dir="data/h5_files"
)

# Convert all files
converted_files = converter.convert_all()

# Verify converted files
for h5_file in converted_files:
    is_valid = converter.verify_h5_file(h5_file)
    print(f"{h5_file}: {'Valid' if is_valid else 'Invalid'}")
```

**Chunked mode (for large datasets):**
```python
from convert_bed_to_h5 import BedToH5Converter

# Initialize converter with chunked processing
converter = BedToH5Converter(
    input_dir="data/simulations",
    output_dir="data/h5_files",
    chunk_size=1000,      # Process 1000 samples at a time
    memory_usage=4.0      # Target 4GB memory usage
)

# Convert all files using chunked processing
converted_files = converter.convert_all_chunked()

# Verify converted files
for h5_file in converted_files:
    is_valid = converter.verify_h5_file(h5_file)
    print(f"{h5_file}: {'Valid' if is_valid else 'Invalid'}")
```

## Reading HDF5 Files

### Using the H5PLINKReader Class

```python
from gensim.h5_utils import H5PLINKReader

# Open and read HDF5 file
with H5PLINKReader("data/simulation.h5") as reader:
    # Get file information
    print(reader.info())
    
    # Get genotype data
    genotypes = reader.get_genotypes()
    print(f"Genotype matrix shape: {genotypes.shape}")
    
    # Get sample and variant information
    samples_df = reader.get_samples_df()
    variants_df = reader.get_variants_df()
    
    # Get subset by variant IDs
    subset_geno = reader.subset_by_variant_ids(['rs1234', 'rs5678'])
    
    # Compute allele frequencies
    frequencies = reader.compute_allele_frequencies()
```

### Using h5py Directly

```python
import h5py
import numpy as np

with h5py.File("data/simulation.h5", 'r') as h5f:
    # Print file structure
    print("Datasets:", list(h5f.keys()))
    print("Metadata:", dict(h5f.attrs))
    
    # Load genotype data
    genotypes = h5f['genotypes'][:]
    
    # Load variant information
    variant_ids = h5f['variants/snp'][:]
    chromosomes = h5f['variants/chrom'][:]
    positions = h5f['variants/pos'][:]
    
    # Load sample information
    sample_ids = h5f['samples/iid'][:]
    family_ids = h5f['samples/fid'][:]
```

## File Structure

The converter creates HDF5 files with the following structure:

```
simulation.h5
├── genotypes               # Genotype matrix (n_samples × n_variants)
├── variants/               # Variant information from .bim file
│   ├── chrom              # Chromosome
│   ├── snp                # Variant ID
│   ├── cm                 # Genetic distance
│   ├── pos                # Base pair position
│   ├── a1                 # Allele 1
│   └── a2                 # Allele 2
├── samples/                # Sample information from .fam file
│   ├── fid                # Family ID
│   ├── iid                # Individual ID
│   ├── father             # Father ID
│   ├── mother             # Mother ID
│   ├── sex                # Sex (1=male, 2=female)
│   └── phenotype          # Phenotype
└── attributes              # Metadata
    ├── n_samples          # Number of samples
    ├── n_variants         # Number of variants
    ├── source_bed_file    # Original .bed file path
    └── encoding           # Genotype encoding information
```

## Chunked Processing

For large datasets (>100K samples or >1M variants), chunked processing is recommended to avoid memory issues.

### How Chunked Processing Works

1. **Automatic chunk size calculation**: Based on target memory usage
2. **Progressive processing**: Samples are processed in batches
3. **Memory efficiency**: Only loads one chunk into memory at a time
4. **Progress tracking**: Shows progress through chunks with progress bars

### Chunk Size Calculation

The chunk size is automatically calculated based on:
- Number of variants in the dataset
- Target memory usage (default: 2.0 GB)
- Available system memory

Formula: `chunk_size = target_memory_bytes / (n_variants × 4 bytes_per_genotype)`

### When to Use Chunked Processing

- **Large datasets**: >100K samples or >500K variants
- **Limited memory**: Systems with <16GB RAM
- **Network storage**: When reading from slow file systems
- **Batch processing**: SLURM jobs with memory constraints

### Example Chunk Sizes

| Dataset Size | Memory Usage | Recommended Chunk Size |
|-------------|--------------|----------------------|
| 50K samples × 500K variants | 2 GB | ~1,000 samples |
| 100K samples × 1M variants | 4 GB | ~1,000 samples |
| 500K samples × 500K variants | 4 GB | ~2,000 samples |
| 1M samples × 1M variants | 8 GB | ~2,000 samples |

## Genotype Encoding

The genotype matrix uses the standard PLINK encoding:
- `0`: Missing genotype
- `1`: Heterozygous genotype
- `2`: Homozygous alternate genotype

## Features

### File Discovery
- Recursively searches directories for .bed files
- Automatically excludes files with "temporary" in the name
- Validates that corresponding .bim and .fam files exist

### Efficient Storage
- Uses gzip compression for reduced file size
- Maintains original data structure and metadata
- Preserves directory structure in output

### Memory-Efficient Processing
- **Chunked processing mode** for large datasets
- **Automatic chunk size calculation** based on available memory
- **Progress tracking** with tqdm progress bars
- Processes samples in batches to avoid memory issues

### Data Validation
- Verifies file integrity after conversion
- Checks for required datasets and attributes
- Provides detailed error reporting

### Flexible Usage
- Command-line interface with extensive options
- Programmatic Python API
- Support for both standard and chunked processing modes
- **Automatic processing mode selection** based on system memory and dataset size

## Automatic Processing Mode Selection

The converter can automatically determine the best processing mode based on:

### Memory Detection
- **System memory analysis**: Uses `psutil` to detect available RAM
- **Dataset size calculation**: Estimates memory required for full dataset loading
- **Safety factor**: Only uses 70% of available memory to prevent system issues
- **Fallback handling**: Works even when `psutil` is not available

### Decision Logic
```
Required Memory = n_samples × n_variants × 4 bytes (float32)
Available Safe Memory = Available RAM × 0.7

If Required Memory > Available Safe Memory:
    → Use chunked processing
Else:
    → Use standard processing
```

### Auto-adjustment Features
- **Dynamic memory targeting**: Adjusts chunk size based on available memory
- **Conservative fallback**: Uses chunking for >1GB datasets when memory detection fails
- **Transparent operation**: Logs all decisions for user awareness

### Usage Recommendations
1. **Default choice**: Use `--auto` for most scenarios
2. **Large datasets**: Auto mode will automatically select chunking
3. **Small datasets**: Auto mode will use faster standard processing
4. **Memory-constrained systems**: Auto mode prevents out-of-memory errors

## Examples

### Example 1: Basic Conversion with Auto Mode

```python
# Run the example script
python bed_to_h5_examples.py

# Or use auto mode directly
python convert_bed_to_h5.py --input-dir data --auto --verbose
```

### Example 2: Batch Processing with SLURM

```bash
# Create a SLURM script for chunked conversion
cat > convert_h5.slurm << 'EOF'
#!/bin/bash
#SBATCH --job-name=bed_to_h5
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1

# Activate environment
source activate gensim

# Convert files using auto mode (recommended)
python convert_bed_to_h5.py \
    --input-dir data/ukb_simulation \
    --output-dir data/h5_files \
    --auto \
    --verbose \
    --log-file conversion_${SLURM_JOB_ID}.log
EOF

# Submit job
sbatch convert_h5.slurm
```

### Example 3: Reading and Analysis

```python
from gensim.h5_utils import H5PLINKReader, list_h5_files
import numpy as np

# Find all H5 files
h5_files = list_h5_files("data/h5_files")
print(f"Found {len(h5_files)} H5 files")

# Analyze each file
for h5_file in h5_files:
    with H5PLINKReader(h5_file) as reader:
        print(f"\nFile: {h5_file.name}")
        print(f"Samples: {reader.n_samples:,}")
        print(f"Variants: {reader.n_variants:,}")
        
        # Compute basic statistics
        freqs = reader.compute_allele_frequencies()
        print(f"Mean allele frequency: {np.mean(freqs):.4f}")
        print(f"MAF < 0.05: {np.sum(freqs < 0.05):,} variants")
```

## Performance Considerations

### File Sizes
- HDF5 files are typically 20-30% smaller than uncompressed PLINK files
- Compression level can be adjusted in the converter code
- Reading is significantly faster than PLINK text files

### Memory Usage
- Converter loads one file at a time to minimize memory usage
- Reader supports partial loading of large datasets
- Consider using chunked reading for very large files (>1M variants)

### Speed
- Conversion speed: ~1-2 minutes per 100K variants
- Reading speed: ~10x faster than PLINK text files
- Random access: Much faster than sequential PLINK reading

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Install h5py and pandas-plink
2. **Memory errors**: Reduce batch size or use chunked processing
3. **File permission errors**: Check write permissions in output directory
4. **Corrupted files**: Use `--verify-only` to check file integrity

### Error Messages

- "Missing .bim/.fam file": Ensure complete PLINK filesets
- "pandas-plink import error": Install with `pip install pandas-plink`
- "HDF5 file not valid": File may be corrupted, reconvert from PLINK

### Performance Tips

1. Use SSD storage for better I/O performance
2. Process files in parallel using job arrays
3. Use compression for network file systems
4. Monitor memory usage with large datasets

## Integration with Gensim

The H5 utilities integrate seamlessly with the gensim simulation pipeline:

```python
from gensim import PLINKParameterGridSimulator
from gensim.h5_utils import list_h5_files, H5PLINKReader

# Run simulations
simulator = PLINKParameterGridSimulator(...)
simulator.run_all()

# Convert results to H5
os.system("python convert_bed_to_h5.py --input-dir data/simulations")

# Analyze results
h5_files = list_h5_files("data/simulations")
for h5_file in h5_files:
    with H5PLINKReader(h5_file) as reader:
        # Perform analysis
        pass
```
