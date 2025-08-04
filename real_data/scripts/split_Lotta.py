import argparse
import gzip
import os

def split_Lotta(filename):
    # The Lotta et al. 2021 paper argues that no filtering is needed based on heterogeneity of meta-analysis effect sizes, thus we don't include this information
    names = [b'trait', b'chr', b'pos', b'Allele1', b'Allele2', b'Freq1_MA', b'Weight_MA', b'Zscore_MA']
    traits = {}

    with gzip.open(filename, 'r') as Lotta_file:
        header = {name: ix for ix, name in enumerate(Lotta_file.readline().split()) if name in names}

        for line in Lotta_file:
            fields = line.split()
            trait_data = [fields[header[name]] for name in names]
            
            trait = trait_data[0]
            
            if trait not in traits:
                traits[trait] = open(trait + b'.txt', 'ba')
                traits[trait].write(b'\t'.join(names) + b'\n')
                
            traits[trait].write(b'\t'.join(trait_data) + b'\n')

    for trait in traits:
        traits[trait].close()
        os.system(b'gzip ' + trait + b'.txt')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', type = str, required = True)
    args = parser.parse_args()

    split_Lotta(args.filename)
