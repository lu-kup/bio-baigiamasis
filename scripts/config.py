CHR18_SUBSET_SIZE = 18261711

CHR_SIZES_GRCm38 = {
    "1": 195471971, "2": 182113224, "3": 160039680, "4": 156508116,
    "5": 151834684, "6": 149736546, "7": 145441459, "8": 129401213,
    "9": 124595110, "10": 130694993, "11": 122082543, "12": 120129022,
    "13": 120421639, "14": 124902244, "15": 104043685, "16": 98810272,
    "17": 93617343, "18": 90153841, "19": 61026292, "X": 171031299,
    "Y": 91000000, "MT": 16299
}

# Chromosome
CHR = '18'

# Inputs
SAMPLE_FILEPATH = '../inputs/subset2.rds'
GENCODE_DATA_PATH = "../inputs/gencode_chr18_M25.gtf"

# Signal columns
SIGNAL_COLUMN = 'TT_S0'
OTHER_SIGNALS = ['TT_S1', 'TT_S2']

# Subset size
CHR_SIZE = CHR_SIZES_GRCm38[CHR]
