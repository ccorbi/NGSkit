
# NUCLEOTIDES

BASES = ['A' 'C',  'G', 'T']

dna_letters = "GATC"
dna_extended_letters = "GATCRYWSMKHBVDN"

rna_letters = "GAUC"
rna_extended_letters = "GAUCRYWSMKHBVDN"

dna_ambiguity = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
}

rna_ambiguity = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
}


nucleotide_names = {
    'A': 'Adenosine',
    'C': 'Cytidine',
    'G': 'Guanine',
    'T': 'Thymidine',
    'U': 'Uracil',
    'R': 'G A (puRine)',
    'Y': 'T C (pYrimidine)',
    'K': 'G T (Ketone)',
    'M': 'A C (aMino group)',
    'S': 'G C (Strong interaction)',
    'W': 'A T (Weak interaction)',
    'B': 'G T C (not A) (B comes after A)',
    'D': 'G A T (not C) (D comes after C)',
    'H': 'A C T (not G) (H comes after G)',
    'V': 'G C A (not T, not U) (V comes after U)',
    'N': 'A G C T (aNy)',
    '-': 'gap',
}


# AMINOACIDS

AMINOACIDS = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
    "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W", ]


THREE_LETTERS = {"ARG":"R", "HIS":"H", 
                "LYS":"K", "ASP": "D",
                "GLU": "E","SER": "S",
                "THR": "T", "ASN": "N",
                "GLN" :"Q", "CYS" :"C",
                "GLY": "G", "PRO" : "P",
                "ALA": "A", "VAL": "V", 
                "ILE": "I", "LEU": "L", 
                "MET":  "M", "PHE": "F", 
                "TYR" :"Y", "TRP": "W" }


amino_acid_letters = "ACDEFGHIKLMNPQRSTVWY"

amino_acid_alternative_letters = "ARNDCQEGHILKMFPSTWYV"

amino_acid_extended_letters = "ACDEFGHIKLMNOPQRSTUVWYBJZX*-"



amino_acid_ambiguity = {
    "A": "A",
    "B": "ND",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    "Y": "Y",
    "Z": "QE",
    "J": "IL",
    'U': 'U',
    'O': 'O',
}


extended_three_to_one = {
    '2as': 'D', '3ah': 'H', '5hp': 'E', 'Acl': 'R', 'Agm': 'R', 'Aib': 'A', 'Ala': 'A', 'Alm': 'A',
    'Alo': 'T', 'Aly': 'K', 'Arg': 'R', 'Arm': 'R', 'Asa': 'D', 'Asb': 'D', 'Ask': 'D', 'Asl': 'D',
    'Asn': 'N', 'Asp': 'D', 'Asq': 'D', 'Asx': 'B', 'Aya': 'A', 'Bcs': 'C', 'Bhd': 'D', 'Bmt': 'T',
    'Bnn': 'A', 'Buc': 'C', 'Bug': 'L', 'C5c': 'C', 'C6c': 'C', 'Ccs': 'C', 'Cea': 'C', 'Cgu': 'E',
    'Chg': 'A', 'Cle': 'L', 'Cme': 'C', 'Csd': 'A', 'Cso': 'C', 'Csp': 'C', 'Css': 'C', 'Csw': 'C',
    'Csx': 'C', 'Cxm': 'M', 'Cy1': 'C', 'Cy3': 'C', 'Cyg': 'C', 'Cym': 'C', 'Cyq': 'C', 'Cys': 'C',
    'Dah': 'F', 'Dal': 'A', 'Dar': 'R', 'Das': 'D', 'Dcy': 'C', 'Dgl': 'E', 'Dgn': 'Q', 'Dha': 'A',
    'Dhi': 'H', 'Dil': 'I', 'Div': 'V', 'Dle': 'L', 'Dly': 'K', 'Dnp': 'A', 'Dpn': 'F', 'Dpr': 'P',
    'Dsn': 'S', 'Dsp': 'D', 'Dth': 'T', 'Dtr': 'W', 'Dty': 'Y', 'Dva': 'V', 'Efc': 'C', 'Fla': 'A',
    'Fme': 'M', 'Ggl': 'E', 'Gl3': 'G', 'Gln': 'Q', 'Glu': 'E', 'Glx': 'Z', 'Gly': 'G', 'Glz': 'G',
    'Gma': 'E', 'Gsc': 'G', 'Hac': 'A', 'Har': 'R', 'Hic': 'H', 'Hip': 'H', 'His': 'H', 'Hmr': 'R',
    'Hpq': 'F', 'Htr': 'W', 'Hyp': 'P', 'Iil': 'I', 'Ile': 'I', 'Iyr': 'Y', 'Kcx': 'K', 'Leu': 'L',
    'Llp': 'K', 'Lly': 'K', 'Ltr': 'W', 'Lym': 'K', 'Lys': 'K', 'Lyz': 'K', 'Maa': 'A', 'Men': 'N',
    'Met': 'M', 'Mhs': 'H', 'Mis': 'S', 'Mle': 'L', 'Mpq': 'G', 'Msa': 'G', 'Mse': 'M', 'Mva': 'V',
    'Nem': 'H', 'Nep': 'H', 'Nle': 'L', 'Nln': 'L', 'Nlp': 'L', 'Nmc': 'G', 'Oas': 'S', 'Ocs': 'C',
    'Omt': 'M', 'Paq': 'Y', 'Pca': 'E', 'Pec': 'C', 'Phe': 'F', 'Phi': 'F', 'Phl': 'F', 'Pr3': 'C',
    'Pro': 'P', 'Prr': 'A', 'Ptr': 'Y', 'Pyl': 'O', 'Sac': 'S', 'Sar': 'G', 'Sch': 'C', 'Scs': 'C',
    'Scy': 'C', 'Sec': 'U', 'Sel': 'U', 'Sep': 'S', 'Ser': 'S', 'Set': 'S', 'Shc': 'C', 'Shr': 'K',
    'Smc': 'C', 'Soc': 'C', 'Sty': 'Y', 'Sva': 'S', 'Ter': '*', 'Thr': 'T', 'Tih': 'A', 'Tpl': 'W',
    'Tpo': 'T', 'Tpq': 'A', 'Trg': 'K', 'Tro': 'W', 'Trp': 'W', 'Tyb': 'Y', 'Tyq': 'Y', 'Tyr': 'Y',
    'Tys': 'Y', 'Tyy': 'Y', 'Unk': 'X', 'Val': 'V', 'Xaa': 'X', 'Xer': 'X', 'Xle': 'J'}



amino_acid_names = {
    'A': 'alanine',
    'M': 'methionine',
    'C': 'cysteine',
    'N': 'asparagine',
    'D': 'aspartic acid',
    'P': 'proline',
    'E': 'glutamic acid',
    'Q': 'glutamine',
    'F': 'phenylalanine',
    'R': 'arginine',
    'G': 'glycine',
    'S': 'serine',
    'H': 'histidine',
    'T': 'threonine',
    'I': 'isoleucine',
    'V': 'valine',
    'K': 'lysine',
    'W': 'tryptophan',
    'L': 'leucine',
    'Y': 'tyrosine',
    'B': 'aspartic acid or asparagine',
    'J': 'leucine or isoleucine',
    'X': 'unknown',
    'Z': 'glutamic acid or glutamine',
    'U': 'selenocysteine',
    'O': 'pyrrolysine',
    '*': 'translation stop',
    '-': 'gap'
}
