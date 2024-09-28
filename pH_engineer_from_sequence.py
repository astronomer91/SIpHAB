#!/Users/wanleiwei/miniconda3/bin/python3.8
import sys
import os
import requests
from pathlib import Path
from abnumber import Chain

# Must install abnumber using conda: instructions: https://abnumber.readthedocs.io/en/latest/


# exclusion data points - if less than 200 data points, then insufficient to make any conclusion
excl_threshold = 60

# data of enrichments (in the scriptdir
his_data_pos_TME = Path(__file__).with_name("His_pos_preference_TME.txt")
his_data_res_TME = Path(__file__).with_name("His_res_preference_TME.txt")

his_data_pos_recycl = Path(__file__).with_name("His_pos_preference_recycling.txt")
his_data_res_recycl = Path(__file__).with_name("His_res_preference_recycling.txt")


amino_acid_dict = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

amino_acid_dict_rev = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}



def main():

    if (not (len(sys.argv) == 3)):
        print("usage: " + sys.argv[0] + "  [sequence.txt] [direction]")
        print("direction could be TME or recycling") 
        print("sequence file should be antibody sequences as one-letter codes, separated by a newline")
        exit()

    # input native sequence of heavy and light chain
    input_seq_name = sys.argv[1]
    input_seq = open(input_seq_name, "r")

    mutation_direction=(sys.argv[2]).lower()
    sequences = []
    seq_kabat = []

    for line in input_seq:
        size = len(line)
        if size > 40:
            sequences += [line.strip("\n")]

    # turn it into Kabat antibody sequence using abnum
    for seq in sequences:

        curr = ""
        try:
            numbered = Chain(seq, scheme="kabat")
            for pos, aa in numbered:
                curr = curr + str(pos) + " " + aa + "\n"

        except:
            print(seq + " is not a variable chain and was discarded") 


        seq_kabat.append(curr)

    seq_dict = {}
    for chain in (seq_kabat):
        for aa in chain.split("\n"):
            if (aa.split() != []):
                ab_seq = aa.split()

                res_pos = ab_seq[0]
                aa_type = ab_seq[1]

                # seq_dict contains kabat-format sequence of input antibody
                seq_dict[res_pos] = aa_type


    # load in POSITIONAL and AMINO ACID preferences of toHIS

    out_app = ""

    if(mutation_direction[0] == "t"):
       his_data_pos = his_data_pos_TME
       his_data_res = his_data_res_TME
       out_app = "tme"

    elif(mutation_direction[0] == "r"):
       his_data_pos = his_data_pos_recycl
       his_data_res = his_data_res_recycl
       out_app = "recycl"
        
    # toHIS: his_data_pos, his_data_res
    his_pos = read_pos_data(his_data_pos)
    his_res = read_res_data(his_data_res)

    # if any of the sequences are enriched (ie. above 1)
    his_EF = synthesis_EF(his_pos, his_res)

    # output sequence with promising pH mutations- sequence matters: asp, glu, his
    EFs = identify_enrichments(seq_dict, his_EF)


    output_file = open(input_seq_name.split(".")[0] + "_enrichment_" + out_app  + ".out", "w")
    output_file.writelines(input_seq_name + "\n")
    output_file.writelines("---------------------------------\n")
    output_file.writelines("Position  Mutation  pH Enrichment\n")
    output_file.writelines("---------------------------------\n")
    rank = 0
    for mutation in EFs:

        mut_vec = mutation.split()

        if len(sys.argv) == 4:
            aa_incl = sys.argv[3].upper()
            if (aa_incl in amino_acid_dict):            
                incl = "->" + amino_acid_dict[aa_incl]
            else:
                print("Error: " + aa_incl + " is a valid amino acid for restriction")
                exit()

            if incl in mutation:
                rank += 1
           

                output_file.writelines(mut_vec[0] + (10 - len(mut_vec[0]))*" " + mut_vec[1] + (10 - len(mut_vec[1]))*" " + str(EFs[mutation]) + (15 - len(str(EFs[mutation]) + str(rank)))*" " + str(rank)  + "\n")


        
        else:
            rank += 1
            output_file.writelines(mut_vec[0] + (10 - len(mut_vec[0]))*" " + mut_vec[1] + (10 - len(mut_vec[1]))*" " + str(EFs[mutation]) + (15 - len(str(EFs[mutation]) + str(rank)))*" " + str(rank)  + "\n")
    output_file.close()



def read_pos_data(data_file_name):

    pos_dict = {}
    data_file = open(data_file_name,"r")

    for line in data_file.readlines():

        line_vec = (line.strip("\n")).split()

        if len(line_vec) < 7 or line_vec[0] == "Position":
            continue

        res_pos = line_vec[0]
        EF = float(line_vec[6])
        sample = float(line_vec[5])

        if (sample > excl_threshold):
            pos_dict[res_pos] = EF

    return pos_dict


def read_res_data(data_file_name):

    res_dict = {}
    data_file = open(data_file_name,"r")

    for line in data_file.readlines():
        line_vec = (line.strip("\n")).split()

        if len(line_vec) < 7 or line_vec[0] == "Position":
            continue

        res_pos = line_vec[0]
        EF = float(line_vec[6])
        sample = float(line_vec[5])

        if (sample > excl_threshold):
            res_dict[res_pos] = EF

    return res_dict


def synthesis_EF(pos_enrichment, res_enrichment):
    EF_dict = {}

    for pos in pos_enrichment:
        for res in res_enrichment:
            pos_res = pos + " " + res
            EF = round(pos_enrichment[pos] * res_enrichment[res], 4)

            EF_dict[pos_res] = EF
    return EF_dict




def identify_enrichments(seq_dict, his_EF):
    all_EF = {}


    for pos in seq_dict:
        res = amino_acid_dict_rev[seq_dict[pos]]
        pos_res = pos + " " + res

        if (pos_res in his_EF):
            mut = pos + " " + seq_dict[pos] + "->H"
            EF = his_EF[pos_res]
            all_EF[mut] = EF


    all_EF = dict(sorted(all_EF.items(), key=lambda item: item[1], reverse=True))
    return all_EF



if __name__ == '__main__':
    main()





