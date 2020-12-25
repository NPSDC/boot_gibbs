import gzip
import pickle
import argparse
import numpy as np
import sys

def read_fastq(file, pi_dump = None, read_count = True):
    trans_dict={}
    with gzip.open(file) as ifile:
        for line in ifile:
            line = line.decode("utf-8").rstrip()
            if(line.startswith(">")):
                req_line = line.split(";")[0].split("/")
                trans = req_line[1]
                read = int(req_line[0][5:])
                if(not trans in trans_dict.keys()):
                    trans_dict[trans] = []
                trans_dict[trans].append(read)
    if not pi_dump is None:
        pickle.dump(trans_dict, open(pi_dump+"_dict.pi", "wb"))
    return trans_dict
    # if(read_count):
    #     read_trans = [len(trans_dict[eq]) for eq in trans_dict.keys()]
    #     if(not pi_dump is None):
    #         pickle.dump(read_trans, open(pi_dump+"_reads.pi", "wb"))
    
def write_read_counts(trans_dict, file_name):
    read_counts = {}
    with open(file_name, "w") as f:
        f.write("Transcript\tCounts\n")
        for trans in trans_dict.keys():
            read_counts[trans] = len(trans_dict[trans])
            f.write(trans+"\t"+str(len(trans_dict[trans]))+"\n")
    return(read_counts)

"""Extracts all the transcripts that occur with the given transcript"""
def extract_eq_trans(trans_list, eqcl_name, transcript):
#    samp_dict = pickle.load(open(pi_file, "rb"))
#    trans = samp_dict["tnames"]
#    eqcl_name = samp_dict["eqclass_name"]

    eqcl = [] ##Storing equivalence class indexes
    eqtrans = []
    if(transcript not in trans_list):
        sys.exit("Invalid transcript")

    trans_ind = np.where(np.array(trans_list) == transcript)[0][0]
    print(trans_ind)    
    for eq_ind, eq_trans in eqcl_name.items():
        #print(eq_trans)
        if(trans_ind in eq_trans):
            eqcl.append(eq_ind)
            eqtrans.extend(list(eq_trans))
    
    return eqcl, [trans[t] for t in set(eqtrans)]

def main():
    parser = argparse.ArgumentParser(
        description="reading fastq files to read counts"
    )
    parser.add_argument(
        '-f', '--file', required=True, type=str, help='fastq file')
    parser.add_argument(
        '-d', '--dump', required=True, type=str, help='Dump name')
    parser.add_argument(
        '-o', '--output', required=True, type=str, help='Output file for storing contents')
    args = parser.parse_args()
    trans_dict = read_fastq(args.file, args.dump)
    read_counts = write_read_counts(trans_dict, args.output)

if __name__ == "__main__":
    main()