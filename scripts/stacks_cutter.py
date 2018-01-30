from __future__ import print_function
import argparse
import sys
from Bio import SeqIO
from Bio.Restriction import *


def multicut(args):
    try:
        stacks_re = args["enzymes"]
    except KeyError:
        stacks_re = ['ApeKI', 'ApoI', 'BamHI', 'BglII', 'BstYI', 'ClaI',
                     'DpnII', 'EaeI', 'EcoRI', 'EcoRV', 'EcoT22I', 'HindIII', 'KpnI',
                     'MluCI', 'MseI', 'MspI', 'NdeI', 'NheI', 'NlaIII', 'NotI', 'NsiI',
                     'PstI', 'SacI', 'Sau3AI', 'SbfI', 'SexAI', 'SgrAI', 'SpeI',
                     'SphI', 'TaqI', 'XbaI', 'XhoI']

    print("! Trying the following enzymes:\n! {}".format(",".join(stacks_re)), file=sys.stderr)
    stacks_re = map(lambda x: x[:1].capitalize() + x[1:], stacks_re)

    enzymes = __import__("Bio.Restriction", fromlist=stacks_re)
    for enzyme in stacks_re:
        try:
            getattr(enzymes, enzyme)
        except AttributeError:
            raise ValueError("No such enzyme available: '{}'".format(enzyme))

    print("! Reading input from fasta file '{}'...".format(args["fasta"]), file=sys.stderr)

    with open(args["fasta"], 'rU') as f:
        seq_records = list(SeqIO.parse(f, "fasta"))
    print("! Calculating length of input fasta sequences...", file=sys.stderr)
    total_length_genome = sum(map(len, seq_records))
    print("! Total length fasta sequences: {:,} bp".format(total_length_genome), file=sys.stderr)
    print("! Sequencing set-up is {}x{} bp".format(args["reads"], args["readlen"]))
    print("! Selecting restriction fragments between {} and {} bp".format(args["lower_bound"], args["upper_bound"]))
    print("! Gathering GBS stats by enzyme\n", file=sys.stderr)
    res_map_full = {}
    for record in seq_records:
        # Use Biopythons digestion analysis function to get the coordinates of the cut size for all enzymes
        anl = Analysis(RestrictionBatch([getattr(enzymes, enz) for enz in stacks_re]), record.seq, linear=True)
        res_map = anl.full()
        len_map = {}
        seq_len = len(record.seq)
        for key, val in res_map.iteritems():
            # Get the fragment sizes from cut site positions
            res_len = map(lambda f: f[1]-f[0], zip([0]+val, val+[seq_len]))
            res_frags = map(lambda f: [(record.id, f[0], f[0]+args["readlen"]), 
                (record.id, f[1]-args["readlen"], f[1])], zip([0]+val, val+[seq_len]))
            res_filtered = []
            for i in range(0,len(res_len)):
                if res_len[i] >= args["lower_bound"] and res_len[i] <= args["upper_bound"]:
                    res_filtered.append(res_frags[i])

            len_map[key.__name__] = res_filtered

        for key, lens in len_map.iteritems():
            try:
                res_map_full[key].extend(lens)
            except KeyError:
                res_map_full[key] = lens

    print("enzyme\t#size_selected_frags\trad_sites_length_(bp)\t#rad_sites\t%genome_in_rad_sites")
    for key, lens in res_map_full.iteritems():
        nfrags = len(lens)
        sumfrags = nfrags * int(args["reads"]) * int(args["readlen"])
        captured = (nfrags * float(args["reads"]) * float(args["readlen"]) / total_length_genome) * 100.0
        print("{}\t{}\t{}\t{}\t{:.2f}".format(
            key, nfrags, sumfrags, nfrags * int(args["reads"]), captured))

    if args["bed"] is not None:
        for key, sites in res_map_full.iteritems():
            fname = "{}_{}_{}-{}.bed".format(args["bed"], key, args["lower_bound"], args["upper_bound"])
            print("! writing RADsites to file {}".format(fname))
            with open(fname, "w") as f:
                for site2 in sites:
                    for s in site2:
                        line = "{}\t{}\t{}\trad\t999\t.\n".format(s[0], s[1], s[2])
                        f.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Script to in silico digest genomes for planning GBS / RADseq experiments at NGI. The idea is to use
    assemblies of closely related species to inform how many RAD sites / % of the genome you would expect
    to sequence, given a choice of sequencing setup and a size selection interval.

    Outputs a tab-separated table to stdOut.
    """)
    parser.add_argument("-f", "--fasta", required=True, help="Fasta file to be digested")
    parser.add_argument("-l", "--lower-bound", required=True, type=int, help="Lower bound for in silico size selection")
    parser.add_argument("-u", "--upper-bound", required=True, type=int, help="Upper bound for in silico size selection")
    parser.add_argument("-r", "--reads", default=2, help="Number of reads, ie. 1 for single-end sequencing and 2 for paired-end")
    parser.add_argument("-n", "--readlen", type=int, default=126, help="Read length")
    parser.add_argument("-e", "--enzymes", nargs="*", metavar=("EnzI", "EnzII"), default=argparse.SUPPRESS,
                        help="Give a list of enzymes to use. If left out, I will use all the ones supported by Stacks.")
    parser.add_argument("-b", "--bed", type=str, help="(optional) file prefix of output BED6 file(s)")
    args = vars(parser.parse_args())
    multicut(args)
