#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A concise Python script to convert multiple-sequence alignment (MSA) formats:
fas, a2m, a3m, sto, psi, clu, ufas.

Essential transformations:
 - Insert→match logic (with -M 'first' or -M <int>)
 - Remove inserts/columns (-r), remove sequences with only gaps
 - Convert gap symbols ('.' vs '-'), uppercase/lowercase, etc.

Example:
    python reformat_msa.py fas a3m input.fas output.a3m -num -r 90 -g '-'
"""

import sys
import os
import re
import glob
import argparse

FORMATS_IN = ["fas","a2m","a3m","sto","psi","clu"]
FORMATS_OUT = ["fas","a2m","a3m","sto","psi","clu","ufas"]

def read_fas_a2m_a3m(path, remove_sa, remove_ss):
    with open(path) as f:
        data = f.read()
    blocks = data.split('>')[1:]  # ignore any leading stuff before first '>'
    names, seqs = [], []
    for block in blocks:
        lines = block.strip().splitlines()
        if not lines: 
            continue
        name = lines[0].strip()
        s = "".join(lines[1:]).replace(" ","")
        if name.startswith("aa_"): continue
        if remove_sa and name.startswith("sa_"): continue
        if remove_ss and name.startswith("ss_"): continue
        names.append(name)
        seqs.append(s)
    return names, seqs

def read_stockholm(path, remove_sa, remove_ss):
    names, seqs, name2idx = [], [], {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): 
                if line.startswith('//'): break
                continue
            parts = line.split()
            if len(parts)<2: 
                continue
            name, aln = parts[0], parts[1]
            if remove_sa and name.startswith("sa_"): continue
            if remove_ss and name.startswith("ss_"): continue
            if name.startswith("aa_"): continue
            if name not in name2idx:
                name2idx[name] = len(names)
                names.append(name)
                seqs.append(aln)
            else:
                seqs[name2idx[name]] += aln
    return names, seqs

def read_clustal(path, remove_sa, remove_ss):
    names, seqs, name2idx = [], [], {}
    with open(path) as f:
        block = []
        for line in f:
            line = line.rstrip()
            if not line or line.upper().startswith("CLUSTAL") or line.startswith('#'):
                _process_clu_block(block, names, seqs, name2idx, remove_sa, remove_ss)
                block = []
            elif line.startswith('//'):
                _process_clu_block(block, names, seqs, name2idx, remove_sa, remove_ss)
                break
            else:
                block.append(line)
        else:
            # leftover
            _process_clu_block(block, names, seqs, name2idx, remove_sa, remove_ss)
    return names, seqs

def _process_clu_block(lines, names, seqs, name2idx, remove_sa, remove_ss):
    """
    Minimal parse of Clustal block lines.
    """
    for line in lines:
        # skip lines that are all alignment markers: * : .
        if re.match(r'^[.*:\s]+$', line): 
            continue
        m = re.match(r'^(\S+)\s+([\w.\-]+)', line)
        if not m: 
            continue
        name, chunk = m.group(1), m.group(2).replace(" ","")
        if name.startswith("aa_"): continue
        if remove_sa and name.startswith("sa_"): continue
        if remove_ss and name.startswith("ss_"): continue
        if name not in name2idx:
            name2idx[name] = len(names)
            names.append(name)
            seqs.append(chunk)
        else:
            idx = name2idx[name]
            seqs[idx] += chunk

def read_psi(path, remove_sa, remove_ss):
    names, seqs, name2idx = [], [], {}
    with open(path) as f:
        blocks = []
        block = []
        for line in f:
            line = line.rstrip()
            if not line:
                if block: blocks.append(block); block=[]
            else:
                block.append(line)
        if block: blocks.append(block)
    for bl in blocks:
        for line in bl:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts)<2: continue
            name, aligned = parts[0], "".join(parts[1:])
            if name.startswith("aa_"): continue
            if remove_sa and name.startswith("sa_"): continue
            if remove_ss and name.startswith("ss_"): continue
            if name not in name2idx:
                name2idx[name] = len(names)
                names.append(name)
                seqs.append(aligned)
            else:
                seqs[name2idx[name]] += aligned
    return names, seqs

READERS = {
    "fas":  read_fas_a2m_a3m,
    "a2m":  read_fas_a2m_a3m,
    "a3m":  read_fas_a2m_a3m,
    "sto":  read_stockholm,
    "psi":  read_psi,
    "clu":  read_clustal
}

def write_alignment(outfile, names, seqs, outformat, opts, titleline=""):
    numres = opts.numres
    if outformat == "sto":
        with open(outfile, "w") as f:
            f.write("# STOCKHOLM 1.0\n\n")
            # Minimal #=GC RF line from first real sequence
            refidx = next((i for i,nm in enumerate(names) if not re.match(r'^(ss_|aa_|sa_)', nm)), 0)
            refline = list(seqs[refidx])
            for i,c in enumerate(refline):
                if c.islower(): refline[i] = '-'
            f.write("#=GC RF {}\n".format("".join(refline)))
            for nm, s in zip(names, seqs):
                nmprint = nm[:opts.lname]
                f.write("{} {}\n".format(nmprint.ljust(opts.lname), s))
            f.write("//\n")
    elif outformat == "psi":
        # blocky output
        length = max(len(s) for s in seqs) if seqs else 0
        with open(outfile, "w") as f:
            start = 0
            while start < length:
                end = start + numres
                for nm,s in zip(names, seqs):
                    nmprint = nm[:opts.lname]
                    chunk = s[start:end]
                    f.write(f"{nmprint.ljust(opts.lname)} {chunk}\n")
                f.write("\n")
                start = end
    elif outformat == "clu":
        with open(outfile, "w") as f:
            f.write("CLUSTAL\n\n\n")
            length = len(seqs[0]) if seqs else 0
            start = 0
            while start < length:
                end = start + numres
                for nm,s in zip(names, seqs):
                    shortnm = nm[:opts.lname]
                    chunk = s[start:end]
                    f.write(f"{shortnm.ljust(opts.lname)} {chunk}\n")
                f.write("\n")
                start = end
    else:
        # fas, a2m, a3m, ufas
        with open(outfile,"w") as f:
            if outformat=="a3m" and titleline.strip():
                f.write(titleline.strip() + "\n")
            for nm,s in zip(names,seqs):
                f.write(f">{nm}\n")
                for i in range(0,len(s), numres):
                    f.write(s[i:i+numres] + "\n")

def apply_transformations(names, seqs, informat, outformat, opts, titleline):
    # If input not a2m/a3m, uppercase everything
    if informat not in ("a2m","a3m"):
        seqs = [s.upper() for s in seqs]
    # Keep only [A-Za-z0-9.\-~], map '~'->'-'
    clean = []
    for s in seqs:
        s = re.sub(r'[^A-Za-z0-9.\-~]', '', s)
        s = s.replace('~','-')
        clean.append(s)
    seqs = clean

    # If a3m input and not removing inserts (or we need match mode), fill up '.' 
    # to keep all sequences same length at each insert. 
    # (We simplify that logic or skip it if not needed.)
    # The original code does a big tokenization to align inserts. We'll do a simpler approach:
    # if needed, do the expansion:
    need_insert_align = (informat=="a3m") and (not opts.remove_inserts or opts.matchmode)
    if need_insert_align:
        seqs = expand_a3m_inserts(seqs)

    # -M first or -M int => set match states
    if opts.matchmode=="first":
        seqs = matchmode_first(names, seqs)
    elif opts.matchmode=="gaprule":
        seqs = matchmode_gaprule(seqs, opts.match_gaprule)

    # -r <int> => remove columns with too many gaps
    if opts.remove_gapped>0:
        seqs = remove_gapped_cols(seqs, opts.remove_gapped)

    # -r => remove all lowercase and '.'
    if opts.remove_inserts:
        newseqs = []
        for s in seqs:
            # remove all a-z and '.' => replicate s.tr/a-z.// in Perl
            newseqs.append("".join(ch for ch in s if not(ch.islower() or ch==".")))
        seqs = newseqs

    # remove sequences that are all gaps
    final_names, final_seqs = [], []
    for nm,s in zip(names, seqs):
        if re.search(r'[A-Za-z0-9]', s):
            final_names.append(nm)
            final_seqs.append(s)
    names, seqs = final_names, final_seqs

    # -d => limit name length
    names = [n[:opts.desclen] for n in names]

    # if outformat=a3m => remove '.'  ( =>  '-')
    # if outformat=fas/clu/sto/psi => '.' => '-'
    # if outformat=ufas => remove all gaps
    if outformat=="a3m":
        seqs = [s.replace('.','') for s in seqs]
    elif outformat in ("fas","clu","sto","psi"):
        seqs = [s.replace('.', '-') for s in seqs]
    elif outformat=="ufas":
        seqs = [s.replace('.','').replace('-','') for s in seqs]

    # -g => override gap chars
    if opts.gapchar!="default":
        def gsub(ch): return opts.gapchar if ch in ('.','-') else ch
        seqs = ["".join(gsub(c) for c in s) for s in seqs]

    # -uc / -lc
    if opts.case=="uc":
        seqs = [s.upper() for s in seqs]
    elif opts.case=="lc":
        seqs = [s.lower() for s in seqs]

    # ensure all sequences same length if outformat not in (a3m, ufas)
    if outformat not in ("a3m","ufas") and seqs:
        length0 = len(seqs[0])
        for i,s in enumerate(seqs[1:],1):
            if len(s)!=length0:
                sys.exit("ERROR: sequences differ in length after processing.")

    return names, seqs, titleline

def expand_a3m_inserts(seqs):
    # Align all sequences' inserts by adding '.' up to max needed for each position
    # Short version of the original approach
    import re
    pat = re.compile(r'([A-Z]|-|~|\d)')
    # For each seq, we do re.split to keep tokens, track max length per token index
    splitted = []
    for s in seqs:
        splitted.append(pat.split(s))  # keep delimiters in the result
    max_tokens = max(len(x) for x in splitted)
    # find max len at each token index
    col_lens = [0]*max_tokens
    for arr in splitted:
        for i,token in enumerate(arr):
            col_lens[i] = max(col_lens[i], len(token))
    # pad each token
    new_seqs = []
    for arr in splitted:
        out = []
        for i,token in enumerate(arr):
            pad_len = col_lens[i] - len(token)
            if pad_len>0: token+=('.'*pad_len)
            out.append(token)
        new_seqs.append("".join(out))
    return new_seqs

def matchmode_gaprule(seqs, threshold):
    if not seqs:
        return seqs
    length = len(seqs[0])
    gap_count = [0]*length
    nseq = len(seqs)
    for s in seqs:
        for i,ch in enumerate(s):
            if ch in ('.','-'):
                gap_count[i]+=1
    new_seqs = []
    for s in seqs:
        out = []
        for i,ch in enumerate(s):
            if gap_count[i]<0.01*threshold*nseq:
                # match => uppercase, '.' => '-'
                out.append('-' if ch=='.' else ch.upper())
            else:
                # insert => lowercase, '-' => '.'
                out.append('.' if ch=='-' else ch.lower())
        new_seqs.append("".join(out))
    return new_seqs

def matchmode_first(names, seqs):
    # pick first real (non-ss_...) sequence
    idx0 = 0
    for i,nm in enumerate(names):
        if not re.match(r'^(ss_|aa_|sa_)', nm):
            idx0 = i
            break
    ref = seqs[idx0]
    new_seqs = []
    for s in seqs:
        out = []
        for c_ref, c in zip(ref, s):
            if c_ref in ('.','-'):  # => insert
                out.append('.' if c=='-' else c.lower())
            else:
                out.append('-' if c=='.' else c.upper())
        new_seqs.append("".join(out))
    return new_seqs

def remove_gapped_cols(seqs, threshold):
    # Remove columns with >= threshold% gaps
    if not seqs: return seqs
    length = len(seqs[0])
    gap_count = [0]*length
    nseq = len(seqs)
    for s in seqs:
        for i,ch in enumerate(s):
            if ch in ('.','-'): gap_count[i]+=1
    keep = [ (gap_count[i]<0.01*threshold*nseq) for i in range(length) ]
    new_seqs = []
    for s in seqs:
        new_s = []
        for i,ch in enumerate(s):
            if keep[i]:
                new_s.append(ch)
        new_seqs.append("".join(new_s))
    return new_seqs

def guess_format(filename, is_input=False):
    ext = os.path.splitext(filename)[1].lower().lstrip('.')
    mapping = {"aln":"clu","fa":"fas","fasta":"fas","afa":"fas","afas":"fas","afasta":"fas"}
    if ext in mapping: 
        return mapping[ext]
    if ext in FORMATS_IN+FORMATS_OUT:
        return ext
    return "fas" if is_input else "fas"

def main():
    parser = argparse.ArgumentParser(
        description="Convert and reformat MSAs among fas, a2m, a3m, sto, psi, clu, ufas."
    )
    parser.add_argument("informat", type=str,
                        help="Input format or file extension guess (fas,a2m,a3m,sto,psi,clu).")
    parser.add_argument("outformat", type=str,
                        help="Output format or file extension guess (fas,a2m,a3m,sto,psi,clu,ufas).")
    parser.add_argument("infile", type=str, help="Input file or fileglob (use quotes).")
    parser.add_argument("outfile", type=str, help="Output file or extension if infile is a fileglob.")
    parser.add_argument("-v", type=int, default=2, help="Verbose level (0/1/2/3)")
    parser.add_argument("-num", action="store_true", help="Add numeric prefix to sequence names")
    parser.add_argument("-noss", action="store_true", help="Remove sequences with names starting ss_")
    parser.add_argument("-sa", action="store_true", help="Do NOT remove sequences with names starting sa_ (default removes them)")
    parser.add_argument("-M", metavar="MATCH", default="", 
                        help="Use 'first' or an integer for match-state assignment.")
    parser.add_argument("-r", nargs="?", const=True, help="Remove inserts/lowercase, or if an int is given, remove columns > X%% gaps.")
    parser.add_argument("-g", metavar="GAP", default="default", help="Set gap char (e.g. '' or '-')")
    parser.add_argument("-uc", action="store_true", help="Convert to uppercase at end")
    parser.add_argument("-lc", action="store_true", help="Convert to lowercase at end")
    parser.add_argument("-l", type=int, default=100, help="Number of residues per line in output")
    parser.add_argument("-d", type=int, default=1000, help="Max description length")
    parser.add_argument("-lname", type=int, default=32, help="Name padding length for sto/psi/clu")
    parser.add_argument("-u", action="store_true", help="Update mode: skip if outfile exists")

    args = parser.parse_args()

    # final options
    class Opts: pass
    opts = Opts()
    opts.verbose = args.v
    opts.numres = args.l
    opts.desclen = args.d
    opts.lname = args.lname
    opts.noss = args.noss
    opts.remove_sa = not args.sa   # default remove sa_ unless user says -sa
    opts.gapchar = args.g
    opts.case = "uc" if args.uc else ("lc" if args.lc else "default")
    opts.matchmode = ""
    opts.match_gaprule = 0
    opts.remove_inserts = False
    opts.remove_gapped = 0

    # interpret -M
    if args.M == "first":
        opts.matchmode = "first"
    else:
        # maybe int
        if args.M.isdigit():
            opts.matchmode = "gaprule"
            opts.match_gaprule = int(args.M)

    # interpret -r
    if args.r is not None:
        if args.r is True:
            # user typed just -r
            opts.remove_inserts = True
        else:
            # user gave an int
            try:
                val = int(args.r)
                opts.remove_gapped = val
            except ValueError:
                opts.remove_inserts = True

    # parse in/out format if not recognized
    informat = args.informat.lower()
    outformat = args.outformat.lower()
    if informat not in FORMATS_IN:
        # guess from infile
        informat = guess_format(args.infile, is_input=True)
    if outformat not in FORMATS_OUT:
        outformat = guess_format(args.outfile, is_input=False)

    # handle update mode
    if args.u and os.path.exists(args.outfile) and not ("*" in args.infile or "?" in args.infile):
        if opts.verbose>=1:
            print(f"Skipping {args.infile} => {args.outfile} because it already exists and -u used.")
        return

    # if glob in infile => process multiple
    if any(sym in args.infile for sym in ("*","?")):
        in_files = glob.glob(args.infile)
        if opts.verbose>=1:
            print(f"{len(in_files)} files to process.")
        # if outfile is extension
        if not args.outfile.startswith("."):
            sys.exit("ERROR: When using fileglob for infile, outfile must be an extension, e.g. .a3m")
        out_ext = args.outfile
        for f_in in in_files:
            base = os.path.splitext(f_in)[0]
            f_out = base + out_ext
            if args.u and os.path.exists(f_out):
                if opts.verbose>=1:
                    print(f"Skipping {f_in} => {f_out}, already exists.")
                continue
            run_conversion(f_in, f_out, informat, outformat, opts)
    else:
        run_conversion(args.infile, args.outfile, informat, outformat, opts)

def run_conversion(infile, outfile, informat, outformat, opts):
    if opts.verbose>=3:
        print(f"Reading {infile} as {informat} -> {outformat}.")
    reader = READERS.get(informat, None)
    if not reader:
        sys.exit(f"Unknown input format {informat}")
    names, seqs = reader(infile, opts.remove_sa, opts.noss)
    if not names:
        sys.exit(f"No sequences read from {infile}.")

    titleline = ""
    # We used read_fas_a2m_a3m for 'fas','a2m','a3m', which might capture a leading # line if needed.
    # For brevity, we skip storing that here. If you want "titleline", parse it from the file.
    names, seqs, titleline = apply_transformations(names, seqs, informat, outformat, opts, titleline)

    # -num => rename sequences with #n
    if opts.num:
        counter = 2
        for i,nm in enumerate(names):
            if re.match(r'^(ss_|aa_|sa_|ss_conf)', nm):
                continue
            base = nm.split('#')[0]
            base = base[:25]
            names[i] = f"{base}#{counter}"
            counter+=1

    if opts.verbose>=2:
        print(f"Writing {len(names)} seqs to {outfile} ({outformat}).")
    write_alignment(outfile, names, seqs, outformat, opts, titleline)

if __name__=="__main__":
    main()
