from variables import address_dict, subfolders

class Options:
    """
    Holds configuration for the run_conversion() function,
    including all flags that control reading, writing, and
    transformations of the MSA.
    """
    def __init__(
        self,
        verbose=2,
        remove_sa=True,
        noss=False,
        num=False,
        remove_inserts=False,
        remove_gapped=0,
        matchmode="",
        match_gaprule=0,
        case="default",
        gapchar="default",
        desclen=1000,
        lname=32,
        numres=100,
    ):
        """
        :param verbose: Verbosity level (int), e.g. 0=quiet, 2=normal, 3=debug
        :param remove_sa: If True, remove 'sa_' sequences
        :param noss: If True, remove 'ss_' sequences
        :param num: If True, add numeric prefixes to sequence names
        :param remove_inserts: If True, remove all lowercase inserts
        :param remove_gapped: If >0, remove columns with more than X% gaps
        :param matchmode: "" (none), "first", or "gaprule"
        :param match_gaprule: When matchmode="gaprule", integer for % gap threshold
        :param case: "uc"=uppercase final, "lc"=lowercase final, "default"=no change
        :param gapchar: Gap character override, e.g. "-", "", or "default" (no override)
        :param desclen: Maximum number of characters for sequence name
        :param lname: Used to pad/align names in certain formats (sto, psi, clu)
        :param numres: Number of residues per line in the output
        """
        self.verbose = verbose
        self.remove_sa = remove_sa
        self.noss = noss
        self.num = num
        self.remove_inserts = remove_inserts
        self.remove_gapped = remove_gapped
        self.matchmode = matchmode
        self.match_gaprule = match_gaprule
        self.case = case
        self.gapchar = gapchar
        self.desclen = desclen
        self.lname = lname
        self.numres = numres

def convert_msa_format(input_file, input_fmt='fas', output_fmt='a3m', msa_dir='', from_perl_or_python='perl', options=None):
    input_fpath = msa_dir+input_file
    output_file = input_file.split('.')[0] + '.' + output_fmt
    output_fpath = msa_dir+output_file
    # use Perl script
    if from_perl_or_python=='perl':
        import subprocess
        cmd = [
            "perl",
            "msa/reformat_msa.pl",
            "fas",
            output_fmt,
            "-i", input_fpath,
            "-o", output_fpath,
        ]
        try:
            subprocess.call(cmd, shell=False)
        except Exception as e:
            print(e)
    # use Python script
    elif from_perl_or_python == 'python':
        from msa.reformat_msa import run_conversion
        run_conversion(input_fpath, output_fpath, input_fmt, output_fmt, options)
    print(f'Converted MSA: {input_file} > {output_file}')

if __name__=='__main__':
    data_folder = address_dict['ECOHARVEST'] # address_dict['PON-Sol2']
    msa_dir = data_folder + subfolders['msa']
    input_file = 'CALA_blastp_uniprot_sprot_E1e-03_mafft.fasta'
    from_perl_or_python = 'perl'
    convert_msa_format(input_file, input_fmt='fas', output_fmt='a3m', msa_dir=msa_dir, from_perl_or_python=from_perl_or_python, options=Options())