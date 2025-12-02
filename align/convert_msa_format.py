import os
import shutil
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
    input_temp_fpath = msa_dir+f'TEMP.{input_fmt}'
    output_temp_fpath = msa_dir+f'TEMP.{output_fmt}'
    input_fpath = msa_dir+input_file
    output_file = input_file.split('.')[0] + '.' + output_fmt
    output_fpath = msa_dir+output_file
    # use Perl script
    if from_perl_or_python=='perl':
        import subprocess
        cmd = [
            "perl",
            "msa/reformat_msa.pl",
            input_fmt,
            output_fmt,
            "-i", input_temp_fpath,
            "-o", output_temp_fpath,
        ]
        try:
            shutil.copyfile(input_fpath, input_temp_fpath)
            print(f'Copied input file: {input_fpath} > {input_temp_fpath}')
            subprocess.call(cmd, shell=False)
            os.rename(output_temp_fpath, output_fpath)
            print(f'Renamed output file: {output_temp_fpath} > {output_fpath}')
            os.remove(input_temp_fpath)
            print(f'Converted MSA: {input_file} > {output_file}')
        except Exception as e:
            print(e)
    # use Python script
    elif from_perl_or_python == 'python':
        from msa.reformat_msa import run_conversion
        run_conversion(input_fpath, output_fpath, input_fmt, output_fmt, options)
        print(f'Converted MSA: {input_file} > {output_file}')
    print()

if __name__=='__main__':
    os.chdir('../')
    data_folder = address_dict['PON-Sol2'] # address_dict['ECOHARVEST'] #
    msa_dir = data_folder + subfolders['msa']
    from_perl_or_python = 'perl'
    input_fmt = 'a3m' # 'fas'
    output_fmt = 'fas' # 'a3m'
    input_file_list = [f for f in os.listdir(msa_dir) if f.find('.'+input_fmt)>-1]
    # input_file = '0_RLBP1_HUMAN__blastp_uniprot_trembl_E1e-03_mafft.a3m' # 'CALA_blastp_uniprot_sprot_E1e-03_mafft.fasta'
    for input_file in input_file_list:
        print(input_file)
        convert_msa_format(input_file, input_fmt=input_fmt, output_fmt=output_fmt, msa_dir=msa_dir, from_perl_or_python=from_perl_or_python, options=Options())