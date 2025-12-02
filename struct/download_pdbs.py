import os
from Bio.PDB import PDBList
from variables import address_dict, subfolders

def download_pdb(pdb_list, pdb_dir):
    # Initialize downloader
    pdbl = PDBList()

    # Download a single PDB entry by 4-character ID
    for pdb in pdb_list:
        if pdb.find('|')>-1:
            pdb_id, pdb_chain = pdb.split('|')
        else:
            pdb_id = pdb
        pdb_fpath = pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format="pdb")
        # rename file
        new_pdb_fpath = os.path.join(pdb_dir, f"{pdb_id.lower()}.pdb")
        os.rename(pdb_fpath, new_pdb_fpath)
        print("Downloaded:", new_pdb_fpath)

if __name__=='__main__':
    os.chdir('../')
    data_folder = address_dict['SoluProtMut'] # address_dict['PON-Sol2']
    ponsol2_pdb_list = ['2cy8|A', '2f2v|A', '1ax8|A', '1i1b|A', '1jy7|S', '3g4o|B', '1lax|C', '3mzs|D', '3qqd|A', '3pg0|A', '1q31|B', '1xww|A', '2d26|A']
    soluprotmut_pdb_list = ['4zfv|A', '1itg|A', '5dmq|A', '1c6v|A', '2l7b|A', '1s5e|A', '1kl9|A', '2hzq|A', '6fwt|A', '2w9s|A', '1ax8|A', '3upa|A', '3zhk|A', '4hqa|A', '1ifr|A', '2xwr|A', '1cv2|A', '2bne|A', '1lni|A', '3ne4|A', '5r8q|A', '1yue|A', '1ccz|A', '1anf|A', '4rov|A', '1z0q|A', '1lvm', '1col|A', '1bff|A', '1cw3|A', '3hx3|A', '2yvq', '2yvq|A', '1oki|A', '2wl1', '2wl1|A', '1xse|A', '3ebb', '3ibd|A', '5u9a|A', '1xu7|A', '1hvy', '1j42|A', '1xqi|A', '1df8|A', '1j1d', '1j1d|A', '4tsx|A', '1xww|A', '4v06|A', '3hrn|A', '2cy8|A', '4jqv', '4jqv|A', '2fjt|A', '1pre|A', '3mzs|A', '3cki', '3k77|A', '1azv|A', '1kg2|A', '4ah6|A', '3vub|A', '4ye6|A', '2x6g|A', '1b55', '1pwt|A', '1dwo|A', '5afh|A', '1h4a|A', '4zlp', '3g7v|A', '5f59', '1ayi|A', '1f0y|A', '3gxp|A', '1ege|A']
    pdb_dir = data_folder + subfolders['pdb']
    print('PDB list:', len(soluprotmut_pdb_list), soluprotmut_pdb_list)
    download_pdb(soluprotmut_pdb_list, pdb_dir)