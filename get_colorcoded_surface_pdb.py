import pandas as pd
import yasara
from variables import address_dict, subfolders
from variables import aa_taylor_colorcode_yasara

data_folder = address_dict['ECOHARVEST']
data_subfolder = 'lipases'
pdb_fname = 'RML_4tgl_open.pdb'
pdb_fpath = data_folder + subfolders['pdb'] + data_subfolder + '/' + pdb_fname
res_to_colorcode = {'green':[181, 220], 'blue':[88,150,156,269]} # {'polarity': [181, 220]}

base_color = 'grey'
base_style = 'ribbon' # 'surface'
label_res = 'RESNUM RESNAME1' # None

# initialize yasara
yasara.info.mode = 'txt'
yasara.Console('Off')
yasara.Clear()
yasara.DelWater()

# load PDB
yasara.LoadPDB(pdb_fpath)

# set simulation cell around protein

# show surface

#

yasara.SaveSce(pdb_fpath.replace('pdb', 'sce').replace('.sce', '_patentmuts.sce'))