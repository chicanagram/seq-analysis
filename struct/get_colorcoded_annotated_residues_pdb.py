import numpy as np
import pandas as pd
import yasara
from variables import address_dict, subfolders
from variables import aa_taylor_colorcode_yasara

data_folder = address_dict['ECOHARVEST']
data_subfolder = 'lipases'
struct_fname = 'Docked_RML_OleicAcid_preOpt.sce' # 'RML_4tgl_open.pdb'
if struct_fname.find('.sce')>-1: struct_fmt = 'sce'
elif struct_fname.find('.pdb')>-1: struct_fmt = 'pdb'
struct_fpath = data_folder + subfolders[struct_fmt] + data_subfolder + '/' + struct_fname
output_sce_fsuffix = '_sift'# '_patentmuts'
colorcode = (0, 180) # 'green' # 'ResType'
get_res_from_csv = {'fpath': f'{data_folder}{subfolders["conservation_analysis"]}{data_subfolder}/RML_blastp_nr_E1e-05_mafft_sift_selected.csv',
                    'pos': 'RealPos',
                    'color_by':'entropy',
                    'res_offset': 94
                    }
res_to_colorcode = list(range(5,271)) # [88, 150, 156, 269, 181, 220]
colors = None # ['green']*2 ['blue']*4
colorby_vals = None
base_color = 'grey'
base_style = 'ribbon' # 'surface'
label_res = 'RESNUM' # None # 'RESNUM RESNAME1' #
label_thres = 0.1 # fraction of 'top' residues to label / show atoms for

# get res to colorcode from CSV
res_offset = 0
if get_res_from_csv is not None:
    pos_col = get_res_from_csv['pos']
    color_by_col = get_res_from_csv['color_by']
    res_offset = get_res_from_csv['res_offset']
    df = pd.read_csv(get_res_from_csv['fpath'])
    df_filt = df.loc[df[pos_col].isin([res+res_offset for res in res_to_colorcode])]
    res_to_colorcode = df_filt[pos_col].tolist()
    colorby_vals = df_filt[color_by_col].tolist()
    print('res_to_colorcode:', res_to_colorcode)
    print('colorby_vals:', colorby_vals)


# get colors
show_res_threshold = None
if colors is None:
    numres = len(res_to_colorcode)
    if isinstance(colorcode, str):
        colors = [colorcode]*numres
    elif isinstance(colorcode, tuple) and colorby_vals is not None:
        colorcode_int = colorcode[1]-colorcode[0]
        val_min = min(colorby_vals)
        val_max = max(colorby_vals)
        val_int = val_max - val_min
        colors =  [int(colorcode[0] + (val-val_min)/val_int*colorcode_int) for val in colorby_vals]
        # get show_res_threshold
        colorby_vals_sorted = np.sort(np.array(colorby_vals))
        thres_idx = int(label_thres * len(colorby_vals))
        thres_val = colorby_vals_sorted[thres_idx]
        print('thres index:', thres_idx, '/', len(colorby_vals), '; threshold:', thres_val)


# initialize yasara
yasara.info.mode = 'txt'
yasara.Console('Off')
yasara.Clear()
yasara.DelWater()

# load PDB
if struct_fmt=='pdb':
    yasara.LoadPDB(struct_fpath)
elif struct_fmt=='sce':
    yasara.LoadSce(struct_fpath)

# set style and color of backbone
yasara.ColorRes('protein', base_color)
yasara.Style(backbone=base_style)

# color specific residues
for res, val, color in zip(res_to_colorcode, colorby_vals, colors):
    res_pdb = res - res_offset
    print(f'{res}({res_pdb})', end=' ')
    # get color to annotate
    if color=='Restype':
        yasara.ColorRes(res_pdb, 'ResType')
    else:
        yasara.ColorRes(res_pdb, color)
    if label_res is not None:
        if val <= thres_val:
            print(res_pdb, val, thres_val)
            yasara.ShowRes(res_pdb)
            yasara.StickRes(res_pdb)
            yasara.LabelRes(res_pdb, label_res)

# color ligand by element
yasara.ColorRes('not protein', 'Element')
# save scene
yasara.SaveSce(struct_fpath.replace(struct_fmt, 'sce').replace('.sce', output_sce_fsuffix + '.sce'))