import os
import numpy as np
import pandas as pd
import yasara
from variables import address_dict, subfolders
from variables import aa_taylor_colorcode_yasara

if __name__=='__main__':
    os.chdir('../')
    data_folder = address_dict['PIPS2']# address_dict['ECOHARVEST']
    data_subfolder = '' # 'lipases'
    struct_fname = 'ET096.pdb' # 'Docked_RML_OleicAcid_preOpt.sce' # 'RML_4tgl_open.pdb'
    if struct_fname.find('.sce')>-1: struct_fmt = 'sce'
    elif struct_fname.find('.pdb')>-1: struct_fmt = 'pdb'
    struct_fpath = data_folder + subfolders[struct_fmt] + data_subfolder + '/' + struct_fname
    output_sce_fsuffix = '_zeroshot' # '_sift'# '_patentmuts'
    colorcode = (0, 180) # 'green' # 'ResType'
    # get_res_from_csv = {'fpath': f'{data_folder}{subfolders["conservation_analysis"]}{data_subfolder}/RML_blastp_nr_E1e-05_mafft_sift_selected.csv',
    #                     'pos': 'RealPos',
    #                     'color_by':'entropy',
    #                     'res_offset': 94
    #                     }
    get_res_from_csv = None
    res_to_colorcode = [82, 163, 87, 195, 171, 29, 174, 182, 78, 223, 230, 70, 79, 20, 64, 121, 169, 83, 190, 165, 143, 57, 236, 102, 67, 120, 239, 212, 207, 234, 214, 84, 217, 238, 119, 138, 208, 150, 38, 86, 32, 220, 114, 90] + \
        [93, 61, 200, 41, 88, 75, 236, 222, 143, 140, 76, 212, 174, 65] + [174, 143, 236, 212]
    # res_to_colorcode = list(range(5,271)) # [88, 150, 156, 269, 181, 220]
    colors = ['green']*44 + [150]*14 + ['red']*4 # None # ['green']*2 ['blue']*4
    colorby_vals = None
    base_color = 'grey'
    base_style = 'ribbon' # 'surface'
    label_res = 'RESNUM RESNAME1' # 'RESNUM' # None #
    label_thres = 0.1 # fraction of 'top' residues to label / show atoms for
    show_residue = True

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

    if colorby_vals is None:
        colorby_vals = [None]*len(colors)

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
        if show_residue:
            yasara.ShowRes(res_pdb)
        if label_res is not None:
            if val is not None and val <= thres_val:
                print(res_pdb, val, thres_val)
                yasara.ShowRes(res_pdb)
                yasara.StickRes(res_pdb)
                yasara.LabelRes(res_pdb, label_res)

    for res in [143,174,212,236]:
        yasara.LabelRes(res,'RESNUM RESNAME1', height=1.2)
    # color ligand by element
    yasara.ColorRes('not protein', 'Element')
    yasara.DelMol('not A')
    # save scene
    yasara.SaveSce(struct_fpath.replace(struct_fmt, 'sce').replace('.sce', output_sce_fsuffix + '.sce'))