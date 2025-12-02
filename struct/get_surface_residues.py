from variables import address_dict, subfolders
import numpy as np
import pandas as pd

def yasara_get_surface_residues(pdb_fpath, surft=2.55, distt=12, minepisize=13):

    import yasara

    surf_types = ['Accessible', 'VdW', 'Molecular']
    res_dist_to_surf = {f'dist2surf_{surf_type}':[] for surf_type in surf_types}
    res2res_dist = {}
    res_surfpatches = {}

    # load pdb
    yasara.LoadPDB(pdb_fpath)
    # Add protein residues to the environment
    yasara.AddEnvRes('Protein')
    numres = yasara.ListRes('Protein')
    # turn off console
    yasara.Console('off')

    # get all residue distances to surface
    for i in range(1, len(numres) + 1):
        for surf_type in surf_types:
            res_dist_to_surf[f'dist2surf_{surf_type}'].append(round(yasara.SurfDisRes(f'Res {i}', surf_type)[0],3))

    # Iterate through each residue and 1) check if it's a surface residue; 2) look for surrounding surface residues
    for i in range(1, len(numres) + 1):
        surfy = res_dist_to_surf['dist2surf_Accessible'][i-1]
        # color all residues grey
        yasara.HideRes('all')
        yasara.ColorRes('all', 'gray')
        # if less than threshold distance (surft), color residue red and add to newlist
        if surfy < surft:
            yasara.ShowRes(i)
            yasara.BallRes(i)
            yasara.ColorRes(i, 'red')
            newlist = [i]

            # iterate over residues and identify those that are in the vicinity of Residue i
            for j in range(1, len(numres) + 1):
                if i != j:
                    # get tuple of residue pair
                    res_pair = [i, j]
                    res_pair.sort()
                    res_pair = tuple(res_pair)
                    if res_pair not in res2res_dist:
                        disty = yasara.Distance(f'res {i}', f'res {j}')[0]
                        res2res_dist[res_pair] = disty
                    else:
                        disty = res2res_dist[res_pair]
                    # if other residue is within 12 angstroms of original residue, check if it's also solvent accessible
                    if disty < distt:
                        surfz = res_dist_to_surf['dist2surf_Accessible'][j-1]
                        if surfz < surft:
                            yasara.ShowRes(j)
                            yasara.BallRes(j)
                            yasara.ColorRes(j, 'yellow')
                            newlist.append(j)
            # get size of surface patch and log patch if it's above a threshold size
            episize = len(newlist)
            if episize >= minepisize:
                res_surfpatches[i] = ','.join([str(res) for res in newlist])
                print(f"Residue: {i}, Epitope size: {episize}, Residues: {newlist}")
    res_dist_to_surf = pd.DataFrame(res_dist_to_surf)
    res_dist_to_surf['RealPos'] = list(range(1,len(numres)+1))
    res_surfpatches = pd.DataFrame.from_dict(res_surfpatches, orient='index', dtype=object, columns=['surface_residue_patch'])
    res_surfpatches['RealPos'] = list(res_surfpatches.index)
    res = res_dist_to_surf.merge(res_surfpatches, on='RealPos', how='left')
    return res


if __name__ == '__main__':
    data_folder = address_dict['PIPS2']
    data_fbase = 'ET096'
    pdb_dir = data_folder + subfolders['pdb']
    pdb_fname = 'ET096.pdb'
    log_fname =  pdb_fname.replace('.pdb', '_surfRes.csv')
    pdb_fpath = f'{pdb_dir}{data_fbase}/{pdb_fname}'
    log_fpath = f'{pdb_dir}{data_fbase}/{log_fname}'
    res = yasara_get_surface_residues(pdb_fpath, surft=2.55, distt=12, minepisize=13)
    res.to_csv(log_fpath)
    print(res)


