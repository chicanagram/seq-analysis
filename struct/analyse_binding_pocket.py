from __future__ import annotations
import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib.pyplot as plt
from typing import Dict, Literal
import yasara
from variables import address_dict, subfolders

def showres_bindingpocket_struct(pdb_fpath, binding_pocket_residues):
    yasara.LoadPdb(pdb_fpath)
    for resnum in binding_pocket_residues:
        yasara.HideRes('protein')
        yasara.ShowRes(f'protein and Res {resnum}')
        yasara.LabelRes(f'protein and Res {resnum}', 'RESNUM')
    # save as scene
    yasara.SaveSce(pdb_fpath.replace('pdb', 'sce'))

# calculate distances between residue and binding pocket center
def get_distances_residues_bindingpocket_centroid(df_bindingpocket, centroid, get_residue_min_distance=False):
    df_bindingpocket['distance_to_centroid'] = np.linalg.norm(df_bindingpocket[['x','y','z']] - centroid, axis=1)
    resnum_list = df_bindingpocket['res_num'].drop_duplicates().tolist()
    if get_residue_min_distance:
        df_bindingpocket.loc[:, 'min_distance_to_centroid'] = None
        for resnum in resnum_list:
            min_dist_res = df_bindingpocket.loc[df_bindingpocket['res_num']==resnum, 'distance_to_centroid'].min()
            df_bindingpocket.loc[(df_bindingpocket['res_num']==resnum) & (df_bindingpocket['atom_name']=='CA'), 'min_distance_to_centroid'] = min_dist_res

    # get stats
    mean_dist_to_centroid = round(df_bindingpocket['distance_to_centroid'].mean(),4)
    mean_min_dist_to_centroid = round(df_bindingpocket['min_distance_to_centroid'].mean(), 4) if get_residue_min_distance else None
    print('Mean distance to binding pocket centroid:', mean_dist_to_centroid)
    print('Mean MIN distance to binding pocket centroid:', mean_min_dist_to_centroid)

    return df_bindingpocket, mean_dist_to_centroid, mean_min_dist_to_centroid

def compute_index_summary(
    df_pocket: pd.DataFrame,
    value_col: Literal["kd_hydro", "hw_polarity"],
    weight_mode: Literal["mean", "weighted"] = "mean",
    dist_col: str = "dist_to_centroid",
    eps: float = 1e-6,
) -> float:
    """
    Compute either an unweighted mean or a distance-weighted mean for a polarity/hydropathy index.

    Parameters
    ----------
    df_pocket : pd.DataFrame
        Must contain `value_col` and `dist_col`.
    value_col : {"kd_hydro", "hw_polarity"}
        Column containing numeric index values per residue.
    weight_mode : {"mean", "weighted"}
        "mean" -> simple mean of value_col
        "weighted" -> weights = 1 / (dist_to_centroid + eps)
    dist_col : str
        Column with distance-to-centroid values (smaller means closer to centroid).
    eps : float
        Small constant to avoid divide-by-zero.

    Returns
    -------
    float
        Summary value (NaNs are ignored).
    """
    if value_col not in df_pocket.columns:
        raise KeyError(f"Missing required column: {value_col}")
    if weight_mode == "weighted" and dist_col not in df_pocket.columns:
        raise KeyError(f"Missing required column: {dist_col}")

    vals = pd.to_numeric(df_pocket[value_col], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(vals)

    if not mask.any():
        return float("nan")

    if weight_mode == "mean":
        return float(np.nanmean(vals))

    # weighted
    d = pd.to_numeric(df_pocket.loc[:, dist_col], errors="coerce").to_numpy(dtype=float)
    w = 1.0 / (d + eps)
    w[~np.isfinite(w)] = np.nan

    # Only keep rows where both value and weight are finite
    m = mask & np.isfinite(w)
    if not m.any():
        return float("nan")

    return float(np.average(vals[m], weights=w[m]))

def pocket_polarity_report(
    df_pocket: pd.DataFrame,
    *,
    aa_col: str = "res",
    aa_polarity_col: str = "aa_polarity",
    dist_col: str = "dist_to_centroid",
    kd_col: str = "kd_hydro",
    hw_col: str = "hw_polarity",
    eps: float = 1e-6,
) -> Dict[str, float]:
    """
    Compute a compact "pocket polarity report" with 6 metrics:
      1) KD_mean
      2) KD_weighted (1/(dist+eps))
      3) HW_mean
      4) HW_weighted (1/(dist+eps))
      5) charged_fraction
      6) polar_fraction

    Assumptions
    -----------
    df_pocket already contains:
      - 'res' (1-letter AA)
      - 'aa_polarity' categorical (e.g. np/p~/p-/p+)
      - 'kd_hydro' numeric
      - 'hw_polarity' numeric
      - 'dist_to_centroid' numeric
      - 'res_num' (not used here but expected to exist upstream)

    Returns
    -------
    Dict[str, float]
    """
    # --- KD / HW summaries (re-using the same function) ---
    kd_mean = compute_index_summary(
        df_pocket, value_col=kd_col, weight_mode="mean", dist_col=dist_col, eps=eps
    )
    kd_weighted = compute_index_summary(
        df_pocket, value_col=kd_col, weight_mode="weighted", dist_col=dist_col, eps=eps
    )
    hw_mean = compute_index_summary(
        df_pocket, value_col=hw_col, weight_mode="mean", dist_col=dist_col, eps=eps
    )
    hw_weighted = compute_index_summary(
        df_pocket, value_col=hw_col, weight_mode="weighted", dist_col=dist_col, eps=eps
    )

    # --- Composition metrics ---
    if aa_col not in df_pocket.columns:
        raise KeyError(f"Missing required column: {aa_col}")
    if aa_polarity_col not in df_pocket.columns:
        raise KeyError(f"Missing required column: {aa_polarity_col}")

    aa = df_pocket[aa_col].astype(str).str.strip().str.upper()
    cat = df_pocket[aa_polarity_col].astype(str).str.strip()

    # Keep only standard 20 AA rows (optional but helps avoid weird residues)
    std_mask = aa.isin(list("ACDEFGHIKLMNPQRSTVWY"))
    cat = cat[std_mask]

    n = int(cat.notna().sum())
    if n == 0:
        charged_fraction = float("nan")
        polar_fraction = float("nan")
    else:
        charged_fraction = float(cat.isin(["p-", "p+"]).sum() / n)
        polar_fraction = float(cat.isin(["p~", "p-", "p+"]).sum() / n)

    return {
        "kd_mean": kd_mean,
        "kd_weighted": kd_weighted,
        "hw_mean": hw_mean,
        "hw_weighted": hw_weighted,
        "charged_fraction": charged_fraction,
        "polar_fraction": polar_fraction,
    }

def pocket_sterics_report(
    df_pocket: pd.DataFrame,
    *,
    dist_col: str = "distance_to_centroid",
    aa_col: str = "res",
    vol_col: str = "aa_vol",
    eps: float = 1e-6,
) -> Dict[str, float]:
    """
    Compute a compact 6-metric sterics report for a binding pocket.

    Required columns (defaults):
      - distance_to_centroid (float)
      - res (1-letter AA)
      - aa_vol (float; side-chain volume in Å^3 or other consistent units)

    Metrics returned
    ----------------
    1) mean_volume
    2) weighted_mean_volume            (weights = 1/(distance_to_centroid + eps))
    3) volume_variance                 (unweighted variance)
    4) small_residue_frac               ( fraction of {G, A, S})
    5) small_residue_frac_weighted     (weighted fraction of {G, A, S})
    6) bulky_residue_frac               (fraction of {G, A, S})
    7) bulky_residue_frac_weighted     (fraction of {F, Y, W, R, K, L, I, M})

    Notes
    -----
    - NaNs are ignored where possible.
    - Weighted metrics use only rows where both distance and volume are finite.
    """
    for c in (dist_col, aa_col, vol_col):
        if c not in df_pocket.columns:
            raise KeyError(f"Missing required column: {c}")

    aa = df_pocket[aa_col].astype(str).str.strip().str.upper()

    dist = pd.to_numeric(df_pocket[dist_col], errors="coerce").to_numpy(dtype=float)
    vol = pd.to_numeric(df_pocket[vol_col], errors="coerce").to_numpy(dtype=float)

    # ---- Unweighted stats (volume only) ----
    if np.isfinite(vol).any():
        mean_volume = float(np.nanmean(vol))
        volume_variance = float(np.nanvar(vol, ddof=1)) if np.isfinite(vol).sum() > 1 else float("nan")
    else:
        mean_volume = float("nan")
        volume_variance = float("nan")

    # ---- Weighted stats ----
    w = 1.0 / (dist + eps)
    w[~np.isfinite(w)] = np.nan

    valid_w = np.isfinite(w) & np.isfinite(vol)
    if valid_w.any():
        weighted_mean_volume = float(np.average(vol[valid_w], weights=w[valid_w]))
        crowding_score = float(np.nansum(vol[valid_w] / (dist[valid_w] + eps)))
    else:
        weighted_mean_volume = float("nan")
        crowding_score = float("nan")

    # ---- composition fractions (unweighted & weighted) ----
    small = {"G", "A", "S"}
    bulky = {"F", "Y", "W", "R", "K", "L", "I", "M"}

    aa_arr = aa.to_numpy()
    valid_comp = np.isfinite(w) & (aa_arr != "")  # use weights even if vol missing

    if valid_comp.any():

        small_residue_frac = float(
            np.average(np.isin(aa_arr[valid_comp], list(small)).astype(float))
        )
        small_residue_frac_weighted = float(
            np.average(np.isin(aa_arr[valid_comp], list(small)).astype(float), weights=w[valid_comp])
        )
        bulky_residue_frac = float(
            np.average(np.isin(aa_arr[valid_comp], list(bulky)).astype(float))
        )
        bulky_residue_frac_weighted = float(
            np.average(np.isin(aa_arr[valid_comp], list(bulky)).astype(float), weights=w[valid_comp])
        )
    else:
        small_residue_frac = float("nan")
        bulky_residue_frac = float("nan")
        small_residue_frac_weighted = float("nan")
        bulky_residue_frac_weighted = float("nan")

    return {
        "mean_volume": mean_volume,
        "weighted_mean_volume": weighted_mean_volume,
        "volume_variance": volume_variance,
        "small_residue_frac": small_residue_frac,
        "small_residue_frac_weighted": small_residue_frac_weighted,
        "bulky_residue_frac": bulky_residue_frac,
        "bulky_residue_frac_weighted": bulky_residue_frac_weighted,
    }


class PocketAnalysis:

    def __init__(
            self,
            pdb_dir,
            struct_csv_dir,
    ):
        self.pdb_dir = pdb_dir
        self.struct_csv_dir = struct_csv_dir

    def plot_pocket_properties(self, bindingpocket_analysis):
        fig, ax = plt.subplots(1, 3, figsize=(12, 4))
        ax[0].scatter(bindingpocket_analysis['mean_min_dist_to_centroid'],
                      bindingpocket_analysis['mean_dist_to_centroid'])
        ax[0].set_xlabel('mean_min_dist_to_centroid')
        ax[0].set_ylabel('mean_dist_to_centroid')
        ax[1].scatter(bindingpocket_analysis['mean_dist_to_centroid'],
                      bindingpocket_analysis['mean_dist_backbone_to_centroid'])
        ax[1].set_xlabel('mean_dist_to_centroid')
        ax[1].set_ylabel('mean_dist_backbone_to_centroid')
        ax[2].scatter(bindingpocket_analysis['mean_dist_backbone_to_centroid'],
                      bindingpocket_analysis['mean_min_dist_to_centroid'])
        ax[2].set_xlabel('mean_dist_backbone_to_centroid')
        ax[2].set_ylabel('mean_min_dist_to_centroid')
        plt.show()

    def __call__(
            self,
            struct_name_list,
            pocket_residues_dict,
            protein_molname='A',
            plot_properties=False,
    ):
        # initialize dict to store binding pocket analyses
        bindingpocket_analysis = []
        df_bindingpocket_backbone_dict = {}
        df_bindingpocket_dict = {}

        for struct_name in struct_name_list:
            csv_fname = struct_name + '.csv'
            csv_fpath = self.struct_csv_dir + csv_fname
            struct_coords = pd.read_csv(csv_fpath)
            struct_backbone_coords = pd.read_csv(csv_fpath.replace('.csv', '_backbone.csv'))

            # get binding pocket residues
            binding_pocket_residues = pocket_residues_dict[struct_name]
            num_res_binding_pocket_ali = len(binding_pocket_residues)
            print(f'[{struct_name}] Binding pocket residues  ({num_res_binding_pocket_ali}): {binding_pocket_residues}')

            # get binding pocket residue df
            df_bindingpocket = struct_coords[struct_coords['res_num'].isin(binding_pocket_residues)].copy()
            df_backbone_bindingpocket = struct_backbone_coords[struct_backbone_coords['res_num'].isin(binding_pocket_residues)].copy()

            # get binding pocket centroid and other key atoms
            centroid = df_bindingpocket[['x', 'y', 'z']].mean(axis=0).to_numpy()
            print('Centroid:', centroid)

            # get distances to binding pocket centroid
            print('--- All Residue Atoms ---')
            df_bindingpocket, mean_dist_to_centroid, mean_min_dist_to_centroid = get_distances_residues_bindingpocket_centroid(df_bindingpocket, centroid, get_residue_min_distance=True)
            print('--- Backbone Only ---')
            df_backbone_bindingpocket, mean_backbone_dist_to_centroid, _ = get_distances_residues_bindingpocket_centroid(df_backbone_bindingpocket, centroid, get_residue_min_distance=False)
            df_backbone_bindingpocket = df_backbone_bindingpocket.rename(columns={'distance_to_centroid': 'distance_to_centroid_CA'})
            min_dist_by_res = df_bindingpocket[['res_num', 'distance_to_centroid', 'min_distance_to_centroid']].dropna(how='any')
            df_backbone_bindingpocket = df_backbone_bindingpocket.merge(min_dist_by_res, on='res_num', how='left')

            # get pocket polarity report
            pocket_polarity = pocket_polarity_report(
                df_backbone_bindingpocket,
                aa_col="res",
                aa_polarity_col="aa_polarity",
                dist_col="distance_to_centroid",
                kd_col="kd_hydro",
                hw_col="hw_polarity",
            )

            # get volume estimate
            pocket_sterics = pocket_sterics_report(
                df_backbone_bindingpocket,
                dist_col="distance_to_centroid",
                aa_col="res",
                vol_col="aa_vol",
            )

            # update analysis for this struct
            struct_analysis = {
                'struct_name': struct_name,
                'mean_min_dist_to_centroid': mean_min_dist_to_centroid,
                'mean_dist_to_centroid': mean_dist_to_centroid,
                'mean_dist_backbone_to_centroid': mean_backbone_dist_to_centroid,
            }
            struct_analysis.update(pocket_sterics)
            struct_analysis.update(pocket_polarity)
            bindingpocket_analysis.append(struct_analysis)
            df_bindingpocket_backbone_dict[struct_name] = df_backbone_bindingpocket
            df_bindingpocket_dict[struct_name] = df_bindingpocket
            print()

        bindingpocket_analysis = pd.DataFrame(bindingpocket_analysis).round(3)
        bindingpocket_analysis.to_csv(self.pdb_dir + 'bindingpocket_analysis.csv')

        # plot
        if plot_properties:
            self.plot_pocket_properties(bindingpocket_analysis)

        return bindingpocket_analysis, df_bindingpocket_dict, df_bindingpocket_backbone_dict



if __name__ == "__main__":
    os.chdir('../')

    # ---- user input ----
    data_folder = address_dict['PIPS2']
    data_subfolder = 'UPOs_peroxygenation_analysis/' # 'CARs' # 'sidestream_cocktail' #
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder
    struct_csv_dir = data_folder + subfolders['pdb'] + data_subfolder + 'structure_csv/'
    residues_near_ligand_fpath = pdb_dir + 'residues_near_ligand.csv'
    protein_molname = 'A'
    plot_properties = False

    # get binding pocket residues
    residues_near_ligand_df = pd.read_csv(residues_near_ligand_fpath)

    # iterate through structures
    struct_name_list = [
        'ET096',
        'CviUPO',
        'CviUPO-F88L+T158A',
        'DcaUPO',
        'OA167',
        'TE314'
    ]

    analyse_pocket = PocketAnalysis(pdb_dir, struct_csv_dir)
    bindingpocket_analysis, df_bindingpocket_dict, df_bindingpocket_backbone_dict = analyse_pocket(struct_name_list, residues_near_ligand_df, protein_molname, plot_properties)
    print(bindingpocket_analysis)
