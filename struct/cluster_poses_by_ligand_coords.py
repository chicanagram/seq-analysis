from __future__ import annotations
import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib.pyplot as plt
import yasara
from variables import address_dict, subfolders
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

def _pose_csv_to_fingerprint(csv_path: Path) -> np.ndarray:
    """
    Load a ligand-pose CSV and convert its atom coordinates into a 1D fingerprint:
    [x1,y1,z1, x2,y2,z2, ..., xN,yN,zN]
    """
    print(csv_path)
    df = pd.read_csv(csv_path)

    # remove rows with 'H'
    df = df[df['element']!='H']

    # extract fingerprint data
    centroid = [df['x'].mean(), df['y'].mean(), df['z'].mean()]
    df_N = df.loc[df['element']=='N'].iloc[0]
    nitrogen_coords = [df_N['x'], df_N['y'], df_N['z']]
    fingerprint = np.array(centroid + nitrogen_coords)
    print(fingerprint)

    return fingerprint


def cluster_ligand_pose_fingerprints(
    struct_base: str,
    csv_dir: Union[str, Path],
    *,
    file_pattern: Optional[str] = None,
    csv_paths: Optional[Sequence[Union[str, Path]]] = None,
    top_k_clusters: int = 3,
    random_state: int = 0,
) -> List[Dict]:
    """
    Load ligand-pose coordinate CSVs, build concatenated xyz fingerprints, cluster them,
    and return the top clusters (by population), with members ranked by distance to centroid.

    Inputs
    ------
    struct_base:
        Base structure name (e.g., "ET096", "CviUPO"). Used for filtering (if csv_paths not provided).
    csv_dir:
        Directory containing pose CSVs (used if csv_paths not provided).
    file_pattern:
        Optional glob pattern relative to csv_dir. If not set, defaults to f"{struct_base}*.csv".
    csv_paths:
        Optional explicit list of CSV paths. If provided, csv_dir/file_pattern are ignored.
    top_k_clusters:
        Number of clusters to report (default 3). If fewer poses than this, uses n_poses.
    random_state:
        Random seed for clustering.

    Assumptions
    ----------
    - Each CSV represents ONE pose, with one row per ligand atom, and columns x,y,z.
    - All poses must have the SAME number of rows (atoms) and same ordering, so the
      concatenated fingerprints are comparable.

    Returns
    -------
    List[dict], sorted by cluster_id, where cluster_id is assigned by decreasing cluster size.
    Each dict contains:
      - struct_base
      - cluster_id
      - n_structures
      - centroid_mean_xyz               (mean x,y,z of centroid's un-concatenated coords)
      - pose_names_ranked_by_centroid_distance
    """
    # --- Collect CSV files ---
    if csv_paths is not None:
        pose_files = [Path(p) for p in csv_paths]
    else:
        csv_dir = Path(csv_dir)
        patt = file_pattern or f"{struct_base}*.csv"
        pose_files = sorted(csv_dir.glob(patt))

    if not pose_files:
        raise FileNotFoundError("No pose CSVs found for the given inputs.")

    # --- Build fingerprints ---
    fingerprints: List[np.ndarray] = []
    pose_names: List[str] = []
    lengths: List[int] = []

    for p in pose_files:
        fp = _pose_csv_to_fingerprint(p)
        fingerprints.append(fp)
        pose_names.append(p.stem)
        lengths.append(fp.size)

    # Enforce same length (same number of atoms)
    uniq_lengths = sorted(set(lengths))
    if len(uniq_lengths) != 1:
        raise ValueError(
            "Pose CSVs do not all have the same number of coordinates (atoms). "
            f"Fingerprint lengths found: {uniq_lengths}. "
            "Ensure all poses have identical atom counts/order."
        )

    X = np.vstack(fingerprints)  # shape (n_poses, 3*N_atoms)
    n_poses, fp_dim = X.shape

    # --- Choose number of clusters ---
    k = int(min(top_k_clusters, n_poses))
    if k < 1:
        raise ValueError("Need at least 1 pose to cluster.")

    # --- Cluster (KMeans) ---
    try:
        from sklearn.cluster import KMeans
    except ImportError as e:
        raise ImportError(
            "scikit-learn is required for clustering in this function. "
            "Install with: pip install scikit-learn"
        ) from e

    # Standardize features to avoid any axis dominating purely by scale (x,y,z should be comparable,
    # but centering/scaling generally helps clustering stability across complexes)
    mu = X.mean(axis=0, keepdims=True)
    sigma = X.std(axis=0, keepdims=True)
    sigma[sigma == 0] = 1.0
    Xz = (X - mu) / sigma

    km = KMeans(n_clusters=k, n_init=20, random_state=random_state)
    labels = km.fit_predict(Xz)

    # Centroids back in original units (inverse transform)
    centroids_z = km.cluster_centers_                         # (k, fp_dim)
    centroids = centroids_z * sigma + mu                      # (k, fp_dim)

    # --- Rank clusters by population and remap IDs (largest cluster -> ID 0, next -> 1, ...) ---
    cluster_sizes = np.bincount(labels, minlength=k)
    order = np.argsort(-cluster_sizes)  # descending
    old_to_new = {int(old): int(new) for new, old in enumerate(order)}

    # Distances to centroid (in standardized space, consistent for ranking)
    dists = np.linalg.norm(Xz - centroids_z[labels], axis=1)

    # --- Build cluster outputs ---
    outputs: List[Dict] = []
    for old_cluster in order:
        old_cluster = int(old_cluster)
        new_cluster_id = old_to_new[old_cluster]

        members_idx = np.where(labels == old_cluster)[0]
        if members_idx.size == 0:
            continue

        # Rank members by distance to centroid (closest first)
        ranked_idx = members_idx[np.argsort(dists[members_idx])]
        ranked_pose_names = '\n'.join([pose_names[i] for i in ranked_idx])

        # Extract centroid and Nitrogen atom coordinates
        cluster_centroid_vec = centroids[old_cluster]
        centroid_xyz = cluster_centroid_vec[:3]
        nitrogen_xyz = cluster_centroid_vec[3:]

        outputs.append(
            {
                "struct_base": struct_base,
                "cluster_id": new_cluster_id,
                "n_structures": int(members_idx.size),
                "centroid_x": centroid_xyz[0],
                "centroid_y": centroid_xyz[1],
                "centroid_z": centroid_xyz[2],
                "nitrogen_x": nitrogen_xyz[0],
                "nitrogen_y": nitrogen_xyz[1],
                "nitrogen_z": nitrogen_xyz[2],
                "pose_names_ranked": ranked_pose_names.replace('_LigD',''),
            }
        )

    # Ensure outputs sorted by new cluster_id
    outputs.sort(key=lambda d: d["cluster_id"])

    # Only top_k_clusters requested (already limited by k, but keep explicit)
    return outputs[:k]

def visualize_clustered_complexes(clusters, struct_base, out_sce, pdb_dir, sce_dir, cmap, show_protein_residues_close_to_ligand=None, reps_only=False):
    yasara.info.mode = 'txt'
    yasara.Console('Off')
    yasara.Clear()

    obj_num = 0
    for cluster_id, cluster in enumerate(clusters):
        pose_names_ranked = cluster['pose_names_ranked'].split('\n')
        print(pose_names_ranked)
        for k, pose_name in enumerate(pose_names_ranked):
            if reps_only and k!=0:
                continue
            obj_num += 1
            pdb_fpath = pdb_dir + pose_name + '.pdb'
            print(obj_num, pdb_fpath)
            yasara.LoadPDB(pdb_fpath)
            # align
            if obj_num > 1:
                res = yasara.AlignObj(f'{obj_num} and Protein', '1 and Protein', method='MUSTANGPP', results=4)
            # formatting
            yasara.HideRes(f'Obj {obj_num} and protein')
            yasara.StickMol(f'Obj {obj_num} and not protein')
            # color protein residues by residue type
            yasara.ColorRes(f'Obj {obj_num} and protein', 'ResType')
            # show residues close to ligand
            if show_protein_residues_close_to_ligand is not None:
                yasara.ShowRes(f'Obj {obj_num} and protein and with distance<{show_protein_residues_close_to_ligand} from Mol D')
            # color non-protein atoms by element
            yasara.ColorRes(f'Obj {obj_num} and not protein', 'Element')
            # color heme carbons orange
            yasara.ColorAtom(f'Obj {obj_num} and Mol B and Element C', 150)
            # color ligand carbons by cluster
            yasara.ColorAtom(f'Obj {obj_num} and Mol D and Element C', cmap[cluster_id])
    sce_fpath = f'{sce_dir}{out_sce}.sce'
    yasara.SaveSce(sce_fpath)
    print(sce_fpath)



if __name__ == "__main__":
    os.chdir('../')

    # ---- user input ----
    data_folder = address_dict['PIPS2']
    data_subfolder = 'UPOs_peroxygenation_analysis/docked/'
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + 'ALL/'
    sce_dir = data_folder + subfolders['pdb'] + data_subfolder + 'ALL/sce/'
    struct_csv_dir = data_folder + subfolders['pdb'] + data_subfolder + 'ALL/structure_csv/'
    residues_near_ligand_fpath = pdb_dir + 'residues_near_ligand.csv'
    struct_base_list = [
        'ET096',
        'CviUPO',
        'CviUPO-F88L+T158A',
        'DcaUPO',
        'TE314',
        'OA167'
    ]

    cluster_res = []
    for struct_base in struct_base_list:

        # get clusters
        clusters = cluster_ligand_pose_fingerprints(
            struct_base=struct_base,
            csv_dir=struct_csv_dir,
            file_pattern=f"{struct_base}_S82_*_LigD.csv",
            top_k_clusters=3,
        )
        cluster_res += clusters

        # visualize clustered results (save as yasara scene) -- ALL structures
        out_sce = f'{struct_base}_clustered_poses_all'
        visualize_clustered_complexes(clusters, struct_base, out_sce, pdb_dir, sce_dir, cmap=[60, 120, 180], show_protein_residues_close_to_ligand=None, reps_only=False)
        # visualize clustered results (save as yasara scene) -- REPRESENTATIVES structures only
        out_sce = f'{struct_base}_clustered_poses_reps'
        visualize_clustered_complexes(clusters, struct_base, out_sce, pdb_dir, sce_dir, cmap=[60, 120, 180], show_protein_residues_close_to_ligand=5, reps_only=True)

    cluster_res = pd.DataFrame(cluster_res).round(3)
    cluster_res.to_csv(f'{struct_csv_dir}CLUSTER_LIGAND_RESULTS.csv')
    print(cluster_res)