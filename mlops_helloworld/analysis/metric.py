from typing import Any, Tuple


def compute_metrics(molecule, traj) -> Tuple[int, Any, str, str]:
    """

    Parameters
    ----------
    molecule: path
        .pdb file of the molecule
    traj: path
        .npy file of array of coordinates (numAtoms × dims × numConfs)

    Returns
    -------
    numConfs: int
        Number of conformers
    df: pd.DataFrame
        DataFrame of RMSD info

    """
    import numpy as np
    import pandas as pd
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    coordinates = np.load(traj)
    u = mda.Universe(molecule, coordinates)
    numConfs = u.trajectory.n_frames

    R = rms.RMSD(atomgroup=u.atoms).run()
    df = pd.DataFrame(R.rmsd, columns=['Frame', 'time', 'RMSD'])

    traj_xtc = "traj.xtc"
    rmsd_csv = "rmsd.csv"
    df.to_csv(rmsd_csv)
    u.atoms.write(traj_xtc, frames="all")
    return numConfs, df, rmsd_csv, traj_xtc
