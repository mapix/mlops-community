import os
from typing import List, Dict
from pathlib import Path

from dflow import Step, Steps, Workflow, Inputs, InputParameter, Outputs, OutputParameter, OutputArtifact, argo_range
from dflow.python import OP, Artifact, Parameter, PythonOPTemplate, Slices


@OP.function
# New Dflow functionality of defining OPs! So pretty ヾ(◍°∇°◍)ﾉﾞ
def Sampling(smiles: str) -> {'molecule': Artifact(Path), 'traj': Artifact(Path), 'image': Artifact(Path)}:
    """
    Demo Simulation Algorithm Run: Generate conformers with a given molecule SMILES.

    Parameters:
    ----------`
    smiles: str

    Returns:
    -------
    molecule: path
        .pdb file of the molecule
    traj: path
        .npy file of array of coordinates (numAtoms × dims × numConfs)
    image: path
        .png 2D structure image of the molecule

    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.AllChem import EmbedMolecule, EmbedMultipleConfs, Compute2DCoords
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Compute2DCoords(mol)
    img = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(200, 200))
    img.save("molecule.png")

    EmbedMultipleConfs(mol, numConfs=5)
    coordinates = np.stack([conf.GetPositions() for conf in mol.GetConformers()])

    Chem.MolToPDBFile(mol, "molecule.pdb", confId=0)
    np.save("traj.npy", coordinates)
    return {'molecule': Path("molecule.pdb"), 'traj': Path("traj.npy"), 'image': Path("molecule.png")}


@OP.function
def MetricReport(molecule: Artifact(Path), traj: Artifact(Path), image: Artifact(Path), smiles: str, config_logger: Dict) -> {'rmsd': List[float]}:
    # Initialize tracking environment
    import os
    print(">>> Install AIM...")
    os.system('pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple')
    os.system('pip install pandas')
    os.system('pip install -U dp-tracking-sdk dp-metadata-sdk -i https://repo.mlops.dp.tech/repository/pypi-group/simple')

    # Compute Metrics. Now we have molecule, traj, image, rmsd.csv in current path
    from mlops_helloworld.analysis.metric import compute_metrics
    numConfs, RMSD_dataframe, rmsd_csv, traj_xtc = compute_metrics(molecule, traj)

    from mlops_helloworld.reporters import DPTrackingReporter
    DPTrackingReporter.report_sampler(config_logger=config_logger, sysname=smiles,
                                      numConfs=numConfs, image=image, molecule=molecule, traj_xtc=Path(traj_xtc), RMSD_dataframe=RMSD_dataframe)
    return {"rmsd": RMSD_dataframe["RMSD"].to_list()}


@OP.function
def Summary(molecules: List[str], metrics: List[List[float]], config_logger: Dict, config_protocol: Dict) -> {}:
    print(">>> Install AIM...")
    os.system('pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple')
    os.system('pip install numpy pandas scipy plotly kaleido')
    os.system('pip install -U dp-tracking-sdk dp-metadata-sdk -i https://repo.mlops.dp.tech/repository/pypi-group/simple')

    print(metrics)
    from mlops_helloworld.analysis.visualizations import draw_multiple_distribution
    fig, imgfile = draw_multiple_distribution(data=metrics, systems=molecules)

    from mlops_helloworld.reporters import DPTrackingReporter, FeishuReporter
    run_hash = DPTrackingReporter.report_summary(config_logger=config_logger, config_protocol=config_protocol,
                                                 figure_summary=fig, image_summary=imgfile)
    FeishuReporter.report_summary(feishu_groups=config_logger["feishu_groups"],
                                  experiment_group=config_logger["experiment"],
                                  project=config_logger["project"],
                                  success_ratio=f"{len(metrics)} / {len(molecules)}",
                                  config_protocol=config_protocol,
                                  imgfile=imgfile,
                                  track_url=f"https://tracking-demo.mlops.dp.tech/runs/{run_hash}",
                                  workflow_url=f"https://workflows.deepmodeling.com/workflows/argo/"
                                               f"{os.environ['ARGO_NODE_ID'].rsplit('-', maxsplit=1)[0]}?tab=workflow")
    return {}


def singleSystemProtocol(config_logger: Dict):
    single = Steps("run-single-system")
    single.inputs.parameters = {"smiles": InputParameter()}

    sampling = Step("sampling",
                    PythonOPTemplate(Sampling, image="daverona/rdkit", image_pull_policy="IfNotPresent"),
                    parameters={"smiles": single.inputs.parameters["smiles"]})
    single.add(sampling)
    metric_log = Step("metric-report",
                      PythonOPTemplate(MetricReport, image="gturtu21/ubuntu_mda", image_pull_policy="IfNotPresent", command="/usr/bin/python3"),
                      artifacts={"molecule": sampling.outputs.artifacts["molecule"],
                                 "traj": sampling.outputs.artifacts["traj"],
                                 "image": sampling.outputs.artifacts["image"]},
                      parameters={"smiles": single.inputs.parameters["smiles"],
                                  "config_logger": config_logger})
    single.add(metric_log)

    single.outputs.parameters["rmsd"] = OutputParameter(value_from_parameter=metric_log.outputs.parameters["rmsd"])
    return single


def buildProjectWorkflow(molecules: List[str] = ["c1ccccc1", "c1ccccn1"], config_logger: Dict = {}):
    molecules_dataset = molecules

    wf = Workflow("project-mlops-demo")
    single = singleSystemProtocol(config_logger=config_logger)

    all = Step(name="run-single-system", template=single, with_param=argo_range(len(molecules_dataset)),
               slices=Slices("{{item}}", input_parameter=["smiles"], output_parameter=["rmsd"]),
               parameters={"smiles": molecules_dataset},
               parallelism=120)
    wf.add(all)

    summary_step = Step("summary-report", PythonOPTemplate(Summary, image="python:3.8", image_pull_policy="IfNotPresent"),
                        parameters={"molecules": molecules_dataset,
                                    "metrics": all.outputs.parameters["rmsd"],
                                    "config_logger": config_logger,
                                    "config_protocol": {}})
    wf.add(summary_step)
    return wf


if __name__ == "__main__":
    wf = buildProjectWorkflow([], {})
    wf.submit()
