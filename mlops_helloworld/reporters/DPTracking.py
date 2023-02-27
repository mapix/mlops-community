import os
import glob
import json
from pathlib import Path
from typing import Dict, Union, List


class DPTrackingReporter:
    @staticmethod
    def report_sampler(config_logger: Dict, sysname: str,
                       numConfs, image, molecule, traj_xtc, RMSD_dataframe,
                       step: int = -1):
        from aim import Run

        os.environ['AIM_ACCESS_TOKEN'] = config_logger["aim_personal_token"]

        # Create run under Project. Runs and Experiments are organized as:
        # AIM Experiment (i.e. DPLC Project, e.g. Uni-FEP)
        #   ├─ Run-summary (i.e. One submission of multiple systems, e.g. 20230201-test-10systems)
        #   │  ├─ Run (system 1)
        #   │  ├─ Run (system 2)
        #   │  ├─ ...
        #
        # More hierarchy can be implemented with more tags. e.g. job type, ...
        aim_run = Run(repo=config_logger["aim_repo"])
        aim_run.experiment = config_logger["project"]
        aim_run.name = config_logger["experiment"] + "-" + sysname

        tags = [config_logger["experiment"], "JobType-MLOps_demo"]
        for tag in tags:
            aim_run.add_tag(tag)

        # Log anything you want: Scalars, Distributions, Images, Figures, Tables, Molecules.
        from aim import Figure, Image, Distribution, Table, TableImage, Molecule
        tracking_data = {}
        tracking_data["numConfs"] = numConfs
        tracking_data["2D structure"] = Image(str(image))
        tracking_data["conformers"] = Molecule(str(molecule), str(traj_xtc))
        tracking_data["RMSD"] = Table(RMSD_dataframe)
        tracking_data["2D/3D structure"] = Table([
            {
                "col1": 1,
                "col2": TableImage(str(image))
            },
            {
                "col1": 1,
                "col2": Molecule(str(molecule))
            }
        ])

        for key, value in tracking_data.items():
            if step >= 0:
                aim_run.track(value, name=key, step=step, epoch=0, context={"subset": "sample"})
            else:
                aim_run.track(value, name=key, epoch=0, context={"subset": "sample"})
            # "step", "epoch" is typically used in ML training
            # only value, name are necessary

    @staticmethod
    def report_summary(config_logger: Dict, config_protocol: Dict,
                       figure_summary, image_summary):
        from aim import Run, Figure, Image, Distribution, Table, TableImage, Molecule
        os.environ['AIM_ACCESS_TOKEN'] = config_logger["aim_personal_token"]
        aim_run = Run(repo=config_logger["aim_repo"])
        aim_run.experiment = config_logger["project"]
        aim_run.name = config_logger["experiment"] + "-" + "summary"

        tags = [config_logger["experiment"], "JobType-MLOps_demo"]
        for tag in tags:
            aim_run.add_tag(tag)

        if figure_summary is not None:
            aim_run.track(Figure(figure_summary), "summary")
        if image_summary is not None:
            aim_run.track(Image(image_summary), "summary_image")
        return aim_run.hash
