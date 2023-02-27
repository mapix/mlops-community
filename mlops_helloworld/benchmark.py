from .workflow import buildProjectWorkflow
import dflow
from dflow import Workflow, download_artifact
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient
import json
import argparse


def setup(config):
    dflow.config["host"] = config["workflow_host"]
    dflow.config["k8s_api_server"] = config["k8s_api_server"]
    dflow.config["token"] = config["workflow_token"]
    bohrium.config["username"] = config["bohrium"]["username"]
    bohrium.config["password"] = config["bohrium"]["password"]
    bohrium.config["project_id"] = config["bohrium"]["project_id"]
    dflow.s3_config["repo_key"] = "oss-bohrium"
    dflow.s3_config["storage_client"] = TiefblueClient()
    from dflow.python import upload_packages
    import mlops_helloworld
    upload_packages.append(mlops_helloworld.__path__[0])


def submit_bm(args):
    with open(args.config, "r") as f:
        text = "".join([l for l in f])
    config = json.loads(text)
    setup(config)
    config_logger = config.get("config_logger", {})
    dataset_info = config.get("dataset_info", {})

    molecules = dataset_info.get("parameters", {}).get("molecules", [])

    wf = buildProjectWorkflow(molecules=molecules, config_logger=config_logger)
    wf_return = wf.submit()


def status_bm(args):
    with open(args.config, "r") as f:
        text = "".join([l for l in f])
    config = json.loads(text)
    setup(config)
    jobid = args.jobid
    wf = Workflow(id=jobid)
    print(wf.query_status())


def download_bm(args):
    with open(args.config, "r") as f:
        text = "".join([l for l in f])
    config = json.loads(text)
    setup(config)
    jobid = args.jobid
    prefix = args.prefix
    wf = Workflow(id=jobid)
    out_step = wf.query_step(name="hello")[0]
    download_artifact(out_step.outputs.artifacts["foo"], path=prefix)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    submit_parser = subparsers.add_parser("submit")
    submit_parser.add_argument(
        'config', type=str, help='Configure file of benchmark test.')
    submit_parser.set_defaults(func=submit_bm)

    status_parser = subparsers.add_parser("status")
    status_parser.add_argument(
        'config', type=str, help='Configure file of benchmark test.')
    status_parser.add_argument(
        'jobid', type=str, help='Job-id of benchmark test.')
    status_parser.set_defaults(func=status_bm)

    download_parser = subparsers.add_parser("download")
    download_parser.add_argument(
        'config', type=str, help='Configure file of benchmark test.')
    download_parser.add_argument(
        'jobid', type=str, help='Job-id of benchmark test.')
    download_parser.add_argument(
        '--prefix', type=str, help='Prefix of download pathway.')
    download_parser.set_defaults(func=download_bm)

    args = parser.parse_args()
    args.func(args)
