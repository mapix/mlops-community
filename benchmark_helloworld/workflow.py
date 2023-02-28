from typing import List
from pathlib import Path

from dflow import Step, Workflow, argo_range
from dflow.python import OP, OPIO, Artifact, OPIOSign, PythonOPTemplate, Slices


# New Dflow functionality of defining OPs! So pretty ヾ(◍°∇°◍)ﾉﾞ
@OP.function
def Hello(filename: str) -> {'foo': Artifact(Path)}:
    open(filename, "w").write("foo")
    return {'foo': Path(filename)}


@OP.function
def Check(foo: Artifact(List[str])) -> {}:
    print(foo)
    return {}


def buildWorkflow():
    wf = Workflow("benchmark-hello")
    hello = Step("hello",
                 PythonOPTemplate(Hello, image="python:3.8", image_pull_policy="IfNotPresent",
                                  slices=Slices("{{item}}",
                                                input_parameter=["filename"],
                                                output_artifact=["foo"]
                                                )
                                  ),
                 parameters={"filename": ["f1.txt", "f2.txt"]},
                 with_param=argo_range(2))
    wf.add(hello)
    check = Step("check",
                 PythonOPTemplate(Check, image="python:3.8", image_pull_policy="IfNotPresent"),
                 artifacts={"foo": hello.outputs.artifacts["foo"]})
    wf.add(check)
    return wf


if __name__ == "__main__":
    wf = buildWorkflow()
    wf.submit()
