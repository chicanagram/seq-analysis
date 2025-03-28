
"""A module wrapping the programs of CD-HIT."""

import os
import subprocess
import warnings
from itertools import chain
from pathlib import Path

PROGS = [
    "cd-hit",
    "cd-hit-2d",
    "cd-hit-est",
    "cd-hit-est-2d",
    "cd-hit-div",
    "cd-hit-454",
]
functions = [prog.replace("-", "_") for prog in PROGS]

AUXTOOLS = [
    "cd-hit-dup",
    "cd-hit-lap",
]
aux_functions = [prog.replace("-", "_") for prog in AUXTOOLS]

__all__ = functions + aux_functions


def _format_options(kwargs: dict):
    keys = (f"-{k.replace('_', '-')}" for k in kwargs)
    values = map(str, kwargs.values())
    options = chain(*zip(keys, values))
    return options


def _format_program(name: str, path: str = None):
    if path:
        dir_ = Path(path).expanduser().resolve()
        return dir_ / name
    return name


def _run(command):
    try:
        completed_process = subprocess.run(
            command,
            capture_output=True,
            check=True,
            text=True,
        )
    except subprocess.CalledProcessError as err:
        print(err.stdout)
        print(err.stderr)
        raise
    except FileNotFoundError:
        print(command)
        raise
    if "Warning" in completed_process.stdout:
        warnings.warn(completed_process.stdout)
    return completed_process


def _create_function(name: str, func_name: str, env_var: str):
    def function(**kwargs) -> subprocess.CompletedProcess:
        command = [_format_program(name, os.getenv(env_var))]
        command.extend(_format_options(kwargs))
        return _run(command)

    function.__name__ = func_name
    function.__doc__ = f"""Run command {name}.

        If environment variable `{env_var}` exists,
        it will be used as the path of the program.

        Args:
            **kwargs: Parameters and arguments of the command.

        Returns:
            The `~subprocess.CompletedProcess`.

        Raises:
            `~subprocess.CalledProcessError`: If command returns
                non-zero exit status.
            `FileNotFoundError`: If program is not installed.

        """

    return function


for prog, func in zip(PROGS, functions):
    globals()[func] = _create_function(prog, func, "CD_HIT_DIR")

for prog, func in zip(AUXTOOLS, aux_functions):
    globals()[func] = _create_function(prog, func, "CD_HIT_AUXTOOLS_DIR")
