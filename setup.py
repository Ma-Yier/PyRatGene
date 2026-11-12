"""Package configuration for RatGene.

Install into the current Python environment with::

    python setup.py install

or, for development, with::

    pip install -e .
"""

from __future__ import annotations

from pathlib import Path

from setuptools import find_packages, setup

ROOT = Path(__file__).parent

_long_description = ""
_readme = ROOT / "README.md"
if _readme.exists():
    _long_description = _readme.read_text(encoding="utf-8")

setup(
    name="ratgenepy",
    version="1.0.0",
    description=(
        "Growth-to-production ratio-based gene and reaction design "
        "algorithms for constraint-based metabolic networks (RatGene)"
    ),
    long_description=_long_description,
    long_description_content_type="text/markdown",
    author="RatGene contributors",
    license="MIT",
    packages=find_packages(exclude=("tests", "tests.*")),
    py_modules=["RatGene"],
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "scipy",
        "pyscipopt",
    ],
    extras_require={
        "cplex": ["cplex"],
        "gurobi": ["gurobipy"],
    },
    entry_points={
        "console_scripts": [
            "ratgene=ratgenepy.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
