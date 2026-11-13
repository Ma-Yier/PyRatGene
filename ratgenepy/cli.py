"""Command line interface for RatGene."""

from __future__ import annotations

import argparse
import os
from typing import Optional, Sequence

from .core import RatGene
from .model import Model

def _load_model_file(path: str) -> Model:
    candidates = [path]
    if not path.lower().endswith(".mat"):
        candidates.append(path + ".mat")
        candidates.append(os.path.join("data", path + ".mat"))
        candidates.append(os.path.join(os.getcwd(), "data", path + ".mat"))
    for candidate in candidates:
        if os.path.exists(candidate):
            return Model.load_mat(candidate)
    raise FileNotFoundError(path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Hyper parameters for RatGene")
    parser.add_argument(
        "model", type=str,
        help="path to the .mat model file downloaded from BiGG"
    )
    parser.add_argument(
        "targetMet", type=str,
        help="the target metabolite (a string in model.mets)"
    )
    parser.add_argument(
        "--biomass", type=str, default=None,
        help="the biomass reaction in the model (default: objective reaction)"
    )
    parser.add_argument(
        "--carbon", type=str, default="EX_glc__D_e",
        help="the carbon source exchange reaction (default: EX_glc__D_e)"
    )
    parser.add_argument(
        "--oxygen", type=str, default="EX_o2_e",
        help="the oxygen exchange reaction (default: EX_o2_e)"
    )
    parser.add_argument(
        "--LBbiomass", type=float, default=0.05,
        help="lower threshold of the biomass reaction (default: 0.05)"
    )
    parser.add_argument(
        "--LBcarbon", type=float, default=-15.0,
        help="lower threshold of the carbon source (default: -15)"
    )
    parser.add_argument(
        "--LBoxygen", type=float, default=-15.0,
        help="lower threshold of the oxygen source (default: -15)"
    )
    parser.add_argument(
        "--maxLoop", type=int, default=1000,
        help="maximum iterations of the ratio-based procedure (default: 1000)"
    )
    parser.add_argument(
        "-t", "--timeLimit", type=float, default=float("inf"),
        help="maximum computation time in seconds (default: no limit)"
    )
    parser.add_argument(
        "-p", "--pool", type=int, default=10,
        help="number of solutions taken from the MILP solution pool (default: 10)"
    )
    parser.add_argument(
        "--type", type=str, default="gene", choices=["gene", "reaction", "1", "0"],
        help="type of modification strategy: gene or reaction (default: gene)"
    )
    parser.add_argument(
        "-s", "--size", type=str, default="true",
        choices=["true", "false", "True", "False", "1", "0"],
        help="whether to reduce the size of the output strategy (default: true)"
    )
    parser.add_argument(
        "--addition", type=str, default=None,
        help="another .mat model supplying additional information (addition problem)"
    )
    parser.add_argument(
        "--backend", "--solver", dest="backend", type=str, default="scip",
        choices=["scip", "cplex", "gurobi"],
        help="LP/MILP solver backend (default: scip)"
    )
    parser.add_argument(
        "--integerTolerance", type=float, default=None,
        help=("MILP integrality tolerance (default: the smallest value "
              "supported by the selected backend)")
    )

    args = parser.parse_args(argv)

    model = _load_model_file(args.model)
    addition = _load_model_file(args.addition) if args.addition else None

    modification_type = 1 if args.type in ("gene", "1") else 0
    reduce_size = 1 if args.size in ("true", "True", "1") else 0

    result = RatGene(
        model,
        args.targetMet,
        biomass=args.biomass,
        carbon=args.carbon,
        oxygen=args.oxygen,
        LBbiomass=args.LBbiomass,
        LBcarbon=args.LBcarbon,
        LBoxygen=args.LBoxygen,
        maxLoop=args.maxLoop,
        timeLimit=args.timeLimit,
        pool=args.pool,
        type=modification_type,
        size=reduce_size,
        addition=addition,
        backend=args.backend,
        integerTolerance=args.integerTolerance,
    )

    x_target = result[0]
    knockouts = result[1]
    print(f"The number of knockouts is: {len(knockouts)}")
    if len(result) > 2:
        print(f"The number of additions is: {len(result[2])}")
    print(f"The target reaction rate with the strategy applied: {x_target:.4f} mmol/gDW/h")
    print("Knockouts:")
    for name in knockouts:
        print(f"  {name}")
    if len(result) > 2:
        print("Additions:")
        for name in result[2]:
            print(f"  {name}")

    return 0
