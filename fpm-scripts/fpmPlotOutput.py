#!/bin/python3

"""
    TODO: Some Explanation on the use cases of this code
"""

from warnings import warn
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from dataclasses import dataclass
from typing import List
import time

from itertools import cycle

COLORS = [
    # "#ffa000",
    "#ffb300",
    "#ff6f00",
    # "#ffecb3",
    # "#ffc107",
    # "#212121",
    "#607d8b",
    "#212121",
    "#757575",
    "#BDBDBD",
]

INFTY = 1e100
COLOR_CYCLER = cycle(COLORS)

LINES = ["-", "--", "-.", ":"]
LINE_CYCLER = cycle(LINES)

plt.rcParams["mathtext.fontset"] = 'stix'

"""
TODO:
* Fix headerlinestart
"""


def log(*args, **kwargs):
    print(*args, **kwargs)


def error(*args, **kwargs):
    print("ERROR:   ", *args, **kwargs)


def warning(*args, **kwargs):
    print("WARNING: ", *args, **kwargs)

LABELS = {
    'Cl' : r"$C_L$",
    'Ckinematic' : r"$C_\kappa$",
    'Cvorticity' : r"$C_\omega$",
    'Cviscous' : r"$C_\sigma$",
    'Cpotential' : r"$C_\phi$",
    'Cboundary' : r"$C_\Sigma$",
    'Csum' : r"$C_\mathcal{B}^{\ FPM}$",
}


@dataclass
class Output:
    name: str = "Not Found"
    values: np.ndarray = None
    fileName: str = None

    def plot(self, time, args):
        if type(self.values) == type(None):
            return 0

        tMin = args.tmin
        tMax = args.tmax or INFTY
        idxMin = 0
        idxMax = -1

        if np.max(time) > tMin:
            idxMin = np.where(time > tMin)[0][0]
        if np.max(time) > tMax > tMin:
            idxMax = np.where(time > tMax)[0][0]

        log(f"\tPlotting {self.name}")
        plt.plot(
            time[idxMin:idxMax],
            self.values[idxMin:idxMax],
            color=next(COLOR_CYCLER),
            linestyle=next(LINE_CYCLER),
            label=LABELS.get(self.name, self.name),
        )
        return 1


@dataclass
class OutputFile:
    name: str = "Not Found"
    times: np.ndarray = None
    outputNames: List[str] = None


def outputResidual(reference: Output, approx: Output, outFiles: dict):
    times = outFiles[reference.fileName].times

    sizes = [len(o.values) for o in [reference, approx]]
    if abs(sizes[0] - sizes[1]) > 10:
        print("Trying to comparing outputs with different length")
        return

    minSize = min(sizes)

    rVals = reference.values[:minSize]
    aVals = approx.values[:minSize]

    rVals = (rVals - np.mean(rVals)) / (np.std(rVals) * minSize)
    aVals = (aVals - np.mean(aVals)) / (np.std(aVals))
    corr = np.correlate(rVals, aVals, "same")
    corrVal = float(np.correlate(rVals, aVals))

    # corr = np.correlate(reference.values[:minSize], approx.values[:minSize], mode='same')
    # print(np.correlate(reference.values[:minSize], approx.values[:minSize]))

    np.save(f"correlation-{reference.name}-{approx.name}", corr)

    plt.figure("Correlation")
    plt.clf()
    plt.plot(
        corr,
        label=(f"{reference.name} vs. {approx.name} (corr = {corrVal:5f})"),
    )
    log(f"\t{reference.name} vs. {approx.name} (corr = {corrVal})")
    ax = plt.gca()
    ax.set_ylabel(r"Relative error [-]")
    ax.axes.autoscale()
    plt.title(f"Correlation")
    plt.legend(loc="upper left", frameon=False)


def readOutputFiles(fileNames, comments="#", headerLineStart="# Time"):
    outputDict = {}
    outFileDict = {}

    for fileName in fileNames:
        log(f"\nReading output from file {fileName}")

        try:
            data = np.genfromtxt(fileName, comments=comments).T
        except FileNotFoundError:
            error("{fileName} not found")
            return

        outputList = [Output(values=vals, fileName=fileName) for vals in data[1:]]

        with open(fileName, "r") as file:
            for line in file.readlines():
                if line.startswith(headerLineStart):
                    outputNames = list(line.strip("#").split())[1:]
                    break

        if not outputNames:
            error(f"{headerLineStart} is not found in file")
            return

        for outputName, output in zip(outputNames, outputList):
            output.name = outputName
            outputDict[outputName] = output

        outFileDict[fileName] = OutputFile(
            name=fileName, times=data[0], outputNames=outputNames
        )

    return outputDict, outFileDict


def plotOutputs(outputs: dict, files: dict, args):
    numPlots = 0

    plt.figure("Outputs", figsize=(6,3))
    plt.clf()

    for name in args.names or outputs:
        output = outputs[name]
        times = files[output.fileName].times
        numPlots += output.plot(times, args)

    if numPlots > 0:
        ax = plt.gca()
        ax.set_xlabel(r"Time $t$ [s]")
        ax.set_ylabel(r"Coefficient [-]")
        ax.axes.autoscale()
        plt.ylim([args.ymin, args.ymax])
        plt.legend(loc="upper left", frameon=False, bbox_to_anchor=(1.04, 1))

    else:
        log("\tNothing to plot")

    return numPlots


def plotAll(args):

    numPlots = 0

    outputDict, outFileDict = readOutputFiles(args.files)

    numPlots += plotOutputs(outputDict, outFileDict, args)

    if numPlots < len(args.names or [None]):
        log(f"\n" f"[WARNING]: Not all of your specified coefficients were found.")
        log(f"The following coefficients are included in the files:")

        for outputName in outputDict:
            log(f"\t- {outputName}")

    if args.error_names is not None:
        if len(eNames := args.error_names) >= 2:
            for name in eNames[1:]:
                outputResidual(outputDict[eNames[0]], outputDict[name], outFileDict)

    return numPlots


def parse_args():
    parser = argparse.ArgumentParser(description="Plot coefficients")

    # Positional Arguments
    parser.add_argument(
        "files", nargs="*", metavar="files", help="File names of coefficients"
    )

    # Optional Arguments
    parser.add_argument(
        "-p",
        "--postprocess-dirs",
        nargs="+",
        metavar="",
        help="Names of postProcessing Directories for standard plotting",
    )
    parser.add_argument(
        "-n", "--names", nargs="*", metavar="", help="Names of coefficients to plot"
    )
    parser.add_argument(
        "-e",
        "--error-names",
        nargs="*",
        metavar="",
        help="Names of coefficients to find the residual off. Provide reference first, then approximation",
    )
    parser.add_argument(
        "--tmin", metavar="", type=float, default=0, help="Minimum time for plot"
    )
    parser.add_argument(
        "--tmax", metavar="", type=float, default=None, help="Minimum time for plot"
    )
    parser.add_argument(
        "--ymin",
        metavar="",
        type=float,
        default=None,
        help="Minimum value for the coefficients",
    )
    parser.add_argument(
        "--ymax",
        metavar="",
        type=float,
        default=None,
        help="Maximum value for the coefficients",
    )
    parser.add_argument(
        "-o", "--outFile", metavar="", help="File name for output image"
    )
    parser.add_argument(
        "-hls",
        "--header-line-start",
        type=str,
        default="# Time",
        metavar="",
        help="The start of the header line",
    )
    parser.add_argument(
        "-st",
        "--sleep-time",
        type=int,
        default=None,
        metavar="",
        help="Time between refreshes",
    )
    parser.add_argument(
        "-c",
        "--comment",
        type=str,
        default="#",
        metavar="",
        help="Symbol that starts a comment",
    )


    args = parser.parse_args()

    return args


def check_args(args):
    if args.postprocess_dirs is not None:
        for d in args.postprocess_dirs:
            args.files.append(
                os.path.join(d, "fpmForcePart_object/0/fpmForcePartCoefficient.dat")
            )
            args.files.append(os.path.join(d, "forceCoeffs_object/0/coefficient.dat"))

    if not args.files:
        warning("You did not provide any files")
        args.files.extend(
            [
                f"/home/simon/Documents/Delft/HPB-Delft/HPB-VIV/openFoam-VIV/fpm-fpm-lam-Re5000U50/OLD/postProcessing_May20_18h01_MATCH_Ckin0/fpmForcePart_object/0/fpmForcePartCoefficient.dat",
                f"/home/simon/Documents/Delft/HPB-Delft/HPB-VIV/openFoam-VIV/fpm-fpm-lam-Re5000U50/OLD/postProcessing_May20_18h01_MATCH_Ckin0/forceCoeffs_object/0/coefficient.dat",
            ]
        )
    return args


def plotTerminalInput():

    # Parse terminal arguments
    args = check_args(parse_args())

    while True and args.sleep_time is not None:
        try:
            # Plot the required coefficients from the required files
            numPlots = plotAll(args)

            # Show plot if any coefficients present
            if numPlots > 0:

                # plt.tight_layout()
                print("Plotting")
                plt.draw()

                plt.pause(args.sleep_time)

        except KeyboardInterrupt:
            break

    numPlots = plotAll(args)

    if numPlots > 0:
        plt.figure("Outputs")
        plt.tight_layout()
        if args.outFile:
            plt.savefig(args.outFile, dpi=500)
        else:
            plt.show()

    log("\n\nDone.\n")


if __name__ == "__main__":
    plotTerminalInput()
