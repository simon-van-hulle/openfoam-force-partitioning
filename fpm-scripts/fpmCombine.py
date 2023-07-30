#!/bin/python3

import numpy as np
import fpmPlotOutput as vplt
from fpmPlotOutput import log, error, warning
import os




def merge(inFileNames, outFileName="combined.dat", commentSym="#"):
    
    if not inFileNames:
        return


    header = ""
    outStr = ""

    for inFileName in inFileNames:
        localHeader = ""
        with open(inFileName, "r") as inFile:
            for line in inFile.readlines():
                if line.startswith(commentSym):
                    localHeader += line
                else:
                    outStr += line
        if not header:
            header += localHeader

    with open(outFileName, "w") as outFile:
        outFile.write(header)
        outFile.write(outStr)



def combineTerminalInput():
    args = vplt.parse_args()

    merge(args.files, args.outFile)

    if args.postprocess_dirs is not None:
        pass

    


if __name__ == "__main__":
    combineTerminalInput()
