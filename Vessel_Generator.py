"""
MIT License

Copyright (c) 2021 Riccardo Ferrero

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import random
import math
from Lsystem_class import *
import VesselInterpreter


def main():

    """ Use this script file to generate your L-system """
    import turtle
    #from Lsystem_class import Lsystem, rule, ruleOutput
    from PIL import Image
    from PIL import EpsImagePlugin
    EpsImagePlugin.gs_windows_binary =  r'C:\Program Files\gs\gs9.55.0\bin\gswin64c'
    # The alphabet used for the L-system
    alphabet = 'ABCE'
    # Symbols used in L-system
    symbols = 'f+-*+[]'
#
    """ L-system definition """
    # Rules definition
    ruleA = [rule('fwturn', 0.2),rule('bif', 0.8)]
    ruleB = [rule('fwturn', 0.3), rule('bif', 0.7)]
    ruleC = [rule('fwturn', 0.5), rule('bif', 0.2),rule('end',0.3)]
    ruleE = [ rule('end', 1.0)]
    # Ruleset definition
    ruleset = {'A': ruleA, 'B': ruleB, 'C': ruleC,'E': ruleE}
    # Lsystem definition (initial state, ruleset)
    ls = Lsystem([['A',1,1]], ruleset,alphabet)
    # generate the string of turtle instructions
    print("Generating L-System Instruction string set")
    instruction_string = ls.processGen(6)
    print("Drawing the Blood Vessel network")
    #print("Drawing the following L-system :\n",instruction_string)
    VesselInterpreter.createPolyline(instruction_string)
    instruction_string=[['E',0,0] if instruction[0]=='A' else instruction for instruction in instruction_string ]
    ls = Lsystem(instruction_string, ruleset,alphabet)
    instruction_string_a = ls.processGen(6)
    VesselInterpreter.createPolyline(instruction_string_a,fileOut="Vessel06")
    #instruction_string_b = ls.processGen(10)
    #VesselInterpreter.createPolyline(instruction_string_b,startingPos=np.array([0,0.2,0.4]),fileOut="vein")
    VesselInterpreter.CreateVoronoiDiagram("vtkVessel06Trunc.vtp",199,ofile="voronoiDiagram06.vtp")
    #VesselInterpreter.visualizePair()
if __name__=='__main__':
    main()

