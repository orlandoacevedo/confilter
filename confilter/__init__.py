FEATURES = [
    'version 0.10 : start',
    'version 0.20 : add bond connection',
    'version 0.30 : add angles, work on square cos',
    'version 0.40 : change angle to sin',
    'version 0.50 : split input as fragments, same on bonds & angles',
    'version 0.60 : change angle to itself, unit on degree',
    'version 0.70 : add more type for bond perception'
    'version 0.80 : refine filter function',
    'version 0.90 : codes execution efficiency',
    'version 1.00 : printout more info, RELEASE',
    'version 1.10 : add new class CrossFiltration',
    'version 1.20 : make class Filtration more robust & general',
    'version 1.30 : add filtration results save to image',
    'version 1.40 : add bulk process',
    'version 1.50 : bulk process is perfectly done, RELEASE',
    'version 1.60 : add argparse, RELEASE',
    'version 1.70 : printout more info, RELEASE',
    'version 1.80 : add argparse features, RELEASE',
    'version 1.90 : save probability bins',
    'version 2.00 : add Boolean for images output',
    'version 2.10 : add supprot for Backup xsf file',
    'version 2.20 : fix ReadFile to read the last molecule',
    'version 2.21 : fix error when check user input connections',
    'version 2.22 : print out final length & ratio for BulkProcess',
    'version 2.30 : preserve energy for xsf file',
    'version 2.31 : fix BulkProcess fname',
    'version 2.40 : add energy for ReadFile read_txt',
    'version 2.50 : add SaveFile xsf',
    'version 2.60 : add FILEFORMAT explanation',
    'version 2.70 : fix CrossFiltration reflist',
    'version 2.80 : fix ReadFile xsf',
    'version 2.90 : make ReadFile xsf more flexible',
    'version 3.00 : add save to xyz file',
    'version 3.10 : make fname and ftype more compatible in save to file',
    'version 3.20 : add read xyz file',
    'version 3.21 : update FILEFORMAT for xyz',
    'version 3.30 : make ReadFile more compatible',
    'version 3.3.0  : BC, add class AtomInfo',
    'version 3.3.1  : refine ReadFile',
    'version 3.3.2  : refine SaveFile',
    'version 3.3.3  : refine BondPerception',
    'version 3.3.4  : add AnglePerception',
    'version 3.3.5  : class perceptions are done',
    'version 3.3.6  : refine Filtration',
    'version 3.3.7  : refine BulkProcess',
    'version 3.3.8  : refine parsecmd',
    'version 3.4.0  : RELEASE',
    'version 3.5.0  : deep refine Filtration, RELEASE',
    'version 3.6.0  : add class PlotSamples',
    'version 3.7.0  : add subcommand plot, RELEASE',
    'version 3.7.1  : small fix on userinputs',
    'version 3.8.0  : add plot random seed',
    'version 3.9.0  : add dynamic and static method in Filtration',
    'version 4.0.0  : RELEASE',
    'version 4.0.1  : add more calculation type in Check',
    'version 4.1.0  : add report of execution time',
    'version 4.2.0  : prompt check on memory usage overhead',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()

FILEFORMAT = """
Input File Format for BOSS Output Filtration

Currently, only following file formats are supported.

For txt file:

    molecules are separated by new line, line starts with char '#' will
    be ignored, if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the first line,
          equal sign, =, can be used.

            [   # x x x [=] energy
            |   atom-1   x    y    z
    mol     |   atom-2   x    y    z
            |   atom-3   x    y    z
            |   ...
            [   <new line>


For xsf file:

    molecules are separated by '#', keyword 'ATOMS' is important,
    case sensitive, if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the first line,
          equal sign, =, can be used.

            [   # x x x [=] energy
            |
            |   ATOMS
    mol     |   atom-1   x    y    z    else
            |   atom-2   x    y    z    else
            |   atom-3   x    y    z    else
            [   ...


For xyz file:

    molecules are separated by new line, line starts with char '#' will
    be ignored, if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the first line,
          equal sign, =, can be used.

            [   #  x x x [=] energy
            |   atom-1   x    y    z    else
    mol     |   atom-2   x    y    z    else
            |   atom-3   x    y    z    else
            [   ...
"""


