# Conformer Filtration (confilter)

Conformer Filtration (confilter) is a package for unique molecule filtration in machine learning databases.


Xiang Zhong and [Orlando Acevedo*](https://web.as.miami.edu/chemistrylabs/acevedogroup/research.html), University of Miami


## Download

```
git clone git://github.com/orlandoacevedo/confilter.git
```

## Installation

```
pip[3] install confilter
```

Or using source codes

```
python[3] setup.py install
```

## Usage

For helpful information, use

```
confilter
```

Or

```
confilter -h
```

Or, for sub-command

```
confilter plot -h
```

```
usage:

    confilter [-h] [-v] [-f file [file ...]] [-d file [file ...]] [--static] [--separat]
    [--vndx VNDX] [--borandom] [--seed SEED] [-btol BTOL] [-atol ATOL]
    [-bcon B1 B2, B1-B2 [B1 B2, B1-B2 ..]] [-acon A1 A2 A3, A1-A2-A3 [A1 A2 A3, A1-A2-A3 ...]]
    [-g G [G ...]] [--no-oball] [-obpar] [--no-oaall] [-oapar] [-nc] [--features] [-p] [-nu]
    [-o FNAME] [-ft FTYPE]
    {plot} ...

Conformation Filtration

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -f file [file ...], --datafilelist file [file ...]
                        Data files, separate by space or comma
  -d file [file ...], --indexfilelist file [file ...]
                        Index files, as the filtration reference
  --static              turn on static mode calculation, default is in dynamic/all mode
  --separate            valid in dynamic mode, change to dynamic/separate mode
  --vndx VNDX           valid in static mode, set filtration index value
  --borandom            valid in static mode, rather than by lowest-bit, filtering out molecules randomly
  --seed SEED           random seed, (optional, highest priority)
  -btol BTOL, --btol BTOL
                        bonds tolerance, any changes in bcon smaller than it will be filtered out
  -atol ATOL, --atol ATOL
                        angles tolerance, any changes in acon smaller than it will be filtered out
  -bcon B1 B2, B1-B2 [B1 B2, B1-B2 ...]
                        Bond connections, in pairs, separate by comma
  -acon A1 A2 A3, A1-A2-A3 [A1 A2 A3, A1-A2-A3 ...]
                        Angle connections, in pairs, separate by comma
  -g G [G ...], --fragments G [G ...]
                        fragments, separate by comma
  --no-oball            turn off bonds probability overall calculation, Boolean
  -obpar, --obpar       turn on bonds probability parameters calculation, Boolean
  --no-oaall            turn off angles probability overall calculation, Boolean
  -oapar, --oapar       turn on angles probability parameters calculation, Boolean
  -nc, --no-force-double-check
                        turn off double check prompt info before execution
  --features            show development features
  -p, --file-format-explanations
                        show input system file format explanations
  -nu, --no-userinputs  specify the indexes of input connections start at 0
  -o FNAME, --fname FNAME
                        output system file name
  -ft FTYPE, --ftype FTYPE
                        output system file type, [txt, xsf, xyz]

continuous subcommand:
  {plot}
    plot                plot cross comparsions based on different chosen samples
```

```
usage: confilter plot [-h] [-bf B [B ...]] [-nl n [n ...]] [-ns n] [-sn n] [-en n] [-inc n] [-nr n]

optional arguments:
  -h, --help            show this help message and exit
  -bf B [B ...], --probdatafilelist B [B ...]
                        bulk process probability files
  -nl n [n ...], --nmlist n [n ...]
                        highest priority, select number of samples for plot, (optional)
  -ns n, --nmsamples n  choose number of samples for plot, default is 5
  -sn n, --startndx n   start index for the inputs, (optional)
  -en n, --endndx n     end index for the inputs, (optional)
  -inc n, --incndx n    increments for choose, (optional)
  -nr n, --nmranges n   random ranges for increments, (optional)
```

Please refer to **doc/** for more information.

## About

**Contributing Authors**: Xiang Zhong and [Orlando Acevedo*](https://web.as.miami.edu/chemistrylabs/acevedogroup/research.html)

**Software License**: Conformer Filtration (confilter) software package.

**Funding**: Gratitude is expressed to the National Science Foundation (CHE-2102038) for the support of this research.

Copyright (C) 2023  Orlando Acevedo


