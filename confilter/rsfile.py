from confilter.atominfo import FAI


class ReadFile:
    """
    Args:
        file (str): input file name
        ext (str): txt | xsf | xyz
        debug (bool): whether printout more info

    Attributes:
        system : 3D List[ List[[atomtype, x,y,z], ...], ...]
        energy : 1D List[float]  :   None means not exist
    """
    def __init__(self,file,ext=None,debug=True,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.file = file
        self.system = []
        self.energy = []
        self.debug = True if debug is True else False

        # decide file format
        if ext is None:
            ndx = file.rfind('.')
            if ndx == -1 or ndx + 1 >= len(file):
                self.ext = 'txt'
            else:
                self.ext = file[ndx+1:].lower()
        else:
            self.ext = ext

        if self.ext not in ['txt','xsf','xyz']:
            self.nice = False
            self.info = 'Fatal: file not support: {:}'.format(self.file)
            return

    def run(self):
        """process input file

        require:
            atomtype entry can be mixed, which means,

            a) 1  x  y  z
            b) H  x  y  z

            they are equivalent
        """
        if self.debug: print('Note: reading data from file: {:}'.format(self.file))
        prolist,enelist,errlist = getattr(self,'read_'+self.ext)()
        if not len(prolist): return

        if self.debug:
            for i in errlist:
                print('Warning: ignoring: {:}: {:}'.format(i[0],i[1]))

        # format: 2D str: [ [sign, number, name], ... ]
        atominfo = [[i[0], str(i[1]), i[4]] for i in FAI.atominfo]

        ndxlist = []
        for i in prolist[0]:
            atype = i[0].capitalize()
            bo = True
            for ndx in atominfo:
                if atype in ndx:
                    bo = False
                    atype = ndx[0]
                    break
            if bo: atype = i[0]
            ndxlist.append(atype)

        nats = len(ndxlist)
        errlist = []
        for n,mol in enumerate(prolist):
            if len(mol) != nats:
                info = 'Warning: ignoring: error: number of atoms'
                errlist.append([info,mol])
                continue

            # mixed atomtype may occur
            # always make sure mol be in real atomtype
            # be aware of how index in mol is changed, python pointer
            bo = True
            for cnt,atom in enumerate(mol):
                if atom[0] != ndxlist[cnt]:
                    bo = False
                    atype = atom[0].capitalize()
                    for ndx in atominfo:
                        if atype in ndx:
                            bo = True
                            atype = ndx[0]
                            break
                    if bo:
                        if atype == ndxlist[cnt]:
                            atom[0] = atype
                        else:
                            bo = False
                    if not bo:
                        break
            if bo:
                self.energy.append(enelist[n])
                self.system.append(mol)
            else:
                info = 'Warning: ignoring: error: not cooresponded'
                errlist.append([info,mol])

        if self.debug:
            for i in errlist:
                print(i[0])
                for j in i[1]: print(j)
                print()

    def read_xsf(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()
        
        promol = []
        i = 0
        while i < len(profile):
            sub = profile[i].strip()
            if not len(sub):
                i += 1
                continue
            if sub[0] == '#':
                ls = [[sub,i], ]
                j = i + 1
                while j < len(profile):
                    sub = profile[j].strip()
                    if not len(sub):
                        j += 1
                        continue
                    if sub[0] == '#': break
                    ls.append([sub,j])
                    j += 1
                promol.append(ls)
                i = j
            else:
                i += 1

        prolist = []
        enelist = []
        errlist = []
        for mol in promol:
            bo = False
            if len(mol) <= 2:
                bo = True
                errnum = mol[0][1] + 1
                errline = 'Wrong format'
            else:
                if mol[1][0] != 'ATOMS':
                    bo = True
                    errnum = mol[1][1] + 1
                    errline = 'Wrong format'

            if not bo:
                ene = None
                ltmp = mol[0][0].replace('=',' ').split()
                if len(ltmp) >= 2:
                    try:
                        ene = float(ltmp[-1])
                    except ValueError:
                        pass

                ls = []
                for t in mol[2:]:
                    # atom info
                    atom = t[0].split()
                    if len(atom) >= 4:
                        try:
                            x = float(atom[1])
                            y = float(atom[2])
                            z = float(atom[3])
                        except ValueError:
                            bo = True
                    else:
                        bo = True
                    if bo:
                        errnum = t[1] + 1
                        errline = t[0]
                        break
                    else:
                        ls.append([atom[0],x,y,z])
            if bo:
                errlist.append([errnum,errline])
            else:
                prolist.append(ls)
                enelist.append(ene)

        return prolist, enelist, errlist

    def read_txt(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()

        # List[List[[atomtype, x, y, z], ...]]
        promol = []
        mol = []
        for cnt,line in enumerate(profile):
            sub = line.strip()
            if not len(sub):
                if not len(mol): continue
                promol.append(mol)
                # initialize
                mol = []
            else:
                mol.append([sub,cnt])
        # last mol
        if len(mol): promol.append(mol)

        prolist = []
        enelist = []
        errlist = []
        for mol in promol:
            # check whether energy exist or not
            ene = None
            if mol[0][0][0] == '#':
                ltmp = mol[0][0].replace('=',' ').split()
                if len(ltmp) >= 2:
                    try:
                        ene = float(ltmp[-1])
                    except ValueError:
                        pass
                mol = mol[1:]

            bo = False
            ls = []
            for t in mol:
                # atom info
                atom = t[0].split()
                if len(atom) >= 4:
                    try:
                        x = float(atom[1])
                        y = float(atom[2])
                        z = float(atom[3])
                    except ValueError:
                        bo = True
                else:
                    bo = True
                if bo:
                    errnum = t[1] + 1
                    errline = t[0]
                    break
                ls.append([atom[0],x,y,z])
            if bo:
                errlist.append([errnum,errline])
            else:
                prolist.append(ls)
                enelist.append(ene)

        return prolist, enelist, errlist

    def read_xyz(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()

        # List[List[[atomtype, x, y, z], ...]]
        promol = []
        mol = []
        for cnt,line in enumerate(profile):
            sub = line.strip()
            if len(sub):
                mol.append([sub,cnt])
            else:
                if not len(mol): continue
                promol.append(mol)
                # initialize
                mol = []
        # last mol
        if len(mol): promol.append(mol)

        prolist = []
        enelist = []
        errlist = []
        for mol in promol:
            bo = False
            if len(mol) <= 2: bo = True
            if not bo:
                try:
                    atomnum = int(mol[0][0])
                    if atomnum + 2 != len(mol):
                        raise ValueError
                except ValueError:
                    bo = True

            if bo:
                errnum = mol[0][1] + 1
                errline = mol[0][0]
            else:
                ene = None
                # check whether energy exist or not
                ltmp = mol[1][0].replace('=',' ').split()
                if len(ltmp) > 1:
                    try:
                        ene = float(ltmp[-1])
                    except ValueError:
                        pass

                ls = []
                for t in mol[2:]:
                    # atom info
                    atom = t[0].split()
                    if len(atom) >= 4:
                        try:
                            x = float(atom[1])
                            y = float(atom[2])
                            z = float(atom[3])
                        except ValueError:
                            bo = True
                    else:
                        bo = True
                    if bo:
                        errnum = t[1] + 1
                        errline = t[0]
                        break
                    ls.append([atom[0],x,y,z])

            if bo:
                errlist.append([errnum,errline])
            else:
                prolist.append(ls)
                enelist.append(ene)

        return prolist, enelist, errlist


class SaveFile:
    """opposite operation to ReadFile

    Args:
        system : 3D List[ List[[atomtype, x,y,z], ...], ...]
        energy : 1D List[float]  :   None means not exist

        ftype (str): Output file type: txt | xsf | xyz
        fname (str): file to be saved, warning, overwritten may happen

    Note:
        1) number of atoms in system will not be crossly checked, which means
        for: system[ mol(8),  mol(5),  mol(20), ...],
        it can be successfully saved

        2) saved atomtype will always be real atomtype, if it is valid
    """
    def __init__(self,system,energy=None,ftype=None,fname=None,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.system = system if isinstance(system,list) else []
        self.energy = energy if energy is not None and isinstance(energy,list) else []

        if not len(self.system):
            self.nice = False
            self.info = 'Fatal: no inputs'
            return

        # Rule
        # self.ftype always takes the precedence
        # a) if self.fname has the period
        #       1) if its file extension matchs with self.ftype
        #          everything is fine
        #       2) if not, append real ftype on it
        # b) append real ftype
        self.ftype = None if ftype is None else ftype.lower()

        self.fname = None
        if fname is not None and isinstance(fname,str) and len(fname.split()) != 0:
            self.fname = fname.strip()

        if self.fname is not None and self.ftype is None:
            # guess ftype from fname
            ndx = self.fname.rfind('.')
            if ndx != -1:
                ext = self.fname[ndx:]
                if ext == '.':
                    self.fname = self.fname[:-1]
                else:
                    self.ftype = ext.replace('.','')

        if self.ftype is None: self.ftype = 'txt'
        if self.ftype not in ['txt','xsf','xyz']:
            self.nice = False
            self.info = 'Fatal: not support: {:}'.format(self.ftype)
            return

        if self.fname is None: self.fname = 'system.txt'
        # keep dot conversion in original fname, if it has
        ndx = self.fname.rfind('.')
        if ndx == -1:
            self.fname = self.fname + '.' + self.ftype
        else:
            if self.fname[ndx:] != '.' + self.ftype:
                self.fname = self.fname[:ndx] + '.' + self.ftype

        # get real atomtype list
        # format: 2D: [ [sign, number-int, number-str, name], ... ]
        atominfo = [[i[0], i[1], str(i[1]), i[4]] for i in FAI.atominfo]
        self.atypelist = []
        for mol in self.system:
            ls = []
            for at in mol:
                atype = at[0].capitalize() if isinstance(at,str) else at[0]
                for ndx in atominfo:
                    if atype in ndx:
                        atype = ndx[0]
                        break
                ls.append(atype)
            self.atypelist.append(ls)

    def run(self):
        fout = getattr(self,'save_'+self.ftype)()
        with open(self.fname,'wt') as f: f.write(fout)

    def save_xsf(self):
        fout = ''
        for ndx,mol in enumerate(self.system):
            if len(self.energy) and self.energy[ndx] is not None:
                fout += '# {:}\n\nATOMS\n'.format(self.energy[ndx])
            else:
                fout += '#\n\nATOMS\n'
            for cnt,at in enumerate(mol):
                atype = self.atypelist[ndx][cnt]
                fout += '{:3} {:>10} {:>10} {:>10}   1.0  1.0  1.0\n'.format(atype,*at[1:])
            fout += '\n\n'
        return fout

    def save_txt(self):
        fout = ''
        for ndx,mol in enumerate(self.system):
            if len(self.energy) and self.energy[ndx] is not None:
                fout += '#  {:}\n'.format(self.energy[ndx])
            for cnt,at in enumerate(mol):
                atype = self.atypelist[ndx][cnt]
                fout += '{:<2} {:>15} {:>15} {:>15}\n'.format(atype,*at[1:])
            fout += '\n\n'
        return fout

    def save_xyz(self):
        fout = ''
        for ndx,mol in enumerate(self.system):
            fout += '{:}\n'.format(len(mol))
            if len(self.energy) and self.energy[ndx] is not None:
                fout += 'Properties=species:S:1:pos:R:3 energy={:}\n'.format(self.energy[ndx])
            else:
                fout += 'Properties=species:S:1:pos:R:3 energy=0.0\n'
            for cnt,at in enumerate(mol):
                atype = self.atypelist[ndx][cnt]
                fout += '{:2} {:>12} {:>12} {:>12}\n'.format(atype,*at[1:])
            fout += '\n\n'
        return fout


