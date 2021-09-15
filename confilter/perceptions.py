from confilter.atominfo import FAI


class BondPerception:
    """Bond connecting perception

    Reference:
        Zhang, Q., et al.
        A rule-based algorithm for automatic bond type perception.
        J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26

    Args:
        system : single: 2D List[[atomtype, x,y,z], ...]
        userinputs (bool): if it is True, index in list starts at 1

    Attributes:
        conb        : all bonds connections
        bcon        : perceptive bonds connections
        nconb       : all nonbonds connections
        fconb       : fragments all bonds connections
        fnconb      : fragments nonbonds connections
        cfnconb     : cross fragments nonbonds connections
        fragments   : 2D List[ List[int], ... ]
        atradius    : atom radius   :   1D List[float]
    """
    def __init__(self,system,userinputs=None,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.system = system if isinstance(system,list) else []
        self.userinputs = True if userinputs is True else False

        if len(self.system) <= 2:
            self.nice = False
            self.info = 'Fatal: too less'
            return

        # get real atomtype list
        # format: 2D str: [ [sign, number, name], ... ]
        atominfo = [[i[0], str(i[1]), i[4]] for i in FAI.atominfo]
        self.atradius = []
        for atom in self.system:
            bo = True
            if not isinstance(atom,list) or len(atom) != 4:
                bo = False
            else:
                bo = False
                atype = str(atom[0]) if isinstance(atom[0],int) else atom[0].capitalize()
                for cnt,ndx in enumerate(atominfo):
                    if atype in ndx:
                        bo = True
                        r = FAI.atominfo[cnt][2]
                        break
            if not bo:
                self.nice = False
                self.info = 'Fatal: wrong defined: atom: {:}'.format(atom[0])
                return
            self.atradius.append(r)

    def run(self):
        self.conb, self.bcon, self.nconb = self.calc_bcons(self.system,self.atradius)

        # based on self.bcon, calculate fragments
        # note: return is zero-based
        self.fragments = self.calc_bfs(self.bcon,len(self.system))

        self.fconb = []
        self.fnconb = []
        for ref in self.fragments:
            fc, fbc, fnbc = self.calc_bcons(self.system,self.atradius,reflist=ref)
            self.fconb.extend(fc)
            self.fnconb.extend(fnbc)

        self.cfnconb = []
        for i,ref in enumerate(self.fragments):
            j = i + 1
            while j < len(self.fragments):
                for t in ref:
                    for k in self.fragments[j]:
                        if t > k: t, k = k, t
                        if [t,k] not in self.cfnconb:
                            self.cfnconb.append([t,k])
                j += 1

        if self.userinputs:
            self.conb = [[i+1 for i in t] for t in self.conb]
            self.bcon = [[i+1 for i in t] for t in self.bcon]
            self.nconb = [[i+1 for i in t] for t in self.nconb]
            self.fragments = [[i+1 for i in t] for t in self.fragments]
            self.fconb = [[i+1 for i in t] for t in self.fconb]
            self.fnconb = [[i+1 for i in t] for t in self.fnconb]
            self.cfnconb = [[i+1 for i in t] for t in self.cfnconb]

    def calc_bcons(self,system,atradius,reflist=None):
        con = []
        bcon = []
        nconb = []
        reflist = list(range(len(system))) if reflist is None else reflist
        for i,ref in enumerate(system[:-1]):
            if i not in reflist: continue
            j = i + 1
            while j < len(system):
                if j in reflist:
                    atom = system[j]
                    dx = ref[1] - atom[1]
                    dy = ref[2] - atom[2]
                    dz = ref[3] - atom[3]
                    ds = dx*dx + dy*dy + dz*dz
                    rij = atradius[i] + atradius[j] + 0.4
                    if ds >= 0.64 and ds <= rij*rij:
                        bcon.append([i,j])
                    else:
                        nconb.append([i,j])
                    con.append([i,j])
                j += 1
        return con, bcon, nconb

    def calc_bfs(self,bcon,tot=None):
        """Breadth First Search to calculate 2D collections

        Return:
            bfs : 2D [fragments, fragments, ...]
        """
        if tot is None:
            tot = max([max(i) for i in bcon]) + 1
        visited = [False for i in range(tot)]
        bfs = []
        for ref in bcon:
            if visited[ref[0]]: continue
            queue = [ref[0]]
            visited[ref[0]] = True
            ls = []
            while queue:
                s = queue.pop()
                ls.append(s)
                for xmp in bcon:
                    if s == xmp[0] and not visited[xmp[1]]:
                        visited[xmp[1]] = True
                        queue.append(xmp[1])
                    if s == xmp[1] and not visited[xmp[0]]:
                        visited[xmp[0]] = True
                        queue.append(xmp[0])
            bfs.append(sorted(ls))
        flatten = []
        for i in bfs: flatten.extend(i)
        for i in range(tot):
            if i not in flatten:
                bfs.append([i])
        return bfs


class AnglePerception(BondPerception):
    """
    Attributes: (new)
        cona        : all angles connections
        acon        : perceptive angles connections
        nacon       : all nonangles connections
        fcona       : fragments all angles connections
        fncona      : fragments nonangles connections
        cfnacon     : cross two fragments nonangles connections
        cfncona     : cross three fragments nonangles connections
    """
    def __init__(self,system,*args,**kwargs):
        super().__init__(system,*args,**kwargs)
        if not self.nice: return

    def run(self):
        super().run()
        self.cona = []
        self.acon = []
        self.ncona = []
        self.fcona = []
        self.fncona = []
        self.cfnacon = []
        self.cfncona = []
        if len(self.system) < 3: return

        # take care results from BondPerception
        mybcon = self.bcon
        myfragments = self.fragments
        if self.userinputs:
            mybcon = [[i-1 for i in j] for j in self.bcon]
            myfragments = [[i-1 for i in j] for j in self.fragments]

        self.cona = self.calc_conas(list(range(len(self.system))))
        self.acon = self.calc_acons(mybcon)

        # deep copy
        for ref in self.cona:
            if ref not in self.acon:
                self.ncona.append([i for i in ref])

        for ref in myfragments:
            newfcona = self.calc_conas(ref)
            self.fcona.extend(newfcona)
            gfbcon = []
            for i in ref:
                for zmp in mybcon:
                    if i in zmp:
                        gfbcon.append(zmp)
            # it may have repeats
            newfbcon = []
            for i,ref in enumerate(gfbcon):
                j = i + 1
                bo = True
                while j < len(gfbcon):
                    if ref[0] in gfbcon[j] and ref[1] in gfbcon[j]:
                        bo = False
                        break
                    j += 1
                if bo:
                    newfbcon.append(ref)

            newfacon = self.calc_acons(newfbcon)
            # deep copy
            for ndx in newfcona:
                if ndx not in newfacon:
                    self.fncona.append([i for i in ndx])

        for i,ref in enumerate(myfragments):
            for j in range(i+1, len(myfragments)):
                for k in range(j+1, len(myfragments)):
                    for ai in ref:
                        for aj in myfragments[j]:
                            for ak in myfragments[k]:
                                # center is ai
                                if aj < ak:
                                    self.cfncona.append([aj,ai,ak])
                                else:
                                    self.cfncona.append([ak,ai,aj])
                                # center is aj
                                if ai < ak:
                                    self.cfncona.append([ai,aj,ak])
                                else:
                                    self.cfncona.append([ak,aj,ai])
                                # center is ak
                                if ai < aj:
                                    self.cfncona.append([ai,ak,aj])
                                else:
                                    self.cfncona.append([aj,ak,ai])

        for i,ref in enumerate(myfragments):
            if len(ref) < 2: continue
            for j in range(i+1, len(myfragments)):
                for s in range(len(ref)):
                    ai = ref[s]
                    for t in range(s+1,len(ref)):
                        aj = ref[t]
                        for ak in myfragments[j]:
                            # center is ai
                            if aj < ak:
                                self.cfnacon.append([aj,ai,ak])
                            else:
                                self.cfnacon.append([ak,ai,aj])
                            # center is aj
                            if ai < ak:
                                self.cfnacon.append([ai,aj,ak])
                            else:
                                self.cfnacon.append([ak,aj,ai])
                            # center is ak
                            if ai < aj:
                                self.cfnacon.append([ai,ak,aj])
                            else:
                                self.cfnacon.append([aj,ak,ai])

        if self.userinputs:
            self.cona = [[i+1 for i in j] for j in self.cona]
            self.acon = [[i+1 for i in j] for j in self.acon]
            self.ncona = [[i+1 for i in j] for j in self.ncona]
            self.fcona = [[i+1 for i in j] for j in self.fcona]
            self.fncona = [[i+1 for i in j] for j in self.fncona]
            self.cfnacon = [[i+1 for i in j] for j in self.cfnacon]
            self.cfncona = [[i+1 for i in j] for j in self.cfncona]

    def calc_conas(self,nlist):
        """calc all angles connections based on given list
        Args:
            nlist: 1D List[int]
        """
        if len(nlist) < 3: return []
        cona = []
        i = 0
        while i < len(nlist):
            j = i + 1
            while j < len(nlist):
                k = j + 1
                while k < len(nlist):
                    # center index i
                    if nlist[j] < nlist[k]:
                        cona.append([nlist[j], nlist[i], nlist[k]])
                    else:
                        cona.append([nlist[k], nlist[i], nlist[j]])
                    # center index j
                    if nlist[i] < nlist[k]:
                        cona.append([nlist[i], nlist[j], nlist[k]])
                    else:
                        cona.append([nlist[k], nlist[j], nlist[i]])
                    # center index k
                    if nlist[j] < nlist[i]:
                        cona.append([nlist[j], nlist[k], nlist[i]])
                    else:
                        cona.append([nlist[i], nlist[k], nlist[j]])
                    k += 1
                j += 1
            i += 1
        return cona

    def calc_acons(self,bcon):
        if len(bcon) < 2: return []
        # rule: acon: [i,j,k]   :   i < k
        acon = []
        for cnt,ref in enumerate(bcon):
            # center is ref[1]
            i, j = ref[0], ref[1]
            for num in range(len(bcon)):
                if num == cnt: continue
                if j == bcon[num][0]:
                    k = bcon[num][1]
                elif j == bcon[num][1]:
                    k = bcon[num][0]
                else:
                    k = None
                if k is not None:
                    if k < i: i, k = k, i
                    if [i,j,k] not in acon:
                        acon.append([i,j,k])
            # center is ref[0]
            j, i = ref[0], ref[1]
            for num in range(len(bcon)):
                if num == cnt: continue
                if j == bcon[num][0]:
                    k = bcon[num][1]
                elif j == bcon[num][1]:
                    k = bcon[num][0]
                else:
                    k = None
                if k is not None:
                    if k < i: i, k = k, i
                    if [i,j,k] not in acon:
                        acon.append([i,j,k])
        return acon


