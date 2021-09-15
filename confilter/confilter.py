import math
import os
import random
import sys
import tempfile
import time
import matplotlib.pyplot as plt
from confilter.perceptions import AnglePerception
from confilter.rsfile import ReadFile, SaveFile
from confilter.utils import file_gen_new, getrealsizeof, plot_save_image


class Filtration:
    """Filter molecules based on bonds & angles connections

    Inputs:
        system : 3D List[ List[[atomtype, x,y,z], ...], ...]
        userinputs (bool): if it is True, index in list starts at 1

        bcon : System[[atom-i, atom-j], ...] : List[[int,int], ...]
        acon : System[[atom-i, atom-j, atom-k], ...] : List[[int,int,int], ...]

        btol : tolerance on bond, Angstrom
        atol : tolerance on angle, degree

        obpar  :  Boolean  :  whether calculate bonds prob_bpar  :  default False
        oball  :  Boolean  :  whether calculate bonds prob_ball  :  default True
        oapar  :  Boolean  :  whether calculate angles prob_apar  :  default False
        oaall  :  Boolean  :  whether calculate angles prob_aall  :  default True


    Attributes:
        system  :  good molecules after filtration
        sysbad  :  filtered out molecules

        prob_begin : begin probability  :   dict
            keys:
                bpar    : 3D : List[ [List[int],float] ]
                ball    : 2D : List[ List[int],  float ]
                apar    : 3D : List[ [List[int],float] ]
                aall    : 2D : List[ List[int],  float ]

        prob_final :  same format as prob_begin

        bondlist   :  2D  :  List[ List[float] ]  :  good, correspond to bcon
        anglelist  :  2D  :  List[ List[float] ]  :  good, correspond to acon

        reflist    :  1D  List[int]     : sorted index of for bad molecules


    Caution:
            Because the BOSS solute move is only for the single time move,
            so Sum(difference for untouched atoms) =~ 0.0,
            thus, when doing filtration, please take care of btol & atol.
    """
    def __init__(self,system=None,keepndxlist=None,userinputs=None,
                bcon=None,acon=None,btol=None,atol=None,seed=None,
                mode=None,vndx=None,borandom=None,boall=None,
                obpar=None,oball=None,oapar=None,oaall=None,
                *args,**kwargs):
        self.system = system
        self.keepndxlist = keepndxlist
        self.userinputs = True if userinputs is True else False

        self.bcon = [] if bcon is None else bcon
        self.acon = [] if acon is None else acon
        self.btol = 0.1 if btol is None else btol   # Angstrom
        self.atol = 0.1 if atol is None else atol   # degree

        self.seed = seed if seed else random.randrange(100000000)
        random.seed(self.seed)

        if mode is None or mode.lower() in ['d','dynamic']:
            self.mode = 'dynamic'
        else:
            self.mode = 'static'
        self.vndx = vndx
        self.borandom = borandom
        self.boall = True if boall is None else boall

        self.obpar = True if obpar is True else False
        self.oball = False if oball is False else True
        self.oapar = True if oapar is True else False
        self.oaall = False if oaall is False else True

        if self.userinputs:
            self.userinputs = False
            self.bcon = [[i-1 for i in j] for j in self.bcon]
            self.acon = [[i-1 for i in j] for j in self.acon]
        self.prob_begin = {'bpar':[], 'ball':[], 'apar':[], 'aall':[]}
        self.prob_final = {'bpar':[], 'ball':[], 'apar':[], 'aall':[]}

    def run(self):
        """attemption on filtering
        """
        # to improve efficiency, bondlist only needs to be calculated once
        if len(self.bcon): print('Note: calculating bonds connections ...')
        bondlist = self.calc_square_distance(self.system,self.bcon)
        if len(self.acon): print('Note: calculating angles connections ...')
        anglelist = self.calc_angle_degree(self.system,self.acon)

        # increments
        binc = self.btol * self.btol
        ainc = self.atol

        if self.keepndxlist is None: self.keepndxlist = []
        # important, increase efficiency
        self.keepndxlist = sorted(self.keepndxlist)

        if self.obpar or self.oball:
            print('Note: calculating begin bonds probability ...')
            self.prob_begin['bpar'], self.prob_begin['ball'] = self.calc_probs(bondlist,binc,self.obpar,self.oball)
        if self.oapar or self.oaall:
            print('Note: calculating begin angles probability ...')
            self.prob_begin['apar'], self.prob_begin['aall'] = self.calc_probs(anglelist,ainc,self.oapar,self.oaall)
        print('Note: calculating repeats reference ...')
        self.reflist = self.calc_filterlists(bondlist,anglelist,binc,ainc,mode=self.mode,vndx=self.vndx,
                                            borandom=self.borandom,boall=self.boall,keepndxlist=self.keepndxlist)
        
        print('Note: updating ...')
        tmpsys = []
        self.bondlist = []
        self.anglelist = []
        self.sysbad = []
        cnt = 0
        trn = 0
        self.reflist.append(-1)
        self.keepndxlist.append(-1)
        for ndx in range(len(self.system)):
            if ndx == self.reflist[cnt]:
                cnt += 1
                self.sysbad.append(self.system[ndx])
            else:
                if ndx == self.keepndxlist[trn]:
                    trn += 1
                else:
                    tmpsys.append(self.system[ndx])
                    if len(bondlist): self.bondlist.append(bondlist[ndx])
                    if len(anglelist): self.anglelist.append(anglelist[ndx])
        self.fratio = 1.0 - len(tmpsys)/len(self.system)
        # alias
        self.system = tmpsys
        self.reflist.pop(len(self.reflist)-1)
        self.keepndxlist.pop(len(self.keepndxlist)-1)
        
        if self.obpar or self.oball:
            print('Note: calculating final bonds probability ...')
            self.prob_final['bpar'], self.prob_final['ball'] = self.calc_probs(self.bondlist,binc,self.obpar,self.oball)
        if self.oapar or self.oaall:
            print('Note: calculating final angles probability ...')
            self.prob_final['apar'], self.prob_final['aall'] = self.calc_probs(self.anglelist,ainc,self.oapar,self.oaall)

    def calc_probs(self,datalist,dt,opar=None,oall=None):
        """
        Inputs:
            dt : float

        Return:
            prob_par  : 3D : List[ [List[int],float] ]

                => the probability of any two molecules whose difference fall
                   within defined range

                explanation:
                    assume three molecules, the number of atoms are 4;
                    A = [a1, a2, a3, a4]
                    B = [b1, b2, b3, b4]
                    C = [c1, c2, c3, b4]

                    define its connection index is: [[0,1], [1,2], [1,3]]
                    means bonded atom pairs in A are [a1-a2, a2-a3, a2-a4]

                    now, calculate their bondlist;

                    DA = [Da01, Da12, Da13]
                    DB = [Db01, Db12, Db13]
                    DC = [Dc01, Dc12, Dc13]

                    thus, the difference of bond-length for any two molecules
                    corresponding to connection will be DA-DB, DA-DC, DB-DC

                    therefore a new array will be got:
                    New = [Da01-Db01, Da01-Dc01, Db01-Dc01]

                    finally we can calculate molecule-increments probability


            prob_all  : 2D : List[ List[int],  float ]

                => means probability on overall difference

                explanation:
                    continuing top, we have already calculated DA, DB, DC

                    calculate their average values;

                    TA = Average(DA)
                    TB = Average(DB)
                    TC = Average(DC)

                    following the same rule, calculate differences,
                    New = [TA-TB, TA-TC, TB-TC]

                    finally we can calculate overall probability

        Note:
            results are not normalized
        """
        def calc_list(ls,dt):
            stls = sorted(ls)
            rmin = stls[0]
            steps = (stls[-1] - rmin) / dt
            # for round-of-errors, as well as endpoint
            steps = int(steps) + 2
            # all on same values
            if steps <= 3:
                return [len(ls)],rmin
            prolist = []
            ndx = 0
            # for endpoint
            stls.append(stls[-1]+dt)
            for v in range(1,steps):
                high = rmin + v*dt
                for n,t in enumerate(stls[ndx:]):
                    if t >= high:
                        prolist.append(n)
                        ndx += n
                        break
            return prolist,rmin

        if len(datalist) <= 1: return [],[]

        prob_par = []
        if opar:
            print('    --> computing on par entry ...')
            for cnt in range(len(datalist[0])):
                ls = [t[cnt] for t in datalist]
                p = calc_list(ls,dt)
                prob_par.append(p)
        prob_all = []
        if oall:
            print('    --> computing on system ...')
            # Caution! only one movement
            #sub = [sum(i)/len(i) for i in datalist]
            sub = [sum(i) for i in datalist]
            prob_all = calc_list(sub,dt*len(datalist[0]))
        return prob_par, prob_all

    def calc_filterlists(self,bondlist,anglelist,binc,ainc,mode=None,vndx=None,
                        borandom=None,boall=None,keepndxlist=None):
        """
        Rule:
            Since the variance is only based on the single movement,
            either can be a bond or angle, so the changes on any two adjacent
            molecules are smaller than tolerance will be removed,
            thus, dt should not multiple the total numbers,
            because Sum(mol untouched atoms) =~ 0.0

        Return:
            reflist  :  List[int]  :  index of molecules waiting to be removed
        """
        if len(bondlist) <= 3 and len(anglelist) <= 3: return []
        bl = [sum(i) for i in bondlist]
        al = [sum(i) for i in anglelist]
        if mode is None or mode.lower() in ['dynamic','d']:
            if boall:
                return self._calc_filterlists_all_dynamic(bl,al,binc,ainc,keepndxlist)
            # dynamic-separate on bondlist
            reflist = self._calc_filterlists_sep_dynamic(bl,al,binc,ainc,keepndxlist)
            # dynamic-separate on anglelist
            reflist.append(-1)
            cnt = 0
            mybl = []
            myal = []
            ndxlist = []
            for i in range(len(bl)):
                if i == reflist[cnt]:
                    cnt += 1
                else:
                    mybl.append(bl[i])
                    myal.append(al[i])
                    ndxlist.append(i)
            reflist.pop(len(reflist)-1)
            tmplist = self._calc_filterlists_sep_dynamic(myal,mybl,ainc,binc,keepndxlist)
            reflist.extend([ndxlist[i] for i in tmplist])
            return sorted(reflist)
        return self._calc_filterlists_static(bl,al,binc,ainc,keepndxlist,vndx,borandom)

    def _calc_filterlists_all_dynamic(self,bl,al,binc,ainc,keepndxlist):
        """This part does is to recursively remove index, whose difference
           separately with two adjacent values is smaller than tolerance
        
        for example, we have inputs like;
        
          bal      = [0.0, 0.1, 0.15, 0.18, 0.19, 0.20, 0.3, 0.4, 0.43, 0.5]
        ndxlist    =   0    1    2     3     4     5     6    7    8     9
        
        inc = 0.1
        
        then we will know if we remove index in [2,3,4,8], the new array
        balnew = [0.0, 0.1, 0.20, 0.3, 0.4, 0.5] will meet the requirement

        sort by index
        Caution: for future debug, bal is not sorted"""
        # make a deep copy of keepndxlist
        if keepndxlist is None: keepndxlist = []
        keepndxlist = [i for i in keepndxlist]
        keepndxlist.append(-1)
        
        if not len(bl) or not len(al):
            if not len(bl):
                bal = al
                inc = ainc
            else:
                bal = bl
                inc = binc
        else:
            bal = [v+al[i] for i,v in enumerate(bl)]
            inc = binc + ainc
        nlist = sorted(range(len(bal)),key=lambda k: bal[k])
        mdel = [bal[j] - bal[nlist[i]] for i,j in enumerate(nlist[1:])]

        reflist = []
        ndx = 0
        cnt = 0
        while ndx < len(mdel):
            if mdel[ndx] < inc:
                v1 = nlist.pop(ndx+1)
                if v1 == keepndxlist[cnt]:
                    cnt += 1
                    if ndx >= len(mdel)-1: break
                    # be aware in here, nlist has been processed
                    v2 = nlist[ndx+1]
                    mdel.pop(ndx)
                    mdel[ndx] = bal[v2] - bal[v1]
                else:
                    reflist.append(v1)
                    if ndx >= len(mdel)-1: break
                    dt = mdel.pop(ndx+1)
                    mdel[ndx] += dt
            else:
                ndx += 1
        # reflist has to be sorted from smaller to bigger for following refinement
        return sorted(reflist)
    
    def _calc_filterlists_sep_dynamic(self,bl,al,binc,ainc,keepndxlist):
        if not len(bl) or not len(al):
            return self._calc_filterlists_all_dynamic(bl,al,binc,ainc,keepndxlist)
        
        # make a deep copy of keepndxlist
        if keepndxlist is None: keepndxlist = []
        keepndxlist = [i for i in keepndxlist]
        keepndxlist.append(-1)

        nlist = sorted(range(len(bl)),key=lambda k: bl[k])
        mdel = [bl[j] - bl[nlist[i]] for i,j in enumerate(nlist[1:])]
        reflist = []
        ndx = 0
        cnt = 0
        while ndx < len(mdel):
            if mdel[ndx] < binc:
                v0 = nlist[ndx]
                v1 = nlist[ndx+1]
                # be aware of negative value
                if abs(al[v1]-al[v0]) < ainc:
                    nlist.pop(ndx+1)
                    if v1 == keepndxlist[cnt]:
                        cnt += 1
                        if ndx >= len(mdel)-1: break
                        v2 = nlist[ndx+1]
                        mdel.pop(ndx)
                        mdel[ndx] = bl[v2] - bl[v1]
                    else:
                        reflist.append(v1)
                        if ndx >= len(mdel)-1: break
                        dt = mdel.pop(ndx+1)
                        mdel[ndx] += dt
                else:
                    ndx += 1
            else:
                ndx += 1
        return sorted(reflist)

    def _calc_filterlists_static(self,bl,al,binc,ainc,keepndxlist,vndx=None,borandom=None):
        # make a deep copy of keepndxlist
        if keepndxlist is None: keepndxlist = []
        keepndxlist = [i for i in keepndxlist]

        if not len(bl) or not len(al):
            if not len(bl):
                bal = al
                inc = ainc
            else:
                bal = bl
                inc = binc
        else:
            bal = [v+al[i] for i,v in enumerate(bl)]
            inc = binc + ainc

        nlist = sorted(range(len(bal)),key=lambda k: bal[k])
        vlist = [bal[i] for i in nlist]

        # always make vndx one-inc less than smallest value
        if vndx is None:
            vndx = vlist[0]
        elif vndx > vlist[0]:
            while vndx > vlist[0]:
                vndx -= inc
        else:
            while vndx < vlist[0]:
                vndx += inc
            vndx -= inc
        
        reflist = []
        cnt = 0
        n = int((vlist[-1]-vndx)/inc) + 1
        tot = len(vlist)
        for i in range(1,n):
            t = vndx + i*inc
            ls = []
            while cnt < tot:
                if vlist[cnt] >= t: break
                ls.append(nlist[cnt])
                cnt += 1
            if len(ls):
                bo = False
                # important, increase efficiency
                if len(keepndxlist):
                    for k in ls:
                        if k in keepndxlist:
                            bo = True
                            keepndxlist.remove(k)
                if bo:
                    reflist.extend(ls)
                elif len(ls) >= 2:
                    if borandom:
                        x = random.randrange(len(ls))
                        ls.pop(x)
                        reflist.extend(ls)
                    else:
                        reflist.extend(ls[1:])
        return sorted(reflist)

    def calc_square_distance(self,system,bcon):
        """
        Return:
            bondlist : 2D : System[Mol[l1,l2, ...], ...] : List[List[float]]
        """
        if not len(bcon): return []
        bondlist = []
        for mol in system:
            ls = []
            for ndx in bcon:
                at1 = mol[ndx[0]]
                at2 = mol[ndx[1]]
                dx = at1[1] - at2[1]
                dy = at1[2] - at2[2]
                dz = at1[3] - at2[3]
                tmp = dx*dx + dy*dy + dz*dz
                ls.append(tmp)
            bondlist.append(ls)
        return bondlist

    def calc_angle_degree(self,system,acon):
        """
        Rule:
            assume cooridnates,

            A (ax, ay, az)
            B (bx, by, bz)
            C (cx, cy, cz)

            Vector,

            BA = (ax-bx, ay-by, az-bz)
            BC = (cx-bx, cy-by, cz-bz)

            cos<ABC> = Sum(bai*bci) / dis(BA) * dis(BC)

            <ABC> = math.acos(sigma) * 180.0 / math.pi

        Return:
            anglelist : 2D : System[Mol[l1,l2, ...], ...] : List[List[float]]
        """
        if not len(acon): return []
        anglelist = []
        cvt = 180.0 / math.pi
        for mol in system:
            ls = []
            for ndx in acon:
                a = mol[ndx[0]]
                b = mol[ndx[1]]
                c = mol[ndx[2]]
                ba = [a[1]-b[1],a[2]-b[2],a[3]-b[3]]
                bc = [c[1]-b[1],c[2]-b[2],c[3]-b[3]]
                tot = ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]
                sub = sum([i*i for i in ba]) * sum([i*i for i in bc])
                rst = math.acos(tot/pow(sub,0.5)) * cvt
                ls.append(rst)
            anglelist.append(ls)
        return anglelist


class BulkProcess:
    """bulk process for datafilelist based on indexfilelist

    Inputs:
        bcon (str|list): conb, bcon, ...
        acon (str|list): cona, acon, ...

    Attributes:
        overall_system  :   final good system
        overall_sysbad  :   all bad system

        overall_prob_begin  :   based on read data file
        overall_prob_final  :   based on read data file

    Note:
        fragments is input as human-readable number, starting at 1
    """
    def __init__(self,datafilelist=None,indexfilelist=None,
                bool_force_double_check=None,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.mytime = time.time()
        self.datafilelist = []
        if datafilelist is not None:
            for f in datafilelist:
                if os.path.isfile(f):
                    self.datafilelist.append(f)
                else:
                    print('Warning: not a data file < {:} >, ignoring'.format(f))

        if not len(self.datafilelist):
            self.nice = False
            self.info = 'Fatal: no inputs'
            return

        self.indexfilelist = []
        if indexfilelist is not None:
            for f in indexfilelist:
                if os.path.isfile(f):
                    self.indexfilelist.append(f)
                else:
                    print('Warning: not an index file < {:} >, ignoring'.format(f))

        self.bool_force_double_check = False if bool_force_double_check is False else True
        self.args = args
        self.kwargs = kwargs

    def run(self,debug=None):
        fsize = sum([os.stat(i).st_size for i in self.datafilelist])
        fsize += sum([os.stat(i).st_size for i in self.indexfilelist])
        # set warning of maximum valid file size
        memmax = 500
        if fsize/1024/1024 > memmax:
            rsize = 0
            msize = 0
            fsize = 0
            total = 100000
            cnt = 0
            for file in [*self.datafilelist, *self.indexfilelist]:
                if not (file.endswith('.xyz') or file.endswith('.xsf') or file.endswith('.txt')):
                    continue
                fsize += os.stat(file).st_size
                if cnt > total: continue
                profile = []
                with open(file,mode='rt') as f:
                    while True:
                        line = f.readline()
                        if not len(line): break
                        profile.append(line)
                        if cnt > total: break
                        cnt += 1
                ftmp = tempfile.NamedTemporaryFile()
                ftmp.write(''.join(profile).encode('utf-8'))
                ftmp.flush()
                rf = ReadFile(ftmp.name,ext=file[file.rfind('.')+1:],debug=False)
                if len(profile) < 3000: rf.nice = False
                if rf.nice:
                    rf.run()
                    if len(rf.system):
                        rsize += os.stat(ftmp.name).st_size
                        msize += getrealsizeof(rf.system) + getrealsizeof(rf.energy)
                ftmp.close()
            
            if fsize/1024/1024 > memmax and rsize > 1.0:
                # convert to MB, get the 1.2 times memory
                fsize = fsize / rsize * msize / 1000 / 1000 * 1.2
                print('Warning: your input is super large')
                print('Warning: memory requested roughly will be: {:} MB'.format(round(fsize,2)))
                print('Be sure that you want to continue? y/yes, else not. Input: ',end='')
                if input().lower() not in ['y','yes']:
                    print('Note: you decided to quit, nothing will be processed')
                    return
                print()

        systemlist,energylist = self.get_datalist(self.datafilelist)
        if not sum([len(i) for i in systemlist]):
            self.nice = False
            self.info = 'Fatal: no inputs after process'
            return
        self.molnms = [len(i) for i in systemlist]

        sysndxlist = []
        if len(self.indexfilelist):
            sysndxlist,tmp = self.get_datalist(self.indexfilelist)

        # connections only need to be calculated once
        choose = [i for i in systemlist if len(i)]
        self.get_connections(choose[0][0])
        if not self.nice: return

        # to make cross filtration happen, sysndxlist should be at the first
        allsystem = []
        allenergy = []
        for i in sysndxlist:
            allenergy.extend([None for j in range(len(i))])
            allsystem.extend(i)
        allkeeps = list(range(len(allsystem)))
        for i in systemlist: allsystem.extend(i)
        for i in energylist: allenergy.extend(i)
        mf = Filtration(system=allsystem,keepndxlist=allkeeps,*self.args,**self.kwargs)

        # prompt for double check
        if self.bool_force_double_check:
            print('\nCheck: current work path:')
            print('   => {:}'.format(os.path.abspath('.')))
            print('\nCheck: data files:')
            for cnt,fd in enumerate(self.datafilelist):
                print('   => {:} -- molnms {:}'.format(fd,self.molnms[cnt]))

            if len(self.indexfilelist):
                print('Check: index files:')
                for cnt,fd in enumerate(self.indexfilelist):
                    print('   => {:} -- molnms {:}'.format(fd,len(sysndxlist[cnt])))

            print('Check: molecule fragments:')
            lt = []
            for i in self.kwargs['fragments']: lt.append([j+1 for j in i])
            print('   => {:}'.format(lt))

            print('Check: bond connection:')
            lt = []
            for i in self.kwargs['bcon']: lt.append([j+1 for j in i])
            print('   => {:}'.format(lt))

            print('Check: angle connection:')
            lt = []
            for i in self.kwargs['acon']: lt.append([j+1 for j in i])
            print('   => {:}'.format(lt))

            print('Check: total inputs < {:} >'.format(sum(self.molnms)))

            stmp = 'ON' if mf.oball else 'OFF'
            print('Check: (image) bonds all probability < {:} >'.format(stmp))
            stmp = 'ON' if mf.obpar else 'OFF'
            print('Check: (images) bonds par probability < {:} > (time consuming)'.format(stmp))
            stmp = 'ON' if mf.oaall else 'OFF'
            print('Check: (image) angles all probability < {:} >'.format(stmp))
            stmp = 'ON' if mf.oapar else 'OFF'
            print('Check: (images) angles par probability < {:} > (time consuming)'.format(stmp))
            print('Check: bonds tolerance < {:} Angstrom >'.format(mf.btol))
            print('Check: angles tolerance < {:} degree >'.format(mf.atol))
            if mf.mode == 'dynamic':
                if mf.boall:
                    print('Check: calculation type: < dynamic/all >')
                else:
                    print('Check: calculation type: < dynamic/separate >')
            else:
                if mf.borandom:
                    print('Check: calculation type: < static/random >')
                else:
                    print('Check: calculation type: < static/lowest-bit >')
                if mf.vndx: print('Check: calculation vndx: < {:} >'.format(mf.vndx))
            imtot = 0
            if mf.oball: imtot += 1
            if mf.oaall: imtot += 1
            if mf.obpar: imtot += len(mf.bcon)
            if mf.oapar: imtot += len(mf.acon)
            print('Check: number of images will be generated: < {:} >'.format(imtot))

            print('\nDo you want to continue? y/yes, else not. Input: ',end='')
            if input().lower() not in ['y','yes']:
                print('Note: you decided to quit, nothing will be processed')
                return
            print()

        if debug: return allsystem
        mf.run()
        self.seed = mf.seed
        self.mode = mf.mode
        self.boall = mf.boall
        self.vndx = mf.vndx
        self.borandom = mf.borandom

        # accumulation only on datafilelist
        tot = 0
        acclist = []
        for i in sysndxlist: tot += len(i)
        acclist.append(tot)
        for i in systemlist:
            tot += len(i)
            acclist.append(tot)

        # test
        #mf.reflist = [i for i in acclist[:-1]]
        # expect: rmlist = [1 for i in range(len(acclist)-1)]
        # mf.reflist is the sorted list
        self.rmnmlist = [0 for i in range(len(acclist)-1)]
        ndx = 1
        i = 0
        while ndx < len(acclist):
            cnt = 0
            while i < len(mf.reflist):
                if mf.reflist[i] < acclist[ndx]:
                    cnt += 1
                    i += 1
                else:
                    break
            if cnt == 0:
                ndx += 1
            else:
                self.rmnmlist[ndx-1] = cnt

        self.overall_energy = []
        self.overall_system = []
        totreflist = sorted([*allkeeps,*mf.reflist])
        totreflist.append(-1)
        ndx = 0
        for i,v in enumerate(allenergy):
            if i == totreflist[ndx]:
                ndx += 1
            else:
                self.overall_energy.append(v)
                self.overall_system.append(allsystem[i])
        self.bcon = mf.bcon
        self.acon = mf.acon
        self.btol = mf.btol
        self.atol = mf.atol
        self.boim = True if mf.oball or mf.obpar or mf.oaall or mf.oapar else False
        self.overall_prob_begin = mf.prob_begin
        self.overall_prob_final = mf.prob_final
        self.save_files()

    def save_files(self):
        print('\nNote: saving bulk process results ...')
        tot = len(self.overall_system)
        print('Note: total molnms: < {:} >'.format(sum(self.molnms)))
        print('Note: final molnms: < {:} >'.format(tot))
        ratio = ('%f' % (1-round(tot/sum(self.molnms),6))).rstrip('0').rstrip('.')
        print('Note: filtration ratio: < {:} >'.format(ratio))

        # files
        fd = SaveFile(self.overall_system,*self.args,**self.kwargs)
        self.kwargs['fname'] = file_gen_new(fd.fname,fextend=fd.ftype)
        self.kwargs['energy'] = self.overall_energy
        fd = SaveFile(self.overall_system,*self.args,**self.kwargs)
        fd.run()
        outfile = fd.fname
        print('Note: file is saved to < {:} >'.format(outfile))

        filedict = {}

        # images
        if self.boim:
            filedict['probability data file'] = self.save_probdata()

            if len(self.overall_prob_begin['ball']):
                fgp = file_gen_new('bonds-all',fextend='png',foriginal=False)
                fbo = plot_save_image(
                    self.overall_prob_begin['ball'],
                    self.overall_prob_final['ball'],
                    dt=self.btol,
                    fname=fgp,
                    key='bonds',
                )
                if fbo: filedict['image all bonds filtration file'] = fgp

            if len(self.overall_prob_begin['aall']):
                fgp = file_gen_new('angles-all',fextend='png',foriginal=False)
                fbo = plot_save_image(
                    self.overall_prob_begin['aall'],
                    self.overall_prob_final['aall'],
                    dt=self.atol,
                    fname=fgp,
                    key='angles',
                )
                if fbo: filedict['image all angles filtration file'] = fgp

            filedict['image bonds par filtration file'] = []
            for i,t in enumerate(self.overall_prob_begin['bpar']):
                mark = 'bonds-par-{:}+{:}'.format(*self.bcon[i])
                fgp = file_gen_new(mark,fextend='png',foriginal=False)
                fbo = plot_save_image(
                    t,
                    self.overall_prob_final['bpar'][i],
                    dt=self.btol,
                    fname=fgp,
                    key='bonds',
                )
                if fbo: filedict['image bonds par filtration file'].append(fgp)

            filedict['image angles par filtration file'] = []
            for i,t in enumerate(self.overall_prob_begin['apar']):
                mark = 'angles-par-{:}+{:}+{:}'.format(*self.acon[i])
                fgp = file_gen_new(mark,fextend='png',foriginal=False)
                fbo = plot_save_image(
                    t,
                    self.overall_prob_final['apar'][i],
                    dt=self.btol,
                    fname=fgp,
                    key='angles',
                )
                if fbo: filedict['image angles par filtration file'].append(fgp)

        ftot = file_gen_new('bulk-process-info')
        print('Note: please check summary file for more info: < {:} >'.format(ftot))
        with open(ftot,'wt') as f:
            f.write('Note: starting  time: {:}\n'.format(time.ctime(self.mytime)))
            now = time.time()
            f.write('Note: finishing time: {:}\n'.format(time.ctime(now)))
            minutes = round((now-self.mytime)/60.0,2)
            f.write('Note: execution time: {:} minutes\n'.format(minutes))
            f.write('Note: current work path:\n')
            f.write('  => {:}\n'.format(os.path.abspath('.')))
            f.write('Note: random seed: {:}\n'.format(self.seed))
            f.write('Note: bulk process for input files:\n')
            for i,fd in enumerate(self.datafilelist):
                f.write('  => {:} -- molnms {:} => remove {:}\n'.format(fd,self.molnms[i],self.rmnmlist[i]))
            f.write('Note: total number of inputs: {:}\n'.format(sum(self.molnms)))
            if len(self.indexfilelist):
                f.write('\nNote: index files:\n')
                for fd in self.indexfilelist:
                    f.write('  => {:}\n'.format(fd))
            if self.mode is None or self.mode == 'dynamic':
                f.write('\nNote: filtration mode is: dynamic\n')
                if self.boall:
                    f.write('  => calculation is performed for all entries\n')
                else:
                    f.write('  => calculation is performed separately\n')
            else:
                f.write('\nNote: filtration mode is: static\n')
                if self.vndx: f.write('  => index value is: {:}\n'.format(self.vndx))
                if self.borandom:
                    f.write('  => randomly filtration')
                else:
                    f.write('  => lowest-bit filtration')
            f.write('\nNote: result file:\n')
            f.write('  => {:} -- molnms {:}\n'.format(outfile,len(self.overall_system)))
            f.write('\nNote: filtration ratio: {:}\n'.format(ratio))
            f.write('\nNote: index of connections start at 1\n')
            f.write('\nNote: bonds connections:\n')
            tot = ''
            out = '  => '
            for tmp in self.bcon:
                tmp = [i+1 for i in tmp]
                out += '{:}, '.format(tmp)
                if len(out) >= 80:
                    tot += out.rstrip() + '\n'
                    out = '  => '
            if out != '  => ': tot += out.rstrip() + '\n'
            tot += '\n'
            f.write(tot)

            f.write('Note: angles connections:\n')
            tot = ''
            out = '  => '
            for tmp in self.acon:
                tmp = [i+1 for i in tmp]
                out += '{:}, '.format(tmp)
                if len(out) >= 80:
                    tot += out.rstrip() + '\n'
                    out = '  => '
            if out != '  => ': tot += out.rstrip() + '\n'
            tot += '\n'
            f.write(tot)

            f.write('Note: btol   : {:} Angstrom\n'.format(self.btol))
            f.write('Note: atol   : {:} degree\n'.format(self.atol))

            for k,v in filedict.items():
                if isinstance(v,str):
                    f.write('Note: {:<40}:   {:}\n'.format(k,v))
                elif len(v):
                    f.write('Note: {:}:\n'.format(k))
                    for i,j in enumerate(v):
                        f.write('  ==> {:>3}: {:}\n'.format(i+1,j))

    def save_probdata(self):
        def gen_outputs(prob,bcon,acon,key):
            def fout(pdata):
                txt = ''
                out = '{:}\n'.format(pdata[1])
                for p in pdata[0]:
                    txt += '{:}  '.format(p)
                    if len(txt) >= 80:
                        out += txt.strip() + '\n'
                        txt = ''
                out += txt.strip() + '\n\n'
                return out

            k = key.upper()
            contents = f'@{k}  BALL\n' + fout(prob['ball']) if len(prob['ball']) else ''
            if len(prob['bpar']):
                for n,data in enumerate(prob['bpar']):
                    contents += '@{:}  BPAR   {:}  {:}\n'.format(k,*bcon[n])
                    contents += fout(data)
            if len(prob['aall']): contents += f'@{k}  AALL\n' + fout(prob['aall'])
            if len(prob['apar']):
                for n,data in enumerate(prob['apar']):
                    contents += '@{:}  APAR   {:}  {:}  {:}\n'.format(k,*acon[n])
                    contents += fout(data)
            return contents

        fdata = file_gen_new('bulk-probability-data')
        print('Note: probability data is saved to < {:} >'.format(fdata))
        with open(fdata,'wt') as f:
            f.write('# Note: filtration process probability data\n\n')
            for tmp in self.datafilelist:
                f.write('@FILE {:}\n'.format(tmp))
            f.write('@BTOL   {:}\n@ATOL   {:}\n\n\n'.format(self.btol, self.atol))
            mbcon = [[i+1 for i in j] for j in self.bcon]
            macon = [[i+1 for i in j] for j in self.acon]
            f.write(gen_outputs(self.overall_prob_begin,mbcon,macon,'begin'))
            f.write(gen_outputs(self.overall_prob_final,mbcon,macon,'final'))
        return fdata

    def get_datalist(self,filelist):
        """return 4D list"""
        datalist = []
        energylist = []
        for f in filelist:
            rf = ReadFile(f)
            if rf.nice:
                rf.run()
                print('Note: for file < {:} >, number of inputs < {:} >'.format(f,len(rf.system)))
            else:
                print(rf.info)
            datalist.append(rf.system)
            energylist.append(rf.energy)
        return datalist,energylist

    def get_connections(self,system):
        """
        Rule:
            now, we think user input is always on the first priority,
            and bcon & acon at the latter place

            connection is got in sequence:
            1) if bcon or acon is set
            2) if fragments is set
            3) default: bcon=AnglePerception.nconb, acon=AnglePerception.ncona
        """
        fn = AnglePerception(system)
        if not fn.nice:
            self.nice = False
            self.info = fn.info
            return
        fn.run()

        if 'userinputs' in self.kwargs and self.kwargs['userinputs'] is True:
            self.kwargs['userinputs'] = True
        else:
            self.kwargs['userinputs'] = False

        bog = False
        if 'fragments' in self.kwargs and self.kwargs['fragments'] is not None:
            bog = True
            if not len(self.kwargs['fragments']): self.kwargs['fragments'] = fn.fragments
            if self.kwargs['userinputs']:
                fg = [[j-1 for j in i] for i in self.kwargs['fragments']]
                if min([min(i) for i in fg]) < 0:
                    self.nice = False
                    self.info = 'Fatal: wrong defined: atom index should start at 1'
                    return
                self.kwargs['fragments'] = fg
            else:
                fg = self.kwargs['fragments']
            if max([max(i) for i in fg]) >= len(system):
                self.nice = False
                self.info = 'Fatal: wrong defined: atom index exceeding: < {:} >'.format(len(system))
                return
            if len(set([j for i in fg for j in i])) != sum([len(i) for i in fg]):
                self.nice = False
                self.info = 'Fatal: wrong defined: fragments: has repeats'
                return
            fbcon, facon = self.func_calc_connections(self.kwargs['fragments'])
        else:
            self.kwargs['fragments'] = fn.fragments

        # take care of special case, when input is None
        if 'bcon' in self.kwargs and self.kwargs['bcon'] is not None:
            tmpbcon = self.kwargs['bcon']
            if isinstance(self.kwargs['bcon'],list):
                self.kwargs['bcon'] = self.check_user_input_connections(
                    self.kwargs['bcon'],
                    fn.fragments,
                    self.kwargs['userinputs']
                )
            elif isinstance(self.kwargs['bcon'],str):
                contmp = tmpbcon.lower()
                if 'no' in contmp or 'non' in contmp or 'none' in contmp:
                    self.kwargs['bcon'] = []
                elif hasattr(fn,contmp) and 'b' in contmp and 'con' in contmp:
                    self.kwargs['bcon'] = getattr(fn,contmp)
                else:
                    self.kwargs['bcon'] = False
            else:
                self.kwargs['bcon'] = False
            if self.kwargs['bcon'] is False:
                self.nice = False
                self.info += '\nFatal: wrong defined: bcon: {:}'.format(tmpbcon)
                return
        elif bog:
            self.kwargs['bcon'] = fbcon
        else:
            if fn.fragments is None or len(fn.fragments) == 1:
                self.kwargs['bcon'] = fn.nconb
            else:
                self.kwargs['bcon'] = fn.fnconb

        if 'acon' in self.kwargs and self.kwargs['acon'] is not None:
            if isinstance(self.kwargs['acon'],list):
                self.kwargs['acon'] = self.check_user_input_connections(
                    self.kwargs['acon'],
                    fn.fragments,
                    self.kwargs['userinputs']
                )
            elif isinstance(self.kwargs['acon'],str):
                tmpacon = self.kwargs['acon']
                contmp = tmpacon.lower()
                if 'no' in contmp or 'non' in contmp or 'none' in contmp:
                    self.kwargs['acon'] = []
                elif hasattr(fn,contmp) and 'a' in contmp and 'con' in contmp:
                    self.kwargs['acon'] = getattr(fn,contmp)
                else:
                    self.kwargs['acon'] = False
            else:
                tmpacon = self.kwargs['acon']
                self.kwargs['acon'] = False
            if self.kwargs['acon'] is False:
                self.nice = False
                self.info = 'Fatal: wrong defined: acon: {:}'.format(tmpacon)
                return
        elif bog:
            self.kwargs['acon'] = facon
        else:
            if fn.fragments is None or len(fn.fragments) == 1:
                self.kwargs['acon'] = fn.ncona
            else:
                self.kwargs['acon'] = fn.fncona

        # always set final userinputs to False
        self.kwargs['userinputs'] = False

    def check_user_input_connections(self,ul,fl,bo=None):
        if not len(ul): return []
        offset = 1 if bo is True else 0
        for ndx in ul:
            if len(set(ndx)) != len(ndx):
                self.nice = False
                self.info = 'Fatal: wrong defined: repeats at: < {:} >'.format(ndx)
                return False
        ux = max([max(i) for i in ul])
        fx = max([max(i) for i in fl])
        if ux-offset > fx:
            self.nice = False
            self.info = 'Fatal: wrong defined: atom index exceeding: < {:} >'.format(ux)
            return False
        # convert to python readable number, starts at 0
        ptmp = []
        for i in ul: ptmp.append([j-offset for j in i])
        for i in ptmp:
            if min(i) < 0:
                self.nice = False
                self.info = 'Fatal: atom index should start at 1'
                return False
        return ptmp

    def func_calc_connections(self,fragments):
        """Bond & Angle connection based on fragments

        Input:
            fragments  :  2D  :  List[ List[int] ]

        Note:
            Based on first 2 atoms, which is chosen by sequence,
            find bond & angle connection with all atoms in other fragments

        Return:
            bcon  :  2D  :  List[ List[int,int], ...]
            acon  :  2D  :  List[ List[int,int,int], ...]
        """
        if len(fragments) <= 1: return [],[]
        bcon = []
        b1 = fragments[0][0]
        for ndx in fragments[1:]:
            for v in ndx:
                bcon.append([b1,v])
        acon = []
        for ndx,mol in enumerate(fragments):
            if len(mol) >= 2:
                break
        if len(mol) >= 2:
            at1 = mol[0]
            at2 = mol[1]
            for cnt,mol in enumerate(fragments):
                if cnt != ndx:
                    for v in mol:
                        acon.append([at1,at2,v])
        else:
            if len(fragments) >= 3:
                at1 = fragments[0][0]
                at2 = fragments[1][0]
                for v in fragments[2:]:
                    acon.append([at1,at2,v[0]])
        return bcon,acon


class PlotSamples(BulkProcess):
    def __init__(self,probdatafilelist=None,nmsamples=None,nmlist=None,
                startndx=None,endndx=None,incndx=None,nmranges=None,
                *args,**kwargs):
        seed = kwargs['seed'] if 'seed' in kwargs else None
        self.seed = seed if seed else random.randrange(100000000)
        random.seed(self.seed)
        kwargs['seed'] = self.seed
        super().__init__(*args,**kwargs)
        if not len(self.datafilelist):
            self.nice = True
            self.info = ''

        self.probdatafilelist = []
        if probdatafilelist is not None:
            for f in probdatafilelist:
                if os.path.isfile(f):
                    self.probdatafilelist.append(f)
                else:
                    print('Warning: not probability data file < {:} >, ignoring'.format(f))

        if not len(self.datafilelist) and not len(self.probdatafilelist):
            self.nice = True
            self.info = 'Fatal: no valid inputs: datafilelist/probdatafilelist'
            return

        if 'userinputs' in kwargs and kwargs['userinputs']:
            if startndx is not None: startndx -= 1
            if endndx is not None: endndx -= 1
        
        # assume data has been properly processed, index starts at 0
        self.nmsamples = nmsamples
        self.startndx = startndx
        self.endndx = endndx
        self.incndx = incndx
        self.nmranges = nmranges
        self.nmlist = nmlist

        self.choices = []
        # after debug run, everything is ready
        allsystems = []
        if self.nice and len(self.datafilelist):
            allsystems = super().run(debug=True)
            if allsystems is None: allsystems = []
            tot = len(allsystems)
        if not len(self.datafilelist) or not self.nice:
            pass
        elif nmlist is None:
            pass
        elif isinstance(nmlist,list):
            for i in nmlist:
                if not isinstance(i,int):
                    self.nice = False
                    self.info = 'Fatal: wrong defined: not a number: {:}'.format(i)
                    return
                if i > tot:
                    self.nice = False
                    self.info = 'Fatal: wrong defined: too large: {:}'.format(i)
                    return
        else:
            self.nice = False
            if nmlist is None: self.info = 'Fatal: wrong defined: nmlist'
        # nmlist is in the highest priority
        bo = True
        if len(self.datafilelist) and self.nice and nmlist is not None:
            bo = False
            for i in nmlist:
                self.choices.append(random.sample(allsystems,i))
        if len(self.datafilelist) and bo and self.nice and 'datafilelist' in kwargs:
            if tot <= 20:
                self.info = 'Warning: too few inputs: datafilelist'
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if endndx is not None and endndx > tot:
                self.info = 'Fatal: too large: endndx --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if startndx is not None and startndx > tot:
                self.info = 'Fatal: too large: startndx --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if incndx is not None and incndx > tot:
                self.info = 'Fatal: too large: incndx --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if nmranges is not None and nmranges > tot:
                self.info = 'Fatal: too large: nmranges --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if nmsamples is not None and nmsamples*8 > tot:
                self.info = 'Fatal: too large: nmsamples --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            # alias
            samples = allsystems
            if startndx is not None: samples = samples[startndx:]
            if endndx is not None: samples = samples[:endndx]
            if nmranges is None: nmranges = 0
            if nmsamples is None: nmsamples = 5
            tot = len(samples)
            if incndx is None:
                t = tot // nmsamples
                for i in range(nmsamples):
                    self.choices.append(random.sample(samples,t))
            else:
                dt = random.randrange(nmranges)
                while dt < tot:
                    while True:
                        tmp = random.randrange(nmranges) * random.choice([-1,0,1])
                        if incndx + tmp > 0: break
                    if dt + incndx + tmp > 0: break
                    self.choices.append(samples[dt:dt+incndx+tmp])
                    dt += incndx + tmp
            # to avoid waste of time,
            # only number of samples bigger than 10 will be kept
            self.choices = [i for i in self.choices if len(i) > 10]

    def run(self):
        bondsdict = {'all':[], }
        anglesdict = {'all':[], }
        if len(self.probdatafilelist):
            print('Note: generate plots on probability data files ...')
            for f in self.probdatafilelist:
                print('Note: reading data from file: < {:} >'.format(f))
                bonds, angles = self.read_probdatafile(f)
                for k in bonds:
                    if k not in bondsdict: bondsdict[k] = []
                    bondsdict[k].append(bonds[k])
                for k in angles:
                    if k not in anglesdict: anglesdict[k] = []
                    anglesdict[k].append(angles[k])

        if len(self.choices):
            overall_prob_begin = []
            overall_prob_final = []
            overall_dt = []
            bconlist = []
            aconlist = []
            for cnt,ms in enumerate(self.choices):
                print('Note: generate plots on < {:} > sample files ...'.format(cnt+1))
                mf = Filtration(system=ms,*self.args,**self.kwargs)
                mf.run()
                bconlist.append(mf.bcon)
                aconlist.append(mf.acon)
                overall_prob_begin.append(mf.prob_begin)
                overall_prob_final.append(mf.prob_final)
                overall_dt.append([mf.btol, mf.atol])
            
            for i,di in enumerate(overall_prob_begin):
                df = overall_prob_final[i]

                if len(di['ball']):
                    p = {}
                    p['begin'] = di['ball']
                    p['final'] = df['ball']
                    p['dt'] = overall_dt[i][0]
                    bondsdict['all'].append(p)
                
                if len(di['aall']):
                    p = {}
                    p['begin'] = di['aall']
                    p['final'] = df['aall']
                    p['dt'] = overall_dt[i][1]
                    anglesdict['all'].append(p)
                
                if bconlist[i] is not None and not len(bconlist[i]) and not len(di['bpar']):
                    for j,k in enumerate(bconlist[i]):
                        if k[0] > k[1]: k[0],k[1] = k[1],k[0]
                        key = '{:}-{:}'.format(k[0],k[1])
                        if key not in bondsdict: bondsdict[key] = []
                        if len(di['bpar'][j]):
                            p = {}
                            p['begin'] = di['bpar'][j]
                            p['final'] = df['bpar'][j]
                            p['dt'] = overall_dt[i][0]
                            bondsdict[key].append(p)

                if aconlist[i] is not None and not len(aconlist[i]) and not len(di['apar']):
                    for j,k in enumerate(aconlist[i]):
                        if k[0] > k[2]: k[0],k[2] = k[2],k[0]
                        key = '{:}-{:}-{:}'.format(k[0],k[1],k[2])
                        if key not in anglesdict: anglesdict[key] = []
                        if len(di['apar'][j]):
                            p = {}
                            p['begin'] = di['apar'][j]
                            p['final'] = df['apar'][j]
                            p['dt'] = overall_dt[i][1]
                            anglesdict[key].append(p)
        # more info
        print('\nNote: random seed: {:}'.format(self.seed))
        if len(self.choices):
            if 'mode' in self.kwargs and self.kwargs['mode'] is not None:
                mode = self.kwargs['mode'].lower()
                mode = 'dynamic' if mode in ['d','dynamic'] else 'static'
            else:
                mode = 'dynamic'
            print('Note: filtration mode is: {:}'.format(mode))
            if mode == 'dynamic':
                boall = self.kwargs['boall'] if 'boall' in self.kwargs else None
                boall = True if boall or boall is None else False
                if boall:
                    print('  => calculation is performed for all entries')
                else:
                    print('  => calculation is performed separately')
            else:
                vndx = self.kwargs['vndx'] if 'vndx' in self.kwargs else None
                if vndx: print('  => index value is: {:}'.format(vndx))
                borandom = self.kwargs['borandom'] if 'borandom' in self.kwargs else None
                if borandom:
                    print('  => randomly filtration')
                else:
                    print('  => lowest-bit filtration')
        for k in bondsdict:
            self.save_image_samples(bondsdict[k],label='bonds')
        for k in anglesdict:
            self.save_image_samples(anglesdict[k],label='angles')

    class CYCLE:
        """infinite cycle works similar like itertools.cycle()"""
        def __init__(self,data):
            self.data = data
            self.nm = 0
        def __iter__(self):
            return self
        def __next__(self):
            if self.nm >= len(self.data): self.nm = 0
            v = self.data[self.nm]
            self.nm += 1
            return v
        def next(self):
            return self.__next__()

    def save_image_samples(self,datalist,label=None,fname=None):
        """
        Input:
            datalist: List[dict]: key in dict begin|final|dt
        """
        # the number of linestyle should be in ODD number
        linestyle = self.CYCLE(['solid','dotted','dashdot',])
        colors = self.CYCLE(['b','r','m','g','y','brown','palegreen','deepskyblue'])
        if label is None or label.lower() == 'bonds':
            label = 'bonds'
        else:
            label = 'angles'
        if fname is None: fname = 'bulk-image-samples-{:}'.format(label)
        fname = file_gen_new(fname,fextend='png')
        info = 'Conformers filtration on {:}'.format(label.capitalize())

        prodatalist = []
        molnms = []
        for data in datalist:
            if not isinstance(data,dict): continue
            if not ('begin' in data and 'final' in data and 'dt' in data): continue
            if len(data['begin'][0]) < 5 or len(data['final'][0]) < 5: continue
            prodatalist.append(data)
            molnms.append(sum(data['begin'][0]))
        if not len(prodatalist): return

        # sort them to make plot more tight
        reflist = sorted(range(len(molnms)), key=lambda k: molnms[k])
        prodatalist = [prodatalist[i] for i in reflist]

        # space [1,1] for plot, space [1,2] for legend
        fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [4,1]})
        fig.set_figheight(6)
        fig.set_figwidth(10)

        lines = []
        for d in prodatalist:
            color = colors.next()
            ils = linestyle.next()
            fls = linestyle.next()
            itxt = 'Begin molnms = {:}'.format(sum(d['begin'][0]))
            ftxt = 'Final molnms = {:}'.format(sum(d['final'][0]))
            ix = [d['begin'][1]+d['dt']*i for i in range(len(d['begin'][0]))]
            fx = [d['final'][1]+d['dt']*i for i in range(len(d['final'][0]))]
            iln, = ax1.plot(ix,d['begin'][0],color=color,linestyle=ils,label=itxt)
            fln, = ax1.plot(fx,d['final'][0],color=color,linestyle=fls,label=ftxt)
            lines.append(iln)
            lines.append(fln)

        ax1.set_title(info)
        ax2.axis('off')
        ax2.legend(handles=lines,loc='center')
        plt.tight_layout()
        print('Note: image file is saved to < {:} >'.format(fname))
        plt.savefig(fname)
        plt.close()
        return fname

    def read_probdatafile(self,file):
        def getdata(text,key=None):
            ltmp = text.split()
            if not len(ltmp): return []
            if key is not None and key >= len(ltmp): return []
            value = None
            data = []
            for i,v in enumerate(ltmp):
                if i == key:
                    try:
                        value = float(v)
                    except ValueError:
                        return []
                    continue
                else:
                    try:
                        v = int(v)
                    except ValueError:
                        return []
                data.append(v)
            return [data,value]

        atol = None
        btol = None
        ball = {'begin':'', 'final':''}
        bpar = {'begin':[], 'final':[]}
        aall = {'begin':'', 'final':''}
        apar = {'begin':[], 'final':[]}
        with open(file,'rt') as f: profile = f.readlines()
        cnt = -1
        while True:
            cnt += 1
            if cnt >= len(profile): break
            ltmp = profile[cnt].split()
            if len(ltmp) < 2 or ltmp[0][0] == '#': continue
            key = ltmp[0].lower()
            if key == '@btol' or key == '@atol':
                try:
                    if key == '@btol':
                        btol = float(ltmp[1])
                    else:
                        atol = float(ltmp[1])
                except ValueError:
                    continue
            elif key not in ['@begin','@final']:
                continue
            mark = ltmp[1].lower()
            if mark in ['ball','aall','bpar','apar']:
                txt = ''
                if mark == 'bpar':
                    if len(ltmp) < 4: continue
                    txt += ltmp[2] + '  ' + ltmp[3]
                if mark == 'apar':
                    if len(ltmp) < 5: continue
                    txt += ltmp[2] + '  ' + ltmp[3] + '  ' + ltmp[4]
                while True:
                    cnt += 1
                    if cnt >= len(profile): break
                    line = profile[cnt].strip()
                    if not len(line): continue
                    if line[0] == '@':
                        # caution: here is important
                        cnt -= 1
                        break
                    txt += '  ' + line
                if key == '@begin':
                    if mark == 'ball':
                        ball['begin'] = txt
                    elif mark == 'aall':
                        aall['begin'] = txt
                    elif mark == 'bpar':
                        bpar['begin'].append(txt)
                    else:
                        apar['begin'].append(txt)
                else:
                    if mark == 'ball':
                        ball['final'] = txt
                    elif mark == 'aall':
                        aall['final'] = txt
                    elif mark == 'bpar':
                        bpar['final'].append(txt)
                    else:
                        apar['final'].append(txt)
        if btol is None: btol = 0.1
        if atol is None: atol = 0.1
        # process for plot
        # bonds format:
        #   b1-b2 (b1<b2):   begin: List[ List[int], float]
        #                    final: List[ List[int], float]
        #                    dt   : float
        # angles format:
        #   a1-a-a2 (a1<a2): begin: List[ List[int], float]
        #                    final: List[ List[int], float]
        #                    dt   : float
        # special key:  'all'
        bonds = {'all':{}, }
        angles = {'all':{}, }

        bonds['all']['begin'] = getdata(ball['begin'],0)
        bonds['all']['final'] = getdata(ball['final'],0)
        bonds['all']['dt'] = btol
        angles['all']['begin'] = getdata(aall['begin'],0)
        angles['all']['final'] = getdata(aall['final'],0)
        angles['all']['dt'] = atol

        for i in bpar:
            for j in bpar[i]:
                data = getdata(j,2)
                if not len(data): continue
                bi,bj = data[0][0], data[0][1]
                if bi > bj: bj, bi = bi, bj
                key = '{:}-{:}'.format(bi,bj)
                data[0] = data[0][2:]
                if key not in bonds: bonds[key] = {}
                bonds[key][i] = data
                bonds[key]['dt'] = btol

        for i in apar:
            for j in apar[i]:
                data = getdata(j,3)
                if not len(data): continue
                ai,aj = data[0][0], data[0][2]
                if ai > aj: aj, ai = ai, aj
                key = '{:}-{:}-{:}'.format(ai,data[0][1],aj)
                data[0] = data[0][3:]
                if key not in angles: angles[key] = {}
                angles[key][i] = data
                angles[key]['dt'] = atol

        return bonds, angles


