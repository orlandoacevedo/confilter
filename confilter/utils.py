import os
import sys
import matplotlib.pyplot as plt


def file_gen_new(fname,fextend='txt',foriginal=True,bool_dot=True):
    """Generate new file name without overwritings

    Args:
        fname   (str)   :   input file fname
        fextend (str)   :   file extension
        foriginal (bool):   whether keep original
        bool_dot (bool) :   force check dot convention or not

    Returns:
        str     :   new file name
    """
    filename = fname
    pos = filename.rfind('.')
    if bool_dot and pos != -1:
        fname = filename[:pos]
        fextend = filename[pos:]
    else:
        fextend = '.' + fextend

    if foriginal is True:
        if not os.path.isfile(fname+fextend):
            return fname+fextend

    i = 1
    filename = fname
    while True:
        fname = filename + '-' + str(i) + fextend
        if not os.path.isfile(fname): break
        i += 1
    return fname


def plot_save_image(ini,fin=None,dt=None,fname=None,key=None):
    """
    Inputs:
        ini     :   2D  :   List [ List[int, ...],  float]
        fin     :   2D  :   List [ List[int, ...],  float]
        dt      :   float   :   increments, optional
        fname   :   str :   warning, overwritten may happen
        key     :   str :   {bonds, angles}, specify plot type

    Return:
        True if file is successfully generated, otherwise, False
    """
    if fname is None: fname = 'filtration-image.png'

    boi = True
    if ini is None or len(ini) == 0: boi = False

    bof = True
    if fin is None or len(fin) == 0: bof = False

    if (not boi) and (not bof): return False

    # check data size
    if len(ini[0]) <= 3: boi = False
    if len(fin[0]) <= 3: bof = False
    if (not boi) and (not bof): return False

    if boi:
        yini = ini[0]
    if bof:
        yfin = fin[0]

    if dt is None:
        if boi:
            xini = list(range(len(yini)))
        if bof:
            xfin = list(range(len(yfin)))
    else:
        if boi:
            xini = [ini[1]+dt*t for t in range(len(yini))]
        if bof:
            xfin = [fin[1]+dt*t for t in range(len(yfin))]

    if key is None:
        title = 'Filtration'
    else:
        if dt is None:
            if key.lower() == 'bonds':
                title = 'Filtration on Bond (Angstrom)'
            elif key.lower() == 'angles':
                title = 'Filtration on Angle (Degree)'
            else:
                title = 'Filtration'
        else:
            if key.lower() == 'bonds':
                title = 'Filtration on Bond ({:} Angstrom)'.format(dt)
            elif key.lower() == 'angles':
                title = 'Filtration on Angle ({:} Degree)'.format(dt)
            else:
                title = 'Filtration'

    if boi and bof:
        # both exist
        plt.plot(xini,yini,'r-',xfin,yfin,'b-')
    elif boi:
        plt.plot(xini,yini,'r-')
    elif bof:
        plt.plot(xfin,yfin,'r-')
    else:
        print('Fatal: this should be never executed')
    plt.title(title)

    # now save figure
    print('Note: filtration plot is saved to file < {:} >'.format(fname))
    plt.savefig(fname)
    plt.close()
    return True


def getrealsizeof(o):
    """recursively get the real size of built-in objects, unit in bytes
    """
    tot = sys.getsizeof(o)
    if isinstance(o,(int,str,float)):
        return tot
    if isinstance(o,(list,tuple)):
        return tot + sum([getrealsizeof(i) for i in o])
    if isinstance(o,dict):
        return tot + sum(getrealsizeof(k)+getrealsizeof(v) for k,v in o.items())
    return tot


