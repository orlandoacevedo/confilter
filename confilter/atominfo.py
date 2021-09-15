class AtomInfo:
    """
    Method:
        get_atom    : return Atom(s,n,r,m,name,xyz)

    Args:
        s (str)     : sign, atomtype either can be real atom type or its number
        n (int)     : number in periodic table
        r (float)   : radius    (Angstrom)
        m (float)   : mass  (g/mol)
        name (str)  : full name
        xyz (List)  : cartesian coordinates (Angstrom), default [0., 0., 0.]

    Reference:
        Visualizing atomic sizes and molecular shapes with the classical
        turning surface of the Kohn–Sham potential

        Egor Ospadov, Jianmin Tao, Viktor N. Staroverov, John P. Perdew

        Proceedings of the National Academy of Sciences Dec 2018,
        115 (50) E11578-E11585; DOI: 10.1073/pnas.1814300115


        Atomic radius got from:

        Zhang, Q., et al.
        A rule-based algorithm for automatic bond type perception.
        J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26

    Note:
        A special dummy atom X is defined:
            sign=X, number=0, radius=0.0, mass=1.0, name=Dummy

        If conflict happens, Atom is defined in the sequence:
            1) sign
            2) number
            3) name

        If no inputs, dummy atom will be returned

        If any errors happen, dummy atom will be returned
    """
    ATOM_PROPERTY = """sign,number,radius,mass,name
        X,0,0,1,Dummy
        H,1,0.23,1.00794,Hydrogen
        He,2,0.93,4.002602,Helium
        Li,3,0.68,6.941,Lithium
        Be,4,0.35,9.01218,Beryllium
        B,5,0.83,10.811,Boron
        C,6,0.68,12.011,Carbon
        N,7,0.68,14.00674,Nitrogen
        O,8,0.68,15.9994,Oxygen
        F,9,0.64,18.998403,Fluorine
        Ne,10,1.12,20.1797,Neon
        Na,11,0.97,22.989768,Sodium
        Mg,12,1.1,24.305,Magnesium
        Al,13,1.35,26.981539,Aluminum
        Si,14,1.2,28.0855,Silicon
        P,15,1.05,30.973762,Phosphorus
        S,16,1.02,32.066,Sulfur
        Cl,17,0.99,35.4527,Chlorine
        Ar,18,1.57,39.948,Argon
        K,19,1.33,39.0983,Potassium
        Ca,20,0.99,40.078,Calcium
        Sc,21,1.44,44.95591,Scandium
        Ti,22,1.47,47.88,Titanium
        V,23,1.33,50.9415,Vanadium
        Cr,24,1.35,51.9961,Chromium
        Mn,25,1.35,54.93805,Manganese
        Fe,26,1.34,55.847,Iron
        Co,27,1.33,58.9332,Cobalt
        Ni,28,1.5,58.6934,Nickel
        Cu,29,1.52,63.546,Copper
        Zn,30,1.45,65.39,Zinc
        Ga,31,1.22,69.723,Gallium
        Ge,32,1.17,72.61,Germanium
        As,33,1.21,74.92159,Arsenic
        Se,34,1.22,78.96,Selenium
        Br,35,1.21,79.904,Bromine
        Kr,36,1.91,83.8,Krypton
        Rb,37,1.47,85.4678,Rubidium
        Sr,38,1.12,87.62,Strontium
        Y,39,1.78,88.90585,Yttrium
        Zr,40,1.56,91.224,Zirconium
        Nb,41,1.48,92.90638,Niobium
        Mo,42,1.47,95.94,Molybdenum
        Tc,43,1.35,97.9072,Technetium
        Ru,44,1.4,101.07,Ruthenium
        Rh,45,1.45,102.9055,Rhodium
        Pd,46,1.5,106.42,Palladium
        Ag,47,1.59,107.8682,Silver
        Cd,48,1.69,112.411,Cadmium
        In,49,1.63,114.818,Indium
        Sn,50,1.46,118.71,Tin
        Sb,51,,121.76,Antimony
        Te,52,1.47,127.6,Tellurium
        I,53,1.4,126.90447,Iodine
        Xe,54,1.98,131.29,Xenon
        Cs,55,1.67,132.90543,Cesium
        Ba,56,1.34,137.327,Barium
        La,57,1.87,138.9055,Lanthanum
        Ce,58,1.83,140.115,Cerium
        Pr,59,1.82,140.90765,Praseodymium
        Nd,60,1.81,144.24,Neodymium
        Pm,61,1.8,144.9127,Promethium
        Sm,62,1.8,150.36,Samarium
        Eu,63,1.99,151.965,Europium
        Gd,64,1.79,157.25,Gadolinium
        Tb,65,1.76,158.92534,Terbium
        Dy,66,1.75,162.5,Dysprosium
        Ho,67,1.74,164.93032,Holmium
        Er,68,1.73,167.26,Erbium
        Tm,69,1.72,168.93421,Thulium
        Yb,70,,173.04,Ytterbium
        Lu,71,,174.967,Lutetium
        Hf,72,,178.49,Hafnium
        Ta,73,,180.9479,Tantalum
        W,74,,183.84,Tungsten
        Re,75,,186.207,Rhenium
        Os,76,,190.23,Osmium
        Ir,77,,192.22,Iridium
        Pt,78,,195.08,Platinum
        Au,79,,196.96654,Gold
        Hg,80,,200.59,Mercury
        Tl,81,,204.3833,Thallium
        Pb,82,,207.2,Lead
        Bi,83,,208.98037,Bismuth
        Po,84,,208.9824,Polonium
        At,85,,209.9871,Astatine
        Rn,86,,222.0176,Radon
        Fr,87,,223.0197,Francium
        Ra,88,,226.0254,Radium
        Ac,89,,227.0278,Actinium
        Th,90,,232.0381,Thorium
        Pa,91,,231.03588,Protactinium
        U,92,,238.0289,Uranium
        Np,93,,237.048,Neptunium
        Pu,94,,244.0642,Plutonium
        Am,95,,243.0614,Americium
        Cm,96,,247.0703,Curium
        Bk,97,,247.0703,Berkelium
        Cf,98,,251.0796,Californium
        Es,99,,252.083,Einsteinium
        Fm,100,,257.0951,Fermium
        Md,101,,258.1,Mendelevium
        No,102,,259.1009,Nobelium
        Lr,103,,262.11,Lawrencium
        Rf,104,,261,Rutherfordium
        Db,105,,262,Dubnium
        Sg,106,,266,Seaborgium
        Bh,107,,264,Bohrium
        Hs,108,,269,Hassium
        Mt,109,,268,Meitnerium
        Ds,110,,269,Darmstadtium
        Rg,111,,272,Roentgenium
        Cn,112,,277,Copernicium
        Uut,113,,,Ununtrium
        Fl,114,,289,Flerovium
        Uup,115,,,Ununpentium
        Lv,116,,,Livermorium
        Uus,117,,,Ununseptium
        Uuo,118,,,Ununoctium"""

    def __init__(self):
        self.atominfo = []
        for line in self.ATOM_PROPERTY.split('\n')[1:]:
            t = line.split(',')
            sign = t[0].strip()
            number = int(t[1])
            try:
                radius = float(t[2])
            except ValueError:
                radius = 0.0
            try:
                mass = float(t[3])
            except ValueError:
                mass = 1.0
            name = t[4].strip()
            self.atominfo.append([sign,number,radius,mass,name])


FAI = AtomInfo()


