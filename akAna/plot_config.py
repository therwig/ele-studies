from collections import namedtuple

pcfg = namedtuple('pcfg','lo, hi, nbins, log, name')

config = {
    'pt':   pcfg(      0,     20, 40,  True, 'p_T [GeV]'),
    'eta':  pcfg(   -2.8,    2.8, 40, False, 'eta'),
    'phi':  pcfg(-3.1415, 3.1415, 40, False, 'phi'),
    'dxy':  pcfg(    -20,     20, 40,  True, 'd_{xy} [mm]'),
    'dz':   pcfg(    -20,     20, 40, False, 'd_{z} [mm]'),
}

