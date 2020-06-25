from collections import namedtuple

pcfg = namedtuple('pcfg','lo, hi, nbins, log, name')

config = {
    'pt'          : pcfg(      0,     20, 40,  True, 'p_T [GeV]'),
    'eta'         : pcfg(   -2.8,    2.8, 40, False, 'eta'),
    'phi'         : pcfg(-3.1415, 3.1415, 40, False, 'phi'),
    'mass'        : pcfg(   -0.5,    0.5, 60,  True, 'mass'),
    'dxy'         : pcfg(    -20,     20, 40,  True, 'd_{xy} [mm]'),
    'dxyErr'      : pcfg(      0,     10, 60,  True, 'd_{xy} Error [mm]'),
    'dz'          : pcfg(    -20,     20, 40, False, 'd_{z} [mm]'),
    'dzErr'       : pcfg(      0,     10, 40,  True, 'd_{z} Error [mm]'),
    'fBrem'       : pcfg(  -0.05,      0, 60,  True, 'fBrem'),
    'ip3d'        : pcfg(      0,     25, 60,  True, 'ip3d'),
    'sip3d'       : pcfg(      0,   2200, 40,  True, 'sip3d'),
    'pfRelIso'    : pcfg(      0,      2, 40,  True, 'pfRelIso'),
    'pfvaId'      : pcfg(     18,     20, 40,  True, 'pfvaId'),
    'trkRelIso'   : pcfg(      0,    150, 60,  True, 'trkRelIso'),
    'vx'          : pcfg(    -10,     10, 60,  True, 'vx'),
    'vy'          : pcfg(    -10,     10, 60,  True, 'vy'),
    'vz'          : pcfg(    -30,     30, 60,  True, 'vz'),
    'genPartIdx'  : pcfg(     -1,      3, 80,  True, 'genPartIdx'),
    'GenPartFlav' : pcfg(      0,     10, 80,  True, 'GenPartFlav'),
}

