from collections import namedtuple

pcfg = namedtuple('pcfg','lo, hi, nbins, log, name')

config = {
    'pt'          : pcfg(      0,     30, 40,  True, 'p_T [GeV]'),
    'eta'         : pcfg(-3.1415, 3.1415, 40, False, 'eta'),
    'phi'         : pcfg(-3.1415, 3.1415, 40, False, 'phi'),
    'mass'        : pcfg(  -0.05,   0.05, 40,  True, 'mass'),
    'dxy'         : pcfg(    -25,     25, 40,  True, 'd_{xy} [mm]'),
    'dxyErr'      : pcfg(      0,     20, 40,  True, 'd_{xy} Error [mm]'),
    'dz'          : pcfg(    -20,     20, 40, False, 'd_{z} [mm]'),
    'dzErr'       : pcfg(      0,     20, 40,  True, 'd_{z} Error [mm]'),
    'fBrem'       : pcfg(   -0.2,      0, 40, False, 'fBrem'),
    'ip3d'        : pcfg(      0,     35, 40,  True, 'ip3d'),
    'sip3d'       : pcfg(      0,   2750, 40,  True, 'sip3d'),
    'pfRelIso'    : pcfg(      0,    3.5, 40,  True, 'pfRelIso'),
    'pfvaId'      : pcfg(   19.9,     20, 40,  True, 'pfvaId'),
    'trkRelIso'   : pcfg(      0,    200, 40,  True, 'trkRelIso'),
    'vx'          : pcfg(    -20,     20, 40,  True, 'vx'),
    'vy'          : pcfg(    -20,     20, 40,  True, 'vy'),
    'vz'          : pcfg(    -45,     45, 40,  True, 'vz'),
    'charge'      : pcfg(     -1,      1, 40, False, 'charge'),
    'pdgId'       : pcfg(    -11,     11, 40, False, 'pdgId'),
    'genPartIdx'  : pcfg(     -1,   -0.5, 40, False, 'genPartIdx'),
    'GenPartFlav' : pcfg(      0,      6, 40, False, 'GenPartFlav'),
}

