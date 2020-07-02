from collections import namedtuple

pcfg = namedtuple('pcfg','lo, hi, nbins, log, name')

config = {
    'pt'          : pcfg(      0,     20, 100, False, 'pt [GeV]'),
    'eta'         : pcfg(-3.1415, 3.1415, 100, False, 'eta'),
    'phi'         : pcfg(-3.1415, 3.1415, 100, False, 'phi'),
    'mass'        : pcfg( -0.001,  0.001, 100, False, 'mass'),
    'dxy'         : pcfg(     -1,      1, 100, False, 'dxy [mm]'),
    'dxyErr'      : pcfg(      0,   0.25, 100, False, 'dxy Error [mm]'),
    'dz'          : pcfg(     -4,      4, 100, False, 'dz [mm]'),
    'dzErr'       : pcfg(      0,   0.25, 100, False, 'dz Error [mm]'),
    'fBrem'       : pcfg(   -0.5,      0, 100, False, 'fBrem'),
    'ptBiased'    : pcfg(     -5,     -2, 100, False, 'ptBiased'),
    'unBiased'    : pcfg(     -2,      2, 100, False, 'unBiased'),
    'ip3d'        : pcfg(      0,     35, 100, False, 'ip3d'),
    'sip3d'       : pcfg(      0,    500, 100, False, 'sip3d'),
    'mvaId'       : pcfg(    -10,     20, 100, False, 'mvaId'),
    'pfRelIso'    : pcfg(      0, 0.0001, 100, False, 'pfRelIso'),
    'pfvaId'      : pcfg(  19.99,     20, 100, False, 'pfvaId'),
    'trkRelIso'   : pcfg(      0,     10, 100, False, 'trkRelIso'),
    'vx'          : pcfg(     -2,      2, 100, False, 'vx'),
    'vy'          : pcfg(     -2,      2, 100, False, 'vy'),
    'vz'          : pcfg(    -10,     10, 100, False, 'vz'),
    'charge'      : pcfg(     -1,      1,  20, False, 'charge'),
    'pdgId'       : pcfg(    -11,     11,  20, False, 'pdgId'),
    'genPartIdx'  : pcfg(     -1,      0,  20,  True, 'genPartIdx'),
    'GenPartFlav' : pcfg(      0,      6,  20,  True, 'GenPartFlav'),
}

