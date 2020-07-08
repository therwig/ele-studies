from collections import namedtuple

pcfg = namedtuple('pcfg','lo, hi, nbins, log, profile, name')

def notNone(x):
    return not (x is None)

config = {
    'pt'          : pcfg(      0,     30, 100, False, False, 'pt [GeV]'),
    'eta'         : pcfg(-3.1415, 3.1415, 100, False, False, 'eta'),
    'phi'         : pcfg(-3.1415, 3.1415, 100, False, False, 'phi'),
    'mass'        : pcfg( -0.001,  0.001, 100, False, False, 'mass'),
    'dxy'         : pcfg(  -0.25,   0.25, 100, False, False, 'dxy [mm]'),
    'dxyErr'      : pcfg(      0,   0.25, 100, False, False, 'dxy Error [mm]'),
    'dz'          : pcfg(   -0.5,    0.5, 100, False, False, 'dz [mm]'),
    'dzErr'       : pcfg(      0,      1, 100, False, False, 'dz Error [mm]'),
    'fBrem'       : pcfg(      0,      1, 100, False, False, 'fBrem'),
    'ptBiased'    : pcfg(     -5,     20, 100, False, False, 'ptBiased'),
    'unBiased'    : pcfg(     -5,     20, 100, False, False, 'unBiased'),
    'ip3d'        : pcfg(      0,     11, 100, False, False, 'ip3d'),
    'sip3d'       : pcfg(      0,     20, 100, False, False, 'sip3d'),
    'mvaId'       : pcfg(    -10,     20, 100, False, False, 'mvaId'),
    'pfRelIso'    : pcfg(      0,   0.15, 100, False, False, 'pfRelIso'),
    'pfvaId'      : pcfg(   19.9,     20, 100, False, False, 'pfvaId'),
    'trkRelIso'   : pcfg(      0,      2, 100, False, False, 'trkRelIso'),
    'vx'          : pcfg(   -0.2,    0.2, 100, False, False, 'vx'),
    'vy'          : pcfg(   -0.2,    0.2, 100, False, False, 'vy'),
    'vz'          : pcfg(    -10,     10, 100, False, False, 'vz'),
    'charge'      : pcfg(     -1,      1,  20, False, False, 'charge'),
    'pdgId'       : pcfg(    -11,     11,  20, False, False, 'pdgId'),
    'genPartIdx'  : pcfg(     -1,      0,  20,  True, False, 'genPartIdx'),
    'GenPartFlav' : pcfg(      0,      6,  20,  True, False, 'GenPartFlav'),
}

truth_pl = {
            "pt"     : True,
            "eta"    : True,
            "phi"    : True,
            "mass"   : True,
            "status" : True,
            "mother" : True,
            "pdgId"  : True,
}
reco_pl = {
            "pt"          : True,
            "phi"         : True,
            "eta"         : True,
            "dxy"         : True,
            "dxyErr"      : True,
            "dz"          : True,
            "dzErr"       : True,
            "fBrem"       : True,
            "ip3d"        : True,
            "mvaId"       : True,
            "pfRelIso"    : True,
            "pfvaId"      : True,
            "ptBiased"    : True,
            "sip3d"       : True,
            "trkRelIso"   : True,
            "unBiased"    : True,
            "vx"          : True,
            "vy"          : True,
            "vz"          : True,
            "charge"      : False,
            "pdgId"       : False,
            "genPartIdx"  : False,
            "GenPartFlav" : False,
            "sxy"         : (lambda x : x['dxy'] / x['dxyErr']),
            "sz"          : (lambda x : x['dz'] / x['dzErr']),
}
