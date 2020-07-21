import uproot
import numpy as np
import numba
import awkward1 as ak
import os

def getData(opts=None):
    
    if opts and opts.input:
        cms_dict = uproot.open(opts.input)["Events"].arrays()
    else:
        user = os.environ['USER']
        if 'dlehner' in user:
           #cms_dict = uproot.open("/uscms/home/dlehner/nobackup/analysis/data/nanoAOD.root")["Events"].arrays()      # <-- smaller file from earlier
            cms_dict = uproot.open("/uscms/home/dlehner/nobackup/analysis/data/bigNanoAOD.root")["Events"].arrays()
        elif 'herwig' in user:
            if opts and opts.bigInput:
                cms_dict = uproot.open("/Users/herwig/Desktop/dominic/data/nanoAOD.root")["Events"].arrays()
            else:
                cms_dict = uproot.open("/Users/herwig/Desktop/dominic/data/orig_nanoAOD.root")["Events"].arrays()
        else:
            raise Exception("Must set data directory for user: {}!".format(user))
    cms_dict_ak1 = {name.decode(): ak.from_awkward0(array) for name, array in cms_dict.items()}
    
    cms_events = ak.zip({
        "genParticles": ak.zip({
            "pt"     : cms_dict_ak1["GenPart_pt"],
            "eta"    : cms_dict_ak1["GenPart_eta"],
            "phi"    : cms_dict_ak1["GenPart_phi"],
            "mass"   : cms_dict_ak1["GenPart_mass"],
            "status" : cms_dict_ak1["GenPart_status"],
            "mother" : cms_dict_ak1["GenPart_genPartIdxMother"],
            "pdgId"  : cms_dict_ak1["GenPart_pdgId"],
        }),
        "electrons": ak.zip({
            "dxy"        : cms_dict_ak1["Electron_dxy"],
            "dxyErr"     : cms_dict_ak1["Electron_dxyErr"],
            "dz"         : cms_dict_ak1["Electron_dz"],
            "dzErr"      : cms_dict_ak1["Electron_dzErr"],
            "eta"        : cms_dict_ak1["Electron_eta"],
            "ip3d"       : cms_dict_ak1["Electron_ip3d"],
            "mass"       : cms_dict_ak1["Electron_mass"],
            "phi"        : cms_dict_ak1["Electron_phi"],
            "pt"         : cms_dict_ak1["Electron_pt"],
            "charge"     : cms_dict_ak1["Electron_charge"],
            "pdgId"      : cms_dict_ak1["Electron_pdgId"],
            "genPartIdx" : cms_dict_ak1["Electron_genPartIdx"],
            "genPartFlav": cms_dict_ak1["Electron_genPartFlav"],
        }),
        "softElectrons": ak.zip({
            "pt"          : cms_dict_ak1["ElectronBPark_pt"],
            "phi"         : cms_dict_ak1["ElectronBPark_phi"],
            "eta"         : cms_dict_ak1["ElectronBPark_eta"],
            "mass"        : cms_dict_ak1["ElectronBPark_mass"],
            "dxy"         : cms_dict_ak1["ElectronBPark_dxy"],
            "dxyErr"      : cms_dict_ak1["ElectronBPark_dxyErr"],
            "dz"          : cms_dict_ak1["ElectronBPark_dz"],
            "dzErr"       : cms_dict_ak1["ElectronBPark_dzErr"],
            "fBrem"       : cms_dict_ak1["ElectronBPark_fBrem"],
            "ip3d"        : cms_dict_ak1["ElectronBPark_ip3d"],
            "mvaId"       : cms_dict_ak1["ElectronBPark_mvaId"],
            "pfRelIso"    : cms_dict_ak1["ElectronBPark_pfRelIso"],
            "pfvaId"      : cms_dict_ak1["ElectronBPark_pfmvaId"],
            "ptBiased"    : cms_dict_ak1["ElectronBPark_ptBiased"],
            "sip3d"       : cms_dict_ak1["ElectronBPark_sip3d"],
            "trkRelIso"   : cms_dict_ak1["ElectronBPark_trkRelIso"],
            "unBiased"    : cms_dict_ak1["ElectronBPark_unBiased"],
            "vx"          : cms_dict_ak1["ElectronBPark_vx"],
            "vy"          : cms_dict_ak1["ElectronBPark_vy"],
            "vz"          : cms_dict_ak1["ElectronBPark_vz"],
            "charge"      : cms_dict_ak1["ElectronBPark_charge"],
            "pdgId"       : cms_dict_ak1["ElectronBPark_pdgId"],
    #        "convVeto"    : cms_dict_ak1["ElectronBPark_convVeto"],         #these are all boolean values (having errors for now)
    #        "isLowPt"     : cms_dict_ak1["ElectronBPark_isLowPt"],
    #        "isPF"        : cms_dict_ak1["ElectronBPark_isPF"],
    #        "isPFoverlap" : cms_dict_ak1["ElectronBPark_isPFoverlap"],
            "genPartIdx"  : cms_dict_ak1["ElectronBPark_genPartIdx"],
            "GenPartFlav" : cms_dict_ak1["ElectronBPark_genPartFlav"],
        }),
    }, depth_limit=1)
    
    return cms_events
    
