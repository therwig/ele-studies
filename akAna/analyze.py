import uproot
import numpy as np
import awkward1 as ak

from utils import plot

cms_dict = uproot.open("/Users/herwig/Desktop/dominic/data/nanoAOD.root")["Events"].arrays()
cms_dict_ak1 = {name.decode(): ak.from_awkward0(array) for name, array in cms_dict.items()}

cms_events = ak.zip({
    "genParticles": ak.zip({
        "pt":     cms_dict_ak1["GenPart_pt"],
        "eta":    cms_dict_ak1["GenPart_eta"],
        "phi":    cms_dict_ak1["GenPart_phi"],
        "mass":   cms_dict_ak1["GenPart_mass"],
        "status": cms_dict_ak1["GenPart_status"],
        "mother": cms_dict_ak1["GenPart_genPartIdxMother"],
        "pdgId":  cms_dict_ak1["GenPart_pdgId"],
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
    })
}, depth_limit=1)


# histogram all attributes
for collection in cms_events.columns:
    for attr in cms_events[collection].columns:
        flat_values = ak.flatten(ak.to_list(cms_events[collection][attr]))
        plot(flat_values,collection+'_'+attr, xtitle=collection+" "+attr)

# test gen_masks
ele_mask = np.abs(cms_events['genParticles']['pdgId']) == 11
for attr in cms_events['genParticles'].columns:
    flat_values = ak.flatten(ak.to_list(cms_events['genParticles'][attr]))
    plot(flat_values,'genElectron_'+attr, xtitle="Gen Particle "+attr)

# electron_mask = np.abs(cms_events['genParticles']['pdgId']) == 11

# electron_pt_values = ak.to_list( cms_events['genParticles']['pt'][electron_mask])
# flat_electron_pt_values = ak.flatten(electron_pt_values)

# plot(flat_pt_values,'test_all')
# plot(flat_electron_pt_values,'test_ele')

