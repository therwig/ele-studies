import uproot
import numpy as np
import numba
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


@numba.jit
def truth_link(gen_parts_per_event, builder):
    # first fill a dictionary containing every electron like
    # daughter_dict[mother idx] = daughter idx
    # where each are electrons
    # if an ele daughter doesn't exist, give -1 index
    # assumes multiple electrons can't have the same mother
    
    for gen_parts in gen_parts_per_event:
        builder.begin_list()
        
        ele_mothers={} # store mother indices
        for i in range(len(gen_parts)):            
            if np.abs(gen_parts[i].pdgId)==11:
                ele_idx = i
                mother_idx = gen_parts[ele_idx].mother
                ele_mothers[ele_idx] = mother_idx if np.abs(gen_parts[mother_idx].pdgId)==11 else -1
                while(mother_idx>=0 and np.abs(gen_parts[mother_idx].pdgId)==11):
                    # step through mother tree, filling dict along the way
                    ele_idx = mother_idx
                    mother_idx = gen_parts[ele_idx].mother
                    ele_mothers[ele_idx]=mother_idx

        first_ancestor={} # store first ancestor with the same pdgId
        for e_idx in ele_mothers:
            idx = e_idx
            while idx in ele_mothers and ele_mothers[idx]>=0:
                idx = ele_mothers[idx]
            # (idx mother should now be -1)
            first_ancestor[e_idx] = idx
        
        for i in range(len(gen_parts)):            
            isFirst=False
            isLast=False
            motherPdgId=0

            # record mother IDs, including those from the first ele even in the case of subseq eles
            if gen_parts[i].mother>=0:
                if i in first_ancestor:
                    motherPdgId = gen_parts[gen_parts[first_ancestor[i]].mother].pdgId
                else:
                    motherPdgId = gen_parts[gen_parts[i].mother].pdgId
                
            if np.abs(gen_parts[i].pdgId)==11:
                if i in set(ele_mothers.keys()) and ele_mothers[i]<0:
                    #has no mother
                    isFirst=True
                if not (i in set(ele_mothers.values())):
                    #no ones mother
                    isFirst=True
            
            builder.begin_record()
            builder.field("isFirst")
            builder.append(isFirst)
            builder.field("isLast")
            builder.append(isLast)
            builder.field("motherPdgId")
            builder.append(motherPdgId)
            builder.end_record()

        builder.end_list()
    return builder



# histogram all attributes for input quantities
for collection in cms_events.columns:
    n_objs = ak.num(cms_events[collection]['pt'])
    plot(n_objs,"n_"+collection, xtitle=collection+" multiplicity")
    continue
    for attr in cms_events[collection].columns:
        flat_values = ak.flatten(ak.to_list(cms_events[collection][attr]))
        plot(flat_values,collection+'_'+attr, xtitle=collection+" "+attr)


        
# derived array with extra truth information
truth_builder = truth_link(cms_events['genParticles'], ak.ArrayBuilder())
truth_extension = truth_builder.snapshot()

ele_mask = np.abs(cms_events['genParticles']['pdgId']) == 11
z_mask = truth_extension['motherPdgId'] == 23
n2_mask = truth_extension['motherPdgId'] == 1000023
gen_ele_mask = ele_mask & (z_mask | n2_mask)
genEles = cms_events['genParticles'][gen_ele_mask]


# vals = np.abs(truth_extension['motherPdgId'][ele_mask])
# flat_values = ak.flatten(ak.to_list(vals))
# plot(flat_values,'ele_mothers', xtitle="ele mother id ")
# print(flat_values)

n_objs = ak.num(genEles)
print( np.sum(n_objs) )
plot(ak.num(genEles),"n_genElectrons", xtitle="gen electron multiplicity")
for attr in genEles.columns:
    flat_values = ak.flatten(ak.to_list(genEles[attr]))
    plot(flat_values,'genElectron_'+attr, xtitle="Gen Electron "+attr)

