import numpy as np
import numba
import awkward1 as ak

PI=3.1415

@numba.jit
def mymod(a,n):
    return a - np.floor(a/n) * n
@numba.jit
def dPhi(a,b):
    return mymod(a-b+PI,2*PI)-PI
    

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


# truth matching
@numba.jit
def dr_match(evts_ts, evts_rs, builder, doReco=False, debug=False, nDebug=5, fullDebug=False):
    # note: allows multiple truths to match 
    # to same reco if they're within 0.05
    for ei in range(len(evts_ts)):
        if debug and ei<nDebug:
            print("\n>>>Checking event",ei)
        builder.begin_list()

        rs = evts_rs[ei]
        ts = evts_ts[ei]

        truth_to_reco={}
        for ti in range(len(ts)):
            t=ts[ti]
            if debug and ei<nDebug:
                print("Checking truth index",ti, "with eta phi",t.eta,t.phi)
            best_dr=-1
            best_ind=-1
            for ri in range(len(rs)):
                r=rs[ri]
                dr = np.sqrt((t.eta - r.eta)**2 + dPhi(t.phi, r.phi)**2)
                if fullDebug and ei<nDebug:
                    print("  dr(t",ti,", r",ri,") is ",dr," -- from reco eta phi",r.eta,r.phi)
                if dr<0.05 and (best_dr<0 or dr<best_dr):
                    best_dr=dr
                    best_ind=ri
                    if debug and ei<nDebug:
                        print("    updated best_dr is ",best_dr,", best ri is", best_ind)
            truth_to_reco[ti]=best_ind
            if debug and ei<nDebug:
                print("  truth_to_reco[",ti,"]=",best_ind)

            if not doReco:
                builder.begin_record()
                builder.field("truth_to_reco_index")
                builder.append(best_ind)
                builder.end_record()
                if debug and ei<nDebug:
                    print("  NORECO: truth_to_reco_index builder.append(",best_ind,")")

        if debug and ei<nDebug:
            print("set",set(truth_to_reco.values()))
        if doReco:
            for ri in range(len(rs)):
                ind=-1
                if ri in set(truth_to_reco.values()):
                    ind = list(truth_to_reco.keys())[list(truth_to_reco.values()).index(ri)]
                    if debug and ei<nDebug:
                        print("found ri=",ri," in matched values, with truth index",ind)
                builder.begin_record()
                builder.field("reco_to_truth_index")
                builder.append(ind)
                builder.end_record()
                if debug and ei<nDebug:
                    print("  RECO: reco_to_truth_index (ri =",ri,") builder.append(",ind,")")
            
        builder.end_list()
    return builder

# # truth matching
# @numba.jit
# def dr_matchR(evts_ts, evts_rs, builder):
#     # note: allows multiple truths to match 
#     # to same reco if they're within 0.05
#     for ei in range(len(evts_ts)):
#         builder.begin_list()

#         rs = evts_rs[ei]
#         ts = evts_ts[ei]

#         truth_to_reco={}
#         for ti in range(len(ts)):
#             t=ts[ti]
#             best_dr=-1
#             best_ind=-1
#             for ri in range(len(rs)):
#                 r=rs[ri]
#                 dr = np.sqrt((t.eta - r.eta)**2 + (t.phi - r.phi)**2)
#                 if dr<0.05 and (best_dr<0 or dr<best_dr):
#                     best_dr=dr
#                     best_ind=ri
#             truth_to_reco[ti]=ri
            
#             # builder.begin_record()
#             # builder.field("truth_to_reco_index")
#             # builder.append(best_ind)
#             # builder.end_record()

#         for ri in range(len(rs)):
#             #r=rs[ri]
#             ind=-1
#             if ri in set(truth_to_reco.values()):
#                 ind = list(truth_to_reco.keys())[list(truth_to_reco.values()).index(ri)]
#             builder.begin_record()
#             builder.field("reco_to_truth_index")
#             builder.append(ind)
#             builder.end_record()
            
#         builder.end_list()
#     return builder
