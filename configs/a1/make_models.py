# A script to auto-generate all the config json files, so that I can
# manage them as a group.

# Set params once here. Sensitive to ee vs eE, I think
m_g = 0.002  # G
m_j = 0.01   # J
m_i = 0.5    # I 0.04
m_c = 0.05   # C
r_int = 0.05 # common interaction radius

# A Python structure containing the base values that should appear in the output json files
def make_pobj():
    pobj = {
        "num_guiders" : 4,
        "totally_random_init" : True,
        "intro_every" : 0,

        #"desc_arrays" : "Each array has num_guiders entries, in pairs. the order is [ x1, y1, x2, y2, etc ] where x1 means the 1st receptor/ligand that varies in the x direction",
        #"desc_forms" : "0: linear  1: quadratic  2: exponential  3: logarithmic  4: alt exponential",
        "ret_receptor_forms"    : [2,2,2,2],
        "ret_ligand_forms"      : [2,2,2,2],
        "tectum_ligand_forms"   : [4,4,4,4],
        #"desc_unused_arrays" : "The forms for tectal receptors are just used when setting up guidingtissue objects",
        "tectum_receptor_forms" : [2,2,2,2],

        #"desc_directions" : "+x: increases with increasing x, etc",
        "ret_receptor_directions"    : ["-x", "-y", "+x", "+y"],
        "ret_ligand_directions"      : ["+x", "+y", "-x", "-y"],
        "tectum_ligand_directions"   : ["+y", "+x", "-y", "-x"],

        #"desc_interaction" : "retina forward interaction gives the behaviour of retinal receptor expressing projections when the receptors are activated. May be attraction (1), no effect (0), or repulsion (-1)",
        "ret_forward_interactions" : [-1,-1,-1,-1],
        #"desc_rcptrcpt_interactions" : "Do receptor receptor interactions generate repulsion or attraction? May be attraction (1), no effect (0), or repulsion (-1)",
        "ret_rcptrcpt_interactions" : [-1,-1,-1,-1],

        #"desc_non_functional_directions" : "These two direction sets matter for visualisation but aren't currently used by the model",
        "tectum_receptor_directions" : ["x", "+y", "-x", "-y"],

        "rgcside" : 20,

        #"desc_bpa" : "Branches per axon",
        "bpa" : 4,
        #"desc_r" : "Radius of branches/cones (for visualisation)",
        "r" : 0.015, # 0.005 accurate ish, 0.015 good for vis
        #"desc_rc" : "Interaction radius via competition interaction",
        "r_c" : r_int,
        #"desc_r_i" : "Interaction radius via rcpt-rcpt interaction",
        "r_i" : r_int, # 0.015 ok. 0.05 S&G
        #"desc_rrl" : "Interaction radius via rcpt-lgnd interaction",
        "r_j" : r_int,
        #"desc_s" : "Rcpt-rcpt signalling threshold value",
        "s" : 1.1, # sensitive.
        #"desc_m_g" : "chemoaffinity axon-tectum rcpt-lgnd (G)",
        "m_g" : -1,
        #"desc_m_c" : "axon-axon competition (C)",
        "m_c" : -1,
        #"desc_m_i" : "the axon-axon rcpt-rcpt interactions (I)",
        "m_i" : -1,
        #"desc_m_j" : "the axon-axon rcpt-lgnd interactions (J)",
        "m_j" : -1,
        #"desc_mborder" : "border effect",
        "mborder" : 0.5
    }
    return pobj

import json
# Save out one set of models, for a given configuration of rcpt/lgnd expression functions
def save_model_set (modeltag, pobj):
    fpath = 'm_{0}_GIJ.json'.format(modeltag)
    pobj["m_g"] = m_g
    pobj["m_j"] = m_j
    pobj["m_i"] = m_i
    pobj["m_c"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_IJ.json'.format(modeltag)
    pobj["m_g"] = 0.0
    pobj["m_j"] = m_j
    pobj["m_i"] = m_i
    pobj["m_c"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GJ.json'.format(modeltag)
    pobj["m_g"] = m_g
    pobj["m_j"] = m_j
    pobj["m_i"] = 0.0
    pobj["m_c"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GI.json'.format(modeltag)
    pobj["m_g"] = m_g
    pobj["m_j"] = 0.0
    pobj["m_i"] = m_i
    pobj["m_c"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_G.json'.format(modeltag)
    pobj["m_g"] = m_g
    pobj["m_j"] = 0.0
    pobj["m_i"] = 0.0
    pobj["m_c"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GC.json'.format(modeltag)
    pobj["m_g"] = m_g
    pobj["m_j"] = 0.0
    pobj["m_i"] = 0.0
    pobj["m_c"] = m_c;
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GCI.json'.format(modeltag)
    pobj["m_g"] = m_g
    pobj["m_j"] = 0.0
    pobj["m_i"] = m_i
    pobj["m_c"] = m_c;
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

# These are all the model combos I might have:
modeltag = 'eE'
pobj = make_pobj()
pobj["ret_receptor_forms"]    = [2,2,2,2]
pobj["tectum_ligand_forms"]   = [4,4,4,4]
pobj["ret_ligand_forms"]      = [2,2,2,2]
pobj["tectum_receptor_forms"] = [2,2,2,2]
save_model_set (modeltag, pobj)

modeltag = 'ee'
pobj = make_pobj()
pobj["ret_receptor_forms"]    = [2,2,2,2]
pobj["tectum_ligand_forms"]   = [2,2,2,2]
pobj["ret_ligand_forms"]      = [2,2,2,2]
pobj["tectum_receptor_forms"] = [2,2,2,2]
save_model_set (modeltag, pobj)

modeltag = 'lin'
pobj = make_pobj()
pobj["ret_receptor_forms"]    = [0,0,0,0]
pobj["tectum_ligand_forms"]   = [0,0,0,0]
pobj["ret_ligand_forms"]      = [0,0,0,0]
pobj["tectum_receptor_forms"] = [0,0,0,0]
save_model_set (modeltag, pobj)

modeltag = 'eel' # One is EphB and attractive
pobj = make_pobj()
pobj["ret_receptor_forms"]    = [2,2,2,2]
pobj["tectum_ligand_forms"]   = [2,3,2,2]
pobj["ret_ligand_forms"]      = [2,2,2,2]
pobj["tectum_receptor_forms"] = [2,2,2,2]
pobj["tectum_ligand_directions"] = ["+y", "-x", "-y", "-x"]
pobj["ret_forward_interactions"] = [-1, 1,-1,-1]
save_model_set (modeltag, pobj)
