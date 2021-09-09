# A script to auto-generate all the config json files, so that I can
# manage them as a group.
import json

# A Python structure containing the base values that should appear in the output json files
def make_pobj():
    pobj = {
        "num_guiders" : 4,
        "totally_random_init" : True,
        "intro_every" : 0,

        "desc_arrays" : "Each array has num_guiders entries, in pairs. the order is [ x1, y1, x2, y2, etc ] where x1 means the 1st receptor/ligand that varies in the x direction",
        "desc_forms" : "0: linear  1: quadratic  2: exponential  3: logarithmic  4: alt exponential",
        "ret_receptor_forms"    : [2,2,2,2],
        "ret_ligand_forms"      : [2,2,2,2],
        "tectum_ligand_forms"   : [4,4,4,4],
        "desc_unused_arrays" : "The forms for tectal receptors are just used when setting up guidingtissue objects",
        "tectum_receptor_forms" : [2,2,2,2],

        "desc_directions" : "+x: increases with increasing x, etc",
        "ret_receptor_directions"    : ["-x", "-y", "+x", "+y"],
        "ret_ligand_directions"      : ["+x", "+y", "-x", "-y"],
        "tectum_ligand_directions"   : ["+y", "+x", "-y", "-x"],

        "desc_interaction" : "retina forward interaction gives the behaviour of retinal receptor expressing projections when the receptors are activated. May be attraction (1), no effect (0), or repulsion (-1)",
        "ret_forward_interactions" : [-1,-1,-1,-1],
        "desc_rcptrcpt_interactions" : "Do receptor receptor interactions generate repulsion or attraction? May be attraction (1), no effect (0), or repulsion (-1)",
        "ret_rcptrcpt_interactions" : [-1,0,-1,0],

        "desc_non_functional_directions" : "These two direction sets matter for visualisation but aren't currently used by the model",
        "tectum_receptor_directions" : ["x", "+y", "-x", "-y"],

        "rgcside" : 20,

        "desc_bpa" : "Branches per axon",                            "bpa" : 4,

        "desc_r" : "Radius of branches/cones (for visualisation)",     "r" : 0.005,
        "desc_rc" : "Interaction radius via competition interaction", "rc" : 0.005,
        "desc_rrr" : "Interaction radius via rcpt-rcpt interaction", "rrr" : 0.015,
        "desc_rrl" : "Interaction radius via rcpt-lgnd interaction", "rrl" : 0.015,

        "desc_s" : "Rcpt-rcpt signalling threshold value",             "s" : 1.1,

        "desc_m1" : "chemoaffinity axon-tectum rcpt-lgnd",            "m1" : 0.002,
        "desc_m2" : "the axon-axon rcpt-lgnd interactions",           "m2" : 0.0001,
        "desc_m3" : "the axon-axon rcpt-rcpt interactions",           "m3" : 0.002,
        "desc_mborder" : "border effect",                        "mborder" : 0.5,

        "desc_c1" : "axon-axon competition (leave out of paper)",     "c1" : 0
    }
    return pobj

# Save out one set of models, for a given configuration of rcpt/lgnd expression functions
def save_model_set (modeltag, pobj):
    fpath = 'm_{0}_GIJ.json'.format(modeltag)
    pobj["m1"] = m1
    pobj["m2"] = m2
    pobj["m3"] = m3
    pobj["c1"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_IJ.json'.format(modeltag)
    pobj["m1"] = 0.0
    pobj["m2"] = m2
    pobj["m3"] = m3
    pobj["c1"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GJ.json'.format(modeltag)
    pobj["m1"] = m1
    pobj["m2"] = 0.0
    pobj["m3"] = m3
    pobj["c1"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GI.json'.format(modeltag)
    pobj["m1"] = m1
    pobj["m2"] = m2
    pobj["m3"] = 0.0
    pobj["c1"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_G.json'.format(modeltag)
    pobj["m1"] = m1
    pobj["m2"] = 0.0
    pobj["m3"] = 0.0
    pobj["c1"] = 0.0
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GC.json'.format(modeltag)
    pobj["m1"] = m1
    pobj["m2"] = 0.0
    pobj["m3"] = 0.0
    pobj["c1"] = 0.009;
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

    fpath = 'm_{0}_GCI.json'.format(modeltag)
    pobj["m1"] = m1
    pobj["m2"] = 0.0
    pobj["m3"] = m3
    pobj["c1"] = 0.009;
    with open(fpath, 'w') as outfile:
        json.dump(pobj, outfile, indent=4)
    print ('Wrote {0}'.format (fpath))

# Set params once here
m1 = 0.002  # G
m2 = 0.0001 # J
m3 = 0.002  # I
c1 = 0.009  # C

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
