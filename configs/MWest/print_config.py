import numpy as np

# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.00017"
# change to the mass of a single particle in units of Msun/h
mp = "2.81981e5"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 8

# This is a bit wacky because I'm pulling it from the spreadsheet.
halo_strs = np.array([
    "004",
    "113",
    "719",
    "282",
    "788",
    "222",
    "666",
    "747",
    "659",
    "169",
    "975",
    "170",
    "476",
    "407",
    "349",
    "453",
    "229",
    "756",
    "983",
], dtype=str)

halo_ids = np.array([
    7216824,
    12577999,
    16625914,
    18541006,
    19800343,
    22209549,
    24120802,
    25734821,
    27530276,
    35512224,
    37061369,
    41313584,
    45959366,
    47755087,
    59736168,
    119579750,
    121339148,
    153433840,
    164134686,
], dtype=int)

halo_snaps = np.array([
    243,
    243,
    243,
    243,
    235,
    243,
    243,
    243,
    243,
    243,
    235,
    243,
    243,
    243,
    243,
    243,
    240,
    243,
    235,
], dtype=int)


order = np.argsort(halo_strs)
halo_strs, halo_ids, halo_snaps = halo_strs[order], halo_ids[order], halo_snaps[order]
haloes = ["Halo%s" % hs for hs in halo_strs]

# Use one string printf verb for the halo name, then two int printf verbs, one
# for the snapshot and the other for the index. Use double %% instead of single
# % for the ints and a single % for the halo name.
snapshot_format = "/oak/stanford/orgs/kipac/users/enadler/MWest_zoomins/%s/output/snapshot_%%03d.%%d"
# Change this to the directory which contains the consistent trees text files.
# Use a string printf verb for the halo name with a single %.
tree_dir = "/oak/stanford/orgs/kipac/users/enadler/MWest_zoomins/%s/output/rockstar/trees/"
# Directory where the output data products will go. Use a string printf verb
# for the halo name with a single %.
data_product_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/MWest/%s/"

# A string describing the style of the tree files. Currently supported are 
# "ct_rvmax", the forked version of consistent trees used by most of the
# Symphony suite, and "ct_rhapsody", the version used by rhapsody.
tree_style = "ct_rvmax"

fmt_string = ("%%d %%d %s %s %d %s %s %s %s nil" %
              (eps, mp, num_snapshot_files, 
               snapshot_format, tree_dir,
               data_product_dir, tree_style))

for i in range(len(haloes)):
    h = haloes[i]
    #print(fmt_string)
    print(fmt_string % (halo_ids[i], halo_snaps[i], h, h, h))
