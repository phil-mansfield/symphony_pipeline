import numpy as np

# Change to the z=0 Plummer-equivalent force softneing scale of the simulation
# in units of Mpc/h
eps = "0.00017"
# change to the mass of a single particle in units of Msun/h
mp = "2.81981e5"
# Number of particle files per snapshot in the Gadget output
num_snapshot_files = 8

#haloes = ['Halo004', 'Halo113', 'Halo169', 'Halo170', 'Halo222', 'Halo229', 'Halo282', 'Halo327', 'Halo349', 'Halo407', 'Halo453', 'Halo476', 'Halo523', 'Halo625', 'Halo659', 'Halo666', 'Halo719', 'Halo747', 'Halo756', 'Halo788', 'Halo858', 'Halo908', 'Halo953', 'Halo975', 'Halo983']
#halo_ids = [7208101, 12594055, 35694747, 41628142, 22040617, 121130872, 18312877, 8419639, 61419320, 47937421, 120639573, 46437616, 33099515, 110809080, 27827883, 24087537, 16553233, 25940933, 155206896, 19800343, 58886322, 8662410, 49383341, 37178726, 164134686]

# This is a bit wacky because I'm pulling it from the spreadsheet.
halo_strs = np.array([
    "004",
    "327",
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
    8232652,
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
], dtype=int)+1

halo_snaps = np.array([
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
