import symlib
import numpy as np

halo_info_csv = """Halo004,"1.00000",7450336,7448193,4914244
Halo113,"1.01283",12917535,12923750,8736472
Halo169,"0.97483",36691240,36691242,27250046
Halo170,"0.97483",42817810,42822484,24135360
Halo222,"1.05231",22703308,22703326,15089066
Halo229,"1.01283",123401061,123401137,88829369
Halo282,"1.03899",18923242,18925117,13101064
Halo327,"0.93825",8711942,8716183,6394862
Halo349,"0.90305",63114621,63114618,48508675
Halo407,"1.00000",49473371,49473886,35063709
Halo453,"0.96248",120639573,120639561,76461101
Halo476,"0.96248",46437616,46395350,36762450
Halo659,"0.96248",28641820,28652192,17758107
Halo666,"1.00000",24852563,24847328,16913933
Halo719,"1.01283",17060652,17062196,11497567
Halo747,"0.97483",26674150,26661126,20015928
Halo756,"1.00000",157680793,157680804,120198019
Halo788,"1.00000",19800343,19800919,12116619
Halo975,"0.98733",37178726,37178727,23433155
Halo983,"1.00000",164134686,164134692,106437843"""

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "MWest"

lines = halo_info_csv.split("\n")
for line in lines:
    name, a, host, lmc, gse = line.split(",")
    a, lmc, gse = float(a.strip('"')), int(lmc), int(gse)
    sim_dir = symlib.get_host_directory(base_dir, suite, name)
    scale = symlib.scale_factors(sim_dir)
    h, hist = symlib.read_rockstar(sim_dir)

    snap = np.argmin(np.abs(a - scale))

    print(lmc in h["id"], gse in h["id"])
    
    i_lmc, i_gse = -1, -1
    for i in range(len(h)):
        if lmc in h["id"][i]: i_lmc = i
        if gse in h["id"][i]: i_gse = i

    print(name, snap, i_lmc, i_gse)
