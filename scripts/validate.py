import symlib
import numpy as np

base_dir = "/sdf/home/p/phil1/ZoomIns"

symlib.validate_symfind(base_dir, "MWest", -1, "../plots/validate", suffix="fid4", examples_per_bin=1, seed=1337)

