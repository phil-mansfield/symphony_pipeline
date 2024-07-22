import symlib
import numpy as np

base_dir = "/sdf/home/p/phil1/ZoomIns"

symlib.validate_symfind(base_dir, "SymphonyGroup", np.arange(17, dtype=int), "../plots/validate", suffix="fid3")
