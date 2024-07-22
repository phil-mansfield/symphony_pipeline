import read_gadget

file_fmt = "/fs/ddn/sdf/group/kipac/u/ycwang/MW4K_disk/Halo738_8K/output_disk/snapshot_%03d.0"
#"/fs/ddn/sdf/group/kipac/u/ycwang/MW4K_disk/Halo349/output_disk/snapshot_%03d.0"

for snap in range(235):
    print(read_gadget.Gadget2Zoom(file_fmt % snap).n_tot)
