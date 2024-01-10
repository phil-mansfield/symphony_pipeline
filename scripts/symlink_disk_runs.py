import os
import os.path as path

def main():
    hosts = [ 23,  88, 119, 188, 247, 268, 270, 288, 327, 349, 364, 374, 414,
             415, 416, 440, 460, 469, 490, 530, 558, 567, 570, 606, 628, 641,
             675, 718, 738, 749, 797, 800, 825, 829, 852, 878, 881, 925, 926,
             937, 939, 967, 9749, 9829, 990]

    for host in hosts:
        link_host(host)

def link_host(host):
    dmo_fmt, disk_fmt, out_fmt = sim_fmt_strings(host)

    out_dir = path.dirname(out_fmt)
    if not path.exists(out_dir):
        os.mkdir(out_dir)

    link_all_files(dmo_fmt, disk_fmt, out_fmt)

def sim_fmt_strings(host):
    dmo_base_fmt = "/sdf/group/kipac/u/ycwang/MWmass_4K/Halo%03d/output/snapshot_%%03d.%%d"
    disk_base_fmt = "/sdf/group/kipac/u/ycwang/MW4K_disk/Halo%03d/output_disk/snapshot_%%03d.%%d"
    out_base_fmt = "/sdf/group/kipac/g/cosmo/ki21/phil1/data/SymphonyMilkyWayDisk/Halo%03d/snapshot_%%03d.%%d" 
    return dmo_base_fmt % host, disk_base_fmt % host, out_base_fmt % host

def link_all_files(dmo_fmt, disk_fmt, out_fmt):
    for i in range(126):
        for j in range(4):
            link_file(dmo_fmt % (i, j), out_fmt % (i, j))
    for i in range(127, 237):
        for j in range(4):
            link_file(disk_fmt % (i, j), out_fmt % (i-1, j))

def link_file(src, dst):
    #print("ln -s %s %s" % (src, dst))
    os.system("ln -s %s %s" % (src, dst))

if __name__ == "__main__": main()
