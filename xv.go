package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"strconv"
	"github.com/phil-mansfield/symphony_pipeline/lib"
)

var (
	HRLevel = 1
)

func main() {
	haloIndex := -1
    if len(os.Args) > 3 || len(os.Args) <= 1 {
        panic("You must supply a config file and an (optional) index to " +
			"restrict analysis to.")
    } else if len(os.Args) == 3 {
		var err error
		haloIndex, err = strconv.Atoi(os.Args[2])
		if err != nil { panic(err.Error()) }
	}
	
	inputName := os.Args[1]
	
	cfg := lib.ParseConfig(inputName)
	
	for i := range cfg.BaseDir {
		if haloIndex != -1 && i != haloIndex { continue }
		XV(cfg, i)
		runtime.GC()
	}
}

func XV(cfg *lib.Config, cfgi int) {
	log.Printf("Assigning positions and velocities for halo %d (%s)",
		cfgi, cfg.BaseDir[cfgi])

	baseDir, snapFmt := cfg.BaseDir[cfgi], cfg.SnapFormat[cfgi]
	mergers := lib.ReadMergers(lib.MergerFileName(baseDir))

	maxSnap := len(mergers.Rvir[0]) - 1

	fileName := fmt.Sprintf(snapFmt, maxSnap, 0)
	header := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})

	for snap := 0; snap <= maxSnap; snap++ {
		fName := fmt.Sprintf(snapFmt, snap, 0)
		f := lib.OpenGadget2Zoom(fName, []string{"x", "v", "id32"})
		fmt.Printf("%3d -> %d\n", snap, f.NTot)
		fmt.Println(fName)
	}

	fmt.Println(header.NTot[HRLevel])

	x := make([][3]float32, header.NTot[HRLevel])
	v := make([][3]float32, header.NTot[HRLevel])

	pHeader := lib.ReadParticleHeader(cfg.BaseDir[cfgi])

	tags := lib.ReadTags(cfg.BaseDir[cfgi], pHeader)

	//fmt.Printf("%d\n", len(tags.ID[10]))
	//panic("")

	xh := NewHaloVector(tags)
	vh := NewHaloVector(tags)

	for snap := 0; snap <= maxSnap; snap++ {
		log.Printf("Snapshot %3d", snap)
		dirName := lib.SnapDirName(cfg.BaseDir[cfgi], snap)
		lib.MaybeMkdir(dirName)

		ClearVectors(x)
		ClearVectors(v)

		fileNames := make([]string, cfg.Blocks[cfgi])
		for b := range fileNames {
			fileNames[b] = fmt.Sprintf(snapFmt, snap, b)
		}
		//CheckIDs(fileNames)

		for b := 0; b < cfg.Blocks[cfgi]; b++ {
			fileName := fmt.Sprintf(snapFmt, snap, b)
			ReadToIDGrid(fileName, x, v)
			runtime.GC()
		}

		CheckVectors(x)
		CheckVectors(v)

		FillHaloVector(tags, x, xh, snap)
		FillHaloVector(tags, v, vh, snap)
	
		lib.WriteVector(pHeader, cfg.BaseDir[cfgi], "x", snap, xh)
		lib.WriteVector(pHeader, cfg.BaseDir[cfgi], "v", snap, vh)

		runtime.GC()
	}
}

func CountFiles(pHeader *lib.ParticleHeader) int {
	max := pHeader.FileIdxs[0]
	for i := range pHeader.FileIdxs {
		if pHeader.FileIdxs[i] > max {
			max = pHeader.FileIdxs[i]
		}
	}
	return int(max)
}

func CheckIDs(fileNames []string) {
	f := lib.OpenGadget2Zoom(fileNames[0], []string{"x", "v", "id32"})
	
	found := make([]int, f.NTot[HRLevel])
	extras := []int32{ }
	blocks := []int{ }
	idxs := []int{ }

	fmt.Printf("NTot = %d\n", f.NTot[HRLevel])

	for b := range fileNames {
		f = lib.OpenGadget2Zoom(fileNames[b], []string{"x", "v", "id32"})
		fmt.Printf("%d\n", f.NTot[HRLevel])
		idp := make([]int32, f.N[HRLevel])
		f.Read("id32", HRLevel, idp)


		min, max := idp[0], idp[0]
		for i := range idp {
			if idp[i] >= 0 && idp[i] < int32(len(found)) {
				found[idp[i]]++
			} else {
				extras = append(extras, idp[i])
				blocks = append(blocks, b)
				idxs = append(idxs, i)
			}

			if idp[i] < min { min = idp[i] }
			if idp[i] > max { max = idp[i] }
		}
		fmt.Printf("%d) min = %8d max = %8d\n", b, min, max)
	}


	j := 0
	for i := range found {
		if found[i] == 0 {
			j++
			fmt.Printf("Missing: %d\n", i)
		} else if found[i] > 1 {
			j++
			fmt.Printf("Multi: %d (%d)\n", i, found[i])
		}
	}
	if j != 0 {
		fmt.Printf("%d errors found\n", j)
	} else {
		fmt.Printf("All IDs found once\n")
	}

	fmt.Printf("Extra IDs: %d\n", extras)
	fmt.Printf("Blocks: %d\n", blocks)
	fmt.Printf("Indices: %d\n", idxs)
}

func ReadToIDGrid(fileName string, x, v [][3]float32) {
	f := lib.OpenGadget2Zoom(fileName, []string{"x", "v", "id32"})
	
	xp := make([][3]float32, f.N[HRLevel])
	idp := make([]int32, f.N[HRLevel])

	f.Read("x", HRLevel, xp)
	f.Read("id32", HRLevel, idp)

	//fmt.Println(fileName)
	for i := range xp {
		//if idp[i] >= int32(len(x)) {
		//	fmt.Printf("|x| = %d, idp[%d] = %d\n", len(x), i, idp[i])
		//}
		x[idp[i]] = xp[i]
	}

	vp := xp
	f.Read("v", HRLevel, vp)
	
	for i := range vp {
		v[idp[i]] = vp[i]
	}
}

func ClearVectors(x [][3]float32) {
	for i := range x {
		for dim := 0; dim < 3; dim++ {
			x[i][dim] = float32(math.NaN())
		}
	}
}

func CheckVectors(x [][3]float32) {
	nUnset := 0
	for i := range x {
		for dim := 0; dim < 3; dim++ {
			if math.IsNaN(float64(x[i][dim])) {
				nUnset++
			}
		}
	}
	
	for i := range x {
		for dim := 0; dim < 3; dim++ {
			if math.IsNaN(float64(x[i][dim])) {
				panic(fmt.Sprintf("Vector with ID %d not set. %d total " +
					"vectors unset.", i + 1, nUnset))
			}
		}
	}
}

func NewHaloVector(tags *lib.Tags) [][][3]float32 {
	vec := make([][][3]float32, len(tags.N0))
	for i := range vec {
		vec[i] = make([][3]float32, tags.N0[i])
	}
	return vec
}

func FillHaloVector(tags *lib.Tags, x [][3]float32,
	xh [][][3]float32, snap int) {

	for i := range xh {
		xh[i] = xh[i][:cap(xh[i])]
		for j := range xh[i] {
			xh[i][j] = x[tags.ID[i][j]]
		}

		xh[i] = tags.TrimVector(i, xh[i], snap)
	}
}
