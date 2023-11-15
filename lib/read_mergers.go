package lib

import (
	"fmt"
	"os"
	"encoding/binary"
	"math"
	"path"
)

type Mergers struct {
	Haloes, Snaps int
	Index []int32
	// left index goes over haloes, right goes over snapshots.
	Mvir, Rvir, Vmax, RVmax [][]float32
	ID [][]int32
	X, V [][][3]float32
}

func getOmegaM(suiteName string) float64 {
	omegaM, ok := map[string]float64{
		"MWest": 0.286,
		"SymphonyLMC": 0.286,
		"SymphonyMilkyWay": 0.286,
		"EDEN_MilkyWay_8K": 0.286,
		"SymphonyMilkyWayDisk": 0.286,
		"SymphonyMilkyWayDiskDMO": 0.286,
		"SymphonyMilkyWayHR": 0.286,
		"SymphonyGroup": 0.286,
		"SymphonyLCluster": 0.3,
		"SymphonyCluster": 0.25,
	}[suiteName]
	if !ok {
		panic(fmt.Sprintf("Unrecognized suite, %s", suiteName))
	}
	return omegaM
}

func getScaleFactors(mergerName string) []float64 {
	haloDir := path.Dir(mergerName)
	fname := path.Join(haloDir, "snap_scale.dat")

    f, err := os.Open(fname)
    if err != nil { panic(err.Error()) }
    defer f.Close()

	nSnap := int64(0)
    err = binary.Read(f, binary.LittleEndian, &nSnap)
    if err != nil { panic(err.Error()) }
	scales := make([]float64, nSnap)
    err = binary.Read(f, binary.LittleEndian, scales)
    if err != nil { panic(err.Error()) }

	return scales
}

func suiteHaloName(mergerName string) (suite, halo string) {
	haloDir := path.Dir(path.Dir(mergerName))
	suiteDir := path.Dir(haloDir)
	return path.Base(suiteDir), path.Base(haloDir)
}

func ScaleFactors(min, max float64, n int) []float64 {
	logMin, logMax := math.Log10(min), math.Log10(max)
	dlog := (logMax - logMin) / float64(n - 1)
	out := make([]float64, n)
	for i := range out {
		out[i] = math.Pow(10, dlog*float64(i) + logMin)
	}
	return out
}

func MvirToRvir(mvir, a, omegaM float64) float64 {
	omegaL := 1 - omegaM
	Ez := math.Sqrt(omegaM/(a*a*a) + omegaL)
	rhoCrit := 2.77519737e11*(Ez*Ez)
	omegaMz := omegaM/(a*a*a)/(Ez*Ez)

	x := omegaMz - 1
	deltaVir := 18*math.Pi*math.Pi + 82*x - 39*x*x
	rhoVir := rhoCrit*deltaVir

	rPhys := math.Pow(mvir/(rhoVir*math.Pi*4/3), 1.0/3)
	rComv := rPhys/a

	return rComv
}

func binaryRead(f *os.File, x interface{}) {
	err := binary.Read(f, binary.LittleEndian, x)
	if err != nil { panic(err.Error()) }
}

func binaryReadVector(f *os.File, n int) [][3]float32 {
	tmp := make([]float32, 3*n)
	binaryRead(f, tmp)
	out := make([][3]float32, n)
	for i := range out {
		for k := 0; k < 3; k++ {
			out[i][k] = tmp[3*i + k]
		}
	}
	return out
}

func ReadMergers(fname string,) *Mergers {
	f, err := os.Open(fname)
	if err != nil { panic(err.Error()) }
	defer f.Close()

	m := &Mergers{ }
	nh32, ns32 := int32(0), int32(0)
	binaryRead(f, &ns32)
	binaryRead(f, &nh32)
	nh, ns := int(nh32), int(ns32)
	m.Haloes, m.Snaps = nh, ns

	a := getScaleFactors(fname)
	suite, _ := suiteHaloName(fname)
	omegaM := getOmegaM(suite)

	m.Index = make([]int32, nh)
	binaryRead(f, m.Index)

	m.Mvir = make([][]float32, nh)
	m.Rvir = make([][]float32, nh)
	for i := 0; i < nh; i++ {
		m.Mvir[i] = make([]float32, ns)
		m.Rvir[i] = make([]float32, ns)
		binaryRead(f, m.Mvir[i])
		for j := range m.Mvir[i] {
			if m.Mvir[i][j] == -1 {
				m.Rvir[i][j] = -1
			} else {
				m.Rvir[i][j] = float32(
					MvirToRvir(float64(m.Mvir[i][j]), a[j], omegaM))
			}
		}
	}

	m.Vmax = make([][]float32, nh)
	for i := 0; i < nh; i++ {
		m.Vmax[i] = make([]float32, ns)
		binaryRead(f, m.Vmax[i])
	}

	m.RVmax = make([][]float32, nh)
	for i := 0; i < nh; i++ {
		m.RVmax[i] = make([]float32, ns)
		binaryRead(f, m.RVmax[i])
	}

	m.ID = make([][]int32, nh)
	for i := 0; i < nh; i++ {
		m.ID[i] = make([]int32, ns)
		binaryRead(f, m.ID[i])
	}

	m.X = make([][][3]float32, nh)
	for i := 0; i < nh; i++ {
		m.X[i] = binaryReadVector(f, ns)
	}

	m.V = make([][][3]float32, nh)
	for i := 0; i < nh; i++ {
		m.V[i] = binaryReadVector(f, ns)
	}

	return m
}
