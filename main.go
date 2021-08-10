package main

import (
	"fmt"
	"os"

	common "git.oa.com/hengrushang/RampCfd/Common"
	"git.oa.com/hengrushang/RampCfd/Geometry"
	"git.oa.com/hengrushang/RampCfd/Solver"
)

// 主程序
func main() {
	// 几何体参数
	const Xl float64 = 0.0
	const Xr float64 = 3.0
	const Yb float64 = 0.0
	const Yt float64 = 1.0
	const rb float64 = 0.6
	const re float64 = 1.2
	const rh float64 = 0.2
	// 网格参数
	const Nx int = 400
	const Ny int = 100
	const Ng int = 2
	// 求解过程控制参数
	const Niter int = 500           // 最大迭代次数
	const gamma float64 = 1.4       // 气体常数
	const kappa float64 = -1        // MUSCL重构全上风格式
	const converCirt float64 = 1e-6 // 残差收敛准则
	const Ma float64 = 0.0          // 马赫数
	const Density float64 = 1.4     // 无量纲密度
	const Pressure float64 = 1.0    // 无量纲压力
	const cfl float64 = 0.1         // cfl数
	const order float64 = 2.0       // 重构器参数
	Ramp := Geometry.InitGeometry(Xl, Xr, Yb, Yt, rb, re, rh)
	SpaceMesh := common.InitSpaceMesh(Nx, Ny, Ng, &Ramp)
	ProtoSolver := Solver.SolverInitlize(SpaceMesh, gamma)
	Roesolver := &Solver.Roe{Solver: *ProtoSolver}
	Roesolver.InitalizeField(Ma, Density, Pressure, gamma)
	iter := 1
	loop_ctrl := true
	for loop_ctrl {
		err1 := Roesolver.UpdateTimeStep(cfl)
		if err1 != "" {
			fmt.Println(err1, iter)
			break
		}
		err1 = Roesolver.UpdateBC()
		if err1 != "" {
			fmt.Println(err1, iter)
			break
		}
		Roesolver.UpdateLimiter()
		err1 = Roesolver.UpdateFlux(kappa, order)
		if err1 != "" {
			fmt.Println(err1, iter)
			break
		}
		err1 = Roesolver.UpdateResidual()
		if err1 != "" {
			fmt.Println(err1, iter)
			break
		}
		err1 = Roesolver.UpdateSolution(gamma)
		if err1 != "" {
			fmt.Println(err1, iter)
			break
		}
		iter++
		fmt.Println("迭代次数+", iter)
		loop_ctrl = iter < Niter
	}
	// 结果文件存储
	pressure := []float64{}
	density := []float64{}
	xlocation := []float64{}
	ylocation := []float64{}
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			xlocation = append(xlocation, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Xlocation)
			ylocation = append(ylocation, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Yloaction)
			density = append(density, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Density)
			pressure = append(pressure, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Pressure)
		}
	}
	f, err := os.Create("results.txt")
	if err != nil {
		fmt.Println("os create error", err)
		return
	}
	defer f.Close()
	for i := range xlocation {
		fmt.Fprintf(f, "%.*e, %.*e, %.*e, %.*e\n", 4, xlocation[i], 4, ylocation[i], 5, density[i], 5, pressure[i])
	}

}
