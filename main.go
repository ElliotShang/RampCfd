package main

import (
	"fmt"
	"math"
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
	const Ny int = 120
	const Ng int = 1
	// 求解过程控制参数
	const Niter int = 60000           // 最大迭代次数
	const gamma float64 = 1.4         // 气体常数
	const kappa float64 = -1          // MUSCL重构全上风格式
	const converCirt float64 = 1e-6   // 残差收敛准则
	const Ma float64 = 2.0            // 马赫数
	const Density float64 = 0.53      // 无量纲密度
	const Pressure float64 = 31836.15 // 无量纲压力
	const cfl float64 = 0.4           // cfl数
	const order float64 = 1.0         // 重构器参数
	Ramp := Geometry.InitGeometry(Xl, Xr, Yb, Yt, rb, re, rh)
	SpaceMesh := common.InitSpaceMesh(Nx, Ny, Ng, &Ramp)
	ProtoSolver := Solver.SolverInitlize(SpaceMesh, gamma)
	Roesolver := &Solver.Roe{Solver: *ProtoSolver}
	Roesolver.InitalizeField(Ma, Density, Pressure, gamma)
	iter := 1
	loop_ctrl := true
	// 残差文件
	f1, err := os.Create("residuals.txt")
	if err != nil {
		fmt.Println("os create error", err)
		return
	}
	defer f1.Close()
	for loop_ctrl {
		err1 := Roesolver.UpdateBC()
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
		err1 = Roesolver.UpdateTimeStep(cfl)
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
		// 残差输出
		if iter%10 == 0 {
			ResidualFlux := common.ConvectiveFlux{}
			for i := 0; i < Nx; i++ {
				for j := 0; j < Ny; j++ {
					ResidualFlux.ConvFlux1 += math.Abs(Roesolver.Residual[i][j].ConvFlux1)
					ResidualFlux.ConvFlux2 += math.Abs(Roesolver.Residual[i][j].ConvFlux2)
					ResidualFlux.ConvFlux3 += math.Abs(Roesolver.Residual[i][j].ConvFlux3)
					ResidualFlux.ConvFlux4 += math.Abs(Roesolver.Residual[i][j].ConvFlux4)
				}
			}
			Norm1 := ResidualFlux.ConvFlux1 / float64(Nx) / float64(Ny)
			Norm2 := ResidualFlux.ConvFlux2 / float64(Nx) / float64(Ny)
			Norm3 := ResidualFlux.ConvFlux3 / float64(Nx) / float64(Ny)
			Norm4 := ResidualFlux.ConvFlux4 / float64(Nx) / float64(Ny)
			fmt.Fprintf(f1, "%d, %f, %f, %f, %f\n", iter, Norm1, Norm2, Norm3, Norm4)
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
	u := []float64{}
	v := []float64{}
	mach := []float64{}
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			xlocation = append(xlocation, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Xlocation)
			ylocation = append(ylocation, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Yloaction)
			density = append(density, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Density)
			pressure = append(pressure, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Pressure)
			c := math.Sqrt(gamma * Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Pressure / Roesolver.GetMesh().Mesh[i+Ng][j+Ng].Density)
			ui := Roesolver.GetMesh().Mesh[i+Ng][j+Ng].VelocityX
			vi := Roesolver.GetMesh().Mesh[i+Ng][j+Ng].VelocityY
			u = append(u, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].VelocityX)
			v = append(v, Roesolver.GetMesh().Mesh[i+Ng][j+Ng].VelocityY)
			mach = append(mach, math.Sqrt(ui*ui+vi*vi)/c)
		}
	}
	f, err := os.Create("results.txt")
	if err != nil {
		fmt.Println("os create error", err)
		return
	}
	defer f.Close()
	for i := range xlocation {
		fmt.Fprintf(f, "%.*e, %.*e, %f, %f, %f, %f, %f\n", 4, xlocation[i], 4, ylocation[i], density[i], u[i], v[i], pressure[i], mach[i])
	}

}
