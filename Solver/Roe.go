package Solver

import (
	"fmt"
	"math"

	common "git.oa.com/hengrushang/RampCfd/Common"
	muscl "git.oa.com/hengrushang/RampCfd/MUSCL"
)

// roe平均变量
type RoeAvarageFlux struct {
	Densitybar   float64
	VelocityXbar float64
	VelocityYbar float64
	Hbar         float64
	Cbar         float64
	Unbar        float64
}

// Roe求解器对象 继承通用求解器接口
type Roe struct {
	Solver
}

// 实现通量更新函数
func (solver *Roe) UpdateFlux(kappa float64, order float64) string {
	Nx := solver.mesh.Nx
	Ny := solver.mesh.Ny
	Ng := solver.mesh.NG
	for i := 0; i < Nx+1; i++ {
		for j := 0; j < Ny; j++ {
			// muscl 插值重构Vl,Vr
			Vl, Vr := muscl.Muscl(kappa, order, i, j, solver.mesh, solver.PositiveXlimiters[i][j],
				solver.NegativeXlimiters[i][j], solver.PositiveXlimiters[i+1][j], solver.NegativeXlimiters[i+1][j], 'x', solver.gamma)
			// 均量
			//Vl := common.PrimtiveFlux{Density: solver.mesh.Mesh[i+Ng-1][j+Ng].Density, VelocityX: solver.mesh.Mesh[i+Ng-1][j+Ng].VelocityX,
			//	VelocityY: solver.mesh.Mesh[i+Ng-1][j+Ng].VelocityY, Pressure: solver.mesh.Mesh[i+Ng-1][j+Ng].Pressure}
			//Vr := common.PrimtiveFlux{Density: solver.mesh.Mesh[i+Ng][j+Ng].Density, VelocityX: solver.mesh.Mesh[i+Ng][j+Ng].VelocityX,
			//	VelocityY: solver.mesh.Mesh[i+Ng][j+Ng].VelocityY, Pressure: solver.mesh.Mesh[i+Ng][j+Ng].Pressure}
			CenterFlux := CenterCFlux(Vl, Vr, solver.mesh.VectorX[i+Ng][j+Ng], solver.gamma)
			MatrixFlux := MatrixConFlux(Vl, Vr, solver.mesh.VectorX[i+Ng][j+Ng], solver.gamma)
			/**
			if i == 159 && j == 0 {
				fmt.Println("通量", CenterFlux, MatrixFlux)
				fmt.Println("Roe特征值", RoeValue(Vl, Vr, 1.4))
				fmt.Println("Roe dw", Amplitude(Vl, Vr, RoeValue(Vl, Vr, 1.4), solver.mesh.VectorX[i+Ng][j+Ng]))
				fmt.Println("Eigen value", EignValue(RoeValue(Vl, Vr, 1.4), solver.mesh.VectorX[i+Ng][j+Ng]))
				fmt.Println("Eign vector", RightEignVec(RoeValue(Vl, Vr, 1.4), solver.mesh.VectorX[i+Ng][j+Ng]))
				fmt.Println("vectorX", solver.mesh.VectorX[i+Ng][j+Ng])
			}
			**/
			solver.FluxX[i][j] = CenterFlux.FluxMinus(MatrixFlux)
			if math.IsNaN(solver.FluxX[i][j].ConvFlux1) || math.IsNaN(solver.FluxX[i][j].ConvFlux2) || math.IsNaN(solver.FluxX[i][j].ConvFlux3) || math.IsNaN(solver.FluxX[i][j].ConvFlux4) {
				return fmt.Sprintf("解更新，x方向(%d,%d)", i, j)
			}
		}
	}
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny+1; j++ {
			// muscl 插值重构Vl,Vr
			//Vl, Vr := muscl.Muscl(kappa, order, i, j, solver.mesh, solver.PositiveYlimiters[i][j],
			//	solver.NegativeYlimiters[i][j], solver.PositiveYlimiters[i][j+1], solver.NegativeYlimiters[i][j+1], 'y', solver.gamma)
			Vl := common.PrimtiveFlux{Density: solver.mesh.Mesh[i+Ng][j+Ng-1].Density, VelocityX: solver.mesh.Mesh[i+Ng][j+Ng-1].VelocityX,
				VelocityY: solver.mesh.Mesh[i+Ng][j+Ng-1].VelocityY, Pressure: solver.mesh.Mesh[i+Ng][j+Ng-1].Pressure}
			Vr := common.PrimtiveFlux{Density: solver.mesh.Mesh[i+Ng][j+Ng].Density, VelocityX: solver.mesh.Mesh[i+Ng][j+Ng].VelocityX,
				VelocityY: solver.mesh.Mesh[i+Ng][j+Ng].VelocityY, Pressure: solver.mesh.Mesh[i+Ng][j+Ng].Pressure}
			CenterFlux := CenterCFlux(Vl, Vr, solver.mesh.VectorY[i+Ng][j+Ng], solver.gamma)
			MatrixFlux := MatrixConFlux(Vl, Vr, solver.mesh.VectorY[i+Ng][j+Ng], solver.gamma)
			/**
			if i == 159 && j == 0 {
				fmt.Println("通量", CenterFlux, MatrixFlux)
				fmt.Println("Roe特征值", RoeValue(Vl, Vr, 1.4))
				fmt.Println("Roe dw", Amplitude(Vl, Vr, RoeValue(Vl, Vr, 1.4), solver.mesh.VectorY[i+Ng][j+Ng]))
				fmt.Println("Eigen value", EignValue(RoeValue(Vl, Vr, 1.4), solver.mesh.VectorY[i+Ng][j+Ng]))
				fmt.Println("Eign vector", RightEignVec(RoeValue(Vl, Vr, 1.4), solver.mesh.VectorY[i+Ng][j+Ng]))
				fmt.Println("vectorX", solver.mesh.VectorY[i+Ng][j+Ng])
			}
			**/
			solver.FluxY[i][j] = CenterFlux.FluxMinus(MatrixFlux)
			if math.IsNaN(solver.FluxY[i][j].ConvFlux1) || math.IsNaN(solver.FluxY[i][j].ConvFlux2) || math.IsNaN(solver.FluxY[i][j].ConvFlux3) || math.IsNaN(solver.FluxY[i][j].ConvFlux4) {
				return fmt.Sprintf("解更新，y方向(%d,%d)", i, j)
			}
		}
	}
	//fmt.Println("通量更新", solver.FluxX[50][50])
	//fmt.Println("通量更新", solver.FluxY[50][50])
	return ""
}

// Harten修正函数
func Harten(lambda float64, roe_ss float64) float64 {
	error := 0.25
	if math.Abs(lambda) <= error {
		lamba_mag := (lambda*lambda + error*error) / (2.0 * error)
		return lamba_mag
	}
	return math.Abs(lambda)
}

// 计算Roe均值变量
func RoeValue(Vl common.PrimtiveFlux, Vr common.PrimtiveFlux, gamma float64, n common.Vector) RoeAvarageFlux {
	var RoeAve RoeAvarageFlux
	var vel_mag float64
	ht_l := common.TotalEnthalpy(Vl, gamma)
	ht_r := common.TotalEnthalpy(Vr, gamma)
	r := math.Sqrt(Vr.Density / Vl.Density)
	Un_l := Vl.VelocityX*n.X + Vl.VelocityY*n.Y
	Un_r := Vr.VelocityX*n.X + Vr.VelocityY*n.Y
	RoeAve.Densitybar = r * Vl.Density
	RoeAve.VelocityXbar = (r*Vr.VelocityX + Vl.VelocityX) / (r + 1.0)
	RoeAve.VelocityYbar = (r*Vr.VelocityY + Vl.VelocityY) / (r + 1.0)
	RoeAve.Hbar = (r*ht_r + ht_l) / (r + 1.0)
	vel_mag = (RoeAve.VelocityXbar*RoeAve.VelocityXbar + RoeAve.VelocityYbar*RoeAve.VelocityYbar)
	c := (gamma - 1.0) * (RoeAve.Hbar - 0.5*vel_mag)
	RoeAve.Unbar = (Un_l*math.Sqrt(Vl.Density) + Un_r*math.Sqrt(Vr.Density)) / (math.Sqrt(Vl.Density) + math.Sqrt(Vr.Density))
	RoeAve.Cbar = math.Sqrt(c)
	return RoeAve
}

// 计算Roe平均量特征值
func EignValue(Roeave RoeAvarageFlux, n common.Vector) []float64 {
	lambda := make([]float64, 4)
	roe_c := Roeave.Cbar
	Un := Roeave.Unbar
	lambda[0] = Harten(Un, roe_c)
	lambda[1] = lambda[0]
	lambda[2] = Harten((Un + roe_c), roe_c)
	lambda[3] = Harten((Un - roe_c), roe_c)
	return lambda
}

// 计算 0.5*(F(w_r)+F(w_l))
func CenterCFlux(Vl common.PrimtiveFlux, Vr common.PrimtiveFlux, n common.Vector, gamma float64) common.ConvectiveFlux {
	dens_l := Vl.Density
	dens_r := Vr.Density
	velx_l := Vl.VelocityX
	velx_r := Vr.VelocityX
	vely_l := Vl.VelocityY
	vely_r := Vr.VelocityY
	press_l := Vl.Pressure
	press_r := Vr.Pressure
	Un_l := velx_l*n.X + vely_l*n.Y
	Un_r := velx_r*n.X + vely_r*n.Y
	ht_l := common.TotalEnthalpy(Vl, gamma)
	ht_r := common.TotalEnthalpy(Vr, gamma)
	CenterFlux := common.ConvectiveFlux{ConvFlux1: 0.5 * (dens_l*Un_l + dens_r*Un_r), ConvFlux2: 0.5 * (dens_l*velx_l*Un_l + n.X*press_l + dens_r*velx_r*Un_r + n.X*press_r),
		ConvFlux3: 0.5 * (dens_l*vely_l*Un_l + n.Y*press_l + dens_r*vely_r*Un_r + n.Y*press_r), ConvFlux4: 0.5 * (dens_l*Un_l*ht_l + dens_r*Un_r*ht_r)}
	//fmt.Println("中心通量", CenterFlux)
	return CenterFlux
}

// 计算 A^_roe *(Wr-Wl)
func MatrixConFlux(Vl common.PrimtiveFlux, Vr common.PrimtiveFlux, n common.Vector, gamma float64) common.ConvectiveFlux {
	var WaveFlux common.ConvectiveFlux
	var dw []float64
	Cflux := make([]float64, 4)
	Roevalue := RoeValue(Vl, Vr, gamma, n)
	lambda := EignValue(Roevalue, n)
	RightEign := RightEignVec(Roevalue, n)
	dw = Amplitude(Vl, Vr, Roevalue, n)
	// 矩阵相乘
	for i := 0; i < 4; i++ {
		Cflux[i] = 0
		for j := 0; j < 4; j++ {
			Cflux[i] += 0.5 * math.Abs(lambda[j]) * dw[j] * RightEign[i][j]
		}
	}
	WaveFlux = common.ConvectiveFlux{ConvFlux1: Cflux[0], ConvFlux2: Cflux[1], ConvFlux3: Cflux[2], ConvFlux4: Cflux[3]}
	//fmt.Println("波通量", WaveFlux)
	return WaveFlux
}

// 计算l单元,r单元之间的dw
func Amplitude(Vl common.PrimtiveFlux, Vr common.PrimtiveFlux, Roeave RoeAvarageFlux, n common.Vector) []float64 {
	dw := make([]float64, 4)
	d_dens := Vr.Density - Vl.Density
	d_velx := Vr.VelocityX - Vl.VelocityX
	d_vely := Vr.VelocityY - Vl.VelocityY
	d_press := Vr.Pressure - Vl.Pressure
	roe_dens := Roeave.Densitybar
	roe_c := Roeave.Cbar
	dw[0] = d_dens - d_press/(roe_c*roe_c)
	dw[1] = n.Y*d_velx - n.X*d_vely
	dw[2] = n.X*d_velx + n.Y*d_vely + d_press/(roe_dens*roe_c)
	dw[3] = n.X*d_velx + n.Y*d_vely - d_press/(roe_dens*roe_c)
	return dw
}

// 计算Roe矩阵的特征向量
func RightEignVec(Roeave RoeAvarageFlux, n common.Vector) [][]float64 {
	RightEign := make([][]float64, 4)
	for i := range RightEign {
		RightEign[i] = make([]float64, 4)
	}
	roe_dens := Roeave.Densitybar
	roe_velx := Roeave.VelocityXbar
	roe_vely := Roeave.VelocityYbar
	roe_ht := Roeave.Hbar
	roe_c := Roeave.Cbar
	Un := roe_velx*n.X + roe_vely*n.Y
	F := roe_dens / (2.0 * roe_c)
	RightEign[0][0] = 1.0
	RightEign[1][0] = roe_velx
	RightEign[2][0] = roe_vely
	RightEign[3][0] = 0.5 * (roe_velx*roe_velx + roe_vely*roe_vely)

	RightEign[0][1] = 0.0
	RightEign[1][1] = n.Y * roe_dens
	RightEign[2][1] = -1.0 * n.X * roe_vely
	RightEign[3][1] = roe_dens * (n.Y*roe_velx - n.X*roe_vely)

	RightEign[0][2] = F
	RightEign[1][2] = F * (roe_velx + n.X*roe_c)
	RightEign[2][2] = F * (roe_vely + n.Y*roe_c)
	RightEign[3][2] = F * (roe_ht + Un*roe_c)

	RightEign[0][3] = -F
	RightEign[1][3] = -F * (roe_velx - n.X*roe_c)
	RightEign[2][3] = -F * (roe_vely - n.Y*roe_c)
	RightEign[3][3] = -F * (roe_ht - Un*roe_c)
	return RightEign
}
