package timestep

import (
	"fmt"
	"math"

	common "git.oa.com/hengrushang/RampCfd/Common"
)

// 自适应当地时间步长计算
func CalculateTimestep(mesh *common.SpaceMesh, cfl float64, gamma float64) ([][]float64, string) {
	Nx := mesh.Nx
	Ny := mesh.Ny
	NG := mesh.NG
	var avg_xnormal common.Vector
	var avg_ynormal common.Vector
	var avg_lengthx float64
	var avg_lengthy float64
	dt := make([][]float64, Nx)
	for i := range dt {
		dt[i] = make([]float64, Ny)
	}
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			avg_lengthx = 0.5 * (mesh.LengthX[i+NG][j+NG] + mesh.LengthX[i+1+NG][j+NG])
			avg_lengthy = 0.5 * (mesh.LengthY[i+NG][j+NG] + mesh.LengthY[i+NG][j+1+NG])
			avg_xnormal.X = 0.5 * (mesh.VectorX[i+NG][j+NG].X + mesh.VectorX[i+1+NG][j+NG].X)
			avg_xnormal.Y = 0.5 * (mesh.VectorX[i+NG][j+NG].Y + mesh.VectorX[i+1+NG][j+NG].Y)
			avg_ynormal.X = 0.5 * (mesh.VectorY[i+NG][j+NG].X + mesh.VectorY[i+NG][j+1+NG].X)
			avg_ynormal.Y = 0.5 * (mesh.VectorY[i+NG][j+NG].Y + mesh.VectorY[i+NG][j+1+NG].Y)
			//xlambda := MaxEignvalue(mesh, i+NG, j+NG, avg_xnormal, gamma)
			//ylambda := MaxEignvalue(mesh, i+NG, j+NG, avg_ynormal, gamma)
			// 计算当地时间步长
			un := avg_xnormal.X*mesh.Mesh[i][j].VelocityX + avg_xnormal.Y*mesh.Mesh[i][j].VelocityY
			vn := avg_ynormal.X*mesh.Mesh[i][j].VelocityX + avg_ynormal.Y*mesh.Mesh[i][j].VelocityY
			//dneum := xlambda*avg_lengthx + ylambda*avg_lengthy
			flux := common.PrimtiveFlux{Density: mesh.Mesh[i+NG][j+NG].Density, VelocityX: mesh.Mesh[i+NG][j+NG].VelocityX,
				VelocityY: mesh.Mesh[i+NG][j+NG].VelocityY, Pressure: mesh.Mesh[i+NG][j+NG].Pressure}
			c := common.SoundSpeed(flux, gamma)
			dneum := c*(avg_lengthx+avg_lengthy) + math.Abs(un) + math.Abs(vn)
			numer := cfl
			dt[i][j] = numer / dneum
			/**
			if i == 159 && j == 0 {
				fmt.Println("分子分母", numer, dneum)
				fmt.Println("特征值", xlambda, ylambda)
				fmt.Println("x均值,y均值", avg_lengthx, avg_lengthy)
				fmt.Println("X方向，y方向", avg_xnormal, avg_ynormal)
				fmt.Println("声速相关", mesh.Mesh[i+NG][j+NG].Density, mesh.Mesh[i+NG][j+NG].Pressure)
				fmt.Println("位置信息", mesh.VectorX[i+NG][j+NG], mesh.VectorY[i+NG][j+NG])
			}
			**/
			// fmt.Println(numer, dneum)
			if math.IsNaN(dt[i][j]) {
				return dt, fmt.Sprintf("时间步更新阶段，网格位置为(%d,%d),网格面积%f, x，y长度为%f,%f\n", i, j, mesh.Area[i+NG][j+NG], avg_lengthx, avg_lengthy)
			}
		}
	}
	return dt, ""
}

// 计算相应单元的欧拉方程组特征值 X方向
func MaxEignvalue(mesh *common.SpaceMesh, i int, j int, vector common.Vector, gamma float64) float64 {
	vx := mesh.Mesh[i][j].VelocityX
	vy := mesh.Mesh[i][j].VelocityY
	nx := vector.X
	ny := vector.Y
	flux := common.PrimtiveFlux{Density: mesh.Mesh[i][j].Density, VelocityX: mesh.Mesh[i][j].VelocityX, VelocityY: mesh.Mesh[i][j].VelocityY, Pressure: mesh.Mesh[i][j].Pressure}
	c := common.SoundSpeed(flux, gamma)
	lambda := math.Abs(nx*vx+ny*vy) + c
	return lambda
}

// 计算相应单元的欧拉方程组特征值 Y方向
func MaxEignvalueY(mesh *common.SpaceMesh, i int, j int, vector common.Vector, gamma float64) float64 {
	vx := mesh.Mesh[i][j].VelocityX
	vy := mesh.Mesh[i][j].VelocityY
	nx := vector.X
	ny := vector.Y
	flux := common.PrimtiveFlux{Density: mesh.Mesh[i][j].Density, VelocityX: mesh.Mesh[i][j].VelocityX, VelocityY: mesh.Mesh[i][j].VelocityY, Pressure: mesh.Mesh[i][j].Pressure}
	c := common.SoundSpeed(flux, gamma)
	lambda := math.Abs(nx*vx+ny*vy) + c
	return lambda
}
