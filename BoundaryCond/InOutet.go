package boundarycond

import (
	"fmt"
	"math"

	common "git.oa.com/hengrushang/RampCfd/Common"
)

// 进出口边界条件更新
// 超声速条件下进口边界全指定
// 入参，进口气流参数，网格
func SupersonicInlet(inflow common.PrimtiveFlux, mesh *common.SpaceMesh) string {
	for i := 0; i < mesh.NG; i++ {
		for j := mesh.NG; j < mesh.NG+mesh.Ny; j++ {
			mesh.Mesh[i][j].Density = inflow.Density
			mesh.Mesh[i][j].VelocityX = inflow.VelocityX
			mesh.Mesh[i][j].VelocityY = inflow.VelocityY
			mesh.Mesh[i][j].Pressure = inflow.Pressure
			if math.IsNaN(mesh.Mesh[i][j].Density) || math.IsNaN(mesh.Mesh[i][j].Pressure) || math.IsNaN(mesh.Mesh[i][j].VelocityX) || math.IsNaN(mesh.Mesh[i][j].VelocityY) {
				return fmt.Sprintf("进口边界，位置(%d,%d)", i, j)
			}
		}
	}
	return ""
}

func SupersonicOutlet(mesh *common.SpaceMesh) string {
	// lengthj := len(mesh.Mesh[0])
	Nx := mesh.Nx
	NG := mesh.NG
	for i := 0; i < mesh.NG; i++ {
		for j := mesh.NG; j < mesh.NG+mesh.Ny; j++ {
			mesh.Mesh[Nx+NG+i][j].Density = mesh.Mesh[Nx+NG-1][j].Density
			mesh.Mesh[Nx+NG+i][j].VelocityX = mesh.Mesh[Nx+NG-1][j].VelocityX
			mesh.Mesh[Nx+NG+i][j].VelocityY = mesh.Mesh[Nx+NG-1][j].VelocityY
			mesh.Mesh[Nx+NG+i][j].Pressure = mesh.Mesh[Nx+NG-1][j].Pressure
			if math.IsNaN(mesh.Mesh[i][j].Density) || math.IsNaN(mesh.Mesh[i][j].Pressure) || math.IsNaN(mesh.Mesh[i][j].VelocityX) || math.IsNaN(mesh.Mesh[i][j].VelocityY) {
				return fmt.Sprintf("出口边界，位置(%d,%d)", i, j)
			}
		}
	}
	return ""
}
