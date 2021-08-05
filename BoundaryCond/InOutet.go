package boundarycond

import (
	common "git.oa.com/hengrushang/RampCfd/Common"
)

// 进出口边界条件更新
// 超声速条件下进口边界全指定
// 入参，进口气流参数，网格
func SupersonicInlet(inflow []float64, mesh *common.SpaceMesh) {
	lengthj := len(mesh.Mesh[0])
	for i := 0; i < mesh.NG; i++ {
		for j := 0; j < lengthj; j++ {
			mesh.Mesh[i][j].Density = inflow[0]
			mesh.Mesh[i][j].VelocityX = inflow[1]
			mesh.Mesh[i][j].VelocityY = inflow[2]
			mesh.Mesh[i][j].Pressure = inflow[3]
		}
	}
}

func SupersonicOutlet(mesh *common.SpaceMesh) {
	lengthj := len(mesh.Mesh[0])
	Nx := mesh.Nx
	NG := mesh.NG
	for i := 0; i < mesh.NG; i++ {
		for j := 0; j < lengthj; j++ {
			mesh.Mesh[Nx+NG+i][j] = mesh.Mesh[Nx+NG-1-i][j]
		}
	}
}
