package boundarycond

import (
	"fmt"
	"math"

	common "git.oa.com/hengrushang/RampCfd/Common"
)

// 无滑移壁面边界条件更新
// 入参，网格对象
func SlipWall(mesh *common.SpaceMesh) string {
	// 上壁面
	for i := mesh.NG; i < mesh.NG+mesh.Nx; i++ {
		xn := mesh.VectorY[i][mesh.NG].X
		yn := mesh.VectorY[i][mesh.NG].Y
		a := xn*xn - yn*yn
		b := 2.0 * xn * yn
		for j := 0; j < mesh.NG; j++ {
			mesh.Mesh[i][j].VelocityX = -a*mesh.Mesh[i][mesh.NG].VelocityX - b*mesh.Mesh[i][mesh.NG].VelocityY
			mesh.Mesh[i][j].VelocityY = -b*mesh.Mesh[i][mesh.NG].VelocityX + a*mesh.Mesh[i][mesh.NG].VelocityY
			mesh.Mesh[i][j].Density = mesh.Mesh[i][mesh.NG].Density
			mesh.Mesh[i][j].Pressure = mesh.Mesh[i][mesh.NG].Pressure
			if math.IsNaN(mesh.Mesh[i][j].Density) || math.IsNaN(mesh.Mesh[i][j].Pressure) || math.IsNaN(mesh.Mesh[i][j].VelocityX) || math.IsNaN(mesh.Mesh[i][j].VelocityY) {
				return fmt.Sprintf("壁面边界下壁面，位置(%d,%d)", i, j)
			}
		}
	}
	// 下壁面
	for i := mesh.NG; i < mesh.NG+mesh.Nx; i++ {
		xn := mesh.VectorY[i][mesh.NG+mesh.Ny].X
		yn := mesh.VectorY[i][mesh.NG+mesh.Ny].Y
		a := xn*xn - yn*yn
		b := 2.0 * xn * yn
		for j := 0; j < mesh.NG; j++ {
			mesh.Mesh[i][mesh.NG+mesh.Ny+j].VelocityX = -a*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityX - b*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityY
			mesh.Mesh[i][mesh.NG+mesh.Ny+j].VelocityY = -b*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityX + a*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityY
			mesh.Mesh[i][j+mesh.NG+mesh.Ny].Density = mesh.Mesh[i][mesh.NG+mesh.Ny-1].Density
			mesh.Mesh[i][j+mesh.NG+mesh.Ny].Pressure = mesh.Mesh[i][mesh.NG+mesh.Ny-1].Pressure
			if math.IsNaN(mesh.Mesh[i][j].Density) || math.IsNaN(mesh.Mesh[i][j].Pressure) || math.IsNaN(mesh.Mesh[i][j].VelocityX) || math.IsNaN(mesh.Mesh[i][j].VelocityY) {
				return fmt.Sprintf("壁面边界上壁面，位置(%d,%d)", i, j)
			}
		}
	}
	return ""
}
