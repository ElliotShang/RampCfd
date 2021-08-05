package boundarycond

import common "git.oa.com/hengrushang/RampCfd/Common"

// 无滑移壁面边界条件更新
// 入参，网格对象
func SlipWall(mesh *common.SpaceMesh) {
	// 上壁面
	for i := mesh.NG; i < mesh.NG+mesh.Nx; i++ {
		xn := mesh.VectorY[i][mesh.NG].GetVectorYX()
		yn := mesh.VectorY[i][mesh.NG].GetVectorYY()
		a := xn*xn - yn*yn
		b := 2.0 * xn * yn
		for j := 0; j < mesh.NG; j++ {
			mesh.Mesh[i][j].VelocityX = -a*mesh.Mesh[i][mesh.NG].VelocityX - b*mesh.Mesh[i][mesh.NG].VelocityY
			mesh.Mesh[i][j].VelocityY = -b*mesh.Mesh[i][mesh.NG].VelocityX + a*mesh.Mesh[i][mesh.NG].VelocityY
		}
	}
	// 下壁面
	for i := mesh.NG; i < mesh.NG+mesh.Nx; i++ {
		xn := mesh.VectorY[i][mesh.NG+mesh.Ny].GetVectorYX()
		yn := mesh.VectorY[i][mesh.NG+mesh.Ny].GetVectorYY()
		a := xn*xn - yn*yn
		b := 2.0 * xn * yn
		for j := 0; j < mesh.NG; j++ {
			mesh.Mesh[i][mesh.NG+mesh.Ny+j].VelocityX = -a*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityX - b*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityY
			mesh.Mesh[i][mesh.NG+mesh.Ny+j].VelocityY = -b*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityX + a*mesh.Mesh[i][mesh.NG+mesh.Ny-1].VelocityY
		}
	}
}
