package initialcond

import (
	"math"

	common "git.oa.com/hengrushang/RampCfd/Common"
)

// 根据进口条件生成初场
// 入参，马赫数，气体参数，初始化之后的网格
// 气流垂直进入，无攻角情况
func Initialize(Mach float64, density float64, pressure float64, gamma float64, mesh *common.SpaceMesh) []float64 {
	inflow := make([]float64, 4)
	c := math.Sqrt(gamma * pressure / density)
	velocityX := c * Mach
	velocityY := 0.0
	lenghti := mesh.Nx + 2*mesh.NG
	lengthj := mesh.Ny + 2*mesh.NG
	inflow[0] = density
	inflow[1] = velocityX
	inflow[2] = velocityY
	inflow[3] = pressure
	for i := 1; i < lenghti-1; i++ {
		for j := 1; j < lengthj-1; j++ {
			mesh.Mesh[i][j].VelocityX = velocityX
			mesh.Mesh[i][j].VelocityY = velocityY
			mesh.Mesh[i][j].Density = density
			mesh.Mesh[i][j].Pressure = pressure
		}
	}
	return inflow
}
