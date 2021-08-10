package muscl

import common "git.oa.com/hengrushang/RampCfd/Common"

// muscl重构函数, 与限制器
func Muscl(kappa float64, order float64, i int, j int, mesh *common.SpaceMesh, limiterLP common.PrimtiveFlux, limiterLN common.PrimtiveFlux,
	limiterRP common.PrimtiveFlux, limiterRN common.PrimtiveFlux, sign rune, gamma float64) (Vl common.PrimtiveFlux, Vr common.PrimtiveFlux) {
	NG := mesh.NG
	if sign == 'x' {
		primFluxL2 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG-2][j+NG].Density, VelocityX: mesh.Mesh[i+NG-2][j+NG].VelocityX,
			VelocityY: mesh.Mesh[i+NG-2][j+NG].VelocityY, Pressure: mesh.Mesh[i+NG-2][j+NG].Pressure}
		primFluxL1 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG-1][j+NG].Density, VelocityX: mesh.Mesh[i+NG-1][j+NG].VelocityX,
			VelocityY: mesh.Mesh[i+NG-1][j+NG].VelocityY, Pressure: mesh.Mesh[i+NG-1][j+NG].Pressure}
		primFluxR1 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG][j+NG].Density, VelocityX: mesh.Mesh[i+NG][j+NG].VelocityX,
			VelocityY: mesh.Mesh[i+NG][j+NG].VelocityY, Pressure: mesh.Mesh[i+NG][j+NG].Pressure}
		primFluxR2 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG+1][j+NG].Density, VelocityX: mesh.Mesh[i+NG+1][j+NG].VelocityX,
			VelocityY: mesh.Mesh[i+NG+1][j+NG].VelocityY, Pressure: mesh.Mesh[i+NG+1][j+NG].Pressure}
		central := primFluxR1.FluxMinus(primFluxL1)
		upwindl := primFluxL1.FluxMinus(primFluxL2)
		upwindr := primFluxR2.FluxMinus(primFluxR1)
		minusfluxL := central.LimiterMulti(limiterLN).ScalarMultiFlux(1 + kappa)
		plusfluxL := upwindl.LimiterMulti(limiterLP).ScalarMultiFlux(1 - kappa)
		sumL := minusfluxL.FluxPlus(plusfluxL).ScalarMultiFlux(0.25 * (order - 1))
		Vl = primFluxL1.FluxPlus(sumL)
		minusfluxR := upwindr.LimiterMulti(limiterRN).ScalarMultiFlux(1 - kappa)
		plusfluxR := central.LimiterMulti(limiterRP).ScalarMultiFlux(1 + kappa)
		sumR := minusfluxR.FluxPlus(plusfluxR).ScalarMultiFlux(0.25 * (order - 1))
		Vr = primFluxR1.FluxMinus(sumR)
	} else {
		primFluxL2 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG][j+NG-2].Density, VelocityX: mesh.Mesh[i+NG][j+NG-2].VelocityX,
			VelocityY: mesh.Mesh[i+NG][j+NG-2].VelocityY, Pressure: mesh.Mesh[i+NG][j+NG-2].Pressure}
		primFluxL1 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG][j+NG-1].Density, VelocityX: mesh.Mesh[i+NG][j+NG-1].VelocityX,
			VelocityY: mesh.Mesh[i+NG][j+NG-1].VelocityY, Pressure: mesh.Mesh[i+NG][j+NG-1].Pressure}
		primFluxR1 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG][j+NG].Density, VelocityX: mesh.Mesh[i+NG][j+NG].VelocityX,
			VelocityY: mesh.Mesh[i+NG][j+NG].VelocityY, Pressure: mesh.Mesh[i+NG][j+NG].Pressure}
		primFluxR2 := common.PrimtiveFlux{Density: mesh.Mesh[i+NG][j+NG+1].Density, VelocityX: mesh.Mesh[i+NG][j+NG+1].VelocityX,
			VelocityY: mesh.Mesh[i+NG][j+NG+1].VelocityY, Pressure: mesh.Mesh[i+NG][j+NG+1].Pressure}
		central := primFluxR1.FluxMinus(primFluxL1)
		upwindl := primFluxL1.FluxMinus(primFluxL2)
		upwindr := primFluxR2.FluxMinus(primFluxR1)
		minusfluxL := central.LimiterMulti(limiterLN).ScalarMultiFlux(1 + kappa)
		plusfluxL := upwindl.LimiterMulti(limiterLP).ScalarMultiFlux(1 - kappa)
		sumL := minusfluxL.FluxPlus(plusfluxL).ScalarMultiFlux(0.25 * (order - 1))
		Vl = primFluxL1.FluxPlus(sumL)
		minusfluxR := upwindr.LimiterMulti(limiterRN).ScalarMultiFlux(1 - kappa)
		plusfluxR := central.LimiterMulti(limiterLP).ScalarMultiFlux(1 + kappa)
		sumR := minusfluxR.FluxPlus(plusfluxR).ScalarMultiFlux(0.25 * (order - 1))
		Vr = primFluxR1.FluxMinus(sumR)
	}
	return
}

// 限制器函数
// 入参：网格，网格单元index，加减判定
func LimiterCalcX(mesh *common.SpaceMesh, i int, j int, sign rune) common.PrimtiveFlux {
	var limiterFlux common.PrimtiveFlux
	if sign == '+' {
		primiVarR := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG].Pressure}
		primiVarL1 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].Pressure}
		primiVarL2 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].Pressure}
		numerator := primiVarR.FluxMinus(primiVarL1)
		denumerator := primiVarL1.FluxMinDenominator(primiVarL2)
		limiterFlux = numerator.FluxDivision(denumerator)
		limiterFlux = common.LimiterFun(limiterFlux)
	} else {
		primiVarR := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG].Pressure}
		primiVarL1 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG-1][j+mesh.NG].Pressure}
		primiVarL2 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG-2][j+mesh.NG].Pressure}
		numerator := primiVarL1.FluxMinus(primiVarL2)
		denumerator := primiVarR.FluxMinDenominator(primiVarL1)
		limiterFlux = numerator.FluxDivision(denumerator)
		limiterFlux = common.LimiterFun(limiterFlux)
	}
	return limiterFlux
}

func LimiterCalcY(mesh *common.SpaceMesh, i int, j int, sign rune) common.PrimtiveFlux {
	var limiterFlux common.PrimtiveFlux
	if sign == '+' {
		primiVarR := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG].Pressure}
		primiVarL1 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].Pressure}
		primiVarL2 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].Pressure}
		numerator := primiVarR.FluxMinus(primiVarL1)
		denumerator := primiVarL1.FluxMinDenominator(primiVarL2)
		limiterFlux = numerator.FluxDivision(denumerator)
		limiterFlux = common.LimiterFun(limiterFlux)
	} else {
		primiVarR := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG].Pressure}
		primiVarL1 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG-1].Pressure}
		primiVarL2 := common.PrimtiveFlux{Density: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].Density,
			VelocityX: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].VelocityX, VelocityY: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].VelocityY, Pressure: mesh.Mesh[i+mesh.NG][j+mesh.NG-2].Pressure}
		numerator := primiVarL1.FluxMinus(primiVarL2)
		denumerator := primiVarR.FluxMinDenominator(primiVarL1)
		limiterFlux = numerator.FluxDivision(denumerator)
		limiterFlux = common.LimiterFun(limiterFlux)
	}
	return limiterFlux
}
