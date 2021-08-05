package common

import "math"

// 定义相关辅助函数与数据类型
// 守恒性通量类型
type Flux struct {
	Density float64 // rho
	MomX    float64 // rho*u
	MomY    float64 // rho*v
	Energy  float64 // rho*E
}

type PrimtiveFlux struct {
	Density   float64
	VelocityX float64
	VelocityY float64
	Pressure  float64
}

func (flux Flux) Conv2Prim(gamma float64) PrimtiveFlux {
	density := flux.Density
	momX := flux.MomX
	momY := flux.MomY
	totale := flux.Energy
	pr := PrimtiveFlux{Density: density, VelocityX: momX / density, VelocityY: momY / density, Pressure: (gamma - 1.0) * (totale - 0.5*(momX*momX+momY*momY)/density)}
	return pr
}

func (pr PrimtiveFlux) Prim2Conv(gamma float64) Flux {
	var flux Flux
	density := pr.Density
	velocityX := pr.VelocityX
	velocityY := pr.VelocityY
	totale := TotalEnergy(pr, gamma)
	flux = Flux{Density: density, MomX: density * velocityX, MomY: density * velocityY, Energy: density * totale}
	return flux
}

func (flux PrimtiveFlux) FluxMinus(flux1 PrimtiveFlux) PrimtiveFlux {
	var res PrimtiveFlux
	res = PrimtiveFlux{Density: flux.Density - flux1.Density, VelocityX: flux.VelocityX - flux1.VelocityX, VelocityY: flux.VelocityY - flux1.VelocityY, Pressure: flux.Pressure - flux1.Pressure}
	return res
}

func (flux Flux) FluxMinus(flux1 Flux) Flux {
	var res Flux
	res = Flux{Density: flux.Density - flux1.Density, MomX: flux.MomX - flux1.MomX, MomY: flux.MomY - flux1.MomY, Energy: flux.Energy - flux1.Energy}
	return res
}

func (flux Flux) FluxPlus(flux1 Flux) Flux {
	var res Flux
	res = Flux{Density: flux.Density + flux1.Density, MomX: flux.MomX + flux1.MomX, MomY: flux.MomY + flux1.MomY, Energy: flux.Energy + flux1.Energy}
	return res
}

func (flux PrimtiveFlux) FluxPlus(flux1 PrimtiveFlux) PrimtiveFlux {
	var res PrimtiveFlux
	res = PrimtiveFlux{Density: flux.Density + flux1.Density, VelocityX: flux.VelocityX + flux1.VelocityX, VelocityY: flux.VelocityY + flux1.VelocityY, Pressure: flux.Pressure + flux1.Pressure}
	return res
}

func (flux Flux) ScalarMultiFlux(saclar float64) Flux {
	var res Flux
	res = Flux{Density: saclar * flux.Density, MomX: saclar * flux.MomX, MomY: saclar * flux.MomY, Energy: saclar * flux.Energy}
	return res
}

func (flux PrimtiveFlux) ScalarMultiFlux(saclar float64) PrimtiveFlux {
	var res PrimtiveFlux
	res = PrimtiveFlux{Density: saclar * flux.Density, VelocityX: saclar * flux.VelocityX, VelocityY: saclar * flux.VelocityY, Pressure: saclar * flux.Pressure}
	return res
}

// 作为限制器分母，相减的同时确保分母不为0
func (flux PrimtiveFlux) FluxMinDenominator(flux1 PrimtiveFlux) PrimtiveFlux {
	var res PrimtiveFlux
	epsilon := 1.0e-12
	densitydiff := flux.Density - flux1.Density
	signdensitydiff := epsilon * densitydiff / (math.Abs(densitydiff))
	if math.Abs(densitydiff) < epsilon {
		densitydiff = signdensitydiff
	}
	velocityXdiff := flux.VelocityX - flux1.VelocityX
	signvelodiff := epsilon * velocityXdiff / (math.Abs(velocityXdiff))
	if math.Abs(velocityXdiff) < epsilon {
		velocityXdiff = signvelodiff
	}
	velocityYdiff := flux.VelocityY - flux1.VelocityY
	signvloydiff := epsilon * velocityYdiff / (math.Abs(velocityYdiff))
	if math.Abs(velocityYdiff) < epsilon {
		velocityYdiff = signvloydiff
	}
	pressurediff := flux.Pressure - flux1.Pressure
	signpressureDiff := epsilon * pressurediff / (math.Abs(pressurediff))
	if math.Abs(pressurediff) < epsilon {
		pressurediff = signpressureDiff
	}
	res = PrimtiveFlux{Density: densitydiff, VelocityX: velocityXdiff, VelocityY: velocityYdiff, Pressure: pressurediff}
	return res
}

// 限制器向量与守恒性变量相乘
// 入参，限制器向量
func (flux PrimtiveFlux) LimiterMulti(flux1 PrimtiveFlux) PrimtiveFlux {
	var res PrimtiveFlux
	res = PrimtiveFlux{flux.Density * flux1.Density, flux.VelocityX * flux1.VelocityX, flux.VelocityY * flux1.VelocityY, flux.Pressure * flux1.Pressure}
	return res
}

// 在通量限制器中用于计算r
func (flux PrimtiveFlux) FluxDivision(flux1 PrimtiveFlux) PrimtiveFlux {
	var res PrimtiveFlux
	densityDiv := flux.Density / flux1.Density
	velocityXDiv := flux.VelocityX / flux1.VelocityX
	velocityYDiv := flux.VelocityY / flux1.VelocityY
	pressureDiv := flux.Pressure / flux1.Pressure
	res = PrimtiveFlux{Density: densityDiv, VelocityX: velocityXDiv, VelocityY: velocityYDiv, Pressure: pressureDiv}
	return res
}

// limiter计算函数
// 入参：已经计算得到的r. 返回：phi(r)
func LimiterFun(flux PrimtiveFlux) PrimtiveFlux {
	var res PrimtiveFlux
	mindensityR := flux.Density
	minvelocityxR := flux.VelocityX
	minvelocityyR := flux.VelocityY
	minpressure := flux.Pressure
	res.Density = math.Max(math.Min(mindensityR, 1.0), 0)
	res.VelocityX = math.Max(math.Min(minvelocityxR, 1.0), 0)
	res.VelocityY = math.Max(math.Min(minvelocityyR, 1.0), 0)
	res.Pressure = math.Max(math.Min(minpressure, 1.0), 0.0)
	return res
}

func TotalEnergy(pr PrimtiveFlux, gamma float64) float64 {
	density := pr.Density
	velocityX := pr.VelocityX
	velocityY := pr.VelocityY
	pressure := pr.Pressure
	totale := (gamma/(gamma-1.0))*(pressure/density) + 0.5*(velocityX*velocityX+velocityY*velocityY)
	return totale
}
