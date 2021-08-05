package roesolver

import common "git.oa.com/hengrushang/RampCfd/Common"

// Roe求解器
type Roe struct {
	mesh              *common.SpaceMesh       // 网格
	PositiveXlimiters [][]common.PrimtiveFlux // 限制器通量
	PositiveYlimiters [][]common.PrimtiveFlux // +
	NegativeXlimiters [][]common.PrimtiveFlux // -
	NegativeYlimiters [][]common.PrimtiveFlux // -
	FluxX             [][]common.Flux         // x方向守恒型通量
	FluxY             [][]common.Flux         // y方向守恒性通量
	gamma             float64                 // 气体参数
}

func SolverInitlize(mesh *common.SpaceMesh, gamma float64) *Roe {
	NG := mesh.NG
	Nx := mesh.Nx
	Ny := mesh.Ny
	PositiveXlimiters := make([][]common.PrimtiveFlux, Nx+NG, Ny)
	NegativeXlimiters := make([][]common.PrimtiveFlux, Nx+NG, Ny)
	PositiveYlimiters := make([][]common.PrimtiveFlux, Nx, Ny+NG)
	NegativeYlimiters := make([][]common.PrimtiveFlux, Nx, Ny+NG)
	FluxX := make([][]common.Flux, Nx+1, Ny)
	FluxY := make([][]common.Flux, Nx, Ny+1)
	Roe := Roe{mesh: mesh, PositiveXlimiters: PositiveXlimiters, PositiveYlimiters: PositiveYlimiters,
		NegativeXlimiters: NegativeXlimiters, NegativeYlimiters: NegativeYlimiters, FluxX: FluxX, FluxY: FluxY, gamma: gamma}
	return &Roe
}

func (solver *Roe) UpdateBC() {

}

func (solver *Roe) UpdateFlux() {

}

func (solver *Roe) UpdateLimiter() {

}

func (solver *Roe) UpdateResidual() {

}

func (slover *Roe) UpdateSolution() {

}
