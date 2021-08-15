package Solver

import (
	"fmt"
	"math"

	boundarycond "git.oa.com/hengrushang/RampCfd/BoundaryCond"
	common "git.oa.com/hengrushang/RampCfd/Common"
	initialcond "git.oa.com/hengrushang/RampCfd/InitialCond"
	muscl "git.oa.com/hengrushang/RampCfd/MUSCL"
	timestep "git.oa.com/hengrushang/RampCfd/TimeStep"
)

// 通用化求解器接口，适配不同求解格式,实现复用性
type Solve interface {
	UpdateFlux(kappa float64, order float64)                                       // 计算通量
	UpdateBC()                                                                     // 更新边界条件
	UpdateLimiter()                                                                // 更新限制器
	InitalizeField(Mach float64, density float64, pressure float64, gamma float64) // 初始化流场
	UpdateSolution(gamma float64)                                                  // 更新解
	UpdateResidual()                                                               // 更新残差
	UpdateTimeStep(cfl float64)                                                    //更新时间步长
}

// 求解器对象
type Solver struct {
	mesh              *common.SpaceMesh         // 网格
	InitialField      common.PrimtiveFlux       //初场均匀变量
	PositiveXlimiters [][]common.PrimtiveFlux   // 限制器通量+
	PositiveYlimiters [][]common.PrimtiveFlux   // +
	NegativeXlimiters [][]common.PrimtiveFlux   // -
	NegativeYlimiters [][]common.PrimtiveFlux   // -
	FluxX             [][]common.ConvectiveFlux // x方向守恒型通量
	FluxY             [][]common.ConvectiveFlux // y方向守恒性通量
	Residual          [][]common.ConvectiveFlux // 残差变量
	dt                [][]float64               // 网格时间步长
	gamma             float64                   // 气体参数
}

// 求解器初始化
func SolverInitlize(mesh *common.SpaceMesh, gamma float64) *Solver {
	NG := mesh.NG
	Nx := mesh.Nx
	Ny := mesh.Ny
	PositiveXlimiters := make([][]common.PrimtiveFlux, Nx+NG)
	for i := range PositiveXlimiters {
		PositiveXlimiters[i] = make([]common.PrimtiveFlux, Ny)
	}
	NegativeXlimiters := make([][]common.PrimtiveFlux, Nx+NG)
	for i := range NegativeXlimiters {
		NegativeXlimiters[i] = make([]common.PrimtiveFlux, Ny)
	}
	PositiveYlimiters := make([][]common.PrimtiveFlux, Nx)
	for i := range PositiveYlimiters {
		PositiveYlimiters[i] = make([]common.PrimtiveFlux, Ny+NG)
	}
	NegativeYlimiters := make([][]common.PrimtiveFlux, Nx)
	for i := range NegativeYlimiters {
		NegativeYlimiters[i] = make([]common.PrimtiveFlux, Ny+NG)
	}
	FluxX := make([][]common.ConvectiveFlux, Nx+1)
	for i := range FluxX {
		FluxX[i] = make([]common.ConvectiveFlux, Ny)
	}
	FluxY := make([][]common.ConvectiveFlux, Nx)
	for i := range FluxY {
		FluxY[i] = make([]common.ConvectiveFlux, Ny+1)
	}
	Residual := make([][]common.ConvectiveFlux, Nx)
	for i := range Residual {
		Residual[i] = make([]common.ConvectiveFlux, Ny)
	}
	dt := make([][]float64, Nx)
	for i := range dt {
		dt[i] = make([]float64, Ny)
	}
	solver := Solver{mesh: mesh, PositiveXlimiters: PositiveXlimiters, PositiveYlimiters: PositiveYlimiters,
		NegativeXlimiters: NegativeXlimiters, NegativeYlimiters: NegativeYlimiters, FluxX: FluxX, FluxY: FluxY, Residual: Residual, dt: dt, gamma: gamma}
	return &solver
}

func (solver *Solver) GetMesh() *common.SpaceMesh {
	return solver.mesh
}

func (solver *Solver) InitalizeField(Mach float64, density float64, pressure float64, gamma float64) {
	initfield := initialcond.Initialize(Mach, density, pressure, gamma, solver.mesh)
	solver.InitialField.Density = initfield[0]
	solver.InitialField.VelocityX = initfield[1]
	solver.InitialField.VelocityY = initfield[2]
	solver.InitialField.Pressure = initfield[3]
}

// 更新边界条件
func (solver *Solver) UpdateBC() string {
	// 更新进出口边界
	err1 := boundarycond.SupersonicInlet(solver.InitialField, solver.mesh)
	err2 := boundarycond.SupersonicOutlet(solver.mesh)
	// 更新壁面条件
	err3 := boundarycond.SlipWall(solver.mesh)
	// fmt.Println("边界", solver.mesh.Mesh[0][50], solver.mesh.Mesh[50][0], solver.mesh.Mesh[100][0], solver.mesh.Mesh[200][0])
	return err1 + err2 + err3
}

// 更新限制器
func (solver *Solver) UpdateLimiter() {
	Nx := solver.mesh.Nx
	//Ng := solver.mesh.NG
	Ny := solver.mesh.Ny
	//更新x方向限制器
	for i := 2; i < Nx; i++ {
		for j := 2; j < Ny-1; j++ {
			solver.PositiveXlimiters[i][j] = muscl.LimiterCalcX(solver.mesh, i, j, '+')
			solver.NegativeXlimiters[i][j] = muscl.LimiterCalcX(solver.mesh, i, j, '-')
		}
	}
	//更新y方向限制器
	for i := 2; i < Nx-1; i++ {
		for j := 2; j < Ny; j++ {
			solver.PositiveYlimiters[i][j] = muscl.LimiterCalcY(solver.mesh, i, j, '+')
			solver.NegativeYlimiters[i][j] = muscl.LimiterCalcY(solver.mesh, i, j, '-')
		}
	}
}

// 通用化抽象求解器接口不做具体实现，其他接口继承Solver对象后进行重写实现功能
func (solver *Solver) UpdateFlux(order float64) {
}

// 残差更新
func (solver *Solver) UpdateResidual() string {
	Nx := solver.mesh.Nx
	Ny := solver.mesh.Ny
	NG := solver.mesh.NG
	var xin common.ConvectiveFlux
	var xout common.ConvectiveFlux
	var yin common.ConvectiveFlux
	var yout common.ConvectiveFlux
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			xin = solver.FluxX[i][j].ScalarMultiFlux(solver.mesh.LengthX[i+NG][j+NG])
			xout = solver.FluxX[i+1][j].ScalarMultiFlux(solver.mesh.LengthX[i+NG+1][j+NG])
			yin = solver.FluxY[i][j].ScalarMultiFlux(solver.mesh.LengthY[i+NG][j+NG])
			yout = solver.FluxY[i][j+1].ScalarMultiFlux(solver.mesh.LengthY[i+NG][j+NG+1])
			tempX := xout.FluxMinus(xin)
			tempY := yout.FluxMinus(yin)
			solver.Residual[i][j] = tempX.FluxPlus(tempY)
			if math.IsNaN(solver.Residual[i][j].ConvFlux1) || math.IsNaN(solver.Residual[i][j].ConvFlux2) || math.IsNaN(solver.Residual[i][j].ConvFlux3) || math.IsNaN(solver.Residual[i][j].ConvFlux4) {
				return fmt.Sprintf("残差更新阶段，位置为(%d,%d)", i, j)
			}
		}
	}
	//fmt.Println("残差单元", solver.Residual[163][0])
	return ""
}

// 更新时间步长
func (solver *Solver) UpdateTimeStep(cfl float64) (errReply string) {
	solver.dt, errReply = timestep.CalculateTimestep(solver.mesh, cfl, solver.gamma)
	return errReply
	// fmt.Println("时间", solver.dt[50][50])
}

// 更新解
func (solver *Solver) UpdateSolution(gamma float64) string {
	Nx := solver.mesh.Nx
	Ny := solver.mesh.Ny
	NG := solver.mesh.NG
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			// 读取当前单元流场变量
			tempFlux := common.PrimtiveFlux{Density: solver.mesh.Mesh[i+NG][j+NG].Density, VelocityX: solver.mesh.Mesh[i+NG][j+NG].VelocityX,
				VelocityY: solver.mesh.Mesh[i+NG][j+NG].VelocityY, Pressure: solver.mesh.Mesh[i+NG][j+NG].Pressure}
			convser_temp := tempFlux.Prim2Conv(solver.gamma)
			temp1 := solver.Residual[i][j].ScalarMultiFlux(solver.dt[i][j])
			// temp2 := temp1.ScalarMultiFlux(1.0 / solver.mesh.Area[i+NG][j+NG])
			convser_temp1 := convser_temp.FluxMinusConvc(temp1)
			tempPr := convser_temp1.Conv2Prim(solver.gamma)
			solver.mesh.Mesh[i+NG][j+NG].Density = tempPr.Density
			solver.mesh.Mesh[i+NG][j+NG].VelocityX = tempPr.VelocityX
			solver.mesh.Mesh[i+NG][j+NG].VelocityY = tempPr.VelocityY
			solver.mesh.Mesh[i+NG][j+NG].Pressure = tempPr.Pressure
			if math.IsNaN(solver.mesh.Mesh[i+NG][j+NG].Density) || math.IsNaN(solver.mesh.Mesh[i+NG][j+NG].Pressure) || math.IsNaN(solver.mesh.Mesh[i+NG][j+NG].VelocityX) || math.IsNaN(solver.mesh.Mesh[i+NG][j+NG].VelocityY) {
				return fmt.Sprintf("解更新阶段，网格位置为(%d,%d)", i, j)
			}
			if solver.mesh.Mesh[i+NG][j+NG].Pressure <= 0 {
				return fmt.Sprintf("解更新阶段,压力小于零，网格位置为(%d,%d)", i, j)
			}
			if solver.mesh.Mesh[i+NG][j+NG].Density <= 0 {
				return fmt.Sprintf("解更新阶段,压力小于零，网格位置为(%d,%d)", i, j)
			}
		}

	}
	return ""
}

func ResidualNorm() {

}
