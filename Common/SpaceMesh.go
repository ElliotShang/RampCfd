package common

import (
	"fmt"
	"math"

	"git.oa.com/hengrushang/RampCfd/Geometry"
)

// 网格结构体定义
type SpaceMesh struct {
	Node    [][]Node    // 网格边节点
	Mesh    [][]Cell    // 网格
	Area    [][]float64 // 网格面积
	VectorX [][]Vector  // 流向网格间方向向量
	LengthX [][]float64 // x方向单元长度
	VectorY [][]Vector  // 法向网格间方向向量
	LengthY [][]float64 // y方向单元长度
	Nx      int         // x方向网格数量，不含虚拟点
	Ny      int         // y方向网格数量， 不含虚拟点
	NG      int         // 边界虚拟单元数量
}

// 网格初始化函数，使用几何法生成网格，方向向量，面积
// 入参 Nx:网格x方向数量 Ny: 网格y方向数量，NG: 边界虚拟点数量，geo: 所定义几何体
func InitSpaceMesh(Nx int, Ny int, NG int, geo *Geometry.Ramp) *SpaceMesh {
	node := make([][]Node, Nx+1+2*NG)
	for i := range node {
		node[i] = make([]Node, Ny+1+2*NG)
	}
	mesh := make([][]Cell, Nx+2*NG)
	for i := range mesh {
		mesh[i] = make([]Cell, Ny+2*NG)
	}
	area := make([][]float64, Nx+2*NG)
	for i := range area {
		area[i] = make([]float64, Ny+2*NG)
	}
	vectorX := make([][]Vector, Nx+2*NG+1)
	for i := range vectorX {
		vectorX[i] = make([]Vector, Ny+2*NG)
	}
	LengthX := make([][]float64, Nx+2*NG+1)
	for i := range LengthX {
		LengthX[i] = make([]float64, Ny+2*NG)
	}
	vectorY := make([][]Vector, Nx+2*NG)
	for i := range vectorY {
		vectorY[i] = make([]Vector, Ny+2*NG+1)
	}
	LengthY := make([][]float64, Nx+2*NG)
	for i := range LengthY {
		LengthY[i] = make([]float64, Ny+2*NG+1)
	}
	length := geo.GetXright() - geo.GetXleft()
	horizon_ratio1 := (geo.GetrampBegin() - geo.GetXleft()) / length
	fmt.Println("mesh生成阶段第一段", horizon_ratio1)
	//ramp_ratio := (geo.GetrampEnd() - geo.GetrampBegin()) / length
	horizon_ratio2 := (geo.GetXright() - geo.GetrampEnd()) / length
	fmt.Println("mesh生成阶段第三段", horizon_ratio2)
	slope := (geo.GetrampHeight()) / (geo.GetrampEnd() - geo.GetrampBegin())
	fmt.Println("mesh生成阶段斜率", slope)
	// 生成网格节点位置 进出口边界单元中心定义在边界上，含虚拟单元节点
	for i := 0; i < Nx+2*NG+1; i++ {
		for j := 0; j < Ny+2*NG+1; j++ {
			if i < int(horizon_ratio1*float64(Nx))+NG {
				dx := length / float64(Nx)
				dy := (geo.GetYtop() - geo.GetYbottom()) / float64(Ny)
				xi := geo.GetXleft() + float64(i)*dx - dx*float64(NG)
				yi := geo.GetYbottom() + float64(j)*dy - dy*float64(NG)
				node[i][j] = Node{xNode: xi, yNode: yi}
			} else if i >= int((1-horizon_ratio2)*float64(Nx))+NG {
				dx := length / float64(Nx)
				dy := (geo.GetYtop() - geo.GetrampHeight()) / float64(Ny)
				xi := geo.GetXleft() + float64(i-NG)*dx
				yi := geo.GetrampHeight() + float64(j-NG)*dy
				node[i][j] = Node{xNode: xi, yNode: yi}
			} else {
				dx := length / float64(Nx)
				yb := slope * (float64(i-NG)*dx - geo.GetrampBegin())
				dy := (geo.GetYtop() - yb) / float64(Ny)
				xi := geo.GetXleft() + float64(i-NG)*dx
				yi := yb + float64(j-NG)*dy
				node[i][j] = Node{xNode: xi, yNode: yi}
			}
		}
	}
	// 根据节点初始化网格单元对象
	for i := 0; i < Nx+2*NG; i++ {
		for j := 0; j < Ny+2*NG; j++ {
			xcenter := 0.25 * (node[i][j].xNode + node[i][j+1].xNode + node[i+1][j].xNode + node[i+1][j+1].xNode)
			ycenter := 0.25 * (node[i][j].yNode + node[i+1][j].yNode + node[i][j+1].yNode + node[i+1][j+1].yNode)
			mesh[i][j] = Cell{Xlocation: xcenter, Yloaction: ycenter, Density: 0, VelocityX: 0, VelocityY: 0, Pressure: 0}
		}
	}

	// 生成单元间方向向量，距离
	for i := 0; i < Nx+2*NG+1; i++ {
		for j := 0; j < Ny+2*NG; j++ {
			dx := node[i][j+1].xNode - node[i][j].xNode
			dy := node[i][j+1].yNode - node[i][j].yNode
			vectorX[i][j].X = dy / Distance(node[i][j], node[i][j+1])
			vectorX[i][j].Y = -1.0 * dx / Distance(node[i][j], node[i][j+1])
		}
	}
	// 生成单元间距离
	for i := 0; i < Nx+2*NG+1; i++ {
		for j := 0; j < Ny+2*NG; j++ {
			LengthX[i][j] = Distance(node[i][j], node[i][j+1])
		}
	}
	for i := 0; i < Nx+2*NG; i++ {
		for j := 0; j < Ny+2*NG+1; j++ {
			dx := node[i][j].xNode - node[i+1][j].xNode
			dy := node[i][j].yNode - node[i+1][j].yNode
			vectorY[i][j].X = dy / Distance(node[i][j], node[i+1][j])
			vectorY[i][j].Y = -1.0 * dx / Distance(node[i][j], node[i+1][j])
		}
	}
	fmt.Println("网格生成,边界", vectorY[100][0].X, vectorY[100][0].Y)
	// 生成单元间距离
	for i := 0; i < Nx+2*NG; i++ {
		for j := 0; j < Ny+2*NG+1; j++ {
			LengthY[i][j] = Distance(node[i][j], node[i+1][j])
		}
	}
	// 生成单元面积
	for i := 0; i < Nx+2*NG; i++ {
		for j := 0; j < Ny+2*NG; j++ {
			area[i][j] = Volume(node[i][j], node[i+1][j], node[i+1][j+1], node[i][j+1])
		}
	}
	fmt.Println("网格面积", area[50][50])

	Space_mesh := &SpaceMesh{Node: node, Mesh: mesh, Area: area, VectorX: vectorX, LengthX: LengthX, VectorY: vectorY, LengthY: LengthY, Nx: Nx, Ny: Ny, NG: NG}
	return Space_mesh
}

// 计算节点间距离
func Distance(node1 Node, node2 Node) float64 {
	dx := node2.xNode - node1.xNode
	dy := node2.yNode - node1.yNode
	length := math.Sqrt(dx*dx + dy*dy)
	return length
}

// 计算单元面积
func Volume(node1 Node, node2 Node, node3 Node, node4 Node) float64 {
	xa := node1.xNode
	xb := node2.xNode
	xc := node3.xNode
	xd := node4.xNode
	ya := node1.yNode
	yb := node2.yNode
	yc := node3.yNode
	yd := node4.yNode
	cross := (xc-xa)*(yd-yb) - (yc-ya)*(xd-xb)
	return 0.5 * math.Abs(cross)
}
