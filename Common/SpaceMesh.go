package common

import (
	"math"

	"git.oa.com/hengrushang/RampCfd/Geometry"
)

// 网格结构体定义
type SpaceMesh struct {
	Node    [][]Node    // 网格边节点
	Mesh    [][]Cell    // 网格
	Area    [][]float64 // 网格面积
	VectorX [][]vectorX // 流向网格间方向向量
	VectorY [][]vectorY // 法向网格间方向向量
	Nx      int         // x方向网格数量，不含虚拟点
	Ny      int         // y方向网格数量， 不含虚拟点
	NG      int         // 边界虚拟单元数量
}

// 网格初始化函数，使用几何法生成网格，方向向量，面积
// 入参 Nx:网格x方向数量 Ny: 网格y方向数量，NG: 边界虚拟点数量，geo: 所定义几何体
func InitSpaceMesh(Nx int, Ny int, NG int, geo *Geometry.Ramp) *SpaceMesh {
	node := make([][]Node, Nx+1+2*NG, Ny+1+2*NG)
	mesh := make([][]Cell, Nx+2*NG, Ny+2*NG)
	area := make([][]float64, Nx+2*NG, Ny+2*NG)
	vectorX := make([][]vectorX, Nx+2*NG+1, Ny+2*NG)
	vectorY := make([][]vectorY, Nx+2*NG, Ny+2*NG+1)
	length := geo.GetXright() - geo.GetXleft()
	horizon_ratio1 := (geo.GetrampBegin() - geo.GetXleft()) / length
	//ramp_ratio := (geo.GetrampEnd() - geo.GetrampBegin()) / length
	horizon_ratio2 := (geo.GetXright() - geo.GetrampEnd()) / length
	slope := (geo.GetrampHeight()) / (geo.GetrampEnd() - geo.GetrampBegin())
	// 生成网格节点位置 进出口边界单元中心定义在边界上，含虚拟单元节点
	for i := 0; i < Nx+2*NG+1; i++ {
		for j := 0; j < Ny+2*NG+1; j++ {
			if i < int(horizon_ratio1)*Nx+NG {
				dx := length / float64(Nx)
				dy := (geo.GetYtop() - geo.GetYbottom()) / float64(Ny)
				xi := geo.GetXleft() + float64(i)*dx - dx*float64(NG)
				yi := geo.GetYbottom() + float64(j)*dy - dy*float64(NG)
				node[i][j] = Node{xNode: xi, yNode: yi}
			} else if i > int(horizon_ratio2)*Nx+NG {
				dx := length / float64(Nx)
				dy := (geo.GetYtop() - geo.GetrampHeight()) / float64(Ny)
				xi := geo.GetXleft() + float64(i-NG)*dx
				yi := geo.GetYbottom() + float64(j-NG)*dy
				node[i][j] = Node{xNode: xi, yNode: yi}
			} else {
				dx := float64(Nx) / length
				yb := slope * (float64(i)*dx - geo.GetrampBegin())
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
			vectorX[i][j].x = dy / Distance(node[i][j], node[i][j+1])
			vectorX[i][j].y = -1.0 * dx / Distance(node[i][j], node[i][j+1])
		}
	}
	for i := 0; i < Nx+2*NG; i++ {
		for j := 0; j < Ny+2*NG+1; j++ {
			dx := node[i][j].xNode - node[i+1][j].xNode
			dy := node[i][j].yNode - node[i+1][j].yNode
			vectorY[i][j].x = dy / Distance(node[i][j], node[i+1][j])
			vectorY[i][j].y = -1.0 * dx / Distance(node[i][j], node[i+1][j])
		}
	}

	// 生成单元面积
	for i := 0; i < Nx+2*NG; i++ {
		for j := 0; j < Ny+2*NG; j++ {
			area[i][j] = Volume(node[i][j], node[i+1][j], node[i+1][j+1], node[i][j+1])
		}
	}

	Space_mesh := &SpaceMesh{Node: node, Mesh: mesh, Area: area, VectorX: vectorX, VectorY: vectorY, Nx: Nx, Ny: Ny, NG: NG}
	return Space_mesh
}

// 计算节点间距离
func Distance(node1 Node, node2 Node) float64 {
	dx := node2.xNode - node1.xNode
	dy := node2.yNode - node2.yNode
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
