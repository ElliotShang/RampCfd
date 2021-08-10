package common

// 定义每个网格单元的信息
type Cell struct {
	Xlocation float64 // 网格中心点
	Yloaction float64 // 网格中心点
	Density   float64
	Pressure  float64
	VelocityX float64
	VelocityY float64
}

// 网格节点位置，用于生成各单元中心点坐标
type Node struct {
	xNode float64
	yNode float64
}

// 定义方向向量,单元切向
type Vector struct {
	X float64
	Y float64
}

// 方向向量，单元法向
type VectorY struct {
	X float64
	Y float64
}

func (vec VectorY) GetVectorYX() float64 {
	return vec.X
}

func (vec *VectorY) SetVectorYX(x float64) {
	vec.X = x
}

func (vec VectorY) GetVectorYY() float64 {
	return vec.X
}
