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
type vectorX struct {
	x float64
	y float64
}

func (vec vectorX) GetVectorXX() float64 {
	return vec.x
}

func (vec vectorX) GetVectorXY() float64 {
	return vec.y
}

// 方向向量，单元法向
type vectorY struct {
	x float64
	y float64
}

func (vec vectorY) GetVectorYX() float64 {
	return vec.x
}

func (vec vectorY) GetVectorYY() float64 {
	return vec.y
}
