package Geometry

// 定义Ramp管道
type Ramp struct {
	xleft      float64 //管道左边界
	xright     float64 //管道右边界
	ybottom    float64 //管道下边界
	ytop       float64 //管道上边界
	rampBegin  float64 //斜坡起始点
	rampEnd    float64 //斜坡终止处
	rampHeight float64 //斜坡高度， 默认斜坡起始高度与下边界相同
}

func InitGeometry(Xl float64, Xr float64, Yb float64, Yt float64, rb float64, re float64, rh float64) Ramp {
	Ramp := Ramp{xleft: Xl, xright: Xr, ybottom: Yb, ytop: Yt, rampBegin: rb, rampEnd: re, rampHeight: rh}
	return Ramp
}

// 变量权限限制，通过函数读取
func (geo Ramp) GetXleft() float64 {
	return geo.xleft
}

func (geo Ramp) GetXright() float64 {
	return geo.xright
}

func (geo Ramp) GetYbottom() float64 {
	return geo.ybottom
}

func (geo Ramp) GetYtop() float64 {
	return geo.ytop
}

func (geo Ramp) GetrampBegin() float64 {
	return geo.rampBegin
}

func (geo Ramp) GetrampEnd() float64 {
	return geo.rampEnd
}

func (geo Ramp) GetrampHeight() float64 {
	return geo.rampHeight
}
