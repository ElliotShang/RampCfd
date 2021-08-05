package solver

// 通用化求解器接口，适配不同求解格式
type Solve interface {
	UpdateFlux()     // 计算通量
	UpdateBC()       // 更新边界条件
	UpdateLimiter()  // 更新限制器
	InitalizeField() // 初始条件
	UpdateSolution() // 更新解
	UpdateResidual() // 更新残差
}
