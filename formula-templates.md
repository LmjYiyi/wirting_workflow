---
description: 论文核心公式的LaTeX格式模板
---

# 公式模板

> **完整公式推导**请查看：[详细公式模板](file:///Users/mac/Desktop/lunwx/.agent/workflows/resources/formula-templates.md)
> 
> ⚠️ 所有公式必须与 [核心研究档案](file:///Users/mac/Desktop/lunwx/核心研究档案.txt) 保持一致

---

## 最常用公式速查

### 群时延（含碰撞，完整式）
```latex
$\tau_g(\omega) = \frac{d}{c} \frac{1}{\sqrt{1-\frac{\omega_p^2}{\omega^2(1+\delta)}}} \left[ 1 - \frac{\omega_p^2}{\omega^2} \frac{\delta}{(1+\delta)^2} \right]$
```

### 相对群时延（工程模型）
```latex
$\Delta\tau_g(f) = \frac{d}{c} \left( \frac{1}{\sqrt{1 - (f_p/f)^2}} - 1 \right)$
```

### 非线性度因子 η
```latex
$\eta(f) = \frac{B}{f} \cdot \frac{(f_p/f)^2}{\left[ 1 - (f_p/f)^2 \right]^{3/2}}$
```

### 工程判据
```latex
$B \cdot \eta \cdot \tau_0 \le 1$
```

### 时变时延模型
```latex
$\tau_g(t) = A_0 + A_1 t + A_2 t^2$
```

### 加权代价函数
```latex
$J(\mathbf{\theta}) = \sum_{i=1}^{N} w_i \cdot \left(Y_i - Y_{theory}(X_i; \mathbf{\theta})\right)^2$
```

### LM迭代公式
```latex
$\mathbf{\theta}_{k+1} = \mathbf{\theta}_{k} - \left(\mathbf{J}^T \mathbf{J} + \mu \mathbf{I}\right)^{-1} \mathbf{J}^T \mathbf{r}$
```

---

## 符号速查

| 符号 | 含义 | LaTeX |
|-----|-----|-------|
| $\omega_p$ | 特征角频率 | `$\omega_p$` |
| $f_p$ | 截止频率 | `$f_p$` |
| $n_e$ | 电子密度 | `$n_e$` |
| $\nu_e$ | 碰撞频率 | `$\nu_e$` |
| $\delta$ | 损耗因子 $(\nu_e/\omega)^2$ | `$\delta$` |
| $\eta$ | 非线性度 | `$\eta$` |
