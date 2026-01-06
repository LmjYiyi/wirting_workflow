---
description: 论文核心公式的LaTeX格式模板
---

# 论文公式 LaTeX 模板

> ⚠️ **重要**：所有公式必须与 [核心研究档案](file:///Users/mac/Desktop/lunwx/核心研究档案.txt) 保持一致

---

## 第三章：色散信道物理建模

### 3.1.1 复介电常数（Drude模型）

**完整表达式**：
```latex
$\tilde{\varepsilon}_r(\omega) = \varepsilon_{r'}(\omega) - j\varepsilon_{r''}(\omega) = \left( 1 - \frac{\omega_p^2}{\omega^2 + v_e^2} \right) - j \left( \frac{v_e}{\omega} \frac{\omega_p^2}{\omega^2 + v_e^2} \right)$
```

**等离子体特征频率**：
```latex
$\omega_p = \sqrt{\frac{n_e e^2}{\varepsilon_0 m_e}}$
```

### 3.1.1 复波数与传播常数

```latex
$\tilde{k}(\omega) = \frac{\omega}{c}\sqrt{\tilde{\varepsilon}_r(\omega)} = \beta(\omega) - j\alpha(\omega)$
```

**相位常数近似（高频弱碰撞）**：
```latex
$\beta(\omega) \approx \frac{\omega}{c}\sqrt{1-\frac{\omega_p^2}{\omega^2+v_e^2}}$
```

### 3.1.1 群时延推导

**定义**：
```latex
$\tau_g(\omega) = d \cdot \frac{d\beta(\omega)}{d\omega}$
```

**β对ω求导（第一步）**：
```latex
$\frac{d\beta}{d\omega} = \frac{1}{c} \frac{d}{d\omega} \left[ \omega \left( 1 - \frac{\omega_p^2}{\omega^2+v_e^2} \right)^{1/2} \right]$
```

**β对ω求导（展开）**：
```latex
$\frac{d\beta}{d\omega} = \frac{1}{c} \left[ \sqrt{1-\frac{\omega_p^2}{\omega^2+v_e^2}} + \omega \cdot \frac{1}{2} \left( 1 - \frac{\omega_p^2}{\omega^2+v_e^2} \right)^{-1/2} \cdot \frac{2\omega\omega_p^2}{(\omega^2+v_e^2)^2} \right]$
```

**群时延完整表达式（未简化）**：
```latex
$\tau_g(\omega, \omega_p, v_e) = \frac{d}{c} \cdot \frac{1}{\sqrt{1-\frac{\omega_p^2}{\omega^2+v_e^2}}} \cdot \left[ 1 - \frac{\omega_p^2 v_e^2}{(\omega^2+v_e^2)^2} \right]$
```

**δ定义**：
```latex
$\delta = \left( \frac{v_e}{\omega} \right)^2$
```

**引入δ后的群时延**：
```latex
$\tau_g(\omega) = \frac{d}{c} \frac{1}{\sqrt{1-\frac{\omega_p^2}{\omega^2(1+\delta)}}} \left[ 1 - \frac{\omega_p^2}{\omega^2} \frac{\delta}{(1+\delta)^2} \right]$
```

**最终近似公式（碰撞为二阶微扰）**：
```latex
$\tau_g(\omega) \approx \frac{d}{c\sqrt{1-(\omega_p / \omega)^2}} \left[ 1 - \left( \frac{1}{2} \frac{(\omega_p / \omega)^2}{1-(\omega_p / \omega)^2} + \frac{\omega_p^2}{\omega^2} \right) \delta \right]$
```

---

### 3.1.2 相对群时延观测模型

**定义**：
```latex
$\Delta\tau_g(\omega) = \tau_g(\omega) - \tau_0 = \tau_g(\omega) - \frac{d}{c}$
```

**完整映射模型（含碰撞微扰）**：
```latex
$\Delta\tau_g(f) = \mathcal{M}_{full}(f; f_p, d, \nu_e) \approx \frac{d}{c} \left\{ \frac{1}{\sqrt{1-(f_p/f)^2}} \left[ 1 - \Psi(f) \cdot \delta(f) \right] - 1 \right\}$
```

**二阶修正系数**：
```latex
$\Psi(f) = \frac{1}{2} \frac{(f_p/f)^2}{1-(f_p/f)^2} + \left(\frac{f_p}{f}\right)^2$
```

**工程主导模型（简化）**：
```latex
$\Delta\tau_g(f) = \mathcal{M}(f; f_p, d) = \frac{d}{c} \left( \frac{1}{\sqrt{1 - (f_p/f)^2}} - 1 \right), \quad (f > f_p)$
```

---

### 3.1.3 非线性度因子 η

**定义**：
```latex
$\eta(f) \triangleq \frac{1}{\tau_0} \cdot \left| \frac{d\tau_g(f)}{df} \right| \cdot B$
```

**GDD（群时延色散率）**：
```latex
$\frac{d\tau_g(f)}{df} = -\tau_0 \cdot \frac{1}{f} \cdot \frac{(f_p/f)^2}{\left[ 1 - (f_p/f)^2 \right]^{3/2}}$
```

**η显式表达式**：
```latex
$\eta(f) = \frac{B}{f} \cdot \frac{(f_p/f)^2}{\left[ 1 - (f_p/f)^2 \right]^{3/2}}, \quad (f > f_p)$
```

---

### 3.3 差频信号误差解析

**LFMCW发射信号**：
```latex
$s_{TX}(t) = \exp\left\{ j \left( \omega_0 t + \pi \frac{B}{T_m} t^2 \right) \right\}, \quad 0 \le t \le T_m$
```

**瞬时角频率**：
```latex
$\omega(t) = \frac{d\Phi_{TX}}{dt} = \omega_0 + \mu t, \quad \text{其中 } \mu = \frac{2\pi B}{T_m}$
```

**群时延二阶泰勒展开**：
```latex
$\tau_g(\omega) \approx \tau_g(\omega_0) + \frac{d\tau_g}{d\omega}\bigg|_{\omega_0} (\omega - \omega_0) + \frac{1}{2} \frac{d^2\tau_g}{d\omega^2}\bigg|_{\omega_0} (\omega - \omega_0)^2$
```

**时变时延模型**：
```latex
$\tau_g(t) = A_0 + A_1 t + A_2 t^2$
```

**展开系数 A₀（常数时延）**：
```latex
$A_0 = \tau_g(\omega_0) = \frac{d}{c} \left[ 1 - \left(\frac{\omega_p}{\omega_0}\right)^2 \right]^{-\frac{1}{2}}$
```

**展开系数 A₁（线性时变）**：
```latex
$A_1 = \mu \cdot \tau_1(\omega_0) = \left( \frac{2\pi B}{T_m} \right) \cdot \left\{ -\frac{d}{c} \frac{\omega_p^2}{\omega_0^3} \left[ 1 - \left(\frac{\omega_p}{\omega_0}\right)^2 \right]^{-\frac{3}{2}} \right\}$
```

**展开系数 A₂（二次非线性）**：
```latex
$A_2 = \frac{1}{2} \mu^2 \cdot \tau_2(\omega_0) = \frac{1}{2} \left( \frac{2\pi B}{T_m} \right)^2 \cdot \left\{ \frac{3d}{c} \frac{\omega_p^2}{\omega_0^4} \left[ 1 - \frac{\omega_p^2}{3\omega_0^2} \right] \left[ 1 - \left(\frac{\omega_p}{\omega_0}\right)^2 \right]^{-\frac{5}{2}} \right\}$
```

---

### 3.4.2 工程判据

```latex
$B \cdot \eta \cdot \tau_0 \le 1$
```

---

## 第四章：反演算法

### 4.1 滑动窗口探测频率

```latex
$X_i = f_{probe}(t_i) = f_{start} + K \cdot t_i$
```

### 4.1.2 MDL准则

```latex
$\text{MDL}(k) = -\ln\left( \frac{\prod_{j=k+1}^{M} \lambda_j^{\frac{1}{M-k}}}{\frac{1}{M-k} \sum_{j=k+1}^{M} \lambda_j} \right)^{(M-k)L} + \frac{1}{2} k (2M-k) \ln L$
```

### 4.1.3 测量时延

```latex
$Y_i = \tau_{meas}(t_i) = \frac{f_{beat,i}}{K}$
```

### 4.2.1 理论观测模型（含误差修正）

```latex
$Y_{theory}(f; \omega_p, d) = \tau_{group}(f) + \underbrace{\frac{f}{K} \cdot \frac{\partial \tau_{group}(f)}{\partial f}}_{\text{一阶色散修正项}}$
```

**Drude群时延**：
```latex
$\tau_{group}(f) = \frac{d}{c} \frac{1}{\sqrt{1 - (\omega_p/2\pi f)^2}}$
```

### 4.2.2 加权代价函数

```latex
$J(\mathbf{\theta}) = \sum_{i=1}^{N} w_i \cdot \left(Y_i - Y_{theory}(X_i; \mathbf{\theta})\right)^2$
```

**权重因子**：
```latex
$w_i = \left(\frac{A_i}{\max(\mathbf{A})}\right)^2$
```

### 4.2.3 LM迭代公式

**残差向量**：
```latex
$r_i = \sqrt{w_i}(Y_i - Y_{theory}(X_i, \mathbf{\theta}))$
```

**雅可比矩阵**：
```latex
$J_{ij} = \frac{\partial r_i}{\partial \theta_j}$
```

**迭代更新**：
```latex
$\mathbf{\theta}_{k+1} = \mathbf{\theta}_{k} - \left(\mathbf{J}^T \mathbf{J} + \mu \mathbf{I}\right)^{-1} \mathbf{J}^T \mathbf{r}$
```

### 电子密度计算

```latex
$n_e = \frac{\varepsilon_0 m_e}{e^2} \hat{\omega}_p^2$
```

---

## 二阶微扰证明公式

**介电常数实部展开**：
```latex
$\varepsilon'(\omega) \approx \underbrace{\left(1 - \frac{\omega_p^2}{\omega^2}\right)}_{\text{主导项(无碰撞)}} + \underbrace{\frac{\omega_p^2}{\omega^2}\left(\frac{\nu_e}{\omega}\right)^2}_{\text{碰撞修正项}}$
```

**介电常数虚部（一阶）**：
```latex
$\varepsilon''(\omega) \approx \frac{\omega_p^2}{\omega^2} \left( \frac{\nu_e}{\omega} \right)$
```

**条件数（病态性）**：
```latex
$\text{cond}(\mathbf{H}) = \frac{\lambda_{max}}{\lambda_{min}} \to \infty$
```
