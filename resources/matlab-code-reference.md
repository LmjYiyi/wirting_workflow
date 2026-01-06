---
description: MATLAB仿真代码与论文章节的映射参考
---

# MATLAB 仿真代码参考指南

> 代码位置：`.agent/workflows/thesis-code/`

---

## 代码目录

```
/Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/
├── LM.m         # 核心反演算法 (第四章)
├── initial.m    # 传统诊断方法 (对比基线)
├── nue.m        # 参数敏感性分析 (3.2.3节)
├── lianghua.m   # 工程判据可视化 (3.1.3/3.4.2节)
├── test.m       # 二阶小量验证 (3.3.3节)
└── README.md    # 代码文档说明
```

---

## 代码与论文章节映射

| 代码文件 | 对应论文章节 | 核心功能 | 生成图表 |
|----------|-------------|----------|---------|
| [`nue.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/nue.m) | 3.2.3 参数敏感性分析 | 验证电子密度主导性 + 碰撞频率解耦 | 图 3-3a, 3-3b |
| [`lianghua.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/lianghua.m) | 3.1.3 & 3.4.2 | 工程判据 $\xi = B \cdot \eta \cdot \tau_0$ | 工程边界图 |
| [`test.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/test.m) | 3.3.3 二阶微扰证明 | 验证 $(\nu_e/\omega)^2$ 是二阶小量 | 敏感度验证图 |
| [`initial.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/initial.m) | 2.2 & 5.2 | 传统LFMCW诊断方法仿真 | Figure 1-10 |
| [`LM.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/LM.m) | **第四章全部** | 完整反演算法实现 | Figure 9, 11 |

---

## 快速运行命令

```matlab
% === 第三章图表 ===
>> nue          % 生成图 3-3a, 3-3b (参数敏感性)
>> lianghua     % 生成工程判据图 (3.4.2节)
>> test         % 生成二阶小量验证图 (3.3.3节)

% === 第四章图表 ===
>> LM           % 生成 Figure 9, 11 (反演算法验证)

% === 对比分析 ===
>> initial      % 生成传统方法结果 (作为对比基线)
```

---

## 代码与公式对应

| 论文公式 | 代码位置 | 实现方式 |
|---------|---------|---------|
| (3-1) 复介电常数 | `LM.m` 第3节 | `eps_complex = (1 - wp^2/(omega^2+nu^2)) - 1j*...` |
| (3-7) 群时延解析式 | `nue.m` 函数 | `tau_g = -diff(phi)./diff(omega)` |
| (3-18) 非线性度因子 η | `lianghua.m` | `eta = (B./f) .* (x.^2) ./ ((1-x.^2).^1.5)` |
| 工程判据 ξ | `lianghua.m` | `Xi = B .* eta .* tau0` |
| MDL 准则 | `LM.m` 第7.3节 | `mdl_cost(k+1) = term1 + term2` |
| 加权代价函数 | `LM.m` 第9节 | `Weights = (W_raw / max(W_raw)).^2` |

---

## 论文引用模板

### 描述仿真设置
> 为了验证上述理论推导，本文基于 MATLAB 搭建了数值仿真平台。仿真参数设置如下：等离子体厚度 $d = 150$ mm，截止频率 $f_p = 29$ GHz（对应电子密度 $n_e = 1.04 \times 10^{19}$ m$^{-3}$），碰撞频率 $\nu_e = 1.5$ GHz。

### 描述图表观察
> 从仿真结果（图 X-X）可以观测到，当探测频率逼近截止频率时，群时延呈现非线性的急剧陡升趋势，验证了理论模型中 $\lim_{f \to f_p^+} \tau_g(f) = \infty$ 的预测。

### 描述算法性能
> 为验证所提算法的有效性，进行了 Monte Carlo 仿真。在信噪比为 20 dB 的条件下，电子密度反演误差的均方根值（RMSE）为 X.XX%。
