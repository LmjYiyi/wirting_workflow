---
description: MATLAB仿真代码与论文章节的映射参考
---

# MATLAB 代码参考

> **代码位置**：`/Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/`
>
> 详细文档请查看：[完整代码参考指南](file:///Users/mac/Desktop/lunwx/.agent/workflows/resources/matlab-code-reference.md)

---

## 快速映射

| 代码 | 章节 | 运行命令 |
|-----|-----|---------|
| [`nue.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/nue.m) | 3.2.3 参数敏感性 | `>> nue` |
| [`lianghua.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/lianghua.m) | 3.1.3 & 3.4.2 工程判据 | `>> lianghua` |
| [`test.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/test.m) | 3.3.3 二阶微扰 | `>> test` |
| [`LM.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/LM.m) | 第四章 反演算法 | `>> LM` |
| [`initial.m`](file:///Users/mac/Desktop/lunwx/.agent/workflows/thesis-code/initial.m) | 传统方法对比 | `>> initial` |

---

## 代码目录结构

```
.agent/workflows/thesis-code/
├── LM.m         # 核心反演算法 (第四章)
├── initial.m    # 传统诊断方法 (对比基线)
├── nue.m        # 参数敏感性分析 (3.2.3节)
├── lianghua.m   # 工程判据可视化 (3.1.3/3.4.2节)
├── test.m       # 二阶小量验证 (3.3.3节)
└── README.md    # 代码文档说明
```
