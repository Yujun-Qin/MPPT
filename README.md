# MPPT-Photovoltaic-EnergyStorage-Research

## 📖 项目概述

本项目深入研究光伏储能微电网系统的控制策略，主要创新点包括：

- **提出了IPSO-SMC混合MPPT算法**：结合改进粒子群算法(IPSO)与滑模控制(SMC)，在局部遮阴和动态光照条件下实现99.6%以上的追踪效率
- **开发了完整的系统控制框架**：涵盖光伏发电、锂电池储能、直流母线电压稳定的全系统控制
- **验证了理论与实践的可行性**：通过MATLAB/Simulink仿真和STM32硬件平台验证算法性能

## ✨ 核心特性

- **全面的MPPT算法库**：包含传统算法(P&O、IncCond)和智能算法(PSO、IPSO、模糊控制、神经网络)的完整实现
- **多场景仿真模型**：提供均匀光照、局部遮阴、动态环境等多种场景的仿真模型
- **硬件实践指南**：完整的电路设计、元器件选型和嵌入式代码
- **渐进式学习路径**：从理论基础到高级应用的系统学习材料
- **完整实验数据**：提供仿真和实际测试数据用于算法性能分析

## 🗂 仓库结构

```
MPPT-Photovoltaic-EnergyStorage-Research/
│
├── 01-Theory-Basics/                 # 光伏电池、储能系统、DC/DC变换器理论基础
├── 02-Algorithms-Implementation/     # MPPT算法实现(传统算法+智能算法)
├── 03-Simulation-Models/             # MATLAB/Simulink和Python仿真模型
├── 04-Experimental-Data/            # 仿真和实验结果数据
└── 05-Hardware-Design/              # 硬件设计文档和原理图
```

## 🚀 快速开始

### 环境要求

- MATLAB R2020a或更高版本
- Simulink仿真环境
- STM32CubeIDE (用于嵌入式开发)
- Python 3.8+ (可选，用于辅助分析)

### 基本用法

1. **克隆仓库**
```bash
git clone https://github.com/KindKeeper/MPPT.git
```

2. **运行仿真示例**
```matlab
% 在MATLAB中打开并运行基本仿真
open_system('03-Simulation-Models/MATLAB-Simulink/Uniform_Illumination/PV_MPPT_Uniform.slx');
sim('PV_MPPT_Uniform');
```

3. **算法测试**
```matlab
% 测试IPSO-SMC算法性能
cd('02-Algorithms-Implementation/Intelligent-MPPT/IPSO-SMC_MPPT/');
run_IPSO_SMC_simulation();
```

## 📊 性能结果

### 算法对比（仿真环境）
| 算法                 | 追踪效率  | 响应时间  | 震荡幅度  | 局部遮阴处理 |
| -------------------- | --------- | --------- | --------- | ------------ |
| P&O                  | 95.5%     | 0.35s     | ±1.5V     | 差           |
| 标准PSO              | 97.8%     | 0.25s     | ±1.2V     | 良           |
| **IPSO-SMC(本论文)** | **99.6%** | **0.15s** | **±1.0V** | **优**       |

### 典型结果展示

**均匀光照条件下MPPT追踪效果**：
- 电压稳定时间：<0.2s
- 功率波动：<2%
- 直流母线电压稳定性：±0.5V

**局部遮阴条件下性能**：
- 全局最大功率点识别准确率：100%
- 避免局部最优解：有效
- 动态响应能力：优秀


## 🔧 硬件实现

本项目提供完整的硬件设计文档：

- **主控制器**：STM32F407ZGT6
- **光伏模拟器**：100W光伏板+1000W碘钨灯光源
- **储能系统**：18650锂电池组
- **采样电路**：基于NSi1311运算放大器的电压采样
- **驱动电路**：IGP06N60T IGBT驱动模块
- **通信接口**：CAN总线通信模块

详细连接方式和配置请参考 `05-Hardware-Design/` 中的文档。

