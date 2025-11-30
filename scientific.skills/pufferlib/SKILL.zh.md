<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 河豚
描述：在处理强化学习任务时应使用此技能，包括高性能 RL 训练、自定义环境开发、矢量化并行模拟、多智能体系统或与现有 RL 环境（Gymnasium、PettingZoo、Atari、Procgen 等）集成。使用此技能来实施 PPO 训练、创建 PufferEnv 环境、优化 RL 性能或使用 CNN/LSTM 开发策略。
---

# PufferLib - 高性能强化学习

## 概述

PufferLib 是一个高性能强化学习库，专为快速并行环境模拟和训练而设计。它通过优化的矢量化、本机多代理支持和高效的 PPO 实现 (PuffeRL) 实现每秒数百万步的训练。该库提供了包含 20 多个环境的 Ocean 套件，并与 Gymnasium、PettingZoo 和专门的 RL 框架无缝集成。

## 何时使用此技能

在以下情况下使用此技能：
- **在任何环境（单或多智能体）上使用 PPO 训练 RL 智能体**
- **使用 PufferEnv API 创建自定义环境**
- **优化并行环境模拟的性能**（矢量化）
- **集成来自 Gymnasium、PettingZoo、Atari、Procgen 等的现有环境**。
- **使用 CNN、LSTM 或自定义架构制定策略**
- **将 RL** 扩展到每秒数百万步，以加快实验速度
- **多代理 RL** 具有本机多代理环境支持

## 核心能力

### 1.高性能训练（PuffeRL）

PuffeRL是PufferLib优化的PPO+LSTM训练算法，实现1M-4M步/秒。

**快速入门培训：**
```bash
# CLI training
puffer train procgen-coinrun --train.device cuda --train.learning-rate 3e-4

# Distributed training
torchrun --nproc_per_node=4 train.py
```

**Python 训练循环：**
<<<代码块_1>>>

**有关全面的培训指导**，请阅读 `references/training.md` 以了解：
- 完整的培训工作流程和 CLI 选项
- 使用蛋白质进行超参数调整
- 分布式多GPU/多节点训练
- 记录器集成（权重和偏差、Neptune）
- 检查点和恢复训练
- 性能优化技巧
- 课程学习模式

### 2.环境开发（PufferEnv）

使用 PufferEnv API 创建自定义高性能环境。

**基本环境结构：**
<<<代码块_2>>>

**使用模板脚本：** `scripts/env_template.py` 提供完整的单代理和多代理环境模板，示例如下：
- 不同的观察空间类型（矢量、图像、字典）
- 动作空间变化（离散、连续、多离散）
- 多智能体环境结构
- 测试实用程序

**有关完整的环境开发**，请阅读 `references/environments.md` 了解：
- PufferEnv API 详细信息和就地操作模式
- 观察和行动空间定义
- 多代理环境创建
- 海洋套件（20 多个预建环境）
- 性能优化（Python 到 C 工作流程）
- 环境包装和最佳实践
- 调试和验证技术

### 3.矢量化和性能

通过优化的并行仿真实现最大吞吐量。

**矢量化设置：**
<<<代码块_3>>>

**关键优化：**
- 用于零拷贝观察传递的共享内存缓冲区
- 忙等待标志而不是管道/队列
- 异步返回的剩余环境
- 每个工人有多个环境

**对于矢量化优化**，请阅读 `references/vectorization.md` 以了解：
- 架构和性能特征
- 工人和批量大小配置
- 串行模式、多处理模式、异步模式
- 共享内存和零复制模式
- 大规模分层矢量化
- 多智能体矢量化策略
- 性能分析和故障排除

### 4. 政策制定

使用可选实用程序将策略构建为标准 PyTorch 模块。

**基本政策结构：**
<<<代码块_4>>>

**如需完整的政策制定**，请阅读 `references/policies.md` 了解：
- CNN 图像观察策略
- 具有优化 LSTM 的循环策略（推理速度提高 3 倍）
- 复杂观测的多输入策略
- 持续行动政策
- 多代理策略（共享参数与独立参数）
- 高级架构（注意力、残差）
- 观察归一化和梯度裁剪
- 策略调试和测试

### 5.环境整合

无缝集成流行的 RL 框架的环境。

**体育馆整合：**
<<<代码块_5>>>

**PettingZoo 多代理：**
<<<代码块_6>>>

**支持的框架：**
- 健身房 / OpenAI 健身房
- PettingZoo（并行和AEC）
- 雅达利 (ALE)
- 普罗克根
- 网络黑客 / 迷你黑客
- 迷你电网
- 神经网络MMO
- 工匠
- GPU驱动
- 微RTS
- 网格状
- 还有更多...

**有关集成详细信息**，请阅读 `references/integration.md` 了解：
- 每个框架的完整集成示例
- 自定义包装（观察、奖励、框架堆叠、动作重复）
- 空间扁平化和反扁平化
- 环境登记
- 兼容模式
- 性能考虑
- 集成调试

## 快速启动工作流程

### 用于训练现有环境

1.从Ocean套件或兼容框架中选择环境
2.使用`scripts/train_template.py`作为起点
3. 为您的任务配置超参数
4. 使用 CLI 或 Python 脚本运行训练
5. 使用权重和偏差或 Neptune 进行监控
6. 参考`references/training.md`进行优化

### 用于创建自定义环境

1. 以 `scripts/env_template.py` 开头
2. 定义观察和行动空间
3. 实现`reset()`和`step()`方法
4.本地测试环境
5. 使用 `pufferlib.emulate()` 或 `make()` 进行矢量化
6. 高级模式请参阅`references/environments.md`
7. 如果需要，使用 `references/vectorization.md` 进行优化

### 用于政策制定

1. 根据观察选择架构：
   - 向量观测 → MLP 策略
   - 图像观察 → CNN 政策
   - 顺序任务 → LSTM 策略
   - 复杂的观察→多输入策略
2. 使用 `layer_init` 进行适当的权重初始化
3.遵循`references/policies.md`中的模式
4. 全面训练前进行环境测试

### 用于性能优化

1. 分析当前吞吐量（每秒步数）
2.检查矢量化配置（num_envs、num_workers）
3.优化环境代码（in-place ops、numpy向量化）
4. 考虑关键路径的 C 实现
5、使用`references/vectorization.md`进行系统优化

## 资源

### 脚本/

**train_template.py** - 完整的训练脚本模板：
- 环境创建和配置
- 策略初始化
- 记录器集成（WandB、Neptune）
- 带检查点的训练循环
- 命令行参数解析
- 多GPU分布式训练设置

**env_template.py** - 环境实现模板：
- 单代理 PufferEnv 示例（网格世界）
- 多代理PufferEnv示例（协作导航）
- 多种观察/行动空间模式
- 测试实用程序

###参考资料/

**training.md** - 综合培训指南：
- 培训工作流程和 CLI 选项
- 超参数配置
- 分布式训练（多GPU、多节点）
- 监控和记录
- 检查点
- 蛋白质超参数调整
- 性能优化
- 常见的训练模式
- 故障排除

**environments.md** - 环境开发指南：
- PufferEnv API 和特性
- 观察和行动空间
- 多代理环境
- 海洋套房环境
- 自定义环境开发工作流程
- Python到C的优化路径
- 第三方环境集成
- 包装和最佳实践
- 调试

**向量化.md** - 矢量化优化：
- 架构及关键优化
- 矢量化模式（串行、多处理、异步）
- 工作人员和批次配置
- 共享内存和零复制模式
- 高级矢量化（分层、自定义）
- 多智能体矢量化
- 性能监控和分析
- 故障排除和最佳实践

**policies.md** - 策略架构指南：
- 基本政策结构
- CNN 图像政策
- 具有优化的 LSTM 策略
- 多投入政策
- 持续行动政策
- 多代理策略
- 高级架构（注意力、残差）
- 观察处理和展平
- 初始化和标准化
- 调试和测试

**integration.md** - 框架集成指南：
- 体育馆一体化
- PettingZoo 集成（并行和 AEC）
- 第三方环境（Procgen、NetHack、Minigrid 等）
- 自定义包装（观察、奖励、框架堆叠等）
- 空间转换和平整
- 环境登记
- 兼容模式
- 性能考虑
- 调试集成

## 成功秘诀

1. **从简单开始**：在创建自定义环境之前从海洋环境或体育馆集成开始

2. **早期分析**：从一开始就测量每秒的步数以识别瓶颈
3. **使用模板**：`scripts/train_template.py` 和 `scripts/env_template.py` 提供可靠的起点

4. **根据需要阅读参考资料**：每个参考文件都是独立的，并且专注于特定的功能

5. **逐步优化**：从 Python 开始，分析，然后根据需要使用 C 优化关键路径

6. **利用矢量化**：PufferLib 的矢量化是实现高吞吐量的关键

7. **监控训练**：使用 WandB 或 Neptune 跟踪实验并及早发现问题

8. **测试环境**：在扩大训练规模之前验证环境逻辑

9. **检查现有环境**：Ocean suite提供20+预建环境

10. **使用正确的初始化**：始终使用 `pufferlib.pytorch` 中的 `layer_init` 作为策略

## 常见用例

### 标准基准培训
```python
# Atari
env = pufferlib.make('atari-pong', num_envs=256)

# Procgen
env = pufferlib.make('procgen-coinrun', num_envs=256)

# Minigrid
env = pufferlib.make('minigrid-empty-8x8', num_envs=256)
```

### 多代理学习
```python
# PettingZoo
env = pufferlib.make('pettingzoo-pistonball', num_envs=128)

# Shared policy for all agents
policy = create_policy(env.observation_space, env.action_space)
trainer = PuffeRL(env=env, policy=policy)
```

### 自定义任务开发
```python
# Create custom environment
class MyTask(PufferEnv):
    # ... implement environment ...

# Vectorize and train
env = pufferlib.emulate(MyTask, num_envs=256)
trainer = PuffeRL(env=env, policy=my_policy)
```

### 高性能优化
```python
# Maximize throughput
env = pufferlib.make(
    'my-env',
    num_envs=1024,      # Large batch
    num_workers=16,     # Many workers
    envs_per_worker=64  # Optimize per worker
)
```

## 安装

```bash
uv pip install pufferlib
```

## 文档

- 官方文档：https://puffer.ai/docs.html
- GitHub：https://github.com/PufferAI/PufferLib
- Discord：提供社区支持