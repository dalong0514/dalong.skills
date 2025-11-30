<!-- 此文件由机器翻译自 policies.md -->

# PufferLib 政策指南

## 概述

PufferLib 策略是标准 PyTorch 模块，具有用于观察处理和 LSTM 集成的可选实用程序。该框架提供默认架构和工具，同时允许策略设计具有充分的灵活性。

## 政策架构

### 基本政策结构

```python
import torch
import torch.nn as nn
from pufferlib.pytorch import layer_init

class BasicPolicy(nn.Module):
    def __init__(self, observation_space, action_space):
        super().__init__()

        self.observation_space = observation_space
        self.action_space = action_space

        # Encoder network
        self.encoder = nn.Sequential(
            layer_init(nn.Linear(observation_space.shape[0], 256)),
            nn.ReLU(),
            layer_init(nn.Linear(256, 256)),
            nn.ReLU()
        )

        # Policy head (actor)
        self.actor = layer_init(nn.Linear(256, action_space.n), std=0.01)

        # Value head (critic)
        self.critic = layer_init(nn.Linear(256, 1), std=1.0)

    def forward(self, observations):
        """Forward pass through policy."""
        # Encode observations
        features = self.encoder(observations)

        # Get action logits and value
        logits = self.actor(features)
        value = self.critic(features)

        return logits, value

    def get_action(self, observations, deterministic=False):
        """Sample action from policy."""
        logits, value = self.forward(observations)

        if deterministic:
            action = logits.argmax(dim=-1)
        else:
            dist = torch.distributions.Categorical(logits=logits)
            action = dist.sample()

        return action, value
```

### 层初始化

PufferLib 提供了 `layer_init` 来进行正确的权重初始化：

<<<代码块_1>>>

## CNN 政策

对于基于图像的观察：

<<<代码块_2>>>

### 高效的 CNN 架构

<<<代码块_3>>>

## 循环策略 (LSTM)

PufferLib 提供了优化的 LSTM 集成与自动递归处理：

<<<代码块_4>>>

### LSTM 优化

PufferLib 的 LSTM 优化在部署期间使用 LSTMCell，在训练期间使用 LSTM，以将推理速度提高 3 倍：

<<<代码块_5>>>

## 多输入策略

对于具有多种观察类型的环境：

<<<代码块_6>>>

## 持续行动政策

对于连续控制任务：

```python
class ContinuousPolicy(nn.Module):
    def __init__(self, observation_space, action_space):
        super().__init__()

        self.encoder = nn.Sequential(
            layer_init(nn.Linear(observation_space.shape[0], 256)),
            nn.ReLU(),
            layer_init(nn.Linear(256, 256)),
            nn.ReLU()
        )

        # Mean of action distribution
        self.actor_mean = layer_init(nn.Linear(256, action_space.shape[0]), std=0.01)

        # Log std of action distribution
        self.actor_logstd = nn.Parameter(torch.zeros(1, action_space.shape[0]))

        # Value head
        self.critic = layer_init(nn.Linear(256, 1), std=1.0)

    def forward(self, observations):
        features = self.encoder(observations)

        action_mean = self.actor_mean(features)
        action_std = torch.exp(self.actor_logstd)

        value = self.critic(features)

        return action_mean, action_std, value

    def get_action(self, observations, deterministic=False):
        action_mean, action_std, value = self.forward(observations)

        if deterministic:
            return action_mean, value
        else:
            dist = torch.distributions.Normal(action_mean, action_std)
            action = dist.sample()
            return torch.tanh(action), value  # Bound actions to [-1, 1]
```

## 观察处理

PufferLib 提供了用于展开观察的实用程序：

```python
from pufferlib.pytorch import unflatten_observations

class PolicyWithUnflatten(nn.Module):
    def __init__(self, observation_space, action_space):
        super().__init__()

        self.observation_space = observation_space

        # Define encoders for each observation component
        self.encoders = nn.ModuleDict({
            'image': self._make_image_encoder(),
            'vector': self._make_vector_encoder()
        })

        # ... rest of policy ...

    def forward(self, flat_observations):
        # Unflatten observations into structured format
        observations = unflatten_observations(
            flat_observations,
            self.observation_space
        )

        # Process each component
        image_features = self.encoders['image'](observations['image'])
        vector_features = self.encoders['vector'](observations['vector'])

        # Combine and continue...
```

## 多代理策略

### 共享参数

所有代理都使用相同的策略：

```python
class SharedMultiAgentPolicy(nn.Module):
    def __init__(self, observation_space, action_space, num_agents):
        super().__init__()

        self.num_agents = num_agents

        # Single policy shared across all agents
        self.encoder = nn.Sequential(
            layer_init(nn.Linear(observation_space.shape[0], 256)),
            nn.ReLU()
        )

        self.actor = layer_init(nn.Linear(256, action_space.n), std=0.01)
        self.critic = layer_init(nn.Linear(256, 1), std=1.0)

    def forward(self, observations):
        """
        Args:
            observations: (batch * num_agents, obs_dim)
        Returns:
            logits: (batch * num_agents, num_actions)
            values: (batch * num_agents, 1)
        """
        features = self.encoder(observations)
        return self.actor(features), self.critic(features)
```

### 独立参数

每个代理都有自己的政策：

```python
class IndependentMultiAgentPolicy(nn.Module):
    def __init__(self, observation_space, action_space, num_agents):
        super().__init__()

        self.num_agents = num_agents

        # Separate policy for each agent
        self.policies = nn.ModuleList([
            self._make_policy(observation_space, action_space)
            for _ in range(num_agents)
        ])

    def _make_policy(self, observation_space, action_space):
        return nn.Sequential(
            layer_init(nn.Linear(observation_space.shape[0], 256)),
            nn.ReLU(),
            layer_init(nn.Linear(256, 256)),
            nn.ReLU()
        )

    def forward(self, observations, agent_ids):
        """
        Args:
            observations: (batch, obs_dim)
            agent_ids: (batch,) which agent each obs belongs to
        """
        outputs = []
        for agent_id in range(self.num_agents):
            mask = agent_ids == agent_id
            if mask.any():
                agent_obs = observations[mask]
                agent_out = self.policies[agent_id](agent_obs)
                outputs.append(agent_out)

        return torch.cat(outputs, dim=0)
```

## 高级架构

### 基于注意力的政策

```python
class AttentionPolicy(nn.Module):
    def __init__(self, observation_space, action_space, d_model=256, nhead=8):
        super().__init__()

        self.encoder = layer_init(nn.Linear(observation_space.shape[0], d_model))

        self.attention = nn.MultiheadAttention(d_model, nhead, batch_first=True)

        self.actor = layer_init(nn.Linear(d_model, action_space.n), std=0.01)
        self.critic = layer_init(nn.Linear(d_model, 1), std=1.0)

    def forward(self, observations):
        # Encode
        features = self.encoder(observations)

        # Self-attention
        features = features.unsqueeze(1)  # Add sequence dimension
        attn_out, _ = self.attention(features, features, features)
        attn_out = attn_out.squeeze(1)

        return self.actor(attn_out), self.critic(attn_out)
```

### 剩余政策

```python
class ResidualBlock(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.block = nn.Sequential(
            layer_init(nn.Linear(dim, dim)),
            nn.ReLU(),
            layer_init(nn.Linear(dim, dim))
        )

    def forward(self, x):
        return x + self.block(x)

class ResidualPolicy(nn.Module):
    def __init__(self, observation_space, action_space, num_blocks=4):
        super().__init__()

        dim = 256

        self.encoder = layer_init(nn.Linear(observation_space.shape[0], dim))

        self.blocks = nn.Sequential(
            *[ResidualBlock(dim) for _ in range(num_blocks)]
        )

        self.actor = layer_init(nn.Linear(dim, action_space.n), std=0.01)
        self.critic = layer_init(nn.Linear(dim, 1), std=1.0)

    def forward(self, observations):
        x = torch.relu(self.encoder(observations))
        x = self.blocks(x)
        return self.actor(x), self.critic(x)
```

## 政策最佳实践

### 初始化

```python
# Always use layer_init for proper initialization
good_layer = layer_init(nn.Linear(256, 256))

# Use small std for actor head (more stable early training)
actor = layer_init(nn.Linear(256, num_actions), std=0.01)

# Use std=1.0 for critic head
critic = layer_init(nn.Linear(256, 1), std=1.0)
```

### 观察标准化

```python
class NormalizedPolicy(nn.Module):
    def __init__(self, observation_space, action_space):
        super().__init__()

        # Running statistics for normalization
        self.obs_mean = nn.Parameter(torch.zeros(observation_space.shape[0]), requires_grad=False)
        self.obs_std = nn.Parameter(torch.ones(observation_space.shape[0]), requires_grad=False)

        # ... rest of policy ...

    def forward(self, observations):
        # Normalize observations
        normalized_obs = (observations - self.obs_mean) / (self.obs_std + 1e-8)

        # Continue with normalized observations
        return self.policy(normalized_obs)

    def update_normalization(self, observations):
        """Update running statistics."""
        self.obs_mean.data = observations.mean(dim=0)
        self.obs_std.data = observations.std(dim=0)
```

### 渐变裁剪

```python
# PufferLib trainer handles gradient clipping automatically
trainer = PuffeRL(
    env=env,
    policy=policy,
    max_grad_norm=0.5  # Clip gradients to this norm
)
```

### 模型编译

```python
# Enable torch.compile for faster training (PyTorch 2.0+)
policy = MyPolicy(observation_space, action_space)

# Compile the model
policy = torch.compile(policy, mode='reduce-overhead')

# Use with trainer
trainer = PuffeRL(env=env, policy=policy, compile=True)
```

## 调试策略

### 检查输出形状

```python
def test_policy_shapes(policy, observation_space, batch_size=32):
    """Verify policy output shapes."""
    # Create dummy observations
    obs = torch.randn(batch_size, *observation_space.shape)

    # Forward pass
    logits, value = policy(obs)

    # Check shapes
    assert logits.shape == (batch_size, policy.action_space.n)
    assert value.shape == (batch_size, 1)

    print("✓ Policy shapes correct")
```

### 验证梯度

```python
def check_gradients(policy, observation_space):
    """Check that gradients flow properly."""
    obs = torch.randn(1, *observation_space.shape, requires_grad=True)

    logits, value = policy(obs)

    # Backward pass
    loss = logits.sum() + value.sum()
    loss.backward()

    # Check gradients exist
    for name, param in policy.named_parameters():
        if param.grad is None:
            print(f"⚠ No gradient for {name}")
        elif torch.isnan(param.grad).any():
            print(f"⚠ NaN gradient for {name}")
        else:
            print(f"✓ Gradient OK for {name}")
```