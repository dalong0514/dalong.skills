<!-- 此文件由机器翻译自 esm-c-api.md -->

# ESM C API 参考

## 概述

ESM C (Cambrian) 是一系列蛋白质语言模型，针对表示学习和高效嵌入生成进行了优化。 ESM C 旨在作为 ESM2 的直接替代品，在所有模型尺寸上显着提高速度和质量。

## 模型架构

**ESM C 系列型号：**

|型号 ID |参数|层 |最适合 |
|----------|------------|--------|----------|
| `esmc-300m` | 300M | 30|快速推理，轻量级应用 |
| `esmc-600m` | 600M | 36 | 36平衡的性能和质量|
| `esmc-6b` | 6B| 80|最高的表现质量 |

**主要特点：**
- 推理速度比 ESM2 快 3 倍
- 提高了困惑度和嵌入质量
- 高效的生产部署架构
- 与 ESM2 工作流程兼容（直接替换）
- 支持长序列（高效可达 1024 个残基）

**相对于 ESM2 的架构改进：**
- 优化注意力机制
- 更好的代币表示
- 强化培训程序
- 减少内存占用

## 核心 API 组件

### ESMC 类

ESM C 型号的主界面。

**模型加载：**

```python
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein

# Load model with automatic device placement
model = ESMC.from_pretrained("esmc-300m").to("cuda")

# Or specify device explicitly
model = ESMC.from_pretrained("esmc-600m").to("cpu")

# For maximum quality
model = ESMC.from_pretrained("esmc-6b").to("cuda")
```

**型号选择标准：**

- **esmc-300m**：开发、实时应用、多个序列的批处理
- **esmc-600m**：生产部署，良好的质量/速度平衡
- **esmc-6b**：研究，下游任务的最大准确性

### 基本嵌入生成

**单一序列：**

<<<代码块_1>>>

**输出形状：**

对于长度为 L 的序列：
- `embeddings.shape`: `(1, L, hidden_dim)` 其中hidden_dim 取决于模型
  - esmc-300m：hidden_dim = 960
  - esmc-600m：hidden_dim = 1152
  - esmc-6b：hidden_dim = 2560
- `logits.shape`: `(1, L, 64)` - 每个位置的氨基酸预测

### 批处理

高效处理多个序列：

<<<代码块_2>>>

**可变长度的高效批处理：**

<<<代码块_3>>>

## 常见用例

### 1. 序列相似度分析

使用嵌入计算蛋白质之间的相似性：

<<<代码块_4>>>

### 2. 蛋白质分类

使用嵌入作为分类特征：

<<<代码块_5>>>

### 3. 蛋白质聚类

基于序列相似性的蛋白质聚类：

<<<代码块_6>>>

### 4. 序列搜索和检索

在数据库中查找相似序列：

```python
import torch
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

def build_sequence_index(model, database_sequences):
    """Build searchable index of sequence embeddings."""
    embeddings = []
    for seq in database_sequences:
        emb = get_sequence_embedding(model, seq)
        embeddings.append(emb.cpu().detach().numpy().flatten())
    return np.array(embeddings)

def search_similar_sequences(model, query_seq, database_embeddings,
                            database_sequences, top_k=10):
    """Find top-k most similar sequences."""
    query_emb = get_sequence_embedding(model, query_seq)
    query_emb_np = query_emb.cpu().detach().numpy().flatten().reshape(1, -1)

    # Compute similarities
    similarities = cosine_similarity(query_emb_np, database_embeddings)[0]

    # Get top-k
    top_indices = np.argsort(similarities)[-top_k:][::-1]

    results = [
        (database_sequences[idx], similarities[idx])
        for idx in top_indices
    ]
    return results

# Example usage
database_seqs = [...]  # Large sequence database
index = build_sequence_index(model, database_seqs)

query = "MPRTKEINDAGLIVHSP"
similar = search_similar_sequences(model, query, index, database_seqs, top_k=5)

for seq, score in similar:
    print(f"Score: {score:.4f} - {seq[:30]}...")
```

### 5.下游模型的特征提取

使用 ESM C 嵌入作为自定义神经网络的输入：

```python
import torch.nn as nn

class ProteinPropertyPredictor(nn.Module):
    """Example: Predict protein properties from ESM C embeddings."""

    def __init__(self, embedding_dim, hidden_dim, output_dim):
        super().__init__()
        self.fc1 = nn.Linear(embedding_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.fc3 = nn.Linear(hidden_dim, output_dim)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(0.3)

    def forward(self, embeddings):
        # embeddings: (batch, seq_len, embedding_dim)
        # Mean pool over sequence
        x = embeddings.mean(dim=1)

        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.fc3(x)
        return x

# Use ESM C as frozen feature extractor
esm_model = ESMC.from_pretrained("esmc-600m").to("cuda")
esm_model.eval()  # Freeze

# Create task-specific model
predictor = ProteinPropertyPredictor(
    embedding_dim=1152,  # esmc-600m dimension
    hidden_dim=512,
    output_dim=1  # e.g., stability score
).to("cuda")

# Training loop
for sequence, target in dataloader:
    protein = ESMProtein(sequence=sequence)
    with torch.no_grad():
        embeddings = esm_model.forward(esm_model.encode(protein))

    prediction = predictor(embeddings)
    loss = criterion(prediction, target)
    # ... backprop through predictor only
```

### 6. 每个残留物分析

提取每个残基的表示以进行详细分析：

```python
def get_per_residue_embeddings(model, sequence):
    """Get embedding for each residue."""
    protein = ESMProtein(sequence=sequence)
    tensor = model.encode(protein)
    embeddings = model.forward(tensor)

    # embeddings shape: (1, seq_len, hidden_dim)
    return embeddings.squeeze(0)  # (seq_len, hidden_dim)

# Analyze specific positions
sequence = "MPRTKEINDAGLIVHSPQWFYK"
residue_embeddings = get_per_residue_embeddings(model, sequence)

# Extract features for position 10
position_10_features = residue_embeddings[10]
print(f"Features for residue {sequence[10]} at position 10:")
print(f"Shape: {position_10_features.shape}")

# Compare residue representations
pos_5 = residue_embeddings[5]
pos_15 = residue_embeddings[15]
similarity = F.cosine_similarity(pos_5, pos_15, dim=0)
print(f"Residue similarity: {similarity.item():.4f}")
```

## 性能优化

### 内存管理

```python
import torch

# Use half precision for memory efficiency
model = ESMC.from_pretrained("esmc-600m").to("cuda").half()

# Process with mixed precision
with torch.cuda.amp.autocast():
    embeddings = model.forward(model.encode(protein))

# Clear cache between batches
torch.cuda.empty_cache()
```

### 批处理最佳实践

```python
def efficient_batch_processing(model, sequences, batch_size=32):
    """Process sequences in optimized batches."""
    results = []

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]

        # Process batch
        batch_embeddings = []
        for seq in batch:
            protein = ESMProtein(sequence=seq)
            emb = model.forward(model.encode(protein))
            batch_embeddings.append(emb)

        results.extend(batch_embeddings)

        # Periodically clear cache
        if i % (batch_size * 10) == 0:
            torch.cuda.empty_cache()

    return results
```

### 缓存嵌入

```python
import pickle
import hashlib

def get_cache_key(sequence):
    """Generate cache key for sequence."""
    return hashlib.md5(sequence.encode()).hexdigest()

class EmbeddingCache:
    """Cache for protein embeddings."""

    def __init__(self, cache_file="embeddings_cache.pkl"):
        self.cache_file = cache_file
        try:
            with open(cache_file, 'rb') as f:
                self.cache = pickle.load(f)
        except FileNotFoundError:
            self.cache = {}

    def get(self, sequence):
        key = get_cache_key(sequence)
        return self.cache.get(key)

    def set(self, sequence, embedding):
        key = get_cache_key(sequence)
        self.cache[key] = embedding

    def save(self):
        with open(self.cache_file, 'wb') as f:
            pickle.dump(self.cache, f)

# Usage
cache = EmbeddingCache()

def get_embedding_cached(model, sequence):
    cached = cache.get(sequence)
    if cached is not None:
        return cached

    # Compute
    protein = ESMProtein(sequence=sequence)
    embedding = model.forward(model.encode(protein))
    cache.set(sequence, embedding)

    return embedding

# Don't forget to save cache
cache.save()
```

## 与 ESM2 的比较

**性能改进：**

|公制| ESM2-650M | ESM C-600M |改进|
|--------|---------|------------|------------|
|推理速度 | 1.0 倍 | 3.0 倍 |快 3 倍 |
|困惑|更高 |降低|更好 |
|内存使用情况| 1.0 倍 | 0.8 倍 |减少 20% |
|嵌入质量|基线|改进| +5-10% |

**从 ESM2 迁移：**

ESM C 被设计为直接替代品：

```python
# Old ESM2 code
from esm import pretrained
model, alphabet = pretrained.esm2_t33_650M_UR50D()

# New ESM C code (similar API)
from esm.models.esmc import ESMC
model = ESMC.from_pretrained("esmc-600m")
```

主要区别：
- 推理速度更快，质量相同或更好
- 通过 ESMProtein 简化 API
- 更好地支持长序列
- 更有效的内存使用

## 高级主题

### 微调 ESM C

ESM C 可以针对特定任务进行微调：

```python
import torch.optim as optim

# Load model
model = ESMC.from_pretrained("esmc-300m").to("cuda")

# Unfreeze for fine-tuning
for param in model.parameters():
    param.requires_grad = True

# Define optimizer
optimizer = optim.Adam(model.parameters(), lr=1e-5)

# Training loop
for epoch in range(num_epochs):
    for sequences, labels in dataloader:
        optimizer.zero_grad()

        # Forward pass
        proteins = [ESMProtein(sequence=seq) for seq in sequences]
        embeddings = [model.forward(model.encode(p)) for p in proteins]

        # Your task-specific loss
        loss = compute_loss(embeddings, labels)

        loss.backward()
        optimizer.step()
```

### 注意力可视化

提取注意力权重以实现可解释性：

```python
def get_attention_weights(model, sequence):
    """Extract attention weights from model."""
    protein = ESMProtein(sequence=sequence)
    tensor = model.encode(protein)

    # Forward with attention output
    output = model.forward(tensor, output_attentions=True)

    return output.attentions  # List of attention tensors per layer

# Visualize attention
attentions = get_attention_weights(model, "MPRTKEINDAGLIVHSP")
# Process and visualize attention patterns
```

## 引文

如果在研究中使用 ESM C，请引用：

```
ESM Cambrian: https://www.evolutionaryscale.ai/blog/esm-cambrian
EvolutionaryScale (2024)
```

## 其他资源

- ESM C 博客文章：https://www.evolutionaryscale.ai/blog/esm-cambrian
- 模型权重：HuggingFace EvolutionaryScale 组织
- 比较基准：请参阅博客文章了解详细的性能比较