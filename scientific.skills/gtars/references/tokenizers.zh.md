<!-- 此文件由机器翻译自 tokenizers.md -->

# 基因组分词器

分词器将基因组区域转换为机器学习应用程序的离散标记，对于训练基因组深度学习模型特别有用。

## Python API

### 创建分词器

从各种来源加载分词器配置：

```python
import gtars

# From BED file
tokenizer = gtars.tokenizers.TreeTokenizer.from_bed_file("regions.bed")

# From configuration file
tokenizer = gtars.tokenizers.TreeTokenizer.from_config("tokenizer_config.yaml")

# From region string
tokenizer = gtars.tokenizers.TreeTokenizer.from_region_string("chr1:1000-2000")
```

### 基因组区域标记化

将基因组坐标转换为标记：

<<<代码块_1>>>

### 令牌属性

访问令牌信息：

<<<代码块_2>>>

## 用例

### 机器学习预处理

分词器对于为 ML 模型准备基因组数据至关重要：

1. **序列建模**：将基因组间隔转换为变压器模型的离散标记
2. **位置编码**：跨数据集创建一致的位置编码
3. **数据增强**：生成用于训练的替代标记化

### 与 geniml 集成

tokenizers 模块与基因组 ML 的 geniml 库无缝集成：

<<<代码块_3>>>

## 配置格式

Tokenizer配置文件支持YAML格式：

<<<代码块_4>>>

## 性能考虑因素

- TreeTokenizer 使用高效的数据结构进行快速标记化
- 建议对大型数据集进行批量标记化
- 预加载分词器可减少重复操作的开销