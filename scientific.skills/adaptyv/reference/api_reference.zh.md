<!-- 此文件由机器翻译自 api_reference.md -->

# Adaptyv API 参考

## 基本网址

```
https://kq5jp7qj7wdqklhsxmovkzn4l40obksv.lambda-url.eu-central-1.on.aws
```

## 身份验证

所有 API 请求都需要在请求标头中进行不记名令牌身份验证：

<<<代码块_1>>>

要获取 API 访问权限：
1.联系support@adaptyvbio.com
2. 在 alpha/beta 期间请求 API 访问
3. 接收您的个人访问令牌

安全地存储您的 API 密钥：
- 使用环境变量：`ADAPTYV_API_KEY`
- 切勿将 API 密钥提交给版本控制
- 使用带有 `.gitignore` 的 `.env` 文件进行本地开发

## 端点

### 实验

#### 创建实验

提交蛋白质序列进行实验测试。

**端点：** `POST /experiments`

**请求正文：**
<<<代码块_2>>>

**序列格式：**
- 带标题的 FASTA 格式
- 支持多个序列
- 标准氨基酸代码

**回应：**
<<<代码块_3>>>

#### 获取实验状态

检查实验的当前状态。

**端点：** `GET /experiments/{experiment_id}`

**回应：**
<<<代码块_4>>>

**状态值：**
- `submitted` - 接收实验并排队
- `processing` - 正在进行主动测试
- `completed` - 结果可供下载
- `failed` - 实验遇到错误

#### 列出实验

检索您组织的所有实验。

**端点：** `GET /experiments`

**查询参数：**
- `status` - 按状态过滤（可选）
- `limit` - 每页结果数（默认值：50）
- `offset` - 分页偏移量（默认值：0）

**回应：**
<<<代码块_5>>>

### 结果

#### 获取实验结果

下载已完成实验的结果。

**端点：** `GET /experiments/{experiment_id}/results`

**回应：**
<<<代码块_6>>>

### 目标

#### 搜索目标目录

搜索 ACROBiosystems 抗原目录。

**端点：** `GET /targets`

**查询参数：**
- `search` - 搜索词（蛋白质名称、UniProt ID 等）
- `species` - 按物种过滤
- `category` - 按类别过滤

**回应：**
```json
{
  "targets": [
    {
      "target_id": "tgt_12345",
      "name": "Human PD-L1",
      "species": "Homo sapiens",
      "uniprot_id": "Q9NZQ7",
      "availability": "in_stock|custom_order",
      "price_usd": 450
    }
  ]
}
```

#### 请求自定义目标

请求标准目录中没有的抗原。

**端点：** `POST /targets/request`

**请求正文：**
```json
{
  "target_name": "Custom target name",
  "uniprot_id": "optional_uniprot_id",
  "species": "species_name",
  "notes": "Additional requirements"
}
```

### 组织

#### 获取积分余额

检查您组织的信用余额和使用情况。

**端点：** `GET /organization/credits`

**回应：**
```json
{
  "balance": 10000,
  "currency": "USD",
  "usage_this_month": 2500,
  "experiments_remaining": 22
}
```

## 网络钩子

配置 Webhook URL 以在实验完成时接收通知。

**Webhook 有效负载：**
```json
{
  "event": "experiment.completed",
  "experiment_id": "exp_abc123xyz",
  "status": "completed",
  "timestamp": "2025-12-15T10:00:00Z",
  "results_url": "/experiments/exp_abc123xyz/results"
}
```

**Webhook 事件：**
- `experiment.submitted` - 收到实验
- `experiment.started` - 处理开始
- `experiment.completed` - 结果可用
- `experiment.failed` - 发生错误

**安全：**
- 验证 webhook 签名（在入职期间提供的详细信息）
- 仅使用 HTTPS 端点
- 回复 200 OK 以确认收到

## 错误处理

**错误响应格式：**
```json
{
  "error": {
    "code": "invalid_sequence",
    "message": "Sequence contains invalid amino acid codes",
    "details": {
      "sequence_id": "protein1",
      "position": 45,
      "character": "X"
    }
  }
}
```

**常见错误代码：**
- `authentication_failed` - API 密钥无效或丢失
- `invalid_sequence` - FASTA 格式错误或氨基酸无效
- `insufficient_credits` - 实验积分不足
- `target_not_found` - 指定的目标 ID 不存在
- `rate_limit_exceeded` - 请求过多
- `experiment_not_found` - 实验 ID 无效
- `internal_error` - 服务器端错误

## 速率限制

- 每个 API 密钥每分钟 100 个请求
- 每个组织每天进行 1000 次实验
- 鼓励批量提交以进行大规模测试

当速率受限时，响应包括：
```
HTTP 429 Too Many Requests
Retry-After: 60
```

## 最佳实践

1. **使用 webhooks** 进行长时间运行的实验而不是轮询
2. **提交多个变体时的批量序列**
3. **缓存结果**以避免冗余的API调用
4. **使用指数退避实现重试逻辑**
5. **监控学分**以避免实验失败
6. **提交前在本地验证序列**
7. **使用描述性元数据**以更好地跟踪实验

## API 版本控制

该 API 目前处于 alpha/beta 阶段。可能会发生重大变化，但会是：
- 通过电子邮件向注册用户公布
- 记录在变更日志中
- 支持迁移指南

当前版本反映在响应标头中：
```
X-API-Version: alpha-2025-11
```

## 支持

对于 API 问题或疑问：
- 电子邮件：support@adaptyvbio.com
- 文档更新：https://docs.adaptyvbio.com
- 使用实验 ID 和请求详细信息报告错误