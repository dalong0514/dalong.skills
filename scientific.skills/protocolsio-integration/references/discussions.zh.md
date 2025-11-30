<!-- 此文件由机器翻译自 discussions.md -->

# 讨论 API

## 概述

讨论 API 支持对协议进行协作评论。可以在协议级别和单个步骤级别添加注释，并支持线程回复、编辑和删除。

## 基本网址

所有讨论端点都使用基本 URL：`https://protocols.io/api/v3`

## 协议级注释

### 列出协议注释

检索协议的所有注释。

**端点：** `GET /protocols/{protocol_id}/comments`

**路径参数：**
- `protocol_id`：协议的唯一标识符

**查询参数：**
- `page_size`：每页结果数（默认：10，最大：50）
- `page_id`：分页的页码（从0开始）

**回复包括：**
- 评论ID和内容
- 作者信息（姓名、单位、头像）
- 时间戳（创建和修改）
- 回复计数和线程结构

### 创建协议注释

向协议添加新注释。

**端点：** `POST /protocols/{protocol_id}/comments`

**请求正文：**
- `body`（必需）：注释文本（支持 HTML 或 Markdown）
- `parent_comment_id`（可选）：线程回复的父评论 ID

**请求示例：**
```bash
curl -X POST \
  -H "Authorization: Bearer YOUR_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "body": "This protocol worked excellently for our CRISPR experiments. We achieved 85% editing efficiency."
  }' \
  "https://protocols.io/api/v3/protocols/12345/comments"
```

### 创建线索回复

要回复现有评论，请包含父评论 ID：

<<<代码块_1>>>

### 更新评论

编辑您自己的评论。

**端点：** `PATCH /protocols/{protocol_id}/comments/{comment_id}`

**请求正文：**
- `body`（必需）：更新的评论文本

**授权**：只有评论作者可以编辑自己的评论

### 删除评论

删除一条评论。

**端点：** `DELETE /protocols/{protocol_id}/comments/{comment_id}`

**授权**：只有评论作者可以删除自己的评论

**注意**：删除父评论可能会影响整个线程，具体取决于 API 实现

## 步骤级注释

### 列出步骤评论

检索特定协议步骤的所有注释。

**端点：** `GET /protocols/{protocol_id}/steps/{step_id}/comments`

**路径参数：**
- `protocol_id`：协议的唯一标识符
- `step_id`：步骤的唯一标识符

**查询参数：**
- `page_size`：每页结果数
- `page_id`：分页的页码

### 创建步骤评论

为特定步骤添加注释。

**端点：** `POST /protocols/{protocol_id}/steps/{step_id}/comments`

**请求正文：**
- `body`（必需）：注释文本
- `parent_comment_id`（可选）：回复的父评论 ID

**请求示例：**
<<<代码块_2>>>

### 更新步骤评论

**端点：** `PATCH /protocols/{protocol_id}/steps/{step_id}/comments/{comment_id}`

**请求正文：**
- `body`（必需）：更新的评论文本

### 删除步骤注释

**端点：** `DELETE /protocols/{protocol_id}/steps/{step_id}/comments/{comment_id}`

## 常见用例

### 1.讨论线索分析

分析围绕协议的讨论：

1. 检索协议注释：`GET /protocols/{id}/comments`
2. 对于每个步骤，检索特定于步骤的注释
3. 使用 `parent_comment_id` 构建讨论线程树
4. 分析反馈模式和常见问题

### 2. 协作协议改进

要收集有关协议的反馈：

1. 发布协议
2. 监控新评论：`GET /protocols/{id}/comments`
3. 通过线索式回复来回答问题
4. 根据反馈更新协议
5. 发布更新版本并附上致谢贡献者的注释

### 3. 社区参与

与协议用户互动：

1. 设置对协议新评论的监控
2. 及时回复问题
3. 使用步骤级注释来提供详细说明
4.为复杂主题创建线程讨论

### 4. 协议故障排除

记录故障排除经验：

1. 确定方案中存在问题的步骤
2.针对遇到的具体问题添加步骤级注释
3. 记录解决方案或解决方法
4. 与遇到类似问题的其他用户创建讨论线程

## 评论格式

评论支持富文本格式：

- **HTML**：使用标准 HTML 标签进行格式化
- **Markdown**：使用 Markdown 语法进行更简单的格式化
- **链接**：包括相关资源或出版物的 URL
- **提及**：引用其他用户（格式可能有所不同）

**Markdown 示例：**
<<<代码块_3>>>

## 最佳实践

1. **具体**：在评论步骤时，参考具体参数或条件
2. **提供背景**：包括相关实验细节（细胞类型、试剂批次、设备）
3. **使用步骤级注释**：在适当的情况下直接反馈到特定步骤而不是协议级
4. **建设性参与**：及时回复问题和反馈
5. **更新协议**：将经过验证的反馈纳入协议更新中
6. **与主题相关的讨论**：使用回复功能将相关评论放在一起
7. **文档变体**：分享您手中有效的协议修改

## 权限和隐私

- **公共协议**：任何人都可以对已发布的公共协议发表评论
- **私有协议**：只有具有访问权限的协作者才能发表评论
- **评论所有权**：只有评论作者可以编辑或删除他们的评论
- **审核**：协议作者可能有额外的审核能力

## 错误处理

常见错误响应：

- `400 Bad Request`：注释格式无效或缺少必填字段
- `401 Unauthorized`：访问令牌丢失或无效
- `403 Forbidden`：权限不足（例如，尝试编辑其他用户的评论）
- `404 Not Found`：未找到协议、步骤或注释
- `429 Too Many Requests`：超出速率限制

## 通知

评论可能会触发通知：

- 协议作者收到新评论的通知
- 评论作者会收到回复通知
- 用户可以在其帐户设置中管理通知首选项