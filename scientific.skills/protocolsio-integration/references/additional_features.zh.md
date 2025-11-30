<!-- 此文件由机器翻译自 additional_features.md -->

# 附加功能

## 概述

本文档涵盖了其他protocols.io API功能，包括用户配置文件、最近发布的协议、实验记录和通知。

## 基本网址

所有端点都使用基本 URL：`https://protocols.io/api/v3`

## 用户档案管理

### 获取用户个人资料

检索经过身份验证的用户的个人资料信息。

**端点：** `GET /profile`

**回复包括：**
- 用户 ID 和用户名
- 全名
- 电子邮件地址
- 隶属关系/机构
- 简介和描述
- 个人资料图片 URL
- 帐户创建日期
- 协议计数和统计

**请求示例：**
```bash
curl -H "Authorization: Bearer YOUR_TOKEN" \
  "https://protocols.io/api/v3/profile"
```

### 更新用户个人资料

更新个人资料信息。

**端点：** `PATCH /profile`

**请求正文：**
- `first_name`：名字
- `last_name`：姓氏
- `email`：电子邮件地址
- `affiliation`：机构或组织
- `bio`：个人资料简介/描述
- `location`：地理位置
- `website`：个人或实验室网站 URL
- `twitter`：Twitter 句柄
- `orcid`：ORCID 标识符

**请求示例：**
<<<代码块_1>>>

### 上传个人资料图片

更新个人资料图片。

**端点：** `POST /profile/image`

**请求格式**：`multipart/form-data`

**表单参数：**
- `image`（必需）：图像文件（JPEG、PNG）

**推荐规格：**
- 最小尺寸：200x200 像素
- 纵横比：正方形 (1:1)
- 格式：JPEG 或 PNG
- 最大文件大小：5 MB

## 最近发布的协议

### 查询已发布的协议

发现最近发布的公共协议。

**端点：** `GET /publications`

**查询参数：**
- `key`：搜索关键字
- `category`：按类别过滤
  - 类别示例：`molecular-biology`、`cell-biology`、`biochemistry` 等。
- `date_from`：开始日期（ISO 8601 格式：YYYY-MM-DD）
- `date_to`：结束日期
- `order_field`：排序字段（`published_on`、`title`、`views`）
- `order_dir`：排序方向（`desc`、`asc`）
- `page_size`：每页结果数（默认：10，最大：50）
- `page_id`：分页的页码

**请求示例：**
<<<代码块_2>>>

**使用案例：**
- 发现趋势协议
- 监控您所在领域的新出版物
- 查找最近发布的特定技术的协议
- 跟踪值得引用的协议

## 实验记录

### 概述

实验记录允许用户记录方案的单独运行或执行，跟踪哪些有效、哪些无效以及所做的任何修改。

### 创建实验记录

记录协议的执行情况。

**端点：** `POST /protocols/{protocol_id}/runs`

**路径参数：**
- `protocol_id`：协议的唯一标识符

**请求正文：**
- `title`（必需）：实验运行标题
- `date`：实验执行日期（ISO 8601 格式）
- `status`：实验结果
  - `success`：实验成功
  - `partial`：部分成功
  - `failed`：实验失败
- `notes`：有关实验运行的详细说明
- `modifications`：协议修改或偏差
- `results`：结果摘要
- `attachments`：数据文件或图像的文件 ID

**请求示例：**
<<<代码块_3>>>

### 列出实验记录

检索方案的所有实验记录。

**端点：** `GET /protocols/{protocol_id}/runs`

**查询参数：**
- `status`：按结果过滤（`success`、`partial`、`failed`）
- `date_from`：开始日期
- `date_to`：结束日期
- `page_size`：每页结果数
- `page_id`：分页的页码

### 更新实验记录

**端点：** `PATCH /protocols/{protocol_id}/runs/{run_id}`

**Request Body**：参数与create相同，均为可选

### 删除实验记录

**端点：** `DELETE /protocols/{protocol_id}/runs/{run_id}`

**使用案例：**
- 跟踪多个实验的再现性
- 文档故障排除和优化
- 与合作者分享成功的修改
- 建立机构知识库
- 支持实验室笔记本要求

## 通知

### 获取用户通知

检索经过身份验证的用户的通知。

**端点：** `GET /notifications`

**查询参数：**
- `type`：按通知类型过滤
  - `comment`：对您的协议的新评论
  - `mention`：评论中提到了您
  - `protocol_update`：您遵循的协议已更新
  - `workspace`：工作区活动
  - `publication`：协议已发布
- `read`：按读取状态过滤
  - `true`：仅读取通知
  - `false`：仅未读通知
  - 忽略所有通知
- `page_size`：每页结果数（默认：20，最大：100）
- `page_id`：分页的页码

**回复包括：**
- 通知ID和类型
- 消息/描述
- 相关协议/评论/工作区
- 时间戳
- 读取状态

**请求示例：**
<<<代码块_4>>>

### 将通知标记为已读

**端点：** `PATCH /notifications/{notification_id}`

**请求正文：**
- `read`：设置为`true`

### 将所有通知标记为已读

**端点：** `POST /notifications/mark-all-read`

### 删除通知

**端点：** `DELETE /notifications/{notification_id}`

## 组织管理

### 导出组织数据

导出组织中的所有协议和工作区数据。

**端点：** `GET /organizations/{organization_id}/export`

**路径参数：**
- `organization_id`：组织的唯一标识符

**查询参数：**
- `format`：导出格式
  - `json`：具有完整元数据的 JSON 格式
  - `csv`：用于电子表格导入的 CSV 格式
  - `xml`：XML 格式
- `include_files`：包括关联文件 (`true`/`false`)
- `include_comments`：包括讨论 (`true`/`false`)

**响应**：导出包的下载 URL

**使用案例：**
- 机构档案
- 合规和审计要求
- 迁移到其他系统
- 备份和灾难恢复
- 数据分析和报告

**请求示例：**
<<<代码块_5>>>

## 常见集成模式

### 1. 协议发现和导入

构建协议发现工作流程：

<<<代码块_6>>>

### 2. 实验跟踪

跟踪所有协议执行：

1. 在实验室执行协议
2、文档执行：`POST /protocols/{id}/runs`
3. 将结果文件上传到工作区
4. 实验记录中的链接文件
5. 分析运行的成功率

### 3.通知系统集成

构建自定义通知系统：

1. 轮询新通知：`GET /notifications?read=false`
2. 处理每种通知类型
3.发送至内部通讯系统
4. 标记为已读：`PATCH /notifications/{id}`

### 4.配置文件同步

保持配置文件跨系统同步：

1. 检索配置文件：`GET /profile`
2、与内部系统比较
3.更新差异
4. 同步个人资料图像和元数据

## API 响应格式

### 标准响应结构

大多数 API 响应都遵循以下结构：

```json
{
  "status_code": 0,
  "status_message": "Success",
  "item": { /* single item data */ },
  "items": [ /* array of items */ ],
  "pagination": {
    "current_page": 0,
    "total_pages": 5,
    "page_size": 10,
    "total_items": 42
  }
}
```

### 错误响应结构

```json
{
  "status_code": 400,
  "status_message": "Bad Request",
  "error_message": "Missing required parameter: title",
  "error_details": {
    "field": "title",
    "issue": "required"
  }
}
```

## 最佳实践

1. **个人资料完整性**
   - 填写所有个人资料字段
   - 添加 ORCID 以进行研究归因
   - 保持隶属关系最新

2. **实验文档**
   - 记录所有协议执行情况
   - 包括成功和失败
   - 记下所有修改
   - 附上相关数据文件

3. **通知管理**
   - 定期查看通知
   - 启用相关通知类型
   - 禁用您不需要的通知类型
   - 及时回复评论

4. **出版物发现**
   - 为您的研究领域设置定期搜索
   - 关注您所在领域的多产作者
   - 为有用的协议添加书签
   - 在出版物中引用协议

5. **数据导出**
   - 定期导出组织数据
   - 测试恢复过程
   - 安全地存储出口
   - 文件导出程序