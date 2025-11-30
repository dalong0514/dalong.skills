<!-- 此文件由机器翻译自 file_manager.md -->

# 文件管理器 API

## 概述

文件管理器 API 支持在 Protocols.io 工作区中进行文件操作，包括上传文件、组织文件夹、搜索内容和管理文件生命周期。这对于将数据文件、图像、文档和其他资源附加到协议非常有用。

## 基本网址

所有文件管理器端点都使用基本 URL：`https://protocols.io/api/v3`

## 搜索和浏览

### 搜索工作区文件

在工作区中搜索文件和文件夹。

**端点：** `GET /workspaces/{workspace_id}/files/search`

**路径参数：**
- `workspace_id`：工作区的唯一标识符

**查询参数：**
- `query`：搜索关键字（搜索文件名和元数据）
- `type`：按类型过滤
  - `file`：仅限文件
  - `folder`：仅限文件夹
  - `all`：文件和文件夹（默认）
- `folder_id`：将搜索限制到特定文件夹
- `page_size`：每页结果数（默认值：20，最大：100）
- `page_id`：分页的页码（从0开始）

**回复包括：**
- 文件/文件夹 ID 和名称
- 文件大小和类型
- 创建和修改日期
- 工作区中的文件路径
- 下载 URL（文件）

**请求示例：**
```bash
curl -H "Authorization: Bearer YOUR_TOKEN" \
  "https://protocols.io/api/v3/workspaces/12345/files/search?query=microscopy&type=file"
```

### 列出文件夹内容

浏览特定文件夹中的文件和文件夹。

**端点：** `GET /workspaces/{workspace_id}/folders/{folder_id}`

**路径参数：**
- `workspace_id`：工作区的唯一标识符
- `folder_id`：文件夹的唯一标识符（使用 `root` 作为工作区根目录）

**查询参数：**
- `order_by`：排序字段（`name`、`size`、`created`、`modified`）
- `order_dir`：排序方向（`asc`、`desc`）
- `page_size`：每页结果数
- `page_id`：分页的页码

**请求示例：**
<<<代码块_1>>>

## 文件上传

### 上传文件

将文件上传到工作区文件夹。

**端点：** `POST /workspaces/{workspace_id}/files/upload`

**请求格式**：`multipart/form-data`

**表单参数：**
- `file`（必需）：要上传的文件
- `folder_id`：目标文件夹 ID（省略或使用 `root` 作为工作区根目录）
- `name`：自定义文件名（可选，如果省略则使用原始文件名）
- `description`：文件描述
- `tags`：逗号分隔的标签

**请求示例：**
<<<代码块_2>>>

### 上传验证

上传后，验证文件是否已正确处理。

**端点：** `GET /workspaces/{workspace_id}/files/{file_id}/status`

**回复包括：**
- 上传状态（`processing`、`complete`、`failed`）
- 文件元数据
- 任何处理错误

## 文件操作

### 下载文件

从工作区下载文件。

**端点：** `GET /workspaces/{workspace_id}/files/{file_id}/download`

**路径参数：**
- `workspace_id`：工作区的唯一标识符
- `file_id`：文件的唯一标识符

**响应**：具有适当 Content-Type 标头的二进制文件数据

**请求示例：**
<<<代码块_3>>>

### 获取文件元数据

无需下载即可检索文件信息。

**端点：** `GET /workspaces/{workspace_id}/files/{file_id}`

**回复包括：**
- 文件名、大小和类型
- 上传日期和作者
- 描述和标签
- 文件路径和位置
- 下载网址
- 共享权限

### 更新文件元数据

更新文件描述、标签或其他元数据。

**端点：** `PATCH /workspaces/{workspace_id}/files/{file_id}`

**请求正文：**
- `name`：新文件名
- `description`：更新了描述
- `tags`：更新的标签（以逗号分隔）
- `folder_id`：移动到不同的文件夹

**请求示例：**
<<<代码块_4>>>

### 删除文件

将文件移至垃圾箱（软删除）。

**端点：** `DELETE /workspaces/{workspace_id}/files/{file_id}`

**注意**：已删除的文件可以在有限的时间内从垃圾箱中恢复

### 恢复文件

从垃圾箱中恢复已删除的文件。

**端点：** `POST /workspaces/{workspace_id}/files/{file_id}/restore`

## 文件夹操作

### 创建文件夹

在工作区中创建一个新文件夹。

**端点：** `POST /workspaces/{workspace_id}/folders`

**请求正文：**
- `name`（必需）：文件夹名称
- `parent_folder_id`：父文件夹 ID（对于工作区根目录省略）
- `description`：文件夹说明

**请求示例：**
<<<代码块_5>>>

### 重命名文件夹

**端点：** `PATCH /workspaces/{workspace_id}/folders/{folder_id}`

**请求正文：**
- `name`：新文件夹名称
- `description`：更新了描述

### 删除文件夹
删除文件夹及其内容（可选）。

**端点：** `DELETE /workspaces/{workspace_id}/folders/{folder_id}`

**查询参数：**
- `recursive`：设置为`true`以删除文件夹和所有内容（默认值：`false`）

**警告**：递归删除无法轻易撤消

## 常见用例

### 1. 协议数据附件

将实验数据文件附加到协议中：

1.上传数据文件：`POST /workspaces/{id}/files/upload`
2.验证上传完成
3. 协议步骤中参考文件 ID
4. 在协议描述中包含下载链接

**工作流程示例：**
<<<代码块_6>>>

### 2. 工作区组织

将文件组织成逻辑文件夹结构：

1. 创建文件夹层次结构：`POST /workspaces/{id}/folders`
2.上传文件到相应的文件夹
3. 使用一致的命名约定
4. 标记文件以便于搜索

**结构示例：**
```
Workspace Root
├── Protocols
│   ├── Published
│   └── Drafts
├── Data
│   ├── Raw
│   └── Processed
├── Images
│   ├── Microscopy
│   └── Gels
└── Documents
    ├── Papers
    └── Presentations
```

### 3. 文件搜索和发现

跨工作区查找文件：

1. 按关键字搜索：`GET /workspaces/{id}/files/search?query=keywords`
2. 按类型和日期过滤
3.下载相关文件
4.更新元数据以更好地组织

### 4.批量文件上传

上传多个相关文件：

1.创建目标文件夹
2. 对于每个文件：
   - 上传文件
   - 验证上传状态
   - 添加一致的标签
3. 创建列出所有上传的索引或清单文件

### 5.数据备份与导出

导出工作区文件进行备份：

1. 列出所有文件夹：`GET /workspaces/{id}/folders/root`
2. 对于每个文件夹，列出文件
3. 下载所有文件：`GET /workspaces/{id}/files/{file_id}/download`
4.本地维护文件夹结构
5.单独存储元数据以供恢复

### 6. 文件版本控制

手动管理文件版本：

1. 上传带有版本名称的新版本（例如，`data_v2.csv`）
2.更新先前版本元数据以指示被取代
3. 在文件夹结构中维护版本历史记录
4. 参考协议中的具体版本

## 支持的文件类型

Protocols.io 支持各种文件类型：

**数据文件：**
- 电子表格：`.xlsx`、`.xls`、`.csv`、`.tsv`
- 统计数据：`.rds`、`.rdata`、`.sav`、`.dta`
- 纯文本：`.txt`、`.log`、`.json`、`.xml`

**图片：**
- 常见格式：`.jpg`、`.jpeg`、`.png`、`.gif`、`.bmp`、`.tif`、`.tiff`
- 科学：`.czi`、`.nd2`、`.lsm`（可能需要特殊处理）

**文件：**
- PDF：`.pdf`
- 字：`.docx`、`.doc`
- PowerPoint：`.pptx`、`.ppt`

**代码和脚本：**
- Python：`.py`、`.ipynb`
- R：`.r`，`.rmd`
- 外壳：`.sh`、`.bash`

**多媒体：**
- 视频：`.mp4`、`.avi`、`.mov`
- 音频：`.mp3`、`.wav`

**档案：**
- 压缩：`.zip`、`.tar.gz`、`.7z`

**文件大小限制：**
- 标准文件：检查工作区限制（通常为 100 MB - 1 GB）
- 大文件：可能需要分块上传或特殊处理

## 最佳实践

1. **文件命名**
   - 使用描述性的、一致的命名约定
   - 包括 ISO 格式的日期 (YYYY-MM-DD)
   - 避免特殊字符和空格（使用下划线）
   - 示例：`experiment_results_2025-10-26.csv`

2. **组织**
   - 创建逻辑文件夹层次结构
   - 将相关文件分组在一起
   - 将原始数据与处理结果分开
   - 将特定于协议的文件保存在专用文件夹中

3. **元数据**
   - 添加详细描述
   - 一致地标记文件
   - 包含版本信息
   - 文件处理步骤

4. **存储管理**
   - 定期审查并归档旧文件
   - 删除不必要的重复项
   - 压缩大型数据集
   - 监控工作区存储限制

5. **合作**
   - 为团队成员使用清晰的文件名
   - 描述中记录文件用途
   - 保持一致的文件夹结构
   - 传达重大组织变革

6. **安全**
   - 避免在没有适当权限的情况下上传敏感数据
   - 注意工作区可见性设置
   - 使用适当的访问控制
   - 定期审核文件访问

## 错误处理

常见错误响应：

- `400 Bad Request`：文件格式或参数无效
- `401 Unauthorized`：访问令牌丢失或无效
- `403 Forbidden`：工作区权限不足
- `404 Not Found`：未找到文件或文件夹
- `413 Payload Too Large`：文件超出大小限制
- `422 Unprocessable Entity`：文件验证失败
- `429 Too Many Requests`：超出速率限制
- `507 Insufficient Storage`：达到工作区存储限制

## 性能考虑因素

1. **大文件**
   - 考虑对 > 100 MB 的文件进行分块上传
   - 对大型数据集使用压缩
   - 如果可能的话，在非高峰时段上传

2. **批量操作**
   - 实现上传失败的重试逻辑
   - 使用指数退避来限制速率
   - 尽可能并行处理上传

3. **下载优化**
   - 在本地缓存经常访问的文件
   - 使用流媒体下载大文件
   - 实现中断下载的恢复功能